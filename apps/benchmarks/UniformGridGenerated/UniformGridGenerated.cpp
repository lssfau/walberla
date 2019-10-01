#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Random.h"
#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"
#include "python_coupling/DictWrapper.h"
#include "blockforest/Initialization.h"
#include "field/vtk/VTKWriter.h"
#include "field/communication/PackInfo.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "timeloop/all.h"
#include "core/timing/TimingPool.h"
#include "core/timing/RemainingTimeLogger.h"
#include "domain_decomposition/SharedSweep.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/vtk/VTKOutput.h"
#include "lbm/gui/Connection.h"
#include "lbm/vtk/Velocity.h"
#include "gui/Gui.h"

#include "UniformGridGenerated_LatticeModel.h"
#include "UniformGridGenerated_Defines.h"


using namespace walberla;

typedef lbm::UniformGridGenerated_LatticeModel LatticeModel_T;
typedef LatticeModel_T::Stencil                Stencil_T;
typedef LatticeModel_T::CommunicationStencil   CommunicationStencil_T;
typedef lbm::PdfField< LatticeModel_T >        PdfField_T;


void initShearVelocity(const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID pdfFieldId,
                       const real_t xMagnitude=0.1, const real_t fluctuationMagnitude=0.05 )
{
    math::seedRandomGenerator(0);
    auto halfZ = blocks->getDomainCellBB().zMax() / 2;
    for( auto & block: *blocks)
    {
        auto pdfField = block.getData<PdfField_T>( pdfFieldId );
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(pdfField,
            Cell globalCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
            real_t randomReal = xMagnitude * math::realRandom<real_t>(-fluctuationMagnitude, fluctuationMagnitude);

            if( globalCell[2] >= halfZ ) {
                pdfField->setDensityAndVelocity(x, y, z, Vector3<real_t>(xMagnitude, 0, randomReal), real_t(1.0));
            } else {
                pdfField->setDensityAndVelocity(x, y, z, Vector3<real_t>(-xMagnitude, 0, randomReal), real_t(1.0));
            }
        );
    }
}


int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );

   for( auto cfg = python_coupling::configBegin( argc, argv ); cfg != python_coupling::configEnd(); ++cfg )
   {
      WALBERLA_MPI_WORLD_BARRIER();

      auto config = *cfg;
      logging::configureLogging( config );
      auto blocks = blockforest::createUniformBlockGridFromConfig( config );

      Vector3<uint_t> cellsPerBlock = config->getBlock( "DomainSetup" ).getParameter<Vector3<uint_t>  >( "cellsPerBlock" );
      // Reading parameters
      auto parameters = config->getOneBlock( "Parameters" );
      const std::string timeStepStrategy = parameters.getParameter<std::string>( "timeStepStrategy", "normal");
      const real_t omega = parameters.getParameter<real_t>( "omega", real_c( 1.4 ));
      const uint_t timesteps = parameters.getParameter<uint_t>( "timesteps", uint_c( 50 ));
      const bool initShearFlow = parameters.getParameter<bool>("initShearFlow", false);

      // Creating fields
      LatticeModel_T latticeModel = LatticeModel_T( omega );
      BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel);

      if( initShearFlow ) {
          initShearVelocity(blocks, pdfFieldId);
      }

      SweepTimeloop timeLoop( blocks->getBlockStorage(), timesteps );
      blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
      communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );

      timeLoop.add() << BeforeFunction( communication, "communication" )
                     << Sweep( LatticeModel_T::Sweep( pdfFieldId ), "LB stream & collide" );

      int warmupSteps = parameters.getParameter<int>( "warmupSteps", 2 );
      int outerIterations = parameters.getParameter<int>( "outerIterations", 1 );
      for(int i=0; i < warmupSteps; ++i )
         timeLoop.singleStep();

      auto remainingTimeLoggerFrequency = parameters.getParameter< double >( "remainingTimeLoggerFrequency", -1.0 ); // in seconds
      if (remainingTimeLoggerFrequency > 0) {
          auto logger = timing::RemainingTimeLogger( timeLoop.getNrOfTimeSteps() * outerIterations, remainingTimeLoggerFrequency );
          timeLoop.addFuncAfterTimeStep( logger, "remaining time logger" );
      }

      // VTK
      uint_t vtkWriteFrequency = parameters.getParameter<uint_t>( "vtkWriteFrequency", 0 );
      if( vtkWriteFrequency > 0 )
      {
          auto vtkOutput = vtk::createVTKOutput_BlockData( *blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                           "simulation_step", false, true, true, false, 0 );
          auto velWriter = make_shared< lbm::VelocityVTKWriter<LatticeModel_T> >(pdfFieldId, "vel");
          vtkOutput->addCellDataWriter(velWriter);
          timeLoop.addFuncAfterTimeStep( vtk::writeFiles( vtkOutput ), "VTK Output" );
      }


      bool useGui = parameters.getParameter<bool>( "useGui", false );
      if( useGui )
      {
          GUI gui( timeLoop, blocks, argc, argv);
          lbm::connectToGui<LatticeModel_T>(gui);
          gui.run();
      }
      else
      {
          for ( int outerIteration = 0; outerIteration < outerIterations; ++outerIteration )
          {
              timeLoop.setCurrentTimeStepToZero();
              WcTimer simTimer;
              WALBERLA_LOG_INFO_ON_ROOT( "Starting simulation with " << timesteps << " time steps" );
              simTimer.start();
              timeLoop.run();
              simTimer.end();
              WALBERLA_LOG_INFO_ON_ROOT( "Simulation finished" );
              auto time = simTimer.last();
              auto nrOfCells = real_c( cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2] );
              auto mlupsPerProcess = nrOfCells * real_c( timesteps ) / time * 1e-6;
              WALBERLA_LOG_RESULT_ON_ROOT( "MLUPS per process " << mlupsPerProcess );
              WALBERLA_LOG_RESULT_ON_ROOT( "Time per time step " << time / real_c( timesteps ));
              WALBERLA_ROOT_SECTION()
              {
                  python_coupling::PythonCallback pythonCallbackResults( "results_callback" );
                  if ( pythonCallbackResults.isCallable())
                  {
                      pythonCallbackResults.data().exposeValue( "mlupsPerProcess", mlupsPerProcess );
                      pythonCallbackResults.data().exposeValue( "stencil", infoStencil );
                      pythonCallbackResults.data().exposeValue( "configName", infoConfigName );
                      pythonCallbackResults.data().exposeValue( "cse_global", infoCseGlobal );
                      pythonCallbackResults.data().exposeValue( "cse_pdfs", infoCsePdfs );
                      // Call Python function to report results
                      pythonCallbackResults();
                  }
              }
          }
      }
   }

   return 0;
}
