#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Random.h"
#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"
#include "python_coupling/DictWrapper.h"
#include "blockforest/Initialization.h"
#include "field/vtk/VTKWriter.h"
#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "timeloop/all.h"
#include "core/timing/TimingPool.h"
#include "core/timing/RemainingTimeLogger.h"
#include "domain_decomposition/SharedSweep.h"
#include "gui/Gui.h"
#include "InitShearVelocity.h"

#include "GenDefines.h"
#include "GenPackInfo.h"
#include "GenLbKernel.h"
#include "GenMacroGetter.h"
#include "GenMacroSetter.h"

using namespace walberla;

using PdfField_T = GhostLayerField< real_t, Stencil_T::Q >;
using VelocityField_T = GhostLayerField< real_t, 3 >;


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
      const real_t shearVelocityMagnitude = parameters.getParameter<real_t>("shearVelocityMagnitude", 0.08);

      // Creating fields
      BlockDataID pdfFieldId = field::addToStorage< PdfField_T >( blocks, "pdfs", real_t( std::nan("") ), field::fzyx );
      BlockDataID velFieldId = field::addToStorage< VelocityField_T >( blocks, "vel", real_t( 0 ), field::fzyx );

      pystencils::GenMacroSetter setterKernel(pdfFieldId, velFieldId);
      pystencils::GenMacroGetter getterKernel(pdfFieldId, velFieldId);

      initShearVelocity(blocks, velFieldId, shearVelocityMagnitude);
      for( auto & b : *blocks)
          setterKernel(&b);

      SweepTimeloop timeLoop( blocks->getBlockStorage(), timesteps );
      blockforest::communication::UniformBufferedScheme< Stencil_T > communication( blocks );
      communication.addPackInfo( make_shared< pystencils::GenPackInfo >( pdfFieldId ) );

      timeLoop.add() << BeforeFunction( communication, "communication" )
                     << Sweep( pystencils::GenLbKernel(pdfFieldId, omega), "LB stream & collide" );

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
          auto velWriter = make_shared< field::VTKWriter< VelocityField_T > >( velFieldId, "vel" );
          vtkOutput->addCellDataWriter( velWriter );
          vtkOutput->addBeforeFunction( [&]()
                                        { for( auto & b : *blocks)
                                            getterKernel(&b);
                                        } );
          timeLoop.addFuncAfterTimeStep( vtk::writeFiles( vtkOutput ), "VTK Output" );
      }


      bool useGui = parameters.getParameter<bool>( "useGui", false );
      if( useGui )
      {
          GUI gui( timeLoop, blocks, argc, argv);
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
