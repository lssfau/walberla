#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/OpenMP.h"
#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"
#include "python_coupling/DictWrapper.h"
#include "blockforest/Initialization.h"
#include "field/vtk/VTKWriter.h"
#include "field/AddToStorage.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/communication/UniformDirectScheme.h"
#include "timeloop/all.h"
#include "core/timing/TimingPool.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/waLBerlaBuildInfo.h"
#include "domain_decomposition/SharedSweep.h"
#include "gui/Gui.h"
#include "InitShearVelocity.h"
#include "ManualKernels.h"
#include "lbm/lattice_model/D3Q19.h"

#include "GenDefines.h"
#include "GenMacroGetter.h"
#include "GenMacroSetter.h"

#include "GenLbKernel.h"
#include "GenLbKernelAAEven.h"
#include "GenLbKernelAAOdd.h"

#include "GenPackInfo.h"
#include "GenPackInfoAAPush.h"
#include "GenPackInfoAAPull.h"
#include "GenMpiDtypeInfo.h"
#include "GenMpiDtypeInfoAAPull.h"
#include "GenMpiDtypeInfoAAPush.h"


#include <iomanip>

using namespace walberla;

using PdfField_T = GhostLayerField< real_t, Stencil_T::Q >;
using VelocityField_T = GhostLayerField< real_t, 3 >;


int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );

   for( auto cfg = python_coupling::configBegin( argc, argv ); cfg != python_coupling::configEnd(); ++cfg )
   {
      WALBERLA_MPI_WORLD_BARRIER()

      auto config = *cfg;
      logging::configureLogging( config );
      auto blocks = blockforest::createUniformBlockGridFromConfig( config );

      Vector3<uint_t> cellsPerBlock = config->getBlock( "DomainSetup" ).getParameter<Vector3<uint_t>  >( "cellsPerBlock" );
      // Reading parameters
      auto parameters = config->getOneBlock( "Parameters" );
      const std::string timeStepMode = parameters.getParameter<std::string>( "timeStepMode", "twoField");
      const real_t omega = parameters.getParameter<real_t>( "omega", real_c( 1.4 ));
            uint_t timesteps = parameters.getParameter<uint_t>( "timesteps", uint_c( 60 ));
      const real_t shearVelocityMagnitude = parameters.getParameter<real_t>("shearVelocityMagnitude", real_t(0.08));
      const bool directComm = parameters.getParameter<bool>("directComm", false);

      auto pdfFieldAdder = [](IBlock* const block, StructuredBlockStorage * const storage) {
          return new PdfField_T(storage->getNumberOfXCells(*block),
                                storage->getNumberOfYCells(*block),
                                storage->getNumberOfZCells(*block),
                                uint_t(1),
                                field::fzyx,
                                make_shared<field::AllocateAligned<real_t, 64>>());
      };

      // Creating fields
      BlockDataID pdfFieldId = blocks->addStructuredBlockData<PdfField_T>(pdfFieldAdder, "pdfs");
      BlockDataID velFieldId = field::addToStorage< VelocityField_T >( blocks, "vel", real_t( 0 ), field::fzyx );

      pystencils::GenMacroSetter setterKernel(pdfFieldId, velFieldId);
      pystencils::GenMacroGetter getterKernel(pdfFieldId, velFieldId);

      if( shearVelocityMagnitude > 0 )
          initShearVelocity(blocks, velFieldId, shearVelocityMagnitude);
      for( auto & b : *blocks)
          setterKernel(&b);

      // Buffered Comm
      blockforest::communication::UniformBufferedScheme< Stencil_T > twoFieldComm(blocks );
      twoFieldComm.addPackInfo(make_shared< pystencils::GenPackInfo >(pdfFieldId ) );

      blockforest::communication::UniformBufferedScheme< Stencil_T > aaPullComm(blocks);
      aaPullComm.addPackInfo(make_shared< pystencils::GenPackInfoAAPull>(pdfFieldId));

      blockforest::communication::UniformBufferedScheme< Stencil_T > aaPushComm(blocks);
      aaPushComm.addPackInfo(make_shared< pystencils::GenPackInfoAAPush>(pdfFieldId));

      // Direct Comm
      blockforest::communication::UniformDirectScheme< Stencil_T > twoFieldCommDirect(blocks);
      twoFieldCommDirect.addDataToCommunicate(make_shared<pystencils::GenMpiDtypeInfo>(pdfFieldId));

      blockforest::communication::UniformDirectScheme< Stencil_T > aaPullCommDirect(blocks);
      aaPullCommDirect.addDataToCommunicate(make_shared<pystencils::GenMpiDtypeInfoAAPull>(pdfFieldId));

      blockforest::communication::UniformDirectScheme< Stencil_T > aaPushCommDirect(blocks);
      aaPushCommDirect.addDataToCommunicate(make_shared<pystencils::GenMpiDtypeInfoAAPush>(pdfFieldId));


      const std::string twoFieldKernelType = parameters.getParameter<std::string>( "twoFieldKernelType", "generated");
      std::function<void(IBlock*)> twoFieldKernel;
      if( twoFieldKernelType == "generated") {
          twoFieldKernel = pystencils::GenLbKernel(pdfFieldId, omega);
      } else if (twoFieldKernelType == "manualGeneric") {
          using MyLM = lbm::D3Q19<lbm::collision_model::SRT>;
          BlockDataID tmpPdfFieldId = blocks->addStructuredBlockData<PdfField_T>(pdfFieldAdder, "pdfs");
          twoFieldKernel = StreamPullCollideGeneric<MyLM>(pdfFieldId, tmpPdfFieldId, omega);
      } else if (twoFieldKernelType == "manualD3Q19") {
          using MyLM = lbm::D3Q19<lbm::collision_model::SRT>;
          BlockDataID tmpPdfFieldId = blocks->addStructuredBlockData<PdfField_T>(pdfFieldAdder, "pdfs");
          twoFieldKernel = StreamPullCollideD3Q19<MyLM>(pdfFieldId, tmpPdfFieldId, omega);
      } else {
          WALBERLA_ABORT_NO_DEBUG_INFO("Invalid option for \"twoFieldKernelType\", "
                                       "valid options are \"generated\", \"manualGeneric\", \"manualD3Q19\"")
      }

      using F = std::function<void()>;
      SweepTimeloop timeLoop( blocks->getBlockStorage(), timesteps / 2 );
      if( timeStepMode == "twoField")
      {
          timeLoop.add() << BeforeFunction(directComm ? F(twoFieldCommDirect) : F(twoFieldComm), "communication" )
                         << Sweep( twoFieldKernel, "LB stream & collide1" );
          timeLoop.add() << BeforeFunction(directComm ? F(twoFieldCommDirect) : F(twoFieldComm), "communication" )
                         << Sweep( twoFieldKernel, "LB stream & collide2" );

      } else if ( timeStepMode == "twoFieldKernelOnly") {
          timeLoop.add() << Sweep( pystencils::GenLbKernel(pdfFieldId, omega), "LB stream & collide1" );
          timeLoop.add() << Sweep( pystencils::GenLbKernel(pdfFieldId, omega), "LB stream & collide2" );
      } else if ( timeStepMode == "aa") {
          timeLoop.add() << Sweep( pystencils::GenLbKernelAAEven(pdfFieldId, omega), "AA Even" );
          timeLoop.add() << BeforeFunction( directComm ? F(aaPullCommDirect) : F(aaPullComm) )
                         << Sweep( pystencils::GenLbKernelAAOdd(pdfFieldId, omega), "AA Odd")
                         << AfterFunction( directComm ? F(aaPushCommDirect) : F(aaPushComm) );
      } else if ( timeStepMode == "aaKernelOnly") {
          timeLoop.add() << Sweep( pystencils::GenLbKernelAAEven(pdfFieldId, omega), "AA Even" );
          timeLoop.add() << Sweep( pystencils::GenLbKernelAAOdd(pdfFieldId, omega), "AA Odd");
      } else {
          WALBERLA_ABORT("Invalid value for timeStepMode")
      }


      int warmupSteps = parameters.getParameter<int>( "warmupSteps", 2 );
      int outerIterations = parameters.getParameter<int>( "outerIterations", 1 );
      for(int i=0; i < warmupSteps; ++i )
         timeLoop.singleStep();

      auto remainingTimeLoggerFrequency = parameters.getParameter< double >( "remainingTimeLoggerFrequency", -1.0 ); // in seconds
      if (remainingTimeLoggerFrequency > 0) {
          auto logger = timing::RemainingTimeLogger( timeLoop.getNrOfTimeSteps() * uint_c( outerIterations ), remainingTimeLoggerFrequency );
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

              auto threads = omp_get_max_threads();

              simTimer.start();
              timeLoop.run();
              simTimer.end();
              auto time = simTimer.last();
              auto nrOfCells = real_c( cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2] );
              auto mlupsPerProcess = nrOfCells * real_c( timesteps ) / time * 1e-6;

              using std::setw;
              WALBERLA_LOG_INFO_ON_ROOT(setw(18) << timeStepMode <<
                                                     "  procs: " << setw(6) << MPIManager::instance()->numProcesses() <<
                                                     "  threads: " << threads <<
                                                     "  direct_comm: " << directComm <<
                                                     "  time steps: " << timesteps <<
                                                     setw(15) << "  block size: " << cellsPerBlock <<
                                                     "  mlups/core:  " << int(mlupsPerProcess/ threads) <<
                                                     "  mlups:  " << int(mlupsPerProcess) *  MPIManager::instance()->numProcesses())

              WALBERLA_ROOT_SECTION()
              {
                  python_coupling::PythonCallback pythonCallbackResults( "results_callback" );
                  if ( pythonCallbackResults.isCallable())
                  {
                      pythonCallbackResults.data().exposeValue( "mlupsPerProcess", mlupsPerProcess );
                      pythonCallbackResults.data().exposeValue( "stencil", infoStencil );
                      pythonCallbackResults.data().exposeValue( "configName", infoConfigName );
                      pythonCallbackResults.data().exposeValue( "timeStepMode", timeStepMode );
                      pythonCallbackResults.data().exposeValue( "twoFieldKernel", twoFieldKernelType );
                      pythonCallbackResults.data().exposeValue( "optimizations", optimizationDict );
                      pythonCallbackResults.data().exposeValue( "githash", core::buildinfo::gitSHA1() );
                      pythonCallbackResults.data().exposeValue( "compilerFlags", core::buildinfo::compilerFlags() );
                      pythonCallbackResults.data().exposeValue( "buildMachine", core::buildinfo::buildMachine() );

                      // Call Python function to report results
                      pythonCallbackResults();
                  }
              }
          }
      }
   }

   return 0;
}
