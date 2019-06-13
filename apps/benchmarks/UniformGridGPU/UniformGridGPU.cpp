#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Random.h"
#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"
#include "python_coupling/DictWrapper.h"
#include "blockforest/Initialization.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/AddToStorage.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/vtk/VTKOutput.h"
#include "lbm/PerformanceLogger.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "timeloop/all.h"
#include "core/math/Random.h"
#include "geometry/all.h"
#include "cuda/HostFieldAllocator.h"
#include "cuda/communication/GPUPackInfo.h"
#include "cuda/ParallelStreams.h"
#include "cuda/NVTX.h"
#include "core/timing/TimingPool.h"
#include "core/timing/RemainingTimeLogger.h"
#include "cuda/AddGPUFieldToStorage.h"
#include "cuda/communication/UniformGPUScheme.h"
#include "cuda/DeviceSelectMPI.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "domain_decomposition/SharedSweep.h"

#include "UniformGridGPU_LatticeModel.h"
#include "UniformGridGPU_LbKernel.h"
#include "UniformGridGPU_PackInfo.h"
#include "UniformGridGPU_UBB.h"
#include "UniformGridGPU_NoSlip.h"
#include "UniformGridGPU_Communication.h"


using namespace walberla;

using LatticeModel_T = lbm::UniformGridGPU_LatticeModel;

using Stencil_T = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;
using CommScheme_T = cuda::communication::UniformGPUScheme<CommunicationStencil_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;


void initShearVelocity(const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID pdfFieldID,
                       const real_t xMagnitude=0.1, const real_t fluctuationMagnitude=0.05 )
{
    math::seedRandomGenerator(0);
    auto halfZ = blocks->getDomainCellBB().zMax() / 2;
    for( auto & block: *blocks)
    {
        PdfField_T * pdfField = block.getData<PdfField_T>( pdfFieldID );
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(pdfField,
            Cell globalCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
            real_t randomReal = xMagnitude * math::realRandom<real_t>(-fluctuationMagnitude, fluctuationMagnitude);
            if( globalCell[2] >= halfZ ) {
                pdfField->setDensityAndVelocity(x, y, z, Vector3<real_t>(xMagnitude, 0, randomReal), real_t(1.0));
            } else {
                pdfField->setDensityAndVelocity(x, y, z, Vector3<real_t>(-xMagnitude, 0,randomReal), real_t(1.0));
            }
        );
    }
}


int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );
   cuda::selectDeviceBasedOnMpiRank();

   for( auto cfg = python_coupling::configBegin( argc, argv ); cfg != python_coupling::configEnd(); ++cfg )
   {
      WALBERLA_MPI_WORLD_BARRIER();

      auto config = *cfg;
      logging::configureLogging( config );
      auto blocks = blockforest::createUniformBlockGridFromConfig( config );

      Vector3<uint_t> cellsPerBlock = config->getBlock( "DomainSetup" ).getParameter<Vector3<uint_t>  >( "cellsPerBlock" );
      // Reading parameters
      auto parameters = config->getOneBlock( "Parameters" );
      const real_t omega = parameters.getParameter<real_t>( "omega", real_c( 1.4 ));
      const uint_t timesteps = parameters.getParameter<uint_t>( "timesteps", uint_c( 50 ));
      const Vector3<real_t> initialVelocity = parameters.getParameter< Vector3<real_t> >( "initialVelocity", Vector3<real_t>() );
      const bool initShearFlow = parameters.getParameter<bool>("initShearFlow", false);

      // Creating fields
      auto latticeModel = LatticeModel_T( omega );
      BlockDataID pdfFieldCpuID = lbm::addPdfFieldToStorage( blocks, "pdfs on CPU", latticeModel, initialVelocity, real_t(1), field::fzyx );
      if( initShearFlow ) {
          WALBERLA_LOG_INFO_ON_ROOT("Initializing shear flow");
          initShearVelocity( blocks, pdfFieldCpuID );
      }

      BlockDataID pdfFieldGpuID = cuda::addGPUFieldToStorage<PdfField_T >( blocks, pdfFieldCpuID, "pdfs on GPU", true );
      BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );


      // Boundaries
      const FlagUID fluidFlagUID( "Fluid" );
      auto boundariesConfig = config->getBlock( "Boundaries" );
      bool disableBoundaries = true;
      if( boundariesConfig )
      {
          disableBoundaries = false;
          geometry::initBoundaryHandling< FlagField_T >( *blocks, flagFieldID, boundariesConfig );
          geometry::setNonBoundaryCellsToDomain< FlagField_T >( *blocks, flagFieldID, fluidFlagUID );
      }

      lbm::UniformGridGPU_UBB ubb(blocks, pdfFieldGpuID);
      lbm::UniformGridGPU_NoSlip noSlip(blocks, pdfFieldGpuID);

      ubb.fillFromFlagField<FlagField_T>( blocks, flagFieldID, FlagUID("UBB"), fluidFlagUID );
      noSlip.fillFromFlagField<FlagField_T>( blocks, flagFieldID, FlagUID("NoSlip"), fluidFlagUID );

       // Communication setup
      bool cudaEnabledMPI = parameters.getParameter<bool>( "cudaEnabledMPI", false );
      Vector3<int32_t> gpuBlockSize = parameters.getParameter<Vector3<int32_t> > ("gpuBlockSize", Vector3<int32_t>(256, 1, 1));
      const std::string communicationSchemeStr = parameters.getParameter<std::string>("communicationScheme", "UniformGPUScheme_Baseline");
      CommunicationSchemeType communicationScheme;
      if( communicationSchemeStr == "GPUPackInfo_Baseline")
          communicationScheme = GPUPackInfo_Baseline;
      else if (communicationSchemeStr == "GPUPackInfo_Streams")
          communicationScheme = GPUPackInfo_Streams;
      else if (communicationSchemeStr == "UniformGPUScheme_Baseline")
          communicationScheme = UniformGPUScheme_Baseline;
      else if (communicationSchemeStr == "UniformGPUScheme_Memcpy")
          communicationScheme = UniformGPUScheme_Memcpy;
      else if (communicationSchemeStr == "MPIDatatypes")
          communicationScheme = MPIDatatypes;
      else {
          WALBERLA_ABORT_NO_DEBUG_INFO("Invalid choice for communicationScheme")
      }

      Vector3<int> innerOuterSplit = parameters.getParameter<Vector3<int> >("innerOuterSplit", Vector3<int>(1, 1, 1));


      int streamHighPriority = 0;
      int streamLowPriority = 0;
      WALBERLA_CUDA_CHECK( cudaDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority) );
      WALBERLA_CHECK(gpuBlockSize[2] == 1);
      pystencils::UniformGridGPU_LbKernel lbKernel( pdfFieldGpuID, omega, gpuBlockSize[0], gpuBlockSize[1],
                                                    Cell(innerOuterSplit[0], innerOuterSplit[1], innerOuterSplit[2]) );
      lbKernel.setOuterPriority( streamHighPriority );
      UniformGridGPU_Communication< CommunicationStencil_T, cuda::GPUField< double > >
         gpuComm( blocks, pdfFieldGpuID, (CommunicationSchemeType) communicationScheme, cudaEnabledMPI );

      auto defaultStream = cuda::StreamRAII::newPriorityStream( streamLowPriority );
      auto innerOuterStreams = cuda::ParallelStreams( streamHighPriority );
      auto boundaryOuterStreams = cuda::ParallelStreams( streamHighPriority );
      auto boundaryInnerStreams = cuda::ParallelStreams( streamHighPriority );

      uint_t currentTimeStep = 0;

      auto simpleOverlapTimeStep = [&] ()
      {
          gpuComm.startCommunication(defaultStream);
          for( auto &block: *blocks )
              lbKernel.inner( &block, defaultStream );
          gpuComm.wait(defaultStream);
          for( auto &block: *blocks )
              lbKernel.outer( &block, defaultStream );
      };

      auto overlapTimeStep = [&]()
      {
         cuda::NvtxRange namedRange("timestep");
         auto innerOuterSection = innerOuterStreams.parallelSection( defaultStream );

         innerOuterSection.run([&]( auto innerStream )
         {
            cuda::nameStream(innerStream, "inner stream");
            for( auto &block: *blocks )
            {
               if(!disableBoundaries)
               {
                  auto p = boundaryInnerStreams.parallelSection( innerStream );
                  p.run( [&block, &ubb]( cudaStream_t s ) { ubb.inner( &block, s ); } );
                  p.run( [&block, &noSlip]( cudaStream_t s ) { noSlip.inner( &block, s ); } );
               }
               lbKernel.inner( &block, innerStream );
            }
         });

         innerOuterSection.run([&]( auto outerStream )
         {
            cuda::nameStream(outerStream, "outer stream");
            gpuComm( outerStream );

            for( auto &block: *blocks )
            {
               if(!disableBoundaries)
               {
                  auto p = boundaryOuterStreams.parallelSection( outerStream );
                  p.run( [&block, &ubb]( cudaStream_t s ) { ubb.outer( &block, s ); } );
                  p.run( [&block, &noSlip]( cudaStream_t s ) { noSlip.outer( &block, s ); } );
               }
               lbKernel.outer( &block, outerStream );
            }
         });
         currentTimeStep += 1;
      };


      auto boundaryStreams = cuda::ParallelStreams( streamHighPriority );
      auto normalTimeStep = [&]()
      {
         gpuComm();
         for( auto &block: *blocks )
         {
            if(!disableBoundaries)
            {
               auto p = boundaryStreams.parallelSection( defaultStream );
               p.run( [&block, &ubb]( cudaStream_t s ) { ubb( &block, s ); } );
               p.run( [&block, &noSlip]( cudaStream_t s ) { noSlip( &block, s ); } );
            }
            lbKernel( &block );
         }
      };

      auto kernelOnlyFunc = [&] ()
      {
          for( auto &block: *blocks )
              lbKernel( &block );
      };

      SweepTimeloop timeLoop( blocks->getBlockStorage(), timesteps );

      const std::string timeStepStrategy = parameters.getParameter<std::string>( "timeStepStrategy", "normal");
      std::function<void()> timeStep;
      if (timeStepStrategy == "noOverlap")
          timeStep = std::function<void()>( normalTimeStep );
      else if (timeStepStrategy == "complexOverlap")
          timeStep = std::function<void()>( overlapTimeStep );
      else if (timeStepStrategy == "simpleOverlap")
          timeStep = simpleOverlapTimeStep;
      else if (timeStepStrategy == "kernelOnly") {
          WALBERLA_LOG_INFO_ON_ROOT("Running only compute kernel without boundary - this makes only sense for benchmarking!")
          timeStep = kernelOnlyFunc;
      }
      else {
          WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'. Allowed values are 'noOverlap', 'complexOverlap', 'simpleOverlap', 'kernelOnly'");
      }

      timeLoop.add() << BeforeFunction( timeStep  )
                     << Sweep( []( IBlock * ) {}, "time step" );

      // VTK
      uint_t vtkWriteFrequency = parameters.getParameter<uint_t>( "vtkWriteFrequency", 0 );
      if( vtkWriteFrequency > 0 )
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData( *blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                          "simulation_step", false, true, true, false, 0 );
         vtkOutput->addCellDataWriter(
                 make_shared<lbm::VelocityVTKWriter<LatticeModel_T> >( pdfFieldCpuID, "Velocity" ));
         vtkOutput->addCellDataWriter( make_shared<lbm::DensityVTKWriter<LatticeModel_T> >( pdfFieldCpuID, "Density" ));
         vtkOutput->addBeforeFunction(
                 cuda::fieldCpyFunctor<PdfField_T, cuda::GPUField<real_t> >( blocks, pdfFieldCpuID, pdfFieldGpuID ));
         timeLoop.addFuncAfterTimeStep( vtk::writeFiles( vtkOutput ), "VTK Output" );
      }



      int warmupSteps = parameters.getParameter<int>( "warmupSteps", 2 );
      int outerIterations = parameters.getParameter<int>( "outerIterations", 1 );
      for(int i=0; i < warmupSteps; ++i )
         timeLoop.singleStep();

      auto remainingTimeLoggerFrequency = parameters.getParameter< double >( "remainingTimeLoggerFrequency", -1.0 ); // in seconds
      if (remainingTimeLoggerFrequency > 0) {
          auto logger = timing::RemainingTimeLogger( timeLoop.getNrOfTimeSteps() * outerIterations, remainingTimeLoggerFrequency );
          timeLoop.addFuncAfterTimeStep( logger, "remaining time logger" );
      }

      for(int outerIteration = 0; outerIteration < outerIterations; ++outerIteration)
      {
          timeLoop.setCurrentTimeStepToZero();
          WcTimer simTimer;
          cudaDeviceSynchronize();
          WALBERLA_LOG_INFO_ON_ROOT( "Starting simulation with " << timesteps << " time steps" );
          simTimer.start();
          timeLoop.run();
          cudaDeviceSynchronize();
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
                  pythonCallbackResults.data().exposeValue( "githash", WALBERLA_GIT_SHA1 );
                  // Call Python function to report results
                  pythonCallbackResults();
              }
          }
      }
   }

   return 0;
}
