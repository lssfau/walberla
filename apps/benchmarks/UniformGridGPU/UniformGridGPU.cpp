#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "cuda/AddGPUFieldToStorage.h"
#include "cuda/DeviceSelectMPI.h"
#include "cuda/ParallelStreams.h"
#include "cuda/communication/UniformGPUScheme.h"
#include "cuda/FieldCopy.h"
#include "cuda/lbm/CombinedInPlaceGpuPackInfo.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "lbm/inplace_streaming/TimestepTracker.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/SweepTimeloop.h"

#include "InitShearVelocity.h"

#include <cmath>

#include "UniformGridGPU_InfoHeader.h"
using namespace walberla;

using FlagField_T            = FlagField<uint8_t>;

int main(int argc, char** argv)
{
   mpi::Environment env(argc, argv);
   cuda::selectDeviceBasedOnMpiRank();

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      WALBERLA_CUDA_CHECK(cudaPeekAtLastError())

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                        SETUP AND CONFIGURATION                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto config = *cfg;
      logging::configureLogging(config);
      auto blocks = blockforest::createUniformBlockGridFromConfig(config);

      Vector3< uint_t > cellsPerBlock =
         config->getBlock("DomainSetup").getParameter< Vector3< uint_t > >("cellsPerBlock");
      // Reading parameters
      auto parameters        = config->getOneBlock("Parameters");
      const real_t omega     = parameters.getParameter< real_t >("omega", real_c(1.4));
      const uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const bool initShearFlow = parameters.getParameter<bool>("initShearFlow", true);

      // Creating fields
      BlockDataID pdfFieldCpuID =
         field::addToStorage< PdfField_T >(blocks, "pdfs cpu", real_t(std::nan("")), field::fzyx);
      BlockDataID velFieldCpuID = field::addToStorage< VelocityField_T >(blocks, "vel", real_t(0), field::fzyx);

      // Initialize velocity on cpu
      if( initShearFlow ){
         WALBERLA_LOG_INFO_ON_ROOT("Initializing shear flow")
         initShearVelocity(blocks, velFieldCpuID);
      }

      BlockDataID pdfFieldGpuID = cuda::addGPUFieldToStorage< PdfField_T >(blocks, pdfFieldCpuID, "pdfs on GPU", true);
      // Velocity field is copied to the GPU
      BlockDataID velFieldGpuID =
         cuda::addGPUFieldToStorage< VelocityField_T >(blocks, velFieldCpuID, "velocity on GPU", true);

      pystencils::UniformGridGPU_MacroSetter setterSweep(pdfFieldGpuID, velFieldGpuID);

      // Set up initial PDF values
      for (auto& block : *blocks)
         setterSweep(&block);

      Vector3< int > innerOuterSplit =
         parameters.getParameter< Vector3< int > >("innerOuterSplit", Vector3< int >(1, 1, 1));

      for (uint_t i = 0; i < 3; ++i)
      {
         if (int_c(cellsPerBlock[i]) <= innerOuterSplit[i] * 2)
         {
            WALBERLA_ABORT_NO_DEBUG_INFO("innerOuterSplit too large - make it smaller or increase cellsPerBlock")
         }
      }

      Cell innerOuterSplitCell(innerOuterSplit[0], innerOuterSplit[1], innerOuterSplit[2]);
      bool cudaEnabledMPI = parameters.getParameter< bool >("cudaEnabledMPI", false);
      Vector3< int32_t > gpuBlockSize =
         parameters.getParameter< Vector3< int32_t > >("gpuBlockSize", Vector3< int32_t >(256, 1, 1));

      int streamHighPriority = 0;
      int streamLowPriority  = 0;
      WALBERLA_CUDA_CHECK(cudaDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority))

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                      LB SWEEPS AND BOUNDARY HANDLING                                       ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      using LbSweep      = lbm::UniformGridGPU_LbKernel;
      using PackInfoEven = lbm::UniformGridGPU_PackInfoEven;
      using PackInfoOdd  = lbm::UniformGridGPU_PackInfoOdd;
      using cuda::communication::UniformGPUScheme;

      LbSweep lbSweep(pdfFieldGpuID, omega, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], innerOuterSplitCell);
      lbSweep.setOuterPriority(streamHighPriority);

      // Boundaries
      const FlagUID fluidFlagUID( "Fluid" );
      BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>(blocks, "Boundary Flag Field");
      auto boundariesConfig = config->getBlock( "Boundaries" );
      bool boundaries = false;
      if( boundariesConfig )
      {
         boundaries = true;
         geometry::initBoundaryHandling< FlagField_T >( *blocks, flagFieldID, boundariesConfig );
         geometry::setNonBoundaryCellsToDomain< FlagField_T >( *blocks, flagFieldID, fluidFlagUID );
      }

      lbm::UniformGridGPU_NoSlip noSlip(blocks, pdfFieldGpuID);
      noSlip.fillFromFlagField<FlagField_T>(blocks, flagFieldID, FlagUID("NoSlip"), fluidFlagUID);

      lbm::UniformGridGPU_UBB ubb(blocks, pdfFieldGpuID);
      ubb.fillFromFlagField<FlagField_T>(blocks, flagFieldID, FlagUID("UBB"), fluidFlagUID);

      // Initial setup is the post-collision state of an even time step
      auto tracker = make_shared< lbm::TimestepTracker >(0);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                           COMMUNICATION SCHEME                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      UniformGPUScheme< Stencil_T > comm(blocks, cudaEnabledMPI);
      auto packInfo = make_shared< lbm::CombinedInPlaceGpuPackInfo< PackInfoEven, PackInfoOdd > >(tracker, pdfFieldGpuID);
      comm.addPackInfo(packInfo);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                          TIME STEP DEFINITIONS                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto defaultStream = cuda::StreamRAII::newPriorityStream(streamLowPriority);

      auto boundarySweep = [&](IBlock * block, uint8_t t, cudaStream_t stream){
         noSlip.run(block, t, stream);
         ubb.run(block, t, stream);
      };

      auto boundaryInner = [&](IBlock * block, uint8_t t, cudaStream_t stream){
         noSlip.inner(block, t, stream);
         ubb.inner(block, t, stream);
      };

      auto boundaryOuter = [&](IBlock * block, uint8_t t, cudaStream_t stream){
         noSlip.outer(block, t, stream);
         ubb.outer(block, t, stream);
      };

      auto simpleOverlapTimeStep = [&]() {
         // Communicate post-collision values of previous timestep...
         comm.startCommunication(defaultStream);
         for (auto& block : *blocks){
            if(boundaries) boundaryInner(&block, tracker->getCounter(), defaultStream);
            lbSweep.inner(&block, tracker->getCounterPlusOne(), defaultStream);
         }
         comm.wait(defaultStream);
         for (auto& block : *blocks){
            if(boundaries) boundaryOuter(&block, tracker->getCounter(), defaultStream);
            lbSweep.outer(&block, tracker->getCounterPlusOne(), defaultStream);
         }

         tracker->advance();
      };

      auto normalTimeStep = [&]() {
         comm.communicate(defaultStream);
         for (auto& block : *blocks){
            if(boundaries) boundarySweep(&block, tracker->getCounter(), defaultStream);
            lbSweep(&block, tracker->getCounterPlusOne(), defaultStream);
         }

         tracker->advance();
      };

      // With two-fields patterns, ghost layer cells act as constant stream-in boundaries;
      // with in-place patterns, ghost layer cells act as wet-node no-slip boundaries.
      auto kernelOnlyFunc = [&]() {
         tracker->advance();
         for (auto& block : *blocks)
            lbSweep(&block, tracker->getCounter(), defaultStream);
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                             TIME LOOP SETUP                                                ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      SweepTimeloop timeLoop(blocks->getBlockStorage(), timesteps);

      const std::string timeStepStrategy = parameters.getParameter< std::string >("timeStepStrategy", "normal");
      std::function< void() > timeStep;
      if (timeStepStrategy == "noOverlap")
         timeStep = std::function< void() >(normalTimeStep);
      else if (timeStepStrategy == "simpleOverlap")
         timeStep = simpleOverlapTimeStep;
      else if (timeStepStrategy == "kernelOnly")
      {
         WALBERLA_LOG_INFO_ON_ROOT(
            "Running only compute kernel without boundary - this makes only sense for benchmarking!")
         // Run initial communication once to provide any missing stream-in populations
         comm.communicate();
         timeStep = kernelOnlyFunc;
      }
      else
      {
         WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'. Allowed values are 'noOverlap', "
                                      "'simpleOverlap', 'kernelOnly'")
      }

      timeLoop.add() << BeforeFunction(timeStep) << Sweep([](IBlock*) {}, "time step");

      // VTK
      uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
         auto velWriter = make_shared< field::VTKWriter< VelocityField_T > >(velFieldCpuID, "vel");
         vtkOutput->addCellDataWriter(velWriter);

         vtkOutput->addBeforeFunction([&]() {
           cuda::fieldCpy< VelocityField_T , cuda::GPUField< real_t > >(blocks, velFieldCpuID, velFieldGpuID);
         });
         timeLoop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                               BENCHMARK                                                    ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      int warmupSteps     = parameters.getParameter< int >("warmupSteps", 2);
      int outerIterations = parameters.getParameter< int >("outerIterations", 1);
      for (int i = 0; i < warmupSteps; ++i)
         timeLoop.singleStep();

      double remainingTimeLoggerFrequency =
         parameters.getParameter< double >("remainingTimeLoggerFrequency", -1.0); // in seconds
      if (remainingTimeLoggerFrequency > 0)
      {
         auto logger = timing::RemainingTimeLogger(timeLoop.getNrOfTimeSteps() * uint_c(outerIterations),
                                                   remainingTimeLoggerFrequency);
         timeLoop.addFuncAfterTimeStep(logger, "remaining time logger");
      }

      for (int outerIteration = 0; outerIteration < outerIterations; ++outerIteration)
      {
         WALBERLA_CUDA_CHECK(cudaPeekAtLastError())

         timeLoop.setCurrentTimeStepToZero();
         WcTimer simTimer;
         cudaDeviceSynchronize();
         WALBERLA_CUDA_CHECK(cudaPeekAtLastError())
         WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
         simTimer.start();
         timeLoop.run();
         cudaDeviceSynchronize();
         simTimer.end();
         WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
         auto time      = simTimer.last();
         auto nrOfCells = real_c(cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2]);

         auto mlupsPerProcess = nrOfCells * real_c(timesteps) / time * 1e-6;
         WALBERLA_LOG_RESULT_ON_ROOT("MLUPS per process " << mlupsPerProcess)
         WALBERLA_LOG_RESULT_ON_ROOT("Time per time step " << time / real_c(timesteps))
         WALBERLA_ROOT_SECTION()
         {
            python_coupling::PythonCallback pythonCallbackResults("results_callback");
            if (pythonCallbackResults.isCallable())
            {
               pythonCallbackResults.data().exposeValue("mlupsPerProcess", mlupsPerProcess);
               pythonCallbackResults.data().exposeValue("stencil", infoStencil);
               pythonCallbackResults.data().exposeValue("streamingPattern", infoStreamingPattern);
               pythonCallbackResults.data().exposeValue("collisionSetup", infoCollisionSetup);
               pythonCallbackResults.data().exposeValue("cse_global", infoCseGlobal);
               pythonCallbackResults.data().exposeValue("cse_pdfs", infoCsePdfs);
               // Call Python function to report results
               pythonCallbackResults();
            }
         }
      }
   }

   return 0;
}
