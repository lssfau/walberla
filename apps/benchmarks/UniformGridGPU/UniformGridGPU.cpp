//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file UniformGridGPU.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/FieldCopy.h"
#include "gpu/GPUWrapper.h"
#include "gpu/HostFieldAllocator.h"
#include "gpu/ParallelStreams.h"
#include "gpu/communication/UniformGPUScheme.h"

#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/gpu/UniformGeneratedGPUPdfPackInfo.h"
#include "lbm_generated/gpu/GPUPdfField.h"
#include "lbm_generated/gpu/AddToStorage.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/SweepTimeloop.h"

#include <cmath>

#include "InitShearVelocity.h"
#include "UniformGridGPU_InfoHeader.h"

using namespace walberla;

using StorageSpecification_T = lbm::UniformGridGPUStorageSpecification;
using Stencil_T = lbm::UniformGridGPUStorageSpecification::Stencil;

using PdfField_T = lbm_generated::PdfField< StorageSpecification_T >;
using GPUPdfField_T = lbm_generated::GPUPdfField< StorageSpecification_T >;
using FlagField_T = FlagField< uint8_t >;
using BoundaryCollection_T = lbm::UniformGridGPUBoundaryCollection< FlagField_T >;

using SweepCollection_T = lbm::UniformGridGPUSweepCollection;

using gpu::communication::UniformGPUScheme;

int main(int argc, char** argv)
{
   mpi::Environment const env(argc, argv);
   gpu::selectDeviceBasedOnMpiRank();

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      WALBERLA_GPU_CHECK(gpuPeekAtLastError())

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                        SETUP AND CONFIGURATION                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto config = *cfg;
      logging::configureLogging(config);
      auto blocks = blockforest::createUniformBlockGridFromConfig(config);

      // Reading parameters
      auto parameters          = config->getOneBlock("Parameters");
      const real_t omega       = parameters.getParameter< real_t >("omega", real_c(1.4));
      const uint_t timesteps   = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const bool initShearFlow = parameters.getParameter< bool >("initShearFlow", true);
      const bool cudaEnabledMPI = parameters.getParameter< bool >("cudaEnabledMPI", false);

      // Creating fields
      const StorageSpecification_T StorageSpec = StorageSpecification_T();
      const BlockDataID pdfFieldCpuID  = lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, uint_c(1), field::fzyx);

      auto allocator = make_shared< gpu::HostFieldAllocator<real_t> >(); // use pinned memory allocator for faster CPU-GPU memory transfers
      const BlockDataID velFieldCpuID = field::addToStorage< VelocityField_T >(blocks, "vel", real_c(0.0), field::fzyx, uint_c(1), allocator);
      const BlockDataID densityFieldCpuID = field::addToStorage< ScalarField_T >(blocks, "density", real_c(1.0), field::fzyx, uint_c(1), allocator);
      const BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field");

      // Initialize velocity on cpu
      if (initShearFlow)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Initializing shear flow")
         initShearVelocity(blocks, velFieldCpuID);
      }

      const BlockDataID pdfFieldGpuID = lbm_generated::addGPUPdfFieldToStorage< PdfField_T >(blocks, pdfFieldCpuID, StorageSpec, "pdfs on GPU", true);
      const BlockDataID velFieldGpuID =
         gpu::addGPUFieldToStorage< VelocityField_T >(blocks, velFieldCpuID, "velocity on GPU", true);
      const BlockDataID densityFieldGpuID =
         gpu::addGPUFieldToStorage< ScalarField_T >(blocks, densityFieldCpuID, "velocity on GPU", true);

      const Cell innerOuterSplit = Cell(parameters.getParameter< Vector3<cell_idx_t> >("innerOuterSplit", Vector3<cell_idx_t>(1, 1, 1)));
      Vector3< int32_t > gpuBlockSize = parameters.getParameter< Vector3< int32_t > >("gpuBlockSize", Vector3< int32_t >(256, 1, 1));
      SweepCollection_T sweepCollection(blocks, pdfFieldGpuID, densityFieldGpuID, velFieldGpuID, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], omega, innerOuterSplit);
      for (auto& block : *blocks)
      {
         sweepCollection.initialise(&block);
      }

      int streamHighPriority = 0;
      int streamLowPriority  = 0;
      WALBERLA_GPU_CHECK(gpuDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority))
      sweepCollection.setOuterPriority(streamHighPriority);
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                      LB SWEEPS AND BOUNDARY HANDLING                                       ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      const pystencils::UniformGridGPU_StreamOnlyKernel StreamOnlyKernel(pdfFieldGpuID, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);

      // Boundaries
      const FlagUID fluidFlagUID("Fluid");
      auto boundariesConfig   = config->getBlock("Boundaries");
      if (boundariesConfig)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Setting boundary conditions")
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldID, boundariesConfig);
      }
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID);
      BoundaryCollection_T boundaryCollection(blocks, flagFieldID, pdfFieldGpuID, fluidFlagUID);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                           COMMUNICATION SCHEME                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      UniformGPUScheme< Stencil_T > communication(blocks, cudaEnabledMPI);
      auto packInfo = std::make_shared<lbm_generated::UniformGeneratedGPUPdfPackInfo< GPUPdfField_T >>(pdfFieldGpuID);
      communication.addPackInfo(packInfo);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                          TIME STEP DEFINITIONS                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      SweepTimeloop timeLoop(blocks->getBlockStorage(), timesteps);
      const std::string timeStepStrategy = parameters.getParameter< std::string >("timeStepStrategy", "normal");

      auto defaultStream = gpu::StreamRAII::newPriorityStream(streamLowPriority);

      if (timeStepStrategy == "noOverlap") {
         if (boundariesConfig){
            timeLoop.add() << BeforeFunction(communication.getCommunicateFunctor(), "communication")
                           << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL, defaultStream), "Boundary Conditions");
            timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL, defaultStream), "LBM StreamCollide");
         }else {
            timeLoop.add() << BeforeFunction(communication.getCommunicateFunctor(), "communication")
                           << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL, defaultStream), "LBM StreamCollide");}

      } else if (timeStepStrategy == "simpleOverlap") {
         if (boundariesConfig){
            timeLoop.add() << BeforeFunction(communication.getStartCommunicateFunctor(), "Start Communication")
                           << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL, defaultStream), "Boundary Conditions");
            timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::INNER, defaultStream), "LBM StreamCollide Inner Frame");
            timeLoop.add() << BeforeFunction(communication.getWaitFunctor(), "Wait for Communication")
                           << Sweep(sweepCollection.streamCollide(SweepCollection_T::OUTER, defaultStream), "LBM StreamCollide Outer Frame");
         }else{
            timeLoop.add() << BeforeFunction(communication.getStartCommunicateFunctor(), "Start Communication")
                           << Sweep(sweepCollection.streamCollide(SweepCollection_T::INNER, defaultStream), "LBM StreamCollide Inner Frame");
            timeLoop.add() << BeforeFunction(communication.getWaitFunctor(), "Wait for Communication")
                           << Sweep(sweepCollection.streamCollide(SweepCollection_T::OUTER,defaultStream), "LBM StreamCollide Outer Frame");}

      } else if (timeStepStrategy == "kernelOnly") {
         timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL, defaultStream), "LBM StreamCollide");
      } else if (timeStepStrategy == "StreamOnly") {
         timeLoop.add() << Sweep(StreamOnlyKernel, "LBM Stream Only");
      } else {
         WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'")
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                             TIME LOOP SETUP                                                ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // VTK
      const uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
         auto velWriter = make_shared< field::VTKWriter< VelocityField_T, float32 > >(velFieldCpuID, "vel");
         vtkOutput->addCellDataWriter(velWriter);

         vtkOutput->addBeforeFunction([&]() {
            for (auto& block : *blocks)
               sweepCollection.calculateMacroscopicParameters(&block);
            gpu::fieldCpy< VelocityField_T, gpu::GPUField< real_t > >(blocks, velFieldCpuID, velFieldGpuID);
         });
         timeLoop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                               BENCHMARK                                                    ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      lbm_generated::PerformanceEvaluation<FlagField_T> const performance(blocks, flagFieldID, fluidFlagUID);
      const uint_t warmupSteps     = parameters.getParameter< uint_t >("warmupSteps", uint_c(2));
      const uint_t outerIterations = parameters.getParameter< uint_t >("outerIterations", uint_c(1));
      for (uint_t i = 0; i < warmupSteps; ++i)
         timeLoop.singleStep();

      auto remainingTimeLoggerFrequency =
         parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(-1.0)); // in seconds
      if (remainingTimeLoggerFrequency > 0)
      {
         auto logger = timing::RemainingTimeLogger(timeLoop.getNrOfTimeSteps() * uint_c(outerIterations),
                                                   remainingTimeLoggerFrequency);
         timeLoop.addFuncAfterTimeStep(logger, "remaining time logger");
      }

      for (uint_t outerIteration = 0; outerIteration < outerIterations; ++outerIteration)
      {
         WALBERLA_GPU_CHECK(gpuPeekAtLastError())

         timeLoop.setCurrentTimeStepToZero();
         WcTimingPool timeloopTiming;
         WcTimer simTimer;

         WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
         WALBERLA_GPU_CHECK( gpuPeekAtLastError() )
         WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
         simTimer.start();
         timeLoop.run(timeloopTiming);
         WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
         simTimer.end();

         WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
         double time = simTimer.max();
         WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
         performance.logResultOnRoot(timesteps, time);

         const auto reducedTimeloopTiming = timeloopTiming.getReduced();
         WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)

         WALBERLA_ROOT_SECTION()
         {
            python_coupling::PythonCallback pythonCallbackResults("results_callback");
            if (pythonCallbackResults.isCallable())
            {
               pythonCallbackResults.data().exposeValue("mlupsPerProcess", performance.mlupsPerProcess(timesteps, time));
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
   return EXIT_SUCCESS;
}
