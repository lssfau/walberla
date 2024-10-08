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
//! \file NonUniformGridGPU.cpp
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
#include "gpu/ErrorChecking.h"
#include "gpu/FieldCopy.h"
#include "gpu/HostFieldAllocator.h"
#include "gpu/ParallelStreams.h"
#include "gpu/communication/NonUniformGPUScheme.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include <cmath>

#include "GridGeneration.h"
#include "LdcSetup.h"
#include "NonUniformGridGPUInfoHeader.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/gpu/AddToStorage.h"
#include "lbm_generated/gpu/BasicRecursiveTimeStepGPU.h"
#include "lbm_generated/gpu/GPUPdfField.h"
#include "lbm_generated/gpu/NonuniformGeneratedGPUPdfPackInfo.h"
using namespace walberla;

using StorageSpecification_T = lbm::NonUniformGridGPUStorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;

using PdfField_T           = lbm_generated::PdfField< StorageSpecification_T >;
using GPUPdfField_T        = lbm_generated::GPUPdfField< StorageSpecification_T >;
using FlagField_T          = FlagField< uint8_t >;
using BoundaryCollection_T = lbm::NonUniformGridGPUBoundaryCollection< FlagField_T >;

using SweepCollection_T = lbm::NonUniformGridGPUSweepCollection;

using gpu::communication::NonUniformGPUScheme;

int main(int argc, char** argv)
{
   const mpi::Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();
   gpu::selectDeviceBasedOnMpiRank();

   const std::string input_filename(argv[1]);
   const bool inputIsPython = string_ends_with(input_filename, ".py");

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                        SETUP AND CONFIGURATION                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto config = *cfg;
      logging::configureLogging(config);
      auto domainSetup        = config->getOneBlock("DomainSetup");
      auto blockForestSetup   = config->getOneBlock("SetupBlockForest");
      const bool writeSetupForestAndReturn = blockForestSetup.getParameter< bool >("writeSetupForestAndReturn", true);

      Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");
      // Reading parameters
      auto parameters   = config->getOneBlock("Parameters");
      const real_t omega             = parameters.getParameter< real_t >("omega", real_c(1.95));
      const uint_t timesteps         = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const bool gpuEnabledMPI       = parameters.getParameter< bool >("gpuEnabledMPI", false);

      shared_ptr< BlockForest > bfs;
      createBlockForest(bfs, domainSetup, blockForestSetup);

      if (writeSetupForestAndReturn && mpi::MPIManager::instance()->numProcesses() == 1)
      {
         WALBERLA_LOG_INFO_ON_ROOT("BlockForest has been created and writen to file. Returning program")
         return EXIT_SUCCESS;
      }

      auto blocks =
         std::make_shared< StructuredBlockForest >(bfs, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
      blocks->createCellBoundingBoxes();

      WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << blocks->getNumberOfBlocks() << " on " << blocks->getNumberOfLevels() << " refinement levels")
      for (uint_t level = 0; level < blocks->getNumberOfLevels(); level++)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << blocks->getNumberOfBlocks(level))
      }

      WALBERLA_LOG_INFO_ON_ROOT("Start field allocation")
      // Creating fields
      const StorageSpecification_T StorageSpec = StorageSpecification_T();
      auto allocator                           = make_shared< gpu::HostFieldAllocator< real_t > >();
      const BlockDataID pdfFieldCpuID =
         lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, uint_c(2), field::fzyx, allocator);
      const BlockDataID velFieldCpuID =
         field::addToStorage< VelocityField_T >(blocks, "vel", real_c(0.0), field::fzyx, uint_c(2), allocator);
      const BlockDataID densityFieldCpuID =
         field::addToStorage< ScalarField_T >(blocks, "density", real_c(1.0), field::fzyx, uint_c(2), allocator);
      const BlockDataID flagFieldID =
         field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field", uint_c(3));

      const BlockDataID pdfFieldGpuID =
         lbm_generated::addGPUPdfFieldToStorage< PdfField_T >(blocks, pdfFieldCpuID, StorageSpec, "pdfs on GPU", true);
      const BlockDataID velFieldGpuID =
         gpu::addGPUFieldToStorage< VelocityField_T >(blocks, velFieldCpuID, "velocity on GPU", true);
      const BlockDataID densityFieldGpuID =
         gpu::addGPUFieldToStorage< ScalarField_T >(blocks, densityFieldCpuID, "velocity on GPU", true);
      WALBERLA_LOG_INFO_ON_ROOT("Finished field allocation")

      const Cell innerOuterSplit =
         Cell(parameters.getParameter< Vector3< cell_idx_t > >("innerOuterSplit", Vector3< cell_idx_t >(1, 1, 1)));
      Vector3< int32_t > gpuBlockSize =
         parameters.getParameter< Vector3< int32_t > >("gpuBlockSize", Vector3< int32_t >(256, 1, 1));
      SweepCollection_T sweepCollection(blocks, pdfFieldGpuID, densityFieldGpuID, velFieldGpuID, gpuBlockSize[0],
                                        gpuBlockSize[1], gpuBlockSize[2], omega, innerOuterSplit);

      for (auto& iBlock : *blocks){
         sweepCollection.initialise(&iBlock, cell_idx_c(1));
      }
      sweepCollection.initialiseBlockPointer();
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                      LB SWEEPS AND BOUNDARY HANDLING                                       ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto ldc = std::make_shared< LDC >(blocks->getDepth());

      const FlagUID fluidFlagUID("Fluid");
      ldc->setupBoundaryFlagField(*blocks, flagFieldID);
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID, 0);
      BoundaryCollection_T boundaryCollection(blocks, flagFieldID, pdfFieldGpuID, fluidFlagUID);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                           COMMUNICATION SCHEME                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      WALBERLA_LOG_INFO_ON_ROOT("Setting up communication...")
      auto communication = std::make_shared< NonUniformGPUScheme< CommunicationStencil_T > >(blocks, gpuEnabledMPI);
      auto packInfo      = lbm_generated::setupNonuniformGPUPdfCommunication< GPUPdfField_T >(blocks, pdfFieldGpuID);
      communication->addPackInfo(packInfo);
      WALBERLA_MPI_BARRIER()

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                          TIME STEP DEFINITIONS                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      int streamHighPriority = 0;
      int streamLowPriority  = 0;
      WALBERLA_GPU_CHECK(gpuDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority))
      sweepCollection.setOuterPriority(streamHighPriority);
      auto defaultStream = gpu::StreamRAII::newPriorityStream(streamLowPriority);

      lbm_generated::BasicRecursiveTimeStepGPU< GPUPdfField_T, SweepCollection_T, BoundaryCollection_T >
         LBMMeshRefinement(blocks, pdfFieldGpuID, sweepCollection, boundaryCollection, communication, packInfo);
      SweepTimeloop timeLoop(blocks->getBlockStorage(), timesteps);

      LBMMeshRefinement.addRefinementToTimeLoop(timeLoop, uint_c(0));

      // VTK
      const uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      const bool useVTKAMRWriter     = parameters.getParameter< bool >("useVTKAMRWriter", false);
      const bool oneFilePerProcess   = parameters.getParameter< bool >("oneFilePerProcess", false);

      auto finalDomain = blocks->getDomain();
      if (vtkWriteFrequency > 0){
         auto vtkOutput =
            vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out", "simulation_step",
                                           false, true, true, false, 0, useVTKAMRWriter, oneFilePerProcess);
         auto velWriter = make_shared< field::VTKWriter< VelocityField_T, float32 > >(velFieldCpuID, "vel");
         vtkOutput->addCellDataWriter(velWriter);

         if (parameters.getParameter< bool >("writeOnlySlice", true)){
            const AABB sliceXY(finalDomain.xMin(), finalDomain.yMin(), finalDomain.center()[2] - blocks->dz(blocks->getDepth()),
                               finalDomain.xMax(), finalDomain.yMax(), finalDomain.center()[2] + blocks->dz(blocks->getDepth()));
            vtkOutput->addCellInclusionFilter(vtk::AABBCellFilter(sliceXY));
         }

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

      auto remainingTimeLoggerFrequency =
         parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(-1.0)); // in seconds
      if (remainingTimeLoggerFrequency > 0){
         auto logger = timing::RemainingTimeLogger(timeLoop.getNrOfTimeSteps(), remainingTimeLoggerFrequency);
         timeLoop.addFuncAfterTimeStep(logger, "remaining time logger");
      }

      lbm_generated::PerformanceEvaluation< FlagField_T > const performance(blocks, flagFieldID, fluidFlagUID);
      field::CellCounter< FlagField_T > fluidCells(blocks, flagFieldID, fluidFlagUID);
      fluidCells();

      WALBERLA_LOG_INFO_ON_ROOT("Non uniform Grid benchmark with " << fluidCells.numberOfCells()
                                                                   << " fluid cells (in total on all levels)")

      WcTimingPool timeloopTiming;
      WcTimer simTimer;

      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
      WALBERLA_MPI_BARRIER()

      simTimer.start();
      timeLoop.run(timeloopTiming);
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
      WALBERLA_MPI_BARRIER()
      simTimer.end();

      WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
      double time = simTimer.max();
      WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
      performance.logResultOnRoot(timesteps, time);

      const auto reducedTimeloopTiming = timeloopTiming.getReduced();
      WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)

      WALBERLA_ROOT_SECTION()
      {
         if (inputIsPython)
         {
            python_coupling::PythonCallback pythonCallbackResults("results_callback");
            if (pythonCallbackResults.isCallable())
            {
               pythonCallbackResults.data().exposeValue("numProcesses", lbm_generated::PerformanceEvaluation< FlagField_T >::processes());
               pythonCallbackResults.data().exposeValue("numThreads", performance.threads());
               pythonCallbackResults.data().exposeValue("numCores", performance.cores());
               pythonCallbackResults.data().exposeValue("numberOfCells", performance.numberOfCells());
               pythonCallbackResults.data().exposeValue("numberOfFluidCells", performance.numberOfFluidCells());
               pythonCallbackResults.data().exposeValue("mlups", performance.mlups(timesteps, time));
               pythonCallbackResults.data().exposeValue("mlupsPerCore", performance.mlupsPerCore(timesteps, time));
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