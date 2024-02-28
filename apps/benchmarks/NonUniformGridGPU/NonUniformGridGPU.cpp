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
#include "blockforest/loadbalancing/StaticCurve.h"

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
#include "gpu/ErrorChecking.h"
#include "gpu/HostFieldAllocator.h"
#include "gpu/ParallelStreams.h"
#include "gpu/communication/NonUniformGPUScheme.h"

#include "lbm_generated/evaluation/PerformanceEvaluation.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/gpu/NonuniformGeneratedGPUPdfPackInfo.h"
#include "lbm_generated/gpu/GPUPdfField.h"
#include "lbm_generated/gpu/AddToStorage.h"
#include "lbm_generated/gpu/BasicRecursiveTimeStepGPU.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"

#include <cmath>

#include "LdcSetup.h"
#include "NonUniformGridGPUInfoHeader.h"
using namespace walberla;

using StorageSpecification_T = lbm::NonUniformGridGPUStorageSpecification;
using Stencil_T = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;

using PdfField_T = lbm_generated::PdfField< StorageSpecification_T >;
using GPUPdfField_T = lbm_generated::GPUPdfField< StorageSpecification_T >;
using FlagField_T = FlagField< uint8_t >;
using BoundaryCollection_T = lbm::NonUniformGridGPUBoundaryCollection< FlagField_T >;

using SweepCollection_T = lbm::NonUniformGridGPUSweepCollection;

using gpu::communication::NonUniformGPUScheme;

namespace {
void createSetupBlockForest(SetupBlockForest& setupBfs, const Config::BlockHandle& domainSetup, LDC& ldcSetup, const uint_t numProcesses=uint_c(MPIManager::instance()->numProcesses())) {
    Vector3<real_t> domainSize = domainSetup.getParameter<Vector3<real_t> >("domainSize");
    Vector3<uint_t> rootBlocks = domainSetup.getParameter<Vector3<uint_t> >("rootBlocks");
    Vector3<bool> periodic = domainSetup.getParameter<Vector3<bool> >("periodic");

    auto refSelection = ldcSetup.refinementSelector();
    setupBfs.addRefinementSelectionFunction(std::function<void(SetupBlockForest &)>(refSelection));
    const AABB domain(real_t(0.0), real_t(0.0), real_t(0.0), domainSize[0], domainSize[1], domainSize[2]);
    setupBfs.addWorkloadMemorySUIDAssignmentFunction(blockforest::uniformWorkloadAndMemoryAssignment);
    setupBfs.init(domain, rootBlocks[0], rootBlocks[1], rootBlocks[2], periodic[0], periodic[1], periodic[2]);
    setupBfs.balanceLoad(blockforest::StaticLevelwiseCurveBalanceWeighted(), numProcesses);
}
}

int main(int argc, char** argv)
{
   const mpi::Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();
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
      auto domainSetup              = config->getOneBlock("DomainSetup");
      Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");
      // Reading parameters
      auto parameters          = config->getOneBlock("Parameters");
      const real_t omega       = parameters.getParameter< real_t >("omega", real_c(1.4));
      const uint_t refinementDepth = parameters.getParameter< uint_t >("refinementDepth", uint_c(1));
      const uint_t timesteps   = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const bool cudaEnabledMPI = parameters.getParameter< bool >("cudaEnabledMPI", false);
      const bool writeSetupForestAndReturn = parameters.getParameter< bool >("writeSetupForestAndReturn", false);
      const bool benchmarkKernelOnly = parameters.getParameter< bool >("benchmarkKernelOnly", false);
      const uint_t numProcesses = parameters.getParameter< uint_t >( "numProcesses");

      auto ldc = std::make_shared< LDC >(refinementDepth );
      SetupBlockForest setupBfs;
      if (writeSetupForestAndReturn)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Creating SetupBlockForest for " << numProcesses << " processes")
         WALBERLA_LOG_INFO_ON_ROOT("Generating SetupBlockForest...")
         createSetupBlockForest(setupBfs, domainSetup, *ldc, numProcesses);

         WALBERLA_ROOT_SECTION() { setupBfs.writeVTKOutput("SetupBlockForest"); }

         WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << setupBfs.getNumberOfBlocks())
         uint_t totalCellUpdates( 0.0 );
         for (uint_t level = 0; level <= refinementDepth; level++)
         {
            const uint_t numberOfBlocks = setupBfs.getNumberOfBlocks(level);
            const uint_t numberOfCells = numberOfBlocks * cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2];
            totalCellUpdates += timesteps * math::uintPow2(level)  * numberOfCells;
            WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << numberOfBlocks)
         }
         cudaDeviceProp prop{};
         WALBERLA_GPU_CHECK(gpuGetDeviceProperties(&prop, 0))

         const uint_t totalNumberCells = setupBfs.getNumberOfBlocks() * cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2];

         const uint_t PDFsPerCell = StorageSpecification_T::inplace ? Stencil_T::Q : 2 * Stencil_T::Q;
         const uint_t valuesPerCell = (PDFsPerCell + VelocityField_T::F_SIZE + ScalarField_T::F_SIZE);
         const uint_t sizePerValue = sizeof(PdfField_T::value_type);
         const double totalGPUMem = double_c(prop.totalGlobalMem) * 1e-9;
         const double expectedMemory = double_c(totalNumberCells * valuesPerCell * sizePerValue) * 1e-9;

         WALBERLA_LOG_INFO_ON_ROOT( "Total number of cells will be " << totalNumberCells << " fluid cells (in total on all levels)")
         WALBERLA_LOG_INFO_ON_ROOT( "Expected total memory demand will be " << expectedMemory << " GB")
         WALBERLA_LOG_INFO_ON_ROOT( "The total cell updates after " << timesteps << " timesteps (on the coarse level) will be " << totalCellUpdates)
         WALBERLA_LOG_INFO_ON_ROOT( "Total GPU memory " << totalGPUMem)

         WALBERLA_LOG_INFO_ON_ROOT("Ending program")
         return EXIT_SUCCESS;
      }

      WALBERLA_LOG_INFO_ON_ROOT("Generating SetupBlockForest...")
      createSetupBlockForest(setupBfs, domainSetup, *ldc);

      // Create structured block forest
      WALBERLA_LOG_INFO_ON_ROOT("Creating structured block forest...")
      auto bfs    = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
      auto blocks = std::make_shared< StructuredBlockForest >(bfs, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
      blocks->createCellBoundingBoxes();

      WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << blocks->getNumberOfBlocks())
      for (uint_t level = 0; level <= refinementDepth; level++)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << blocks->getNumberOfBlocks(level))
      }

      WALBERLA_LOG_INFO_ON_ROOT("Start field allocation")
      // Creating fields
      const StorageSpecification_T StorageSpec = StorageSpecification_T();
      auto allocator = make_shared< gpu::HostFieldAllocator<real_t> >();
      const BlockDataID pdfFieldCpuID  = lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, uint_c(2), field::fzyx, allocator);
      const BlockDataID velFieldCpuID = field::addToStorage< VelocityField_T >(blocks, "vel", real_c(0.0), field::fzyx, uint_c(2), allocator);
      const BlockDataID densityFieldCpuID = field::addToStorage< ScalarField_T >(blocks, "density", real_c(1.0), field::fzyx, uint_c(2), allocator);
      const BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field", uint_c(3));

      const BlockDataID pdfFieldGpuID = lbm_generated::addGPUPdfFieldToStorage< PdfField_T >(blocks, pdfFieldCpuID, StorageSpec, "pdfs on GPU", true);
      const BlockDataID velFieldGpuID =
         gpu::addGPUFieldToStorage< VelocityField_T >(blocks, velFieldCpuID, "velocity on GPU", true);
      const BlockDataID densityFieldGpuID =
         gpu::addGPUFieldToStorage< ScalarField_T >(blocks, densityFieldCpuID, "velocity on GPU", true);
      WALBERLA_LOG_INFO_ON_ROOT("Finished field allocation")

      const Cell innerOuterSplit = Cell(parameters.getParameter< Vector3<cell_idx_t> >("innerOuterSplit", Vector3<cell_idx_t>(1, 1, 1)));
      Vector3< int32_t > gpuBlockSize = parameters.getParameter< Vector3< int32_t > >("gpuBlockSize", Vector3< int32_t >(256, 1, 1));
      SweepCollection_T sweepCollection(blocks, pdfFieldGpuID, densityFieldGpuID, velFieldGpuID, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], omega, innerOuterSplit);
      for (auto& iBlock : *blocks)
      {
         sweepCollection.initialise(&iBlock, 2, nullptr);
      }
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                      LB SWEEPS AND BOUNDARY HANDLING                                       ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      const FlagUID fluidFlagUID("Fluid");
      ldc->setupBoundaryFlagField(*blocks, flagFieldID);
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID, 2);
      BoundaryCollection_T boundaryCollection(blocks, flagFieldID, pdfFieldGpuID, fluidFlagUID);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                           COMMUNICATION SCHEME                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      WALBERLA_LOG_INFO_ON_ROOT("Setting up communication...")
      auto communication = std::make_shared< NonUniformGPUScheme <CommunicationStencil_T>> (blocks, cudaEnabledMPI);
      auto packInfo = lbm_generated::setupNonuniformGPUPdfCommunication<GPUPdfField_T>(blocks, pdfFieldGpuID);
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

      lbm_generated::BasicRecursiveTimeStepGPU< GPUPdfField_T, SweepCollection_T, BoundaryCollection_T > LBMMeshRefinement(blocks, pdfFieldGpuID, sweepCollection, boundaryCollection, communication, packInfo);
      SweepTimeloop timeLoop(blocks->getBlockStorage(), timesteps);

      // LBMMeshRefinement.test(5);
      // return EXIT_SUCCESS;

      if(benchmarkKernelOnly){
         timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");
      }
      else{
         LBMMeshRefinement.addRefinementToTimeLoop(timeLoop);
      }

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

      auto remainingTimeLoggerFrequency =
         parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(-1.0)); // in seconds
      if (remainingTimeLoggerFrequency > 0)
      {
         auto logger = timing::RemainingTimeLogger(timeLoop.getNrOfTimeSteps(), remainingTimeLoggerFrequency);
         timeLoop.addFuncAfterTimeStep(logger, "remaining time logger");
      }

      lbm_generated::PerformanceEvaluation<FlagField_T> const performance(blocks, flagFieldID, fluidFlagUID);
      field::CellCounter< FlagField_T > fluidCells( blocks, flagFieldID, fluidFlagUID );
      fluidCells();

      WALBERLA_LOG_INFO_ON_ROOT( "Non uniform Grid benchmark with " << fluidCells.numberOfCells() << " fluid cells (in total on all levels)")

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
   }
   return EXIT_SUCCESS;
}