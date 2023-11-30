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
//! \file NonUniformGridCPU.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/SweepTimeloop.h"

#include <cmath>

#include "LdcSetup.h"
#include "NonUniformGridCPUInfoHeader.h"
#include "lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/refinement/BasicRecursiveTimeStep.h"

using namespace walberla;

using StorageSpecification_T = lbm::NonUniformGridCPUStorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;

using PdfField_T           = lbm_generated::PdfField< StorageSpecification_T >;
using BoundaryCollection_T = lbm::NonUniformGridCPUBoundaryCollection< FlagField_T >;

using SweepCollection_T = lbm::NonUniformGridCPUSweepCollection;

using blockforest::communication::NonUniformBufferedScheme;

int main(int argc, char** argv)
{
   const mpi::Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   const std::string input_filename(argv[1]);
   const bool inputIsPython = string_ends_with(input_filename, ".py");

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_BARRIER()

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                            BLOCK FOREST SETUP                                              ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto config = *cfg;
      logging::configureLogging(config);

      auto blockForestSetup = config->getOneBlock("SetupBlockForest");
      const std::string blockForestFilestem =
         blockForestSetup.getParameter< std::string >("blockForestFilestem", "blockforest");
      const uint_t refinementDepth = blockForestSetup.getParameter< uint_t >("refinementDepth", uint_c(1));

      auto domainSetup                = config->getOneBlock("DomainSetup");
      Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");

      // Load structured block forest from file
      std::ostringstream oss;
      oss << blockForestFilestem << ".bfs";
      const std::string setupBlockForestFilepath = oss.str();

      WALBERLA_LOG_INFO_ON_ROOT("Creating structured block forest...")
      auto bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()),
                                                 setupBlockForestFilepath.c_str(), false);
      auto blocks =
         std::make_shared< StructuredBlockForest >(bfs, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
      blocks->createCellBoundingBoxes();

      WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << blocks->getNumberOfBlocks())
      for (uint_t level = 0; level <= refinementDepth; level++)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << blocks->getNumberOfBlocks(level))
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                    DISTRIBUTED DATA STRUCTURES                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto ldc                     = std::make_shared< LDC >(refinementDepth);

      // Reading parameters
      auto parameters                = config->getOneBlock("Parameters");
      const real_t omega             = parameters.getParameter< real_t >("omega", real_c(1.4));
      const uint_t timesteps         = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const bool benchmarkKernelOnly = parameters.getParameter< bool >("benchmarkKernelOnly", false);

      // Creating fields
      const StorageSpecification_T StorageSpec = StorageSpecification_T();
      const BlockDataID pdfFieldID =
         lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, uint_c(2), field::fzyx);
      const BlockDataID velFieldID =
         field::addToStorage< VelocityField_T >(blocks, "vel", real_c(0.0), field::fzyx, uint_c(2));
      const BlockDataID densityFieldID =
         field::addToStorage< ScalarField_T >(blocks, "density", real_c(1.0), field::fzyx, uint_c(2));
      const BlockDataID flagFieldID =
         field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field", uint_c(3));

      const Cell innerOuterSplit =
         Cell(parameters.getParameter< Vector3< cell_idx_t > >("innerOuterSplit", Vector3< cell_idx_t >(1, 1, 1)));
      SweepCollection_T sweepCollection(blocks, pdfFieldID, densityFieldID, velFieldID, omega, innerOuterSplit);
      for (auto& block : *blocks)
      {
         sweepCollection.initialise(&block, 2);
      }
      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                      LB SWEEPS AND BOUNDARY HANDLING                                       ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      const FlagUID fluidFlagUID("Fluid");
      ldc->setupBoundaryFlagField(*blocks, flagFieldID);
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID, 2);
      BoundaryCollection_T boundaryCollection(blocks, flagFieldID, pdfFieldID, fluidFlagUID);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                           COMMUNICATION SCHEME                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      WALBERLA_LOG_INFO_ON_ROOT("Setting up communication...")
      auto communication = std::make_shared< NonUniformBufferedScheme< CommunicationStencil_T > >(blocks);
      auto packInfo      = lbm_generated::setupNonuniformPdfCommunication< PdfField_T >(blocks, pdfFieldID);
      communication->addPackInfo(packInfo);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                          TIME STEP DEFINITIONS                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      lbm_generated::BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T > LBMMeshRefinement(
         blocks, pdfFieldID, sweepCollection, boundaryCollection, communication, packInfo);

      SweepTimeloop timeLoop(blocks->getBlockStorage(), timesteps);

      if (benchmarkKernelOnly)
      {
         timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");
      }
      else { LBMMeshRefinement.addRefinementToTimeLoop(timeLoop); }

      // VTK
      const uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
         auto velWriter = make_shared< field::VTKWriter< VelocityField_T > >(velFieldID, "vel");
         vtkOutput->addCellDataWriter(velWriter);

         vtkOutput->addBeforeFunction([&]() {
            for (auto& block : *blocks)
               sweepCollection.calculateMacroscopicParameters(&block);
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

      lbm_generated::PerformanceEvaluation< FlagField_T > const performance(blocks, flagFieldID, fluidFlagUID);
      field::CellCounter< FlagField_T > fluidCells(blocks, flagFieldID, fluidFlagUID);
      fluidCells();

      WALBERLA_LOG_INFO_ON_ROOT("Non uniform Grid benchmark with " << fluidCells.numberOfCells()
                                                                   << " fluid cells (in total on all levels)")

      WcTimingPool timeloopTiming;
      WcTimer simTimer;

      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Starting benchmark with " << timesteps << " time steps")
      WALBERLA_MPI_BARRIER()

      simTimer.start();
      timeLoop.run(timeloopTiming);
      WALBERLA_MPI_BARRIER()
      simTimer.end();

      WALBERLA_LOG_INFO_ON_ROOT("Benchmark finished")
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
               pythonCallbackResults.data().exposeValue("numProcesses", performance.processes());
               pythonCallbackResults.data().exposeValue("numThreads", performance.threads());
               pythonCallbackResults.data().exposeValue("numCores", performance.cores());
               pythonCallbackResults.data().exposeValue("mlups", performance.mlups(timesteps, time));
               pythonCallbackResults.data().exposeValue("mlupsPerCore", performance.mlupsPerCore(timesteps, time));
               pythonCallbackResults.data().exposeValue("mlupsPerProcess",
                                                        performance.mlupsPerProcess(timesteps, time));
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