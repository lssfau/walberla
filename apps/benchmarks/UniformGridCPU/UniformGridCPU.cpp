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
//! \file UniformGridCPU.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/OpenMP.h"
#include "core/logging/Initialization.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/communication/UniformGeneratedPdfPackInfo.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/all.h"

#include <iomanip>

#include "InitShearVelocity.h"
#include "ManualKernels.h"
#include "UniformGridCPU_InfoHeader.h"

using namespace walberla;

using StorageSpecification_T = lbm::UniformGridCPUStorageSpecification;
using Stencil_T = lbm::UniformGridCPUStorageSpecification::Stencil;

using PdfField_T = lbm_generated::PdfField< StorageSpecification_T >;
using FlagField_T = FlagField< uint8_t >;
using BoundaryCollection_T = lbm::UniformGridCPUBoundaryCollection< FlagField_T >;

using SweepCollection_T = lbm::UniformGridCPUSweepCollection;

using blockforest::communication::UniformBufferedScheme;

int main(int argc, char** argv)
{
   const mpi::Environment env(argc, argv);

   const std::string input_filename(argv[1]);
   const bool inputIsPython = string_ends_with(input_filename, ".py");

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      auto config = *cfg;
      logging::configureLogging(config);
      auto blocks = blockforest::createUniformBlockGridFromConfig(config);

      // Reading parameters
      auto parameters          = config->getOneBlock("Parameters");
      const real_t omega       = parameters.getParameter< real_t >("omega", real_c(1.4));
      const uint_t timesteps   = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const bool initShearFlow = parameters.getParameter< bool >("initShearFlow", true);

      // Creating fields
      const StorageSpecification_T StorageSpec = StorageSpecification_T();
      auto fieldAllocator = make_shared< field::AllocateAligned< real_t, 64 > >();
      const BlockDataID pdfFieldId  = lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, field::fzyx, fieldAllocator);
      const BlockDataID velFieldId = field::addToStorage< VelocityField_T >(blocks, "vel", real_c(0.0), field::fzyx);
      const BlockDataID densityFieldId = field::addToStorage< ScalarField_T >(blocks, "density", real_c(1.0), field::fzyx);
      const BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field");

      // Initialize velocity on cpu
      if (initShearFlow)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Initializing shear flow")
         initShearVelocity(blocks, velFieldId);
      }

      const Cell innerOuterSplit = Cell(parameters.getParameter< Vector3<cell_idx_t> >("innerOuterSplit", Vector3<cell_idx_t>(1, 1, 1)));
      SweepCollection_T sweepCollection(blocks, pdfFieldId, densityFieldId, velFieldId, omega, innerOuterSplit);

      for (auto& block : *blocks)
      {
         sweepCollection.initialise(&block);
      }

      const pystencils::UniformGridCPU_StreamOnlyKernel StreamOnlyKernel(pdfFieldId);

      // Boundaries
      const FlagUID fluidFlagUID("Fluid");
      auto boundariesConfig   = config->getBlock("Boundaries");
      if (boundariesConfig)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Setting boundary conditions")
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldID, boundariesConfig);
      }
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID);
      BoundaryCollection_T boundaryCollection(blocks, flagFieldID, pdfFieldId, fluidFlagUID);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                           COMMUNICATION SCHEME                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto packInfo = std::make_shared<lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T >>(pdfFieldId);
      UniformBufferedScheme< Stencil_T > communication(blocks);
      communication.addPackInfo(packInfo);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                          TIME STEP DEFINITIONS                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      SweepTimeloop timeLoop(blocks->getBlockStorage(), timesteps);
      const std::string timeStepStrategy = parameters.getParameter< std::string >("timeStepStrategy", "normal");

      if (timeStepStrategy == "noOverlap") {
         if (boundariesConfig){
            timeLoop.add() << BeforeFunction(communication, "communication")
                           << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL), "Boundary Conditions");
            timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");
         }else {
            timeLoop.add() << BeforeFunction(communication, "communication")
                           << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");}

      } else if (timeStepStrategy == "simpleOverlap") {
         if (boundariesConfig){
            timeLoop.add() << BeforeFunction(communication.getStartCommunicateFunctor(), "Start Communication")
                           << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL), "Boundary Conditions");
            timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::INNER), "LBM StreamCollide Inner Frame");
            timeLoop.add() << BeforeFunction(communication.getWaitFunctor(), "Wait for Communication")
                           << Sweep(sweepCollection.streamCollide(SweepCollection_T::OUTER), "LBM StreamCollide Outer Frame");
         }else{
            timeLoop.add() << BeforeFunction(communication.getStartCommunicateFunctor(), "Start Communication")
                           << Sweep(sweepCollection.streamCollide(SweepCollection_T::INNER), "LBM StreamCollide Inner Frame");
            timeLoop.add() << BeforeFunction(communication.getWaitFunctor(), "Wait for Communication")
                           << Sweep(sweepCollection.streamCollide(SweepCollection_T::OUTER), "LBM StreamCollide Outer Frame");}

      } else if (timeStepStrategy == "kernelOnly") {
         timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");
      } else if (timeStepStrategy == "StreamOnly") {
         timeLoop.add() << Sweep(sweepCollection.stream(SweepCollection_T::ALL), "LBM Stream-Only");
      } else {
         WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'")
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                             TIME LOOP SETUP                                                ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      const uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
         auto velWriter = make_shared< field::VTKWriter< VelocityField_T, float32 > >(velFieldId, "vel");
         vtkOutput->addCellDataWriter(velWriter);

         vtkOutput->addBeforeFunction([&]() {
           for (auto& block : *blocks){
              sweepCollection.calculateMacroscopicParameters(&block);}
         });

         timeLoop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
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
         auto logger = timing::RemainingTimeLogger(timeLoop.getNrOfTimeSteps() * outerIterations,
                                                   remainingTimeLoggerFrequency);
         timeLoop.addFuncAfterTimeStep(logger, "remaining time logger");
      }

      for (uint_t outerIteration = 0; outerIteration < outerIterations; ++outerIteration)
      {
         timeLoop.setCurrentTimeStepToZero();

         WcTimingPool timeloopTiming;
         WcTimer simTimer;

         WALBERLA_MPI_WORLD_BARRIER()
         WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")

         simTimer.start();
         timeLoop.run(timeloopTiming);
         simTimer.end();

         WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
         double time = simTimer.max();
         WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
         performance.logResultOnRoot(timesteps, time);

         const auto reducedTimeloopTiming = timeloopTiming.getReduced();
         WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)

         WALBERLA_ROOT_SECTION()
         {
            if(inputIsPython){
               python_coupling::PythonCallback pythonCallbackResults("results_callback");
               if (pythonCallbackResults.isCallable())
               {
                  pythonCallbackResults.data().exposeValue("numProcesses", performance.processes());
                  pythonCallbackResults.data().exposeValue("numThreads", performance.threads());
                  pythonCallbackResults.data().exposeValue("numCores", performance.cores());
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
   }
   return EXIT_SUCCESS;
}
