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

#include "lbm/communication/CombinedInPlaceCpuPackInfo.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/all.h"

#include <iomanip>

#include "InitShearVelocity.h"
#include "ManualKernels.h"
#include "UniformGridCPU_InfoHeader.h"

using namespace walberla;

using PackInfoEven_T = lbm::UniformGridCPU_PackInfoEven;
using PackInfoOdd_T = lbm::UniformGridCPU_PackInfoOdd;
using LbSweep = lbm::UniformGridCPU_LbKernel;

using FlagField_T = FlagField< uint8_t >;

auto pdfFieldAdder = [](IBlock* const block, StructuredBlockStorage* const storage) {
   return new PdfField_T(storage->getNumberOfXCells(*block), storage->getNumberOfYCells(*block),
                         storage->getNumberOfZCells(*block), uint_t(1), field::fzyx,
                         make_shared< field::AllocateAligned< real_t, 64 > >());
};

int main(int argc, char** argv)
{
   mpi::Environment env(argc, argv);

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      auto config = *cfg;
      logging::configureLogging(config);
      auto blocks = blockforest::createUniformBlockGridFromConfig(config);

      Vector3< uint_t > cellsPerBlock =
         config->getBlock("DomainSetup").getParameter< Vector3< uint_t > >("cellsPerBlock");
      // Reading parameters
      auto parameters          = config->getOneBlock("Parameters");
      const real_t omega       = parameters.getParameter< real_t >("omega", real_c(1.4));
      const uint_t timesteps   = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const bool initShearFlow = parameters.getParameter< bool >("initShearFlow", true);

      // Creating fields
      BlockDataID pdfFieldId = blocks->addStructuredBlockData< PdfField_T >(pdfFieldAdder, "pdfs");
      BlockDataID velFieldId = field::addToStorage< VelocityField_T >(blocks, "vel", real_t(0), field::fzyx);
      BlockDataID densityFieldId = field::addToStorage< ScalarField_T >(blocks, "density", real_t(1.0), field::fzyx);

      // Initialize velocity on cpu
      if (initShearFlow)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Initializing shear flow")
         initShearVelocity(blocks, velFieldId);
      }

      pystencils::UniformGridCPU_MacroSetter setterSweep(densityFieldId, pdfFieldId, velFieldId);
      pystencils::UniformGridCPU_MacroGetter getterSweep(densityFieldId, pdfFieldId, velFieldId);

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

      LbSweep lbSweep(pdfFieldId, omega, innerOuterSplitCell);
      pystencils::UniformGridCPU_StreamOnlyKernel StreamOnlyKernel(pdfFieldId);

      // Boundaries
      const FlagUID fluidFlagUID("Fluid");
      BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field");
      auto boundariesConfig   = config->getBlock("Boundaries");
      bool boundaries         = false;
      if (boundariesConfig)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Setting boundary conditions")
         boundaries = true;
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldID, boundariesConfig);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID);
      }

      lbm::UniformGridCPU_NoSlip noSlip(blocks, pdfFieldId);
      noSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldID, FlagUID("NoSlip"), fluidFlagUID);

      lbm::UniformGridCPU_UBB ubb(blocks, pdfFieldId);
      ubb.fillFromFlagField< FlagField_T >(blocks, flagFieldID, FlagUID("UBB"), fluidFlagUID);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                           COMMUNICATION SCHEME                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Initial setup is the post-collision state of an even time step
      auto tracker = make_shared< lbm::TimestepTracker >(0);
      auto packInfo =
         make_shared< lbm::CombinedInPlaceCpuPackInfo< PackInfoEven_T , PackInfoOdd_T > >(tracker, pdfFieldId);

      blockforest::communication::UniformBufferedScheme< Stencil_T > communication(blocks);
      communication.addPackInfo(packInfo);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                          TIME STEP DEFINITIONS                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto boundarySweep = [&](IBlock* block, uint8_t t) {
         noSlip.run(block, t);
         ubb.run(block, t);
      };

      auto boundaryInner = [&](IBlock* block, uint8_t t) {
         noSlip.inner(block, t);
         ubb.inner(block, t);
      };

      auto boundaryOuter = [&](IBlock* block, uint8_t t) {
         noSlip.outer(block, t);
         ubb.outer(block, t);
      };

      auto simpleOverlapTimeStep = [&]() {
         // Communicate post-collision values of previous timestep...
         communication.startCommunication();
         for (auto& block : *blocks)
         {
            if (boundaries) boundaryInner(&block, tracker->getCounter());
            lbSweep.inner(&block, tracker->getCounterPlusOne());
         }
         communication.wait();
         for (auto& block : *blocks)
         {
            if (boundaries) boundaryOuter(&block, tracker->getCounter());
            lbSweep.outer(&block, tracker->getCounterPlusOne());
         }

         tracker->advance();
      };

      auto normalTimeStep = [&]() {
         communication.communicate();
         for (auto& block : *blocks)
         {
            if (boundaries) boundarySweep(&block, tracker->getCounter());
            lbSweep(&block, tracker->getCounterPlusOne());
         }

         tracker->advance();
      };

      // With two-fields patterns, ghost layer cells act as constant stream-in boundaries;
      // with in-place patterns, ghost layer cells act as wet-node no-slip boundaries.
      auto kernelOnlyFunc = [&]() {
         tracker->advance();
         for (auto& block : *blocks)
            lbSweep(&block, tracker->getCounter());
      };

      // Stream only function to test a streaming pattern without executing lbm operations inside
      auto StreamOnlyFunc = [&]() {
         for (auto& block : *blocks)
            StreamOnlyKernel(&block);
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
         communication.communicate();
         timeStep = kernelOnlyFunc;
      }
      else if (timeStepStrategy == "StreamOnly")
      {
         WALBERLA_LOG_INFO_ON_ROOT(
            "Running only streaming kernel without LBM - this makes only sense for benchmarking!")
         // Run initial communication once to provide any missing stream-in populations
         timeStep = StreamOnlyFunc;
      }
      else
      {
         WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'. Allowed values are 'noOverlap', "
                                      "'simpleOverlap', 'kernelOnly'")
      }

      timeLoop.add() << BeforeFunction(timeStep) << Sweep([](IBlock*) {}, "time step");

      uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
         auto velWriter = make_shared< field::VTKWriter< VelocityField_T > >(velFieldId, "vel");
         vtkOutput->addCellDataWriter(velWriter);

         vtkOutput->addBeforeFunction([&]() {
           for (auto& block : *blocks){
              getterSweep(&block);}
         });

         timeLoop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                               BENCHMARK                                                    ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      int warmupSteps     = parameters.getParameter< int >("warmupSteps", 2);
      int outerIterations = parameters.getParameter< int >("outerIterations", 1);
      for (int i = 0; i < warmupSteps; ++i)
         timeLoop.singleStep();

      real_t remainingTimeLoggerFrequency =
         parameters.getParameter< real_t >("remainingTimeLoggerFrequency", -1.0); // in seconds
      if (remainingTimeLoggerFrequency > 0)
      {
         auto logger = timing::RemainingTimeLogger(timeLoop.getNrOfTimeSteps() * uint_c(outerIterations),
                                                   remainingTimeLoggerFrequency);
         timeLoop.addFuncAfterTimeStep(logger, "remaining time logger");
      }

      for (int outerIteration = 0; outerIteration < outerIterations; ++outerIteration)
      {
         timeLoop.setCurrentTimeStepToZero();
         WcTimer simTimer;
         WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
         simTimer.start();
         timeLoop.run();
         simTimer.end();
         WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
         real_t time      = simTimer.last();
         WALBERLA_MPI_SECTION()
         {
            walberla::mpi::reduceInplace(time, walberla::mpi::MAX);
         }
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
   return EXIT_SUCCESS;
}
