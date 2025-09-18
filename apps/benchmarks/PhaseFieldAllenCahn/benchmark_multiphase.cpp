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
//! \file benchmark_multiphase.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/SweepTimeloop.h"

#include "InitializerFunctions.h"

//////////////////////////////
// INCLUDE GENERATED FILES //
////////////////////////////

#if defined(WALBERLA_BUILD_WITH_CUDA)
#   include "gpu/AddGPUFieldToStorage.h"
#   include "gpu/DeviceSelectMPI.h"
#   include "gpu/ParallelStreams.h"
#   include "gpu/communication/MemcpyPackInfo.h"
#   include "gpu/communication/UniformGPUScheme.h"
#else
#   include <blockforest/communication/UniformBufferedScheme.h>
#endif

#include "GenDefines.h"

using namespace walberla;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

#if defined(WALBERLA_BUILD_WITH_CUDA)
typedef gpu::GPUField< real_t > GPUField;
#endif

int main(int argc, char** argv)
{
   const mpi::Environment env(argc, argv);
#if defined(WALBERLA_BUILD_WITH_CUDA)
   gpu::selectDeviceBasedOnMpiRank();
#endif

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      auto config = *cfg;
      logging::configureLogging(config);
      shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGridFromConfig(config);

      // Reading parameters
      auto parameters                    = config->getOneBlock("Parameters");
      const std::string timeStepStrategy = parameters.getParameter< std::string >("timeStepStrategy", "normal");
      const uint_t timesteps             = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const uint_t scenario = parameters.getParameter< uint_t >("scenario", uint_c(1));
      const uint_t warmupSteps  = parameters.getParameter< uint_t >("warmupSteps", uint_t(2));

#if defined(WALBERLA_BUILD_WITH_CUDA)
      // CPU fields
      const BlockDataID vel_field   = field::addToStorage< VelocityField_T >(blocks, "vel", real_c(0.0), field::fzyx);
      BlockDataID phase_field = field::addToStorage< PhaseField_T >(blocks, "phase", real_c(0.0), field::fzyx);
      // GPU fields
      const BlockDataID lb_phase_field_gpu = gpu::addGPUFieldToStorage< gpu::GPUField< real_t > >(
         blocks, "lb phase field on GPU", Stencil_phase_T::Size, field::fzyx, 1);
      const BlockDataID lb_velocity_field_gpu = gpu::addGPUFieldToStorage< gpu::GPUField< real_t > >(
         blocks, "lb velocity field on GPU", Stencil_hydro_T::Size, field::fzyx, 1);
      const BlockDataID vel_field_gpu =
         gpu::addGPUFieldToStorage< VelocityField_T >(blocks, vel_field, "velocity field on GPU", true);
      BlockDataID phase_field_gpu =
         gpu::addGPUFieldToStorage< PhaseField_T >(blocks, phase_field, "phase field on GPU", true);
      BlockDataID phase_field_tmp = gpu::addGPUFieldToStorage< PhaseField_T >(blocks, phase_field, "temporary phasefield", true);
#else
      BlockDataID lb_phase_field =
         field::addToStorage< PdfField_phase_T >(blocks, "lb phase field", real_c(0.0), field::fzyx);
      BlockDataID lb_velocity_field =
         field::addToStorage< PdfField_hydro_T >(blocks, "lb velocity field", real_c(0.0), field::fzyx);
      BlockDataID vel_field   = field::addToStorage< VelocityField_T >(blocks, "vel", real_c(0.0), field::fzyx);
      BlockDataID phase_field = field::addToStorage< PhaseField_T >(blocks, "phase", real_c(0.0), field::fzyx);
      BlockDataID phase_field_tmp = field::addToStorage< PhaseField_T >(blocks, "phase tmp", real_c(0.0), field::fzyx);
#endif

      if (timeStepStrategy != "phase_only" && timeStepStrategy != "hydro_only" && timeStepStrategy != "kernel_only")
      {
         WALBERLA_LOG_INFO_ON_ROOT("initialization of the phase field")
         if (scenario == 1)
         {
            auto bubbleParameters = config->getOneBlock("Bubble");
            const Vector3< real_t > bubbleMidPoint =
               bubbleParameters.getParameter< Vector3< real_t > >("bubbleMidPoint");
            const real_t bubbleRadius = bubbleParameters.getParameter< real_t >("bubbleRadius", real_c(20.0));
            initPhaseField_bubble(blocks, phase_field, bubbleRadius, bubbleMidPoint);
         }
         else if (scenario == 2)
         {
            initPhaseField_RTI(blocks, phase_field);
         }
#if defined(WALBERLA_BUILD_WITH_CUDA)
         gpu::fieldCpy< GPUField, PhaseField_T >(blocks, phase_field_gpu, phase_field);
#endif
         WALBERLA_LOG_INFO_ON_ROOT("initialization of the phase field done")
      }

#if defined(WALBERLA_BUILD_WITH_CUDA)
      Vector3< int32_t > gpuBlockSize =
         parameters.getParameter< Vector3< int32_t > >("gpuBlockSize", Vector3< int32_t >(256, 1, 1));
      pystencils::initialize_phase_field_distributions init_h(lb_phase_field_gpu, phase_field_gpu, vel_field_gpu);
      pystencils::initialize_velocity_based_distributions init_g(lb_velocity_field_gpu, vel_field_gpu);

      pystencils::phase_field_LB_step phase_field_LB_step(
         lb_phase_field_gpu, phase_field_gpu, phase_field_tmp, vel_field_gpu, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
      pystencils::hydro_LB_step hydro_LB_step(lb_velocity_field_gpu, phase_field_gpu, vel_field_gpu, gpuBlockSize[0],
                                              gpuBlockSize[1], gpuBlockSize[2]);
#else
      pystencils::initialize_phase_field_distributions init_h(lb_phase_field, phase_field, vel_field);
      pystencils::initialize_velocity_based_distributions init_g(lb_velocity_field, vel_field);
      pystencils::phase_field_LB_step phase_field_LB_step(lb_phase_field, phase_field, phase_field_tmp, vel_field);
      pystencils::hydro_LB_step hydro_LB_step(lb_velocity_field, phase_field, vel_field);
#endif

// add communication
#if defined(WALBERLA_BUILD_WITH_CUDA)
      const bool gpuEnabledMpi = parameters.getParameter< bool >("cudaEnabledMpi", false);
      const int streamLowPriority  = 0;
      const int streamHighPriority = 0;
      auto defaultStream     = gpu::StreamRAII::newPriorityStream(streamLowPriority);
      auto innerOuterStreams = gpu::ParallelStreams(streamHighPriority);

      auto generatedPackInfo_phase_field_distributions = make_shared< lbm::PackInfo_phase_field_distributions>(lb_phase_field_gpu);
      auto generatedPackInfo_velocity_based_distributions = make_shared< lbm::PackInfo_velocity_based_distributions >(lb_velocity_field_gpu);
      auto generatedPackInfo_phase_field = make_shared< pystencils::PackInfo_phase_field >(phase_field_gpu);

      auto UniformGPUSchemeVelocityBasedDistributions = make_shared< gpu::communication::UniformGPUScheme< Stencil_hydro_T > >(blocks, gpuEnabledMpi, false);
      auto UniformGPUSchemePhaseFieldDistributions = make_shared< gpu::communication::UniformGPUScheme< Full_Stencil_T > >(blocks, gpuEnabledMpi, false);
      auto UniformGPUSchemePhaseField = make_shared< gpu::communication::UniformGPUScheme< Stencil_hydro_T > >(blocks, gpuEnabledMpi, false, 65432);

      UniformGPUSchemeVelocityBasedDistributions->addPackInfo(generatedPackInfo_velocity_based_distributions);
      UniformGPUSchemePhaseFieldDistributions->addPackInfo(generatedPackInfo_phase_field_distributions);
      UniformGPUSchemePhaseField->addPackInfo(generatedPackInfo_phase_field);

      auto Comm_velocity_based_distributions_start = std::function< void() >([&]() { UniformGPUSchemeVelocityBasedDistributions->startCommunication(); });
      auto Comm_velocity_based_distributions_wait = std::function< void() >([&]() { UniformGPUSchemeVelocityBasedDistributions->wait(); });

      auto Comm_phase_field_distributions_start = std::function< void() >([&]() { UniformGPUSchemePhaseFieldDistributions->startCommunication(); });
      auto Comm_phase_field_distributions_wait = std::function< void() >([&]() { UniformGPUSchemePhaseFieldDistributions->wait(); });

      auto Comm_phase_field = std::function< void() >([&]() { UniformGPUSchemePhaseField->communicate(); });

      auto swapPhaseField = std::function< void(IBlock *) >([&](IBlock * b)
        {
           auto phaseField    = b->getData< gpu::GPUField<real_t> >(phase_field_gpu);
           auto phaseFieldTMP = b->getData< gpu::GPUField<real_t> >(phase_field_tmp);
           phaseField->swapDataPointers(phaseFieldTMP);
        });

#else

      auto generatedPackInfo_phase_field_distributions = make_shared< lbm::PackInfo_phase_field_distributions>(lb_phase_field);
      auto generatedPackInfo_velocity_based_distributions = make_shared< lbm::PackInfo_velocity_based_distributions >(lb_velocity_field);
      auto generatedPackInfo_phase_field = make_shared< pystencils::PackInfo_phase_field >(phase_field);

      auto UniformGPUSchemeVelocityBasedDistributions = make_shared< blockforest::communication::UniformBufferedScheme< Full_Stencil_T > >(blocks);
      auto UniformGPUSchemePhaseFieldDistributions = make_shared< blockforest::communication::UniformBufferedScheme< Full_Stencil_T > >(blocks);
      auto UniformGPUSchemePhaseField = make_shared< blockforest::communication::UniformBufferedScheme< Full_Stencil_T > >(blocks, 65432);

      UniformGPUSchemeVelocityBasedDistributions->addPackInfo(generatedPackInfo_velocity_based_distributions);
      UniformGPUSchemePhaseFieldDistributions->addPackInfo(generatedPackInfo_phase_field_distributions);
      UniformGPUSchemePhaseField->addPackInfo(generatedPackInfo_phase_field);

      auto Comm_velocity_based_distributions_start = std::function< void() >([&]() { UniformGPUSchemeVelocityBasedDistributions->startCommunication(); });
      auto Comm_velocity_based_distributions_wait = std::function< void() >([&]() { UniformGPUSchemeVelocityBasedDistributions->wait(); });

      auto Comm_phase_field_distributions = std::function< void() >([&]() { UniformGPUSchemePhaseFieldDistributions->communicate(); });
      auto Comm_phase_field_distributions_start = std::function< void() >([&]() { UniformGPUSchemePhaseFieldDistributions->startCommunication(); });
      auto Comm_phase_field_distributions_wait = std::function< void() >([&]() { UniformGPUSchemePhaseFieldDistributions->wait(); });

      auto Comm_phase_field = std::function< void() >([&]() { UniformGPUSchemePhaseField->communicate(); });

      auto swapPhaseField = std::function< void(IBlock *) >([&](IBlock * b)
        {
           auto phaseField    = b->getData< PhaseField_T >(phase_field);
           auto phaseFieldTMP = b->getData< PhaseField_T >(phase_field_tmp);
           phaseField->swapDataPointers(phaseFieldTMP);
        });
#endif

      BlockDataID const flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");
      // Boundaries
      const FlagUID fluidFlagUID("Fluid");
      auto boundariesConfig = config->getBlock("Boundaries_GPU");
      if (boundariesConfig)
      {
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldID, boundariesConfig);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID);
      }

      // initialize the two lattice Boltzmann fields
      if (timeStepStrategy != "phase_only" && timeStepStrategy != "hydro_only" && timeStepStrategy != "kernel_only")
      {
         WALBERLA_LOG_INFO_ON_ROOT("initialization of the distributions")
         for (auto& block : *blocks)
         {
            init_h(&block);
            init_g(&block);
         }
#if defined(WALBERLA_BUILD_WITH_CUDA)
         WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
         WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#endif
         WALBERLA_MPI_BARRIER()
         WALBERLA_LOG_INFO_ON_ROOT("initialization of the distributions done")
      }

      SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);
#if defined(WALBERLA_BUILD_WITH_CUDA)
      timeloop.add() << BeforeFunction(Comm_velocity_based_distributions_start, "Start Hydro PDFs Communication")
                     << Sweep(phase_field_LB_step.getSweep(defaultStream), "Phase LB Step")
                     << AfterFunction(Comm_velocity_based_distributions_wait, "Wait Hydro PDFs Communication");

      timeloop.add() << BeforeFunction(Comm_phase_field_distributions_start, "Start Phase PDFs Communication")
                     << Sweep(hydro_LB_step.getSweep(defaultStream), "Hydro LB Step");
      timeloop.add() << Sweep(swapPhaseField, "Swap PhaseField")
                     << AfterFunction(Comm_phase_field_distributions_wait, "Wait Phase PDFs Communication");

      timeloop.addFuncAfterTimeStep(Comm_phase_field, "Communication Phase field");

#else
      timeloop.add() << BeforeFunction(Comm_velocity_based_distributions_start, "Start Hydro PDFs Communication")
                     << Sweep(phase_field_LB_step.getSweep(), "Phase LB Step")
                     << AfterFunction(Comm_velocity_based_distributions_wait, "Wait Hydro PDFs Communication");

      timeloop.add() << BeforeFunction(Comm_phase_field_distributions_start, "Start Phase PDFs Communication")
                     << Sweep(hydro_LB_step.getSweep(), "Hydro LB Step");
      timeloop.add() << Sweep(swapPhaseField, "Swap PhaseField")
                     << AfterFunction(Comm_phase_field_distributions_wait, "Wait Phase PDFs Communication");

      timeloop.addFuncAfterTimeStep(Comm_phase_field, "Communication Phase field");
#endif

      uint_t const vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 1)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
#if defined(WALBERLA_BUILD_WITH_CUDA)
         vtkOutput->addBeforeFunction(
            [&]() { gpu::fieldCpy< PhaseField_T, GPUField >(blocks, phase_field, phase_field_gpu); });
#endif
         auto phaseWriter = make_shared< field::VTKWriter< PhaseField_T > >(phase_field, "phase");
         vtkOutput->addCellDataWriter(phaseWriter);

         timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      lbm_generated::PerformanceEvaluation< FlagField_T > const performance(blocks, flagFieldID, fluidFlagUID);
      field::CellCounter< FlagField_T > fluidCells(blocks, flagFieldID, fluidFlagUID);
      fluidCells();

      WALBERLA_LOG_INFO_ON_ROOT("Multiphase benchmark with " << fluidCells.numberOfCells() << " fluid cells")
      WALBERLA_LOG_INFO_ON_ROOT("Running " << warmupSteps << " timesteps to warm up the system")

      for (uint_t i = 0; i < warmupSteps; ++i)
         timeloop.singleStep();

#if defined(WALBERLA_BUILD_WITH_CUDA)
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#endif
      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Warmup timesteps done")

      timeloop.setCurrentTimeStepToZero();
      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
      WcTimingPool timeloopTiming;
      WcTimer simTimer;
#if defined(WALBERLA_BUILD_WITH_CUDA)
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
#endif
      simTimer.start();
      timeloop.run(timeloopTiming);
#if defined(WALBERLA_BUILD_WITH_CUDA)
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#endif
      WALBERLA_MPI_BARRIER()
      simTimer.end();
      WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
      double time = simTimer.max();
      WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
      performance.logResultOnRoot(timesteps, time);

      const auto reducedTimeloopTiming = timeloopTiming.getReduced();
      WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)

      WALBERLA_LOG_RESULT_ON_ROOT("MLUPS per process: " << performance.mlupsPerProcess(timesteps, time))
      WALBERLA_LOG_RESULT_ON_ROOT("Time per time step: " << time / real_c(timesteps) << " s")
      WALBERLA_ROOT_SECTION()
      {
         python_coupling::PythonCallback pythonCallbackResults("results_callback");
         if (pythonCallbackResults.isCallable())
         {
            pythonCallbackResults.data().exposeValue("mlupsPerProcess", performance.mlupsPerProcess(timesteps, time));
            pythonCallbackResults.data().exposeValue("stencil_phase", StencilNamePhase);
            pythonCallbackResults.data().exposeValue("stencil_hydro", StencilNameHydro);
            #if defined(WALBERLA_BUILD_WITH_CUDA)
               pythonCallbackResults.data().exposeValue("cuda_enabled_mpi", gpuEnabledMpi);
            #endif
            // Call Python function to report results
            pythonCallbackResults();
         }
      }
   }

   return EXIT_SUCCESS;
}
