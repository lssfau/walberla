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
#include "blockforest/communication/UniformDirectScheme.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/SweepTimeloop.h"

#include <blockforest/communication/UniformBufferedScheme.h>

#include "InitializerFunctions.h"

//////////////////////////////
// INCLUDE GENERATED FILES //
////////////////////////////

#if defined(WALBERLA_BUILD_WITH_CUDA)
#   include "cuda/AddGPUFieldToStorage.h"
#   include "cuda/DeviceSelectMPI.h"
#   include "cuda/HostFieldAllocator.h"
#   include "cuda/ParallelStreams.h"
#   include "cuda/communication/GPUPackInfo.h"
#   include "cuda/communication/MemcpyPackInfo.h"
#   include "cuda/communication/UniformGPUScheme.h"

#   include "GenDefines.h"
#   include "PackInfo_phase_field.h"
#   include "PackInfo_phase_field_distributions.h"
#   include "PackInfo_velocity_based_distributions.h"
#   include "hydro_LB_step.h"
#   include "initialize_phase_field_distributions.h"
#   include "initialize_velocity_based_distributions.h"
#   include "phase_field_LB_step.h"

#else
#   include "GenDefines.h"
#   include "PackInfo_phase_field.h"
#   include "PackInfo_phase_field_distributions.h"
#   include "PackInfo_velocity_based_distributions.h"
#   include "hydro_LB_step.h"
#   include "initialize_phase_field_distributions.h"
#   include "initialize_velocity_based_distributions.h"
#   include "phase_field_LB_step.h"
#endif

using namespace walberla;

using PdfField_phase_T = GhostLayerField< real_t, Stencil_phase_T::Size >;
using PdfField_hydro_T = GhostLayerField< real_t, Stencil_hydro_T::Size >;
using VelocityField_T  = GhostLayerField< real_t, Stencil_hydro_T::Dimension >;
using PhaseField_T     = GhostLayerField< real_t, 1 >;
using flag_t           = walberla::uint8_t;
using FlagField_T      = FlagField< flag_t >;

#if defined(WALBERLA_BUILD_WITH_CUDA)
typedef cuda::GPUField< real_t > GPUField;
#endif

int main(int argc, char** argv)
{
   mpi::Environment env(argc, argv);
#if defined(WALBERLA_BUILD_WITH_CUDA)
   cuda::selectDeviceBasedOnMpiRank();
#endif

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      auto config = *cfg;
      logging::configureLogging(config);
      shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGridFromConfig(config);

      Vector3< uint_t > cellsPerBlock =
         config->getBlock("DomainSetup").getParameter< Vector3< uint_t > >("cellsPerBlock");
      // Reading parameters
      auto parameters                    = config->getOneBlock("Parameters");
      const std::string timeStepStrategy = parameters.getParameter< std::string >("timeStepStrategy", "normal");
      const uint_t timesteps             = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const real_t remainingTimeLoggerFrequency =
         parameters.getParameter< real_t >("remainingTimeLoggerFrequency", 3.0);
      const uint_t scenario = parameters.getParameter< uint_t >("scenario", uint_c(1));
      Vector3< int > overlappingWidth =
         parameters.getParameter< Vector3< int > >("overlappingWidth", Vector3< int >(1, 1, 1));
      const uint_t warmupSteps = parameters.getParameter< uint_t >("warmupSteps", uint_t(2));

#if defined(WALBERLA_BUILD_WITH_CUDA)
      // CPU fields
      BlockDataID vel_field   = field::addToStorage< VelocityField_T >(blocks, "vel", real_t(0), field::fzyx);
      BlockDataID phase_field = field::addToStorage< PhaseField_T >(blocks, "phase", real_t(0), field::fzyx);
      // GPU fields
      BlockDataID lb_phase_field_gpu = cuda::addGPUFieldToStorage< cuda::GPUField< real_t > >(
         blocks, "lb phase field on GPU", Stencil_phase_T::Size, field::fzyx, 1);
      BlockDataID lb_velocity_field_gpu = cuda::addGPUFieldToStorage< cuda::GPUField< real_t > >(
         blocks, "lb velocity field on GPU", Stencil_hydro_T::Size, field::fzyx, 1);
      BlockDataID vel_field_gpu =
         cuda::addGPUFieldToStorage< VelocityField_T >(blocks, vel_field, "velocity field on GPU", true);
      BlockDataID phase_field_gpu =
         cuda::addGPUFieldToStorage< PhaseField_T >(blocks, phase_field, "phase field on GPU", true);
#else
      BlockDataID lb_phase_field =
         field::addToStorage< PdfField_phase_T >(blocks, "lb phase field", real_t(0), field::fzyx);
      BlockDataID lb_velocity_field =
         field::addToStorage< PdfField_hydro_T >(blocks, "lb velocity field", real_t(0), field::fzyx);
      BlockDataID vel_field   = field::addToStorage< VelocityField_T >(blocks, "vel", real_t(0), field::fzyx);
      BlockDataID phase_field = field::addToStorage< PhaseField_T >(blocks, "phase", real_t(0), field::fzyx);
#endif

      if (timeStepStrategy != "phase_only" && timeStepStrategy != "hydro_only" && timeStepStrategy != "kernel_only")
      {
         WALBERLA_LOG_INFO_ON_ROOT("initialization of the phase field")
         if (scenario == 1)
         {
            auto bubbleParameters = config->getOneBlock("Bubble");
            const Vector3< real_t > bubbleMidPoint =
               bubbleParameters.getParameter< Vector3< real_t > >("bubbleMidPoint");
            const real_t bubbleRadius = bubbleParameters.getParameter< real_t >("bubbleRadius", 20.0);
            initPhaseField_bubble(blocks, phase_field, bubbleRadius, bubbleMidPoint);
         }
         else if (scenario == 2)
         {
            initPhaseField_RTI(blocks, phase_field);
         }
#if defined(WALBERLA_BUILD_WITH_CUDA)
         cuda::fieldCpy< GPUField, PhaseField_T >(blocks, phase_field_gpu, phase_field);
#endif
         WALBERLA_LOG_INFO_ON_ROOT("initialization of the phase field done")
      }

#if defined(WALBERLA_BUILD_WITH_CUDA)
      Vector3< int32_t > gpuBlockSize =
         parameters.getParameter< Vector3< int32_t > >("gpuBlockSize", Vector3< int32_t >(256, 1, 1));
      pystencils::initialize_phase_field_distributions init_h(lb_phase_field_gpu, phase_field_gpu, vel_field_gpu);
      pystencils::initialize_velocity_based_distributions init_g(lb_velocity_field_gpu, vel_field_gpu);

      pystencils::phase_field_LB_step phase_field_LB_step(
         lb_phase_field_gpu, phase_field_gpu, vel_field_gpu, gpuBlockSize[0], gpuBlockSize[1],
         Cell(overlappingWidth[0], overlappingWidth[1], overlappingWidth[2]));
      pystencils::hydro_LB_step hydro_LB_step(lb_velocity_field_gpu, phase_field_gpu, vel_field_gpu, gpuBlockSize[0],
                                              gpuBlockSize[1],
                                              Cell(overlappingWidth[0], overlappingWidth[1], overlappingWidth[2]));
#else
      pystencils::initialize_phase_field_distributions init_h(lb_phase_field, phase_field, vel_field);
      pystencils::initialize_velocity_based_distributions init_g(lb_velocity_field, vel_field);
      pystencils::phase_field_LB_step phase_field_LB_step(
         lb_phase_field, phase_field, vel_field, Cell(overlappingWidth[0], overlappingWidth[1], overlappingWidth[2]));
      pystencils::hydro_LB_step hydro_LB_step(lb_velocity_field, phase_field, vel_field,
                                              Cell(overlappingWidth[0], overlappingWidth[1], overlappingWidth[2]));
#endif

// add communication
#if defined(WALBERLA_BUILD_WITH_CUDA)
      auto Comm_velocity_based_distributions =
         make_shared< cuda::communication::UniformGPUScheme< Stencil_hydro_T > >(blocks, 0);
      auto generatedPackInfo_velocity_based_distributions =
         make_shared< lbm::PackInfo_velocity_based_distributions >(lb_velocity_field_gpu);
      Comm_velocity_based_distributions->addPackInfo(generatedPackInfo_velocity_based_distributions);
      auto generatedPackInfo_phase_field = make_shared< pystencils::PackInfo_phase_field >(phase_field_gpu);
      Comm_velocity_based_distributions->addPackInfo(generatedPackInfo_phase_field);

      auto Comm_phase_field_distributions =
         make_shared< cuda::communication::UniformGPUScheme< Stencil_hydro_T > >(blocks, 0);
      auto generatedPackInfo_phase_field_distributions =
         make_shared< lbm::PackInfo_phase_field_distributions >(lb_phase_field_gpu);
      Comm_phase_field_distributions->addPackInfo(generatedPackInfo_phase_field_distributions);
#else


      blockforest::communication::UniformBufferedScheme< Stencil_hydro_T > Comm_velocity_based_distributions(blocks);

      auto generatedPackInfo_phase_field = make_shared< pystencils::PackInfo_phase_field >(phase_field);
      auto generatedPackInfo_velocity_based_distributions =
         make_shared< lbm::PackInfo_velocity_based_distributions >(lb_velocity_field);

      Comm_velocity_based_distributions.addPackInfo(generatedPackInfo_phase_field);
      Comm_velocity_based_distributions.addPackInfo(generatedPackInfo_velocity_based_distributions);

      blockforest::communication::UniformBufferedScheme< Stencil_hydro_T > Comm_phase_field_distributions(blocks);
      auto generatedPackInfo_phase_field_distributions =
         make_shared< lbm::PackInfo_phase_field_distributions >(lb_phase_field);
      Comm_phase_field_distributions.addPackInfo(generatedPackInfo_phase_field_distributions);
#endif

      BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");
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
         WALBERLA_LOG_INFO_ON_ROOT("initialization of the distributions done")
      }

#if defined(WALBERLA_BUILD_WITH_CUDA)
      int streamLowPriority  = 0;
      int streamHighPriority = 0;
      auto defaultStream     = cuda::StreamRAII::newPriorityStream(streamLowPriority);
      auto innerOuterStreams = cuda::ParallelStreams(streamHighPriority);
#endif

      auto timeLoop = make_shared< SweepTimeloop >(blocks->getBlockStorage(), timesteps);
#if defined(WALBERLA_BUILD_WITH_CUDA)
      auto normalTimeStep = [&]() {
         for (auto& block : *blocks)
         {
            Comm_phase_field_distributions->communicate(nullptr);
            phase_field_LB_step(&block);

            Comm_velocity_based_distributions->communicate(nullptr);
            hydro_LB_step(&block);
         }
      };
      auto simpleOverlapTimeStep = [&]() {
         Comm_phase_field_distributions->startCommunication(defaultStream);
         for (auto& block : *blocks)
            phase_field_LB_step.inner(&block, defaultStream);
         Comm_phase_field_distributions->wait(defaultStream);
         for (auto& block : *blocks)
            phase_field_LB_step.outer(&block, defaultStream);

         Comm_velocity_based_distributions->startCommunication(defaultStream);
         for (auto& block : *blocks)
            hydro_LB_step.inner(&block, defaultStream);
         Comm_velocity_based_distributions->wait(defaultStream);
         for (auto& block : *blocks)
            hydro_LB_step.outer(&block, defaultStream);
      };
      auto phase_only = [&]() {
         for (auto& block : *blocks)
            phase_field_LB_step(&block);
      };
      auto hydro_only = [&]() {
         for (auto& block : *blocks)
            hydro_LB_step(&block);
      };
      auto without_comm = [&]() {
         for (auto& block : *blocks)
            phase_field_LB_step(&block);
         for (auto& block : *blocks)
            hydro_LB_step(&block);
      };
#else
      auto normalTimeStep = [&]() {
         for (auto& block : *blocks)
         {
            Comm_phase_field_distributions();
            phase_field_LB_step(&block);

            Comm_velocity_based_distributions();
            hydro_LB_step(&block);
         }
      };
      auto simpleOverlapTimeStep = [&]() {
        Comm_phase_field_distributions.startCommunication();
        for (auto& block : *blocks)
           phase_field_LB_step.inner(&block);
        Comm_phase_field_distributions.wait();
        for (auto& block : *blocks)
           phase_field_LB_step.outer(&block);

        Comm_velocity_based_distributions.startCommunication();
        for (auto& block : *blocks)
           hydro_LB_step.inner(&block);
        Comm_velocity_based_distributions.wait();
        for (auto& block : *blocks)
           hydro_LB_step.outer(&block);
      };
      auto phase_only = [&]() {
         for (auto& block : *blocks)
            phase_field_LB_step(&block);
      };
      auto hydro_only = [&]() {
         for (auto& block : *blocks)
            hydro_LB_step(&block);
      };
      auto without_comm = [&]() {
         for (auto& block : *blocks)
            phase_field_LB_step(&block);
         for (auto& block : *blocks)
            hydro_LB_step(&block);
      };
#endif
      std::function< void() > timeStep;
      if (timeStepStrategy == "overlap")
      {
         timeStep = std::function< void() >(simpleOverlapTimeStep);
         WALBERLA_LOG_INFO_ON_ROOT("overlapping timestep")
      }
      else if (timeStepStrategy == "phase_only")
      {
         timeStep = std::function< void() >(phase_only);
         WALBERLA_LOG_INFO_ON_ROOT("started only phasefield step without communication for benchmarking")
      }
      else if (timeStepStrategy == "hydro_only")
      {
         timeStep = std::function< void() >(hydro_only);
         WALBERLA_LOG_INFO_ON_ROOT("started only hydro step without communication for benchmarking")
      }
      else if (timeStepStrategy == "kernel_only")
      {
         timeStep = std::function< void() >(without_comm);
         WALBERLA_LOG_INFO_ON_ROOT("started complete phasefield model without communication for benchmarking")
      }
      else
      {
         timeStep = std::function< void() >(normalTimeStep);
         WALBERLA_LOG_INFO_ON_ROOT("normal timestep with no overlapping")
      }

      timeLoop->add() << BeforeFunction(timeStep) << Sweep([](IBlock*) {}, "time step");

      // remaining time logger
      if (remainingTimeLoggerFrequency > 0)
         timeLoop->addFuncAfterTimeStep(
            timing::RemainingTimeLogger(timeLoop->getNrOfTimeSteps(), remainingTimeLoggerFrequency),
            "remaining time logger");

      uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 1)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
#if defined(WALBERLA_BUILD_WITH_CUDA)
         vtkOutput->addBeforeFunction([&]() {
            cuda::fieldCpy< PhaseField_T, GPUField >(blocks, phase_field, phase_field_gpu);
            // cuda::fieldCpy<VelocityField_T, GPUField>( blocks, vel_field, vel_field_gpu );
         });
#endif
         auto phaseWriter = make_shared< field::VTKWriter< PhaseField_T > >(phase_field, "phase");
         vtkOutput->addCellDataWriter(phaseWriter);

         timeLoop->addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      for (uint_t i = 0; i < warmupSteps; ++i)
         timeLoop->singleStep();

      timeLoop->setCurrentTimeStepToZero();
      WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
      WcTimer simTimer;
#if defined(WALBERLA_BUILD_WITH_CUDA)
      cudaDeviceSynchronize();
#endif
      simTimer.start();
      timeLoop->run();
#if defined(WALBERLA_BUILD_WITH_CUDA)
      cudaDeviceSynchronize();
#endif
      simTimer.end();
      WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
      auto time            = simTimer.last();
      auto nrOfCells       = real_c(cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2]);
      auto mlupsPerProcess = nrOfCells * real_c(timesteps) / time * 1e-6;
      WALBERLA_LOG_RESULT_ON_ROOT("MLUPS per process: " << mlupsPerProcess)
      WALBERLA_LOG_RESULT_ON_ROOT("Time per time step: " << time / real_c(timesteps) << " s")
      WALBERLA_ROOT_SECTION()
      {
         python_coupling::PythonCallback pythonCallbackResults("results_callback");
         if (pythonCallbackResults.isCallable())
         {
            pythonCallbackResults.data().exposeValue("mlupsPerProcess", mlupsPerProcess);
            // Call Python function to report results
            pythonCallbackResults();
         }
      }
   }

   return EXIT_SUCCESS;
}
