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
//! \file thermocapillary.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"
#include "core/math/IntegerFactorization.h"
#include "core/timing/RemainingTimeLogger.h"

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
#include "gpu/GPURAII.h"
#include "gpu/GPUWrapper.h"
#include "gpu/ErrorChecking.h"
#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/communication/UniformGPUScheme.h"
#endif

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "lbm/PerformanceEvaluation.h"

#include "postprocessing/FieldToSurfaceMesh.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/SweepTimeloop.h"

#include "GenDefines.h"
#include "InitializerFunctions.h"
#include "PythonExports.h"

////////////
// USING //
//////////

using namespace std::placeholders;
using namespace walberla;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

auto DiffusionCallback = [](const Cell& pos, const shared_ptr< StructuredBlockForest >& blocks, IBlock& block,
                            const real_t Th, const real_t T0) {

      auto x_half = real_c(blocks->getDomainCellBB().xMax()) / real_c(2.0);
      Cell global_cell;
      blocks->transformBlockLocalToGlobalCell(global_cell, block, pos);

      return Th + T0 * cos( (math::pi / x_half) * (real_c(global_cell[0]) - x_half) );
};

int main(int argc, char** argv)
{
   mpi::Environment const Env(argc, argv);
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   gpu::selectDeviceBasedOnMpiRank();
#endif
   exportDataStructuresToPython();

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      auto config = *cfg;
      logging::configureLogging(config);

      auto DomainSetup = config->getOneBlock("DomainSetup");
      const Vector3<bool> periodic = DomainSetup.getParameter<Vector3<bool>>("periodic");
      const bool weakScaling = DomainSetup.getParameter<bool>("weakScaling", false); // weak or strong scaling

      const uint_t totalNumProcesses = uint_c(MPIManager::instance()->numProcesses());

      Vector3<uint_t> cellsPerBlock;
      Vector3<uint_t> NumberOfBlocks;
      Vector3<uint_t> numProcesses;

      if (!DomainSetup.isDefined("blocks"))
      {
         if (weakScaling)
         {
            const Vector3<uint_t> cells = DomainSetup.getParameter<Vector3<uint_t> >("cellsPerBlock");
            blockforest::calculateCellDistribution(cells, totalNumProcesses, NumberOfBlocks, cellsPerBlock);
            cellsPerBlock = cells;
            numProcesses = NumberOfBlocks;
         }
         else
         {
            cellsPerBlock = DomainSetup.getParameter<Vector3<uint_t> >("cellsPerBlock");
            const Vector3<uint_t> domainSize = DomainSetup.getParameter<Vector3<uint_t> >("domainSize");
            std::vector< uint_t > tmp = math::getFactors( totalNumProcesses, 3);
            // Round up to be divisible by eight (SIMD)
            for (uint_t i = 0; i < 3; ++i)
            {
               NumberOfBlocks[i] = uint_c(std::ceil(double_c(domainSize[i]) / double_c(cellsPerBlock[i])));
               numProcesses[i] = tmp[i];
            }
         }
      }
      else
      {
         cellsPerBlock = DomainSetup.getParameter<Vector3<uint_t>>("cellsPerBlock");
         NumberOfBlocks = DomainSetup.getParameter<Vector3<uint_t>>("blocks");
         numProcesses = NumberOfBlocks;
      }

      if ((NumberOfBlocks[0] * NumberOfBlocks[1] * NumberOfBlocks[2]) % totalNumProcesses != 0) {
         WALBERLA_ABORT("The total number of blocks is " << NumberOfBlocks[0] * NumberOfBlocks[1] * NumberOfBlocks[2]
                                                         << " they can not be equally distributed on " << totalNumProcesses
                                                         << " Processes")
      }

      auto aabb = AABB(real_c(0), real_c(0), real_c(0), real_c(NumberOfBlocks[0] * cellsPerBlock[0]),
                       real_c(NumberOfBlocks[1] * cellsPerBlock[1]),
                       real_c(NumberOfBlocks[2] * cellsPerBlock[2]));

      auto blocks = blockforest::createUniformBlockGrid( aabb,
                                                        NumberOfBlocks[0], NumberOfBlocks[1], NumberOfBlocks[2],
                                                        cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
                                                        numProcesses[0], numProcesses[1], numProcesses[2],
                                                        periodic[0], periodic[1], periodic[2], //periodicity
                                                        false // keepGlobalBlockInformation
      );
      blocks->createCellBoundingBoxes();


      ///////////////////////////////////
      // GENERAL SIMULATION PARAMETERS //
      //////////////////////////////////

      auto parameters                  = config->getOneBlock("Parameters");
      const uint_t timesteps           = parameters.getParameter< uint_t >("timesteps");
      const uint_t thermal_timesteps   = parameters.getParameter< uint_t >("pre_thermal_timesteps", uint_c(0));
      const bool heat_solver_RK_or_LBM = parameters.getParameter< bool >("HeatSolverRKOrLBM", false);
      const uint_t orderRKSolver = parameters.getParameter< uint_t >("orderRKSolver", 2);
      const real_t remaining_time_logger_frequency = parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(0.0));


      /////////////////////////
      // ADD DATA TO BLOCKS //
      ///////////////////////
      BlockDataID velocity_field_ID = field::addToStorage< VectorField_T >(blocks, "vel", real_c(0.0), field::fzyx);
      BlockDataID phase_field_ID    = field::addToStorage< ScalarField_T >(blocks, "phase", real_c(0.0), field::fzyx);
      BlockDataID temperature_field_ID = field::addToStorage< ScalarField_T >(blocks, "temperature", real_c(0.0), field::fzyx);

      BlockDataID temperature_PDFs_ID;
      BlockDataID rk_field_one_ID;
      BlockDataID rk_field_two_ID;
      BlockDataID rk_field_three_ID;

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      BlockDataID phase_field_gpu = gpu::addGPUFieldToStorage< ScalarField_T >(blocks, phase_field_ID, "phase field on GPU", true);
      BlockDataID phase_field_tmp = gpu::addGPUFieldToStorage< ScalarField_T >(blocks, phase_field_ID, "temporary phasefield", true);

      BlockDataID vel_field_gpu =
         gpu::addGPUFieldToStorage< VectorField_T >(blocks, velocity_field_ID, "velocity field on GPU", true);
      BlockDataID temperature_field_gpu =
         gpu::addGPUFieldToStorage< ScalarField_T >(blocks, temperature_field_ID, "temperature field on GPU", true);

      const BlockDataID allen_cahn_PDFs_ID = gpu::addGPUFieldToStorage< GPUFieldPDFs >(
         blocks, "lb phase field on GPU", Stencil_phase_T::Size, field::fzyx, uint_c(1));
      const BlockDataID hydrodynamic_PDFs_ID = gpu::addGPUFieldToStorage< GPUFieldPDFs >(
         blocks, "lb velocity field on GPU", Stencil_hydro_T::Size, field::fzyx, uint_c(1));

      if(heat_solver_RK_or_LBM)
      {
         rk_field_one_ID = gpu::addGPUFieldToStorage< gpu::GPUField< real_t > >(blocks, "RK1", 1, field::fzyx, 1);
         rk_field_two_ID = gpu::addGPUFieldToStorage< gpu::GPUField< real_t > >(blocks, "RK2", 1, field::fzyx, 1);
         rk_field_three_ID = gpu::addGPUFieldToStorage< gpu::GPUField< real_t > >(blocks, "RK3", 1, field::fzyx, 1);
      }
      else
      {
         temperature_PDFs_ID = gpu::addGPUFieldToStorage< GPUFieldPDFs >(
            blocks, "lb temperature field on GPU", Stencil_thermal_T::Size, field::fzyx, uint_c(1));
      }
#else
      BlockDataID allen_cahn_PDFs_ID =
         field::addToStorage< PdfField_phase_T >(blocks, "lb phase field", PdfField_phase_T::value_type(0.0), field::fzyx);
      BlockDataID hydrodynamic_PDFs_ID =
         field::addToStorage< PdfField_hydro_T >(blocks, "lb velocity field", PdfField_hydro_T::value_type(0.0), field::fzyx);
      BlockDataID phase_field_tmp_ID    = field::addToStorage< ScalarField_T >(blocks, "phase tmp", real_c(0.0), field::fzyx);

      if(heat_solver_RK_or_LBM)
      {
         rk_field_one_ID = field::addToStorage< ScalarField_T >(blocks, "RK1", real_c(0.0), field::fzyx);
         rk_field_two_ID = field::addToStorage< ScalarField_T >(blocks, "RK2", real_c(0.0), field::fzyx);
         rk_field_three_ID = field::addToStorage< ScalarField_T >(blocks, "RK3", real_c(0.0), field::fzyx);
      }
      else
      {
         temperature_PDFs_ID = field::addToStorage< PdfField_thermal_T >(blocks, "lb temperature field", PdfField_thermal_T::value_type(0.0), field::fzyx);
      }

#endif

      const BlockDataID flag_field_allen_cahn_LB_ID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field phase");
      const BlockDataID flag_field_hydrodynamic_LB_ID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field hydro");
      const BlockDataID flag_field_thermal_LB_ID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field thermal");

      auto physical_parameters  = config->getOneBlock("PhysicalParameters");
      const real_t density_liquid           = physical_parameters.getParameter< real_t >("density_liquid");
      const real_t density_gas              = physical_parameters.getParameter< real_t >("density_gas");
      const real_t sigma_ref                = physical_parameters.getParameter< real_t >("sigma_ref");
      const real_t sigma_t                  = physical_parameters.getParameter< real_t >("sigma_t");
      const real_t mobility                 = physical_parameters.getParameter< real_t >("mobility");
      const real_t temperature_ref          = physical_parameters.getParameter< real_t >("temperature_ref");
      const real_t heat_conductivity_liquid = physical_parameters.getParameter< real_t >("heat_conductivity_liquid");
      const real_t heat_conductivity_gas    = physical_parameters.getParameter< real_t >("heat_conductivity_gas");
      const real_t relaxation_time_liquid   = physical_parameters.getParameter< real_t >("relaxation_time_liquid");
      const real_t relaxation_time_gas      = physical_parameters.getParameter< real_t >("relaxation_time_gas");
      const real_t interface_thickness      = physical_parameters.getParameter< real_t >("interface_thickness");
      const real_t velocity_wall            = physical_parameters.getParameter< real_t >("velocity_wall");
      const real_t angle                    = physical_parameters.getParameter< real_t >("contact_angle", real_c(90));
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      Vector3< int32_t > gpuBlockSize = parameters.getParameter< Vector3< int32_t > >("gpuBlockSize", Vector3< int32_t >(128, 1, 1));
#endif

      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the phase-field")
      auto t_h = real_c(0.0);
      auto t_0 = real_c(0.0);

#ifdef GENERATED_HEAT_SOURCE
      auto ws = real_c(0.0);
      auto ds = real_c(0.0);
      auto qsOne = real_c(0.0);
      auto qsTwo = real_c(0.0);
      Vector3< int64_t > heatSourceMidPointOne(0, 0, 0);
      Vector3< int64_t > heatSourceMidPointTwo(0, 0, 0);
#endif

      bool benchmark = false;
      std::string timestep_strategy = "NormalTimestep";
      if (config->getNumBlocks("Droplet") == 1)
      {
         if(!GeneratedHeatSource)
            WALBERLA_ABORT("For the Droplet application 'heat_source' has to be activated in the generation script")
         if(heat_solver_RK_or_LBM)
            WALBERLA_ABORT("Runge Kutta method is only available for the MicroChannel test case so far")
         auto droplet_parameters    = config->getOneBlock("Droplet");
         const real_t droplet_radius = droplet_parameters.getParameter< real_t >("dropletRadius");
         const Vector3< real_t > droplet_midpoint =
            droplet_parameters.getParameter< Vector3< real_t > >("dropletMidPoint");
         initPhaseFieldDroplet(blocks, phase_field_ID, droplet_radius, droplet_midpoint, interface_thickness);

#ifdef GENERATED_HEAT_SOURCE
         auto HeatSource    = config->getOneBlock("HeatSource");
         ws = HeatSource.getParameter< real_t >("ws");
         ds = HeatSource.getParameter< real_t >("sizeDiffusedHotspot");

         qsOne = HeatSource.getParameter< real_t >("maximumHeatFluxOne");
         qsTwo = HeatSource.getParameter< real_t >("maximumHeatFluxTwo");

         heatSourceMidPointOne = HeatSource.getParameter< Vector3< int64_t > >("heatSourceMidPointOne");
         heatSourceMidPointTwo = HeatSource.getParameter< Vector3< int64_t > >("heatSourceMidPointTwo");
#endif
      }
      else if (config->getNumBlocks("MicroChannel") == 1)
      {
         if(GeneratedHeatSource)
            WALBERLA_ABORT("For the MicroChannel application 'heat_source' has to be deactivated in the generation script")
         auto micro_channel_parameters = config->getOneBlock("MicroChannel");
         t_h                         = micro_channel_parameters.getParameter< real_t >("Th");
         t_0                         = micro_channel_parameters.getParameter< real_t >("T0");
         initMicroChannel(blocks, phase_field_ID, temperature_field_ID, t_h, t_0, temperature_ref, interface_thickness);
      }
      else if (config->getNumBlocks("Benchmark") == 1)
      {
         auto benchmark_parameters = config->getOneBlock("Benchmark");
         benchmark = true;
         timestep_strategy = benchmark_parameters.getParameter<std::string>("timeStepStrategy", "NormalTimestep");
         WALBERLA_LOG_INFO_ON_ROOT("Benchmarking Thermocapillary flows with: " << timestep_strategy)
      }
      else
      {
         WALBERLA_ABORT("Only Droplet or MicroChannel or Benchmark Scenario is available. Thus a Parameter Block for these "
                        "Scenarios must exist!")
      }
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the phase-field done")

      /////////////////
      // ADD SWEEPS //
      ///////////////
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)

      int streamHighPriority = 0;
      int streamLowPriority = 0;
      WALBERLA_GPU_CHECK( gpuDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority) );
      auto defaultStream = gpu::StreamRAII::newPriorityStream( streamLowPriority );

      pystencils::initialize_phase_field_distributions init_h(allen_cahn_PDFs_ID, phase_field_gpu, vel_field_gpu,
                                                              interface_thickness);
      pystencils::initialize_velocity_based_distributions init_g(hydrodynamic_PDFs_ID, vel_field_gpu);
      pystencils::initialize_thermal_distributions init_f(temperature_PDFs_ID, temperature_field_gpu, vel_field_gpu);


      pystencils::phase_field_LB_step phase_field_LB_step(allen_cahn_PDFs_ID, phase_field_gpu, phase_field_tmp, vel_field_gpu, mobility,
                                                          interface_thickness, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
      pystencils::hydro_LB_step hydro_LB_step(hydrodynamic_PDFs_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu,
                                              temperature_ref, interface_thickness, density_liquid, density_gas,
                                              sigma_t, sigma_ref, relaxation_time_liquid, relaxation_time_gas, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);

#if defined(GENERATED_HEAT_SOURCE)
      pystencils::thermal_LB_step thermal_LB_step(temperature_PDFs_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu, ds,
                                                  heatSourceMidPointOne[0], heatSourceMidPointOne[1], heatSourceMidPointOne[2],
                                                  heatSourceMidPointTwo[0], heatSourceMidPointTwo[1], heatSourceMidPointTwo[2],
                                                  heat_conductivity_liquid, heat_conductivity_gas,
                                                  qsOne, qsTwo, ws,
                                                  gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
#else
      pystencils::thermal_LB_step thermal_LB_step(temperature_PDFs_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu,
                                                  heat_conductivity_liquid, heat_conductivity_gas, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
#endif
      pystencils::initialize_rk2 initRK2(rk_field_one_ID, temperature_field_gpu, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
      pystencils::rk2_first_stage RK2FirstStage(rk_field_one_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu, heat_conductivity_liquid,
                                                heat_conductivity_gas, density_liquid, density_gas, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
      pystencils::rk2_second_stage RK2SecondStage(rk_field_one_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu,
                                                  heat_conductivity_liquid, heat_conductivity_gas, density_liquid, density_gas,
                                                  gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);

      pystencils::initialize_rk4 initRK4(rk_field_one_ID, rk_field_two_ID, rk_field_three_ID, temperature_field_gpu, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
      pystencils::rk4_first_stage RK4FirstStage(rk_field_one_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu, heat_conductivity_liquid,
                                                heat_conductivity_gas, density_liquid, density_gas,
                                                gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
      pystencils::rk4_second_stage RK4SecondStage(rk_field_one_ID, rk_field_two_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu,
                                                  heat_conductivity_liquid, heat_conductivity_gas, density_liquid, density_gas,
                                                  gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
      pystencils::rk4_third_stage RK4ThirdStage(rk_field_two_ID, rk_field_three_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu,
                                                heat_conductivity_liquid, heat_conductivity_gas, density_liquid, density_gas,
                                                gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
      pystencils::rk4_fourth_stage RK4FourthStage(rk_field_one_ID, rk_field_two_ID, rk_field_three_ID, phase_field_gpu, temperature_field_gpu, vel_field_gpu,
                                                  heat_conductivity_liquid, heat_conductivity_gas, density_liquid, density_gas,
                                                  gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);

      auto swapPhaseField = std::function< void(IBlock *) >([&](IBlock * b)
      {
            auto phaseField = b->getData< gpu::GPUField<real_t> >(phase_field_gpu);
            auto phaseFieldTMP = b->getData< gpu::GPUField<real_t> >(phase_field_tmp);
            phaseField->swapDataPointers(phaseFieldTMP);
      });

#else
      pystencils::initialize_phase_field_distributions init_h(allen_cahn_PDFs_ID, phase_field_ID, velocity_field_ID,
                                                              interface_thickness);
      pystencils::initialize_velocity_based_distributions init_g(hydrodynamic_PDFs_ID, velocity_field_ID);
      pystencils::initialize_thermal_distributions init_f(temperature_PDFs_ID, temperature_field_ID, velocity_field_ID);

      pystencils::phase_field_LB_step phase_field_LB_step(allen_cahn_PDFs_ID, phase_field_ID, phase_field_tmp_ID, velocity_field_ID, mobility,
                                                          interface_thickness);
      pystencils::hydro_LB_step hydro_LB_step(hydrodynamic_PDFs_ID, phase_field_ID, temperature_field_ID, velocity_field_ID,
                                              temperature_ref, interface_thickness, density_liquid, density_gas,
                                              sigma_t, sigma_ref, relaxation_time_liquid, relaxation_time_gas);
#if defined(GENERATED_HEAT_SOURCE)
      pystencils::thermal_LB_step thermal_LB_step(temperature_PDFs_ID, phase_field_ID, temperature_field_ID, velocity_field_ID, ds,
                                                  heatSourceMidPointOne[0], heatSourceMidPointOne[1], heatSourceMidPointOne[2],
                                                  heatSourceMidPointTwo[0], heatSourceMidPointTwo[1], heatSourceMidPointTwo[2],
                                                  heat_conductivity_liquid, heat_conductivity_gas,
                                                  qsOne, qsTwo, ws);
#else
      pystencils::thermal_LB_step thermal_LB_step(temperature_PDFs_ID, phase_field_ID, temperature_field_ID, velocity_field_ID,
                                                  heat_conductivity_liquid, heat_conductivity_gas);
#endif
      pystencils::initialize_rk2 initRK2(rk_field_one_ID, temperature_field_ID);
      pystencils::rk2_first_stage RK2FirstStage(rk_field_one_ID, phase_field_ID, temperature_field_ID, velocity_field_ID, heat_conductivity_liquid,
                                                heat_conductivity_gas, density_liquid, density_gas);
      pystencils::rk2_second_stage RK2SecondStage(rk_field_one_ID, phase_field_ID, temperature_field_ID, velocity_field_ID,
                                                  heat_conductivity_liquid, heat_conductivity_gas, density_liquid, density_gas);

      pystencils::initialize_rk4 initRK4(rk_field_one_ID, rk_field_two_ID, rk_field_three_ID, temperature_field_ID);
      pystencils::rk4_first_stage RK4FirstStage(rk_field_one_ID, phase_field_ID, temperature_field_ID, velocity_field_ID, heat_conductivity_liquid,
                                                heat_conductivity_gas, density_liquid, density_gas);
      pystencils::rk4_second_stage RK4SecondStage(rk_field_one_ID, rk_field_two_ID, phase_field_ID, temperature_field_ID, velocity_field_ID,
                                                  heat_conductivity_liquid, heat_conductivity_gas, density_liquid, density_gas);
      pystencils::rk4_third_stage RK4ThirdStage(rk_field_two_ID, rk_field_three_ID, phase_field_ID, temperature_field_ID, velocity_field_ID,
                                                heat_conductivity_liquid, heat_conductivity_gas, density_liquid, density_gas);
      pystencils::rk4_fourth_stage RK4FourthStage(rk_field_one_ID, rk_field_two_ID, rk_field_three_ID, phase_field_ID, temperature_field_ID, velocity_field_ID,
                                                  heat_conductivity_liquid, heat_conductivity_gas, density_liquid, density_gas);

      auto swapPhaseField = std::function< void(IBlock *) >([&](IBlock * b)
      {
         auto phaseField = b->getData< ScalarField_T >(phase_field_ID);
         auto phaseFieldTMP = b->getData< ScalarField_T >(phase_field_tmp_ID);
         phaseField->swapDataPointers(phaseFieldTMP);
      });
#endif

      for (auto& block : *blocks)
      {
         thermal_LB_step.configure(blocks, &block);
      }

      ////////////////////////
      // ADD COMMUNICATION //
      //////////////////////
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      const bool gpuEnabledMpi = parameters.getParameter< bool >("gpuEnabledMpi", false);
      if(gpuEnabledMpi)
         WALBERLA_LOG_INFO_ON_ROOT("GPU-Direct communication is used for MPI")
      else
         WALBERLA_LOG_INFO_ON_ROOT("No GPU-Direct communication is used for MPI")

      auto generatedPackInfo_phase_field_distributions = make_shared< lbm::PackInfo_phase_field_distributions>(allen_cahn_PDFs_ID);
      auto generatedPackInfo_velocity_based_distributions = make_shared< lbm::PackInfo_velocity_based_distributions >(hydrodynamic_PDFs_ID);
      auto generatedPackInfo_thermal_distributions = make_shared< lbm::PackInfo_thermal_field_distributions>(temperature_PDFs_ID);
      auto generatedPackInfo_phase_field = make_shared< pystencils::PackInfo_phase_field >(phase_field_gpu);
      auto generatedPackInfo_temperature_field = make_shared< pystencils::PackInfo_temperature_field >(temperature_field_gpu);

      auto UniformGPUSchemeVelocityBasedDistributions = make_shared< gpu::communication::UniformGPUScheme< CommStencil_hydro_T > >(blocks, gpuEnabledMpi, false);
      auto UniformGPUSchemePhaseFieldDistributions = make_shared< gpu::communication::UniformGPUScheme< CommStencil_phase_T > >(blocks, gpuEnabledMpi, false);
      auto UniformGPUSchemeThermalDistributions = make_shared< gpu::communication::UniformGPUScheme< CommStencil_phase_T > >(blocks, gpuEnabledMpi, false);
      auto UniformGPUSchemePhaseField = make_shared< gpu::communication::UniformGPUScheme< Full_Stencil_T > >(blocks, gpuEnabledMpi, false, 65432);
      auto UniformGPUSchemeTemperatureField = make_shared< gpu::communication::UniformGPUScheme< Full_Stencil_T > >(blocks, gpuEnabledMpi, false, 65432);

      UniformGPUSchemeVelocityBasedDistributions->addPackInfo(generatedPackInfo_velocity_based_distributions);
      UniformGPUSchemeVelocityBasedDistributions->addPackInfo(generatedPackInfo_phase_field);
      UniformGPUSchemePhaseFieldDistributions->addPackInfo(generatedPackInfo_phase_field_distributions);
      UniformGPUSchemeThermalDistributions->addPackInfo(generatedPackInfo_thermal_distributions);
      UniformGPUSchemeThermalDistributions->addPackInfo(generatedPackInfo_temperature_field);
      UniformGPUSchemePhaseField->addPackInfo(generatedPackInfo_phase_field);
      UniformGPUSchemeTemperatureField->addPackInfo(generatedPackInfo_temperature_field);

      auto Comm_velocity_based_distributions = std::function< void() >([&]() { UniformGPUSchemeVelocityBasedDistributions->communicate(); });
      auto Comm_velocity_based_distributions_start = std::function< void() >([&]() { UniformGPUSchemeVelocityBasedDistributions->startCommunication(); });
      auto Comm_velocity_based_distributions_wait = std::function< void() >([&]() { UniformGPUSchemeVelocityBasedDistributions->wait(); });

      auto Comm_phase_field_distributions = std::function< void() >([&]() { UniformGPUSchemePhaseFieldDistributions->communicate(); });
      auto Comm_phase_field_distributions_start = std::function< void() >([&]() { UniformGPUSchemePhaseFieldDistributions->startCommunication(); });
      auto Comm_phase_field_distributions_wait = std::function< void() >([&]() { UniformGPUSchemePhaseFieldDistributions->wait(); });

      auto Comm_thermal_distributions = std::function< void() >([&]() { UniformGPUSchemeThermalDistributions->communicate(); });
      auto Comm_thermal_distributions_start = std::function< void() >([&]() { UniformGPUSchemeThermalDistributions->startCommunication(); });
      auto Comm_thermal_distributions_wait = std::function< void() >([&]() { UniformGPUSchemeThermalDistributions->wait(); });

      auto Comm_phase_field = std::function< void() >([&]() { UniformGPUSchemePhaseField->communicate(); });
      auto Comm_phase_field_start = std::function< void() >([&]() { UniformGPUSchemePhaseField->startCommunication(); });
      auto Comm_phase_field_wait = std::function< void() >([&]() { UniformGPUSchemePhaseField->wait(); });

      auto Comm_temperature_field = std::function< void() >([&]() { UniformGPUSchemeTemperatureField->communicate(); });
      auto Comm_temperature_field_start = std::function< void() >([&]() { UniformGPUSchemeTemperatureField->startCommunication(); });
      auto Comm_temperature_field_wait = std::function< void() >([&]() { UniformGPUSchemeTemperatureField->wait(); });

#else
      auto UniformBufferedSchemeVelocityDistributions = make_shared<blockforest::communication::UniformBufferedScheme< CommStencil_hydro_T >> (blocks, 1111);
      auto generatedPackInfo_velocity_based_distributions =
         make_shared< lbm::PackInfo_velocity_based_distributions >(hydrodynamic_PDFs_ID);
      UniformBufferedSchemeVelocityDistributions->addPackInfo(generatedPackInfo_velocity_based_distributions);
      auto Comm_velocity_based_distributions = std::function< void() >([&]() { UniformBufferedSchemeVelocityDistributions->communicate(); });
      auto Comm_velocity_based_distributions_start = std::function< void() >([&]() { UniformBufferedSchemeVelocityDistributions->startCommunication(); });
      auto Comm_velocity_based_distributions_wait = std::function< void() >([&]() { UniformBufferedSchemeVelocityDistributions->wait(); });

      auto UniformBufferedSchemePhaseField = make_shared<blockforest::communication::UniformBufferedScheme< Full_Stencil_T >> (blocks, 2222);
      auto generatedPackInfo_phase_field = make_shared< pystencils::PackInfo_phase_field >(phase_field_ID);
      UniformBufferedSchemePhaseField->addPackInfo(generatedPackInfo_phase_field);
      auto Comm_phase_field = std::function< void() >([&]() { UniformBufferedSchemePhaseField->communicate(); });
      auto Comm_phase_field_start = std::function< void() >([&]() { UniformBufferedSchemePhaseField->startCommunication(); });
      auto Comm_phase_field_wait = std::function< void() >([&]() { UniformBufferedSchemePhaseField->wait(); });

      auto UniformBufferedSchemeTemperatureField = make_shared<blockforest::communication::UniformBufferedScheme< Full_Stencil_T >>(blocks, 3333);
      auto generatedPackInfo_temperature_field =
         make_shared< pystencils::PackInfo_temperature_field >(temperature_field_ID);
      UniformBufferedSchemeTemperatureField->addPackInfo(generatedPackInfo_temperature_field);
      auto Comm_temperature_field = std::function< void() >([&]() { UniformBufferedSchemeTemperatureField->communicate(); });
      auto Comm_temperature_field_start = std::function< void() >([&]() { UniformBufferedSchemeTemperatureField->startCommunication(); });
      auto Comm_temperature_field_wait = std::function< void() >([&]() { UniformBufferedSchemeTemperatureField->wait(); });

      auto UniformBufferedSchemePhaseFieldDistributions = make_shared<blockforest::communication::UniformBufferedScheme< CommStencil_phase_T >>(blocks, 4444);
      auto generatedPackInfo_phase_field_distributions =
         make_shared< lbm::PackInfo_phase_field_distributions >(allen_cahn_PDFs_ID);
      UniformBufferedSchemePhaseFieldDistributions->addPackInfo(generatedPackInfo_phase_field_distributions);
      auto Comm_phase_field_distributions = std::function< void() >([&]() { UniformBufferedSchemePhaseFieldDistributions->communicate(); });
      auto Comm_phase_field_distributions_start = std::function< void() >([&]() { UniformBufferedSchemePhaseFieldDistributions->startCommunication(); });
      auto Comm_phase_field_distributions_wait = std::function< void() >([&]() { UniformBufferedSchemePhaseFieldDistributions->wait(); });

      auto UniformBufferedSchemeThermalDistributions = make_shared<blockforest::communication::UniformBufferedScheme< CommStencil_thermal_T >>(blocks, 5555);
      auto generatedPackInfo_thermal_distributions =
         make_shared< lbm::PackInfo_thermal_field_distributions >(temperature_PDFs_ID);
      UniformBufferedSchemeThermalDistributions->addPackInfo(generatedPackInfo_thermal_distributions);
      auto Comm_thermal_distributions = std::function< void() >([&]() { UniformBufferedSchemeThermalDistributions->communicate(); });
      auto Comm_thermal_distributions_start = std::function< void() >([&]() { UniformBufferedSchemeThermalDistributions->startCommunication(); });
      auto Comm_thermal_distributions_wait = std::function< void() >([&]() { UniformBufferedSchemeThermalDistributions->wait(); });
#endif

      ////////////////////////
      // BOUNDARY HANDLING //
      //////////////////////

      const FlagUID fluidFlagUID("Fluid");
      const FlagUID wallFlagUID("NoSlip");
      const FlagUID ubbFlagUID("UBB");
      const FlagUID diffusionStaticFlagUID("DiffusionDirichletStatic");
      const FlagUID diffusionDynamicFlagUID("DiffusionDirichletDynamic");
      const FlagUID neumannFlagUID("NeumannByCopy");

      auto boundariesConfigPhase = config->getBlock("BoundariesAllenCahn");
      if (boundariesConfigPhase)
      {
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flag_field_allen_cahn_LB_ID, boundariesConfigPhase);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flag_field_allen_cahn_LB_ID, fluidFlagUID);
      }

      auto boundariesConfigHydro = config->getBlock("BoundariesHydro");
      if (boundariesConfigHydro)
      {
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flag_field_hydrodynamic_LB_ID, boundariesConfigHydro);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flag_field_hydrodynamic_LB_ID, fluidFlagUID);
      }

      auto boundariesConfigThermal = config->getBlock("BoundariesThermal");
      if (boundariesConfigPhase)
      {
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flag_field_thermal_LB_ID, boundariesConfigThermal);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flag_field_thermal_LB_ID, fluidFlagUID);
      }

      lbm::phase_field_LB_NoSlip phase_field_LB_NoSlip(blocks, allen_cahn_PDFs_ID);
      phase_field_LB_NoSlip.fillFromFlagField< FlagField_T >(blocks, flag_field_allen_cahn_LB_ID, wallFlagUID, fluidFlagUID);

      lbm::hydro_LB_NoSlip hydro_LB_NoSlip(blocks, hydrodynamic_PDFs_ID);
      hydro_LB_NoSlip.fillFromFlagField< FlagField_T >(blocks, flag_field_hydrodynamic_LB_ID, wallFlagUID, fluidFlagUID);
      lbm::hydro_LB_UBB hydro_LB_UBB(blocks, hydrodynamic_PDFs_ID, velocity_wall);
      hydro_LB_UBB.fillFromFlagField< FlagField_T >(blocks, flag_field_hydrodynamic_LB_ID, ubbFlagUID, fluidFlagUID);


      std::function< real_t(const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
         DiffusionInitialisation = std::bind(DiffusionCallback, _1, _2, _3, t_h, t_0);

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      lbm::thermal_LB_DiffusionDirichlet_static thermal_LB_DiffusionStatic(blocks, temperature_PDFs_ID, vel_field_gpu, temperature_ref);
      thermal_LB_DiffusionStatic.fillFromFlagField< FlagField_T >(blocks, flag_field_thermal_LB_ID, diffusionStaticFlagUID, fluidFlagUID);
      lbm::thermal_LB_DiffusionDirichlet_dynamic thermal_LB_DiffusionDynamic(blocks, temperature_PDFs_ID, vel_field_gpu, DiffusionInitialisation);
      thermal_LB_DiffusionDynamic.fillFromFlagField< FlagField_T >(blocks, flag_field_thermal_LB_ID, diffusionDynamicFlagUID, fluidFlagUID);
      lbm::thermal_LB_NeumannByCopy thermal_LB_Neumann(blocks, temperature_PDFs_ID);
      thermal_LB_Neumann.fillFromFlagField< FlagField_T >(blocks, flag_field_thermal_LB_ID, neumannFlagUID, fluidFlagUID);

      pystencils::ContactAngle contact_angle(blocks, phase_field_gpu, interface_thickness, angle);
      contact_angle.fillFromFlagField< FlagField_T >(blocks, flag_field_allen_cahn_LB_ID, wallFlagUID, fluidFlagUID);
#else
      lbm::thermal_LB_DiffusionDirichlet_static thermal_LB_DiffusionStatic(blocks, temperature_PDFs_ID,
                                                                           velocity_field_ID, temperature_ref);
      thermal_LB_DiffusionStatic.fillFromFlagField< FlagField_T >(blocks, flag_field_thermal_LB_ID, diffusionStaticFlagUID, fluidFlagUID);
      lbm::thermal_LB_DiffusionDirichlet_dynamic thermal_LB_DiffusionDynamic(blocks, temperature_PDFs_ID, velocity_field_ID, DiffusionInitialisation);
      thermal_LB_DiffusionDynamic.fillFromFlagField< FlagField_T >(blocks, flag_field_thermal_LB_ID, diffusionDynamicFlagUID, fluidFlagUID);
      lbm::thermal_LB_NeumannByCopy thermal_LB_Neumann(blocks, temperature_PDFs_ID);
      thermal_LB_Neumann.fillFromFlagField< FlagField_T >(blocks, flag_field_thermal_LB_ID, neumannFlagUID, fluidFlagUID);

      pystencils::ContactAngle contact_angle(blocks, phase_field_ID, interface_thickness, angle);
      contact_angle.fillFromFlagField< FlagField_T >(blocks, flag_field_allen_cahn_LB_ID, wallFlagUID, fluidFlagUID);
#endif


#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      ////////////////////////////
      // TIMELOOP THERMAL ONLY //
      //////////////////////////

      SweepTimeloop thermalTimeloop(blocks->getBlockStorage(), thermal_timesteps);

      if(heat_solver_RK_or_LBM)
      {
         if (orderRKSolver == 2)
         {
            thermalTimeloop.add() << Sweep(initRK2.getSweep(defaultStream), "initialise RK");
            thermalTimeloop.add() << Sweep(RK2FirstStage.getSweep(defaultStream), "Runge Kutta step one");
            thermalTimeloop.add() << Sweep(RK2SecondStage.getSweep(defaultStream), "Runge Kutta step two")
                                  << AfterFunction(Comm_temperature_field, "Temperature field Communication");
         }
         else
         {
            thermalTimeloop.add() << Sweep(initRK4.getSweep(defaultStream), "initialise RK");
            thermalTimeloop.add() << Sweep(RK4FirstStage.getSweep(defaultStream), "Runge Kutta step one");
            thermalTimeloop.add() << Sweep(RK4SecondStage.getSweep(defaultStream), "Runge Kutta step two");
            thermalTimeloop.add() << Sweep(RK4ThirdStage.getSweep(defaultStream), "Runge Kutta step three");
            thermalTimeloop.add() << Sweep(RK4FourthStage.getSweep(defaultStream), "Runge Kutta step four")
                                  << AfterFunction(Comm_temperature_field, "Temperature field Communication");
         }
      }
      else
      {
         thermalTimeloop.add() << BeforeFunction(Comm_thermal_distributions, "Thermal PDFs Communication")
                               << Sweep(thermal_LB_Neumann.getSweep(defaultStream), "Neumann Thermal");
         thermalTimeloop.add() << Sweep(thermal_LB_DiffusionStatic.getSweep(defaultStream), "Static Diffusion Thermal");
         thermalTimeloop.add() << Sweep(thermal_LB_DiffusionDynamic.getSweep(defaultStream), "Dynamic Diffusion Thermal");
         thermalTimeloop.add() << Sweep(thermal_LB_step.getSweep(defaultStream), "Thermal LB Step");
      }

      ////////////////
      // TIME LOOP //
      //////////////

      SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);
      if(heat_solver_RK_or_LBM)
      {
         if (timestep_strategy == "NormalTimestep")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Normal time step for production runs and full application benchmarks with " << std::to_string(orderRKSolver) << ". order RK scheme as thermal solver")
            timeloop.add() << BeforeFunction(Comm_velocity_based_distributions_start, "Start Hydro PDFs Communication")
                           << BeforeFunction(Comm_temperature_field_start, "Start Communication temperature field")
                           << Sweep(phase_field_LB_NoSlip.getSweep(defaultStream), "NoSlip Phase");
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(defaultStream), "Phase LB Step")
                           << AfterFunction(Comm_velocity_based_distributions_wait, "Wait Hydro PDFs Communication")
                           << AfterFunction(Comm_temperature_field_wait, "Wait Communication temperature field");

            timeloop.add() << BeforeFunction(Comm_phase_field_distributions_start, "Start Phase PDFs Communication")
                           << Sweep(hydro_LB_UBB.getSweep(defaultStream), "UBB Hydro");
            timeloop.add() << Sweep(hydro_LB_NoSlip.getSweep(defaultStream), "NoSlip Hydro");
            timeloop.add() << Sweep(hydro_LB_step.getSweep(defaultStream), "Hydro LB Step");
            timeloop.add() << Sweep(swapPhaseField, "Swap PhaseField");
            timeloop.add() << Sweep(contact_angle.getSweep(defaultStream), "Contact Angle")
                           << AfterFunction(Comm_phase_field_distributions_wait, "Wait Phase PDFs Communication");

            if(orderRKSolver == 2)
            {
               timeloop.add() << BeforeFunction(Comm_phase_field, "Communication Phase Field")
                              << Sweep(RK2FirstStage.getSweep(defaultStream), "Runge Kutta step one");
               timeloop.add() << Sweep(RK2SecondStage.getSweep(defaultStream), "Runge Kutta step two");
            }
            else
            {
               timeloop.add() << BeforeFunction(Comm_phase_field, "Communication Phase Field")
                              << Sweep(RK4FirstStage.getSweep(defaultStream), "Runge Kutta step one");
               timeloop.add() << Sweep(RK4SecondStage.getSweep(defaultStream), "Runge Kutta step two");
               timeloop.add() << Sweep(RK4ThirdStage.getSweep(defaultStream), "Runge Kutta step three");
               timeloop.add() << Sweep(RK4FourthStage.getSweep(defaultStream), "Runge Kutta step four");

            }
         }
         else if (timestep_strategy == "PhaseOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the Allen Cahn LB kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(defaultStream), "Phase LB Step");
         }
         else if (timestep_strategy == "HydroOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the hydrodynamic LB kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(hydro_LB_step.getSweep(defaultStream), "Hydro LB Step");
         }
         else if (timestep_strategy == "ThermalOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the RK kernels. This makes only sense for benchmarking")
            if(orderRKSolver == 2)
            {
               timeloop.add() << Sweep(RK2FirstStage.getSweep(defaultStream), "Runge Kutta step one");
               timeloop.add() << Sweep(RK2SecondStage.getSweep(defaultStream), "Runge Kutta step two");
            }
            else
            {
               timeloop.add() << Sweep(RK4FirstStage.getSweep(defaultStream), "Runge Kutta step one");
               timeloop.add() << Sweep(RK4SecondStage.getSweep(defaultStream), "Runge Kutta step two");
               timeloop.add() << Sweep(RK4ThirdStage.getSweep(defaultStream), "Runge Kutta step three");
               timeloop.add() << Sweep(RK4FourthStage.getSweep(defaultStream), "Runge Kutta step four");

            }
         }
         else if (timestep_strategy == "KernelOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the compute kernels without boundary conditions and communication. This makes only sense for benchmarking")
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(defaultStream), "Phase LB Step");
            timeloop.add() << Sweep(hydro_LB_step.getSweep(defaultStream), "Hydro LB Step");
            if(orderRKSolver == 2)
            {
               timeloop.add() << Sweep(RK2FirstStage.getSweep(defaultStream), "Runge Kutta step one");
               timeloop.add() << Sweep(RK2SecondStage.getSweep(defaultStream), "Runge Kutta step two");
            }
            else
            {
               timeloop.add() << Sweep(RK4FirstStage.getSweep(defaultStream), "Runge Kutta step one");
               timeloop.add() << Sweep(RK4SecondStage.getSweep(defaultStream), "Runge Kutta step two");
               timeloop.add() << Sweep(RK4ThirdStage.getSweep(defaultStream), "Runge Kutta step three");
               timeloop.add() << Sweep(RK4FourthStage.getSweep(defaultStream), "Runge Kutta step four");

            }
         }
         else
         {
            WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'")
         }
      }
      else
      {
         if (timestep_strategy == "NormalTimestep")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Normal time step for production runs and full application benchmarks with LBM thermal solver")
            timeloop.add() << BeforeFunction(Comm_thermal_distributions_start, "Start Thermal PDFs Communication")
                           << Sweep(phase_field_LB_NoSlip.getSweep(defaultStream), "NoSlip Phase");
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(defaultStream), "Phase LB Step")
                           << AfterFunction(Comm_thermal_distributions_wait, "Wait Thermal PDFs Communication");

            timeloop.add() << BeforeFunction(Comm_phase_field_distributions_start, "Start Phase PDFs Communication")
                           << Sweep(hydro_LB_UBB.getSweep(defaultStream), "UBB Hydro");
            timeloop.add() << Sweep(hydro_LB_NoSlip.getSweep(defaultStream), "NoSlip Hydro");
            timeloop.add() << Sweep(hydro_LB_step.getSweep(defaultStream), "Hydro LB Step");
            timeloop.add() << Sweep(swapPhaseField, "Swap PhaseField");
            timeloop.add() << Sweep(contact_angle.getSweep(defaultStream), "Contact Angle")
                           << AfterFunction(Comm_phase_field_distributions_wait, "Wait Phase PDFs Communication");

            timeloop.add() << BeforeFunction(Comm_velocity_based_distributions_start, "Start Hydro PDFs Communication")
                           << Sweep(thermal_LB_Neumann.getSweep(defaultStream), "Neumann Thermal");
            timeloop.add() << Sweep(thermal_LB_DiffusionStatic.getSweep(defaultStream), "Static Diffusion Thermal");
            timeloop.add() << Sweep(thermal_LB_DiffusionDynamic.getSweep(defaultStream), "Dynamic Diffusion Thermal");
            timeloop.add() << Sweep(thermal_LB_step.getSweep(defaultStream), "Thermal LB Step")
                           << AfterFunction(Comm_velocity_based_distributions_wait, "Wait Hydro PDFs Communication");
         }
         else if (timestep_strategy == "PhaseOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the Allen Cahn LB kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(defaultStream), "Phase LB Step");
         }
         else if (timestep_strategy == "HydroOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the hydrodynamic LB kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(hydro_LB_step.getSweep(defaultStream), "Hydro LB Step");
         }
         else if (timestep_strategy == "ThermalOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the thermal LB step kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(thermal_LB_step.getSweep(defaultStream), "Thermal LB Step");

         }
         else if (timestep_strategy == "KernelOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the compute kernels without boundary conditions and communication. This makes only sense for benchmarking")
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(defaultStream), "Phase LB Step");
            timeloop.add() << Sweep(hydro_LB_step.getSweep(defaultStream), "Hydro LB Step");
            timeloop.add() << Sweep(thermal_LB_step.getSweep(defaultStream), "Thermal LB Step");
         }
         else
         {
            WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'")
         }

      }
      gpu::fieldCpy< GPUField, ScalarField_T >(blocks, phase_field_gpu, phase_field_ID);
      gpu::fieldCpy< GPUField, VectorField_T >(blocks, vel_field_gpu, velocity_field_ID);
      gpu::fieldCpy< GPUField, ScalarField_T >(blocks, temperature_field_gpu, temperature_field_ID);
#else
      ////////////////////////////
      // TIMELOOP THERMAL ONLY //
      //////////////////////////

      SweepTimeloop thermalTimeloop(blocks->getBlockStorage(), thermal_timesteps);
      if(heat_solver_RK_or_LBM)
      {
         if (orderRKSolver == 2)
         {
            thermalTimeloop.add() << Sweep(initRK2.getSweep(), "initialise RK");
            thermalTimeloop.add() << Sweep(RK2FirstStage.getSweep(), "Runge Kutta step one");
            thermalTimeloop.add() << Sweep(RK2SecondStage.getSweep(), "Runge Kutta step two")
                                  << AfterFunction(Comm_temperature_field, "Temperature field Communication");
         }
         else
         {
            thermalTimeloop.add() << Sweep(initRK4.getSweep(), "initialise RK");
            thermalTimeloop.add() << Sweep(RK4FirstStage.getSweep(), "Runge Kutta step one");
            thermalTimeloop.add() << Sweep(RK4SecondStage.getSweep(), "Runge Kutta step two");
            thermalTimeloop.add() << Sweep(RK4ThirdStage.getSweep(), "Runge Kutta step three");
            thermalTimeloop.add() << Sweep(RK4FourthStage.getSweep(), "Runge Kutta step four")
                                  << AfterFunction(Comm_temperature_field, "Temperature field Communication");
         }
      }
      else
      {
         thermalTimeloop.add() << BeforeFunction(Comm_thermal_distributions, "Thermal PDFs Communication")
                               << Sweep(thermal_LB_Neumann.getSweep(), "Neumann Thermal");
         thermalTimeloop.add() << Sweep(thermal_LB_DiffusionStatic.getSweep(), "Static Diffusion Thermal");
         thermalTimeloop.add() << Sweep(thermal_LB_DiffusionDynamic.getSweep(), "Dynamic Diffusion Thermal");
         thermalTimeloop.add() << Sweep(thermal_LB_step.getSweep(), "Thermal LB Step");
      }

      ////////////////
      // TIME LOOP //
      //////////////

      SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);
      if(heat_solver_RK_or_LBM)
      {
         if (timestep_strategy == "NormalTimestep")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Normal time step for production runs and full application benchmarks")
            timeloop.add() << BeforeFunction(Comm_velocity_based_distributions_start, "Start Hydro PDFs Communication")
                           << BeforeFunction(Comm_temperature_field_start, "Start Communication temperature field")
                           << Sweep(phase_field_LB_NoSlip.getSweep(), "NoSlip Phase");
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(), "Phase LB Step")
                           << AfterFunction(Comm_velocity_based_distributions_wait, "Wait Hydro PDFs Communication")
                           << AfterFunction(Comm_temperature_field_wait, "Wait Communication temperature field");

            timeloop.add() << BeforeFunction(Comm_phase_field_distributions_start, "Start Phase PDFs Communication")
                           << Sweep(hydro_LB_UBB.getSweep(), "UBB Hydro");
            timeloop.add() << Sweep(hydro_LB_NoSlip.getSweep(), "NoSlip Hydro");
            timeloop.add() << Sweep(hydro_LB_step.getSweep(), "Hydro LB Step");
            timeloop.add() << Sweep(swapPhaseField, "Swap PhaseField");
            timeloop.add() << Sweep(contact_angle.getSweep(), "Contact Angle")
                           << AfterFunction(Comm_phase_field_distributions_wait, "Wait Phase PDFs Communication")
                           << AfterFunction(Comm_phase_field, "Communication Phase Field");

            if(orderRKSolver == 2)
            {
               timeloop.add() << BeforeFunction(Comm_phase_field, "Communication Phase Field")
                              << Sweep(RK2FirstStage.getSweep(), "Runge Kutta step one");
               timeloop.add() << Sweep(RK2SecondStage.getSweep(), "Runge Kutta step two");
            }
            else
            {
               timeloop.add() << BeforeFunction(Comm_phase_field, "Communication Phase Field")
                              << Sweep(RK4FirstStage.getSweep(), "Runge Kutta step one");
               timeloop.add() << Sweep(RK4SecondStage.getSweep(), "Runge Kutta step two");
               timeloop.add() << Sweep(RK4ThirdStage.getSweep(), "Runge Kutta step three");
               timeloop.add() << Sweep(RK4FourthStage.getSweep(), "Runge Kutta step four");

            }
         }
         else if (timestep_strategy == "PhaseOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the Allen Cahn LB kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(), "Phase LB Step");
         }
         else if (timestep_strategy == "HydroOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the hydrodynamic LB kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(hydro_LB_step.getSweep(), "Hydro LB Step");
         }
         else if (timestep_strategy == "ThermalOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the RK kernels. This makes only sense for benchmarking")
            if(orderRKSolver == 2)
            {
               timeloop.add() << Sweep(RK2FirstStage.getSweep(), "Runge Kutta step one");
               timeloop.add() << Sweep(RK2SecondStage.getSweep(), "Runge Kutta step two");
            }
            else
            {
               timeloop.add() << Sweep(RK4FirstStage.getSweep(), "Runge Kutta step one");
               timeloop.add() << Sweep(RK4SecondStage.getSweep(), "Runge Kutta step two");
               timeloop.add() << Sweep(RK4ThirdStage.getSweep(), "Runge Kutta step three");
               timeloop.add() << Sweep(RK4FourthStage.getSweep(), "Runge Kutta step four");

            }
         }
         else if (timestep_strategy == "KernelOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the compute kernels without boundary conditions and communication. This makes only sense for benchmarking")
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(), "Phase LB Step");
            timeloop.add() << Sweep(hydro_LB_step.getSweep(), "Hydro LB Step");
            if(orderRKSolver == 2)
            {
               timeloop.add() << Sweep(RK2FirstStage.getSweep(), "Runge Kutta step one");
               timeloop.add() << Sweep(RK2SecondStage.getSweep(), "Runge Kutta step two");
            }
            else
            {
               timeloop.add() << Sweep(RK4FirstStage.getSweep(), "Runge Kutta step one");
               timeloop.add() << Sweep(RK4SecondStage.getSweep(), "Runge Kutta step two");
               timeloop.add() << Sweep(RK4ThirdStage.getSweep(), "Runge Kutta step three");
               timeloop.add() << Sweep(RK4FourthStage.getSweep(), "Runge Kutta step four");

            }
         }
         else
         {
            WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'")
         }
      }
      else
      {
         if (timestep_strategy == "NormalTimestep")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Normal time step for production runs and full application benchmarks with LBM thermal solver")
            timeloop.add() << BeforeFunction(Comm_thermal_distributions_start, "Start Thermal PDFs Communication")
                           << BeforeFunction(Comm_temperature_field_start, "Start Communication temperature field")
                           << Sweep(phase_field_LB_NoSlip.getSweep(), "NoSlip Phase");
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(), "Phase LB Step")
                           << AfterFunction(Comm_thermal_distributions_wait, "Wait Thermal PDFs Communication")
                           << AfterFunction(Comm_temperature_field_wait, "Wait Communication temperature field");

            timeloop.add() << BeforeFunction(Comm_phase_field_distributions_start, "Start Phase PDFs Communication")
                           << Sweep(hydro_LB_UBB.getSweep(), "UBB Hydro");
            timeloop.add() << Sweep(hydro_LB_NoSlip.getSweep(), "NoSlip Hydro");
            timeloop.add() << Sweep(hydro_LB_step.getSweep(), "Hydro LB Step");
            timeloop.add() << Sweep(swapPhaseField, "Swap PhaseField");
            timeloop.add() << Sweep(contact_angle.getSweep(), "Contact Angle")
                           << AfterFunction(Comm_phase_field_distributions_wait, "Wait Phase PDFs Communication");

            timeloop.add() << BeforeFunction(Comm_velocity_based_distributions_start, "Start Hydro PDFs Communication")
                           << BeforeFunction(Comm_phase_field_start, "Start Communication Phase Field")
                           << Sweep(thermal_LB_Neumann.getSweep(), "Neumann Thermal");
            timeloop.add() << Sweep(thermal_LB_DiffusionStatic.getSweep(), "Static Diffusion Thermal");
            timeloop.add() << Sweep(thermal_LB_DiffusionDynamic.getSweep(), "Dynamic Diffusion Thermal");
            timeloop.add() << Sweep(thermal_LB_step.getSweep(), "Thermal LB Step")
                           << AfterFunction(Comm_velocity_based_distributions_wait, "Wait Hydro PDFs Communication")
                           << AfterFunction(Comm_phase_field_wait, "Wait Communication Phase Field");
         }
         else if (timestep_strategy == "PhaseOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the Allen Cahn LB kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(), "Phase LB Step");
         }
         else if (timestep_strategy == "HydroOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the hydrodynamic LB kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(hydro_LB_step.getSweep(), "Hydro LB Step");
         }
         else if (timestep_strategy == "ThermalOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the thermal LB step kernel. This makes only sense for benchmarking")
            timeloop.add() << Sweep(thermal_LB_step.getSweep(), "Thermal LB Step");

         }
         else if (timestep_strategy == "KernelOnly")
         {
            WALBERLA_LOG_INFO_ON_ROOT("Running only the compute kernels without boundary conditions and communication. This makes only sense for benchmarking")
            timeloop.add() << Sweep(phase_field_LB_step.getSweep(), "Phase LB Step");
            timeloop.add() << Sweep(hydro_LB_step.getSweep(), "Hydro LB Step");
            timeloop.add() << Sweep(thermal_LB_step.getSweep(), "Thermal LB Step");
         }
         else
         {
            WALBERLA_ABORT_NO_DEBUG_INFO("Invalid value for 'timeStepStrategy'")
         }

      }
#endif

      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the PDFs")
      for (auto& block : *blocks)
      {
         init_h(&block);
         init_g(&block);
         if(!heat_solver_RK_or_LBM)
            init_f(&block);
      }
      Comm_phase_field_distributions();
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the PDFs done")

      // remaining time logger, higher frequency than 0.5 seconds is not allowed
      if (remaining_time_logger_frequency > 0.5)
      {
         timeloop.addFuncAfterTimeStep(
            timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remaining_time_logger_frequency),
            "remaining time logger");

         thermalTimeloop.addFuncAfterTimeStep(
            timing::RemainingTimeLogger(thermalTimeloop.getNrOfTimeSteps(), remaining_time_logger_frequency),
            "remaining time logger");
      }

      const uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      const std::string vtkPath = parameters.getParameter<std::string>("vtkPath", "thermocapillary_vtk");
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput   = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, vtkPath,
                                                           "simulation_step", false, true, true, false, 0);
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
         vtkOutput->addBeforeFunction([&]() {
            gpu::fieldCpy< ScalarField_T , GPUField >(blocks, phase_field_ID, phase_field_gpu);
            gpu::fieldCpy< VectorField_T , GPUField >(blocks, velocity_field_ID, vel_field_gpu);
            gpu::fieldCpy< ScalarField_T , GPUField >(blocks, temperature_field_ID, temperature_field_gpu);
         });
#endif

         auto phaseWriter = make_shared< field::VTKWriter< ScalarField_T > >(phase_field_ID, "PhaseField");
         vtkOutput->addCellDataWriter(phaseWriter);

         auto velWriter = make_shared< field::VTKWriter< VectorField_T > >(velocity_field_ID, "Velocity");
         vtkOutput->addCellDataWriter(velWriter);

         auto tempWriter = make_shared< field::VTKWriter< ScalarField_T > >(temperature_field_ID, "Temperature");
         vtkOutput->addCellDataWriter(tempWriter);

         timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      uint_t meshWriteFrequency = parameters.getParameter< uint_t >("meshWriteFrequency", 0);
      const int targetRank          = 0;
      int counter                   = 0;
      if (meshWriteFrequency > 0)
      {
         timeloop.addFuncAfterTimeStep(
            [&]() {
               if (timeloop.getCurrentTimeStep() % uint_t(meshWriteFrequency) == 0)
               {
                  auto mesh = postprocessing::realFieldToSurfaceMesh< ScalarField_T >(blocks, phase_field_ID, 0.5, 0, true,
                                                                                      targetRank, MPI_COMM_WORLD);
                  WALBERLA_EXCLUSIVE_WORLD_SECTION(targetRank)
                  {
                     const std::string path = "./Meshes/Droplet_contact_angle_" + std::to_string(int32_c(angle));
                     std::ostringstream out;
                     out << std::internal << std::setfill('0') << std::setw(6) << counter;
                     geometry::writeMesh(path + "_" + out.str() + ".obj", *mesh);
                     counter++;
                  }
               }
            },
            "Mesh writer");
      }


      const uint_t dbWriteFrequency = parameters.getParameter< uint_t >("dbWriteFrequency", 0);
      if (dbWriteFrequency > 0)
      {
         timeloop.addFuncAfterTimeStep(
            [&]() {
               if (timeloop.getCurrentTimeStep() % dbWriteFrequency == 0)
               {
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
                  gpu::fieldCpy< ScalarField_T , GPUField >(blocks, phase_field_ID, phase_field_gpu);
                  gpu::fieldCpy< VectorField_T , GPUField >(blocks, velocity_field_ID, vel_field_gpu);
                  gpu::fieldCpy< ScalarField_T , GPUField >(blocks, temperature_field_ID, temperature_field_gpu);
                  WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#endif
                  python_coupling::PythonCallback callback("at_end_of_time_step");
                  if (callback.isCallable())
                  {
                     callback.data().exposeValue("blocks", blocks);
                     callback.data().exposeValue("timeStep", timeloop.getCurrentTimeStep());
                     callback.data().exposeValue("stencil_phase", StencilNamePhase);
                     callback.data().exposeValue("stencil_hydro", StencilNameHydro);
                     callback.data().exposeValue("stencil_thermal", StencilNameThermal);
                     callback.data().exposeValue("collision_space_phase", CollisionSpacePhase);
                     callback.data().exposeValue("collision_space_hydro", CollisionSpaceHydro);
                     callback.data().exposeValue("collision_space_thermal", CollisionSpaceThermal);
                     callback.data().exposeValue("field_type", fieldType);
                     callback.data().exposeValue("field_type_pdfs", fieldTypePDFs);
                     callback();
                  }
               }
            },
            "Python callback");
      }

      const lbm::PerformanceEvaluation< FlagField_T > performance_evaluation( blocks, flag_field_hydrodynamic_LB_ID, fluidFlagUID );

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
#endif
      WALBERLA_LOG_INFO_ON_ROOT("Running "
                                << thermal_timesteps
                                << " timesteps with only the thermal LB step to initialise the temperature field")
      thermalTimeloop.run();
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
#endif
      WALBERLA_MPI_WORLD_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Initialised temperature field")

      WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")

      WcTimingPool timeloopTiming;
      WcTimer simTimer;

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
#endif
      simTimer.start();
      timeloop.run(timeloopTiming);

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
#endif
      simTimer.end();

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
      WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
#endif
      WALBERLA_MPI_WORLD_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")

      real_t time = real_c(simTimer.max());
      WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
      performance_evaluation.logResultOnRoot( timesteps, time );

      const auto reducedTimeloopTiming = timeloopTiming.getReduced();
      WALBERLA_LOG_RESULT_ON_ROOT( "Time loop timing:\n" << *reducedTimeloopTiming )

      if(benchmark)
      {
         WALBERLA_ROOT_SECTION()
         {
            python_coupling::PythonCallback callback("write_benchmark_results");
            if (callback.isCallable())
            {
               callback.data().exposeValue("stencil_phase", StencilNamePhase);
               callback.data().exposeValue("stencil_hydro", StencilNameHydro);
               callback.data().exposeValue("stencil_thermal", StencilNameThermal);
               callback.data().exposeValue("collision_space_phase", CollisionSpacePhase);
               callback.data().exposeValue("collision_space_hydro", CollisionSpaceHydro);
               callback.data().exposeValue("collision_space_thermal", CollisionSpaceThermal);
               callback.data().exposeValue("field_type", fieldType);
               callback.data().exposeValue("field_type_pdfs", fieldTypePDFs);
               callback.data().exposeValue("number_of_processes", totalNumProcesses);
               callback.data().exposeValue("threads", performance_evaluation.threads());
               callback.data().exposeValue("MLUPS", double_c(performance_evaluation.mlups(timesteps, time)));
               callback.data().exposeValue("MLUPS_process",
                                           double_c(performance_evaluation.mlupsPerProcess(timesteps, time)));
               callback.data().exposeValue(
                  "timeStepsPerSecond",
                  double_c(lbm::PerformanceEvaluation< FlagField_T >::timeStepsPerSecond(timesteps, time)));

               callback();
            }
         }
      }
   }
   return EXIT_SUCCESS;
}
