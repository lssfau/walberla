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
//! \file multiphase.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"
#include "core/timing/RemainingTimeLogger.h"

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

using namespace walberla;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

int main(int argc, char** argv)
{
   mpi::Environment Env(argc, argv);
   exportDataStructuresToPython();

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      auto config = *cfg;
      logging::configureLogging(config);
      shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGridFromConfig(config);

      ////////////////////////////
      // ADD DOMAIN PARAMETERS //
      //////////////////////////

      auto domainSetup             = config->getOneBlock("DomainSetup");
      Vector3< uint_t > domainSize = domainSetup.getParameter< Vector3< uint_t > >("domainSize");
      const bool tube              = domainSetup.getParameter< bool >("tube", false);

      ////////////////////////////////////////
      // ADD GENERAL SIMULATION PARAMETERS //
      //////////////////////////////////////

      auto parameters                    = config->getOneBlock("Parameters");
      const std::string timeStepStrategy = parameters.getParameter< std::string >("timeStepStrategy", "normal");
      const uint_t timesteps             = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const real_t remainingTimeLoggerFrequency =
         parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(3.0));
      const uint_t scenario = parameters.getParameter< uint_t >("scenario", uint_c(1));

      /////////////////////////
      // ADD DATA TO BLOCKS //
      ///////////////////////

      BlockDataID allen_cahn_PDFs_ID =
         field::addToStorage< PdfField_phase_T >(blocks, "LB phase field", real_c(0.0), field::fzyx);
      BlockDataID hydrodynamic_PDFs_ID =
         field::addToStorage< PdfField_hydro_T >(blocks, "LB velocity field", real_c(0.0), field::fzyx);
      BlockDataID velocity_field_ID =
         field::addToStorage< VelocityField_T >(blocks, "velocity", real_c(0.0), field::fzyx);
      BlockDataID density_field_ID = field::addToStorage< PhaseField_T >(blocks, "density", real_c(1.0), field::fzyx);
      BlockDataID phase_field_ID   = field::addToStorage< PhaseField_T >(blocks, "phase", real_c(0.0), field::fzyx);

      BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

      auto physical_parameters     = config->getOneBlock("PhysicalParameters");
      const real_t density_liquid  = physical_parameters.getParameter< real_t >("density_liquid", real_c(1.0));
      const real_t density_gas     = physical_parameters.getParameter< real_t >("density_gas");
      const real_t surface_tension = physical_parameters.getParameter< real_t >("surface_tension");
      const real_t mobility        = physical_parameters.getParameter< real_t >("mobility");
      const real_t gravitational_acceleration =
         physical_parameters.getParameter< real_t >("gravitational_acceleration");
      const real_t relaxation_time_liquid = physical_parameters.getParameter< real_t >("relaxation_time_liquid");
      const real_t relaxation_time_gas    = physical_parameters.getParameter< real_t >("relaxation_time_gas");
      const real_t interface_thickness    = physical_parameters.getParameter< real_t >("interface_thickness");

      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the phase-field")
      switch (scenario)
      {
      case 1: {
         auto bubbleParameters                  = config->getOneBlock("Bubble");
         const Vector3< real_t > bubbleMidPoint = bubbleParameters.getParameter< Vector3< real_t > >("bubbleMidPoint");
         const real_t bubbleRadius              = bubbleParameters.getParameter< real_t >("bubbleRadius", real_c(20.0));
         const bool bubble                      = bubbleParameters.getParameter< bool >("bubble", true);
         initPhaseField_sphere(blocks, phase_field_ID, velocity_field_ID, bubbleRadius, bubbleMidPoint, bubble);
         break;
      }
      case 2: {
         initPhaseField_RTI(blocks, phase_field_ID, interface_thickness, tube);
         break;
      }
      case 3: {
         auto dropParameters                     = config->getOneBlock("Drop");
         const Vector3< real_t > dropMidPoint    = dropParameters.getParameter< Vector3< real_t > >("drop_mid_point");
         const real_t dropRadius                 = dropParameters.getParameter< real_t >("drop_radius");
         const real_t poolDepth                  = dropParameters.getParameter< real_t >("pool_depth");
         const Vector3< real_t > impact_velocity = dropParameters.getParameter< Vector3< real_t > >("impact_velocity");
         init_hydrostatic_pressure(blocks, density_field_ID, gravitational_acceleration, poolDepth);

         initPhaseField_sphere(blocks, phase_field_ID, velocity_field_ID, dropRadius, dropMidPoint, false,
                               interface_thickness, impact_velocity);
         initPhaseField_pool(blocks, phase_field_ID, interface_thickness, poolDepth);
         break;
      }
      case 4: {
         auto TaylorBubbleParameters = config->getOneBlock("TaylorBubble");
         const real_t BubbleRadius   = TaylorBubbleParameters.getParameter< real_t >("bubble_radius");
         const real_t InitialHeight  = TaylorBubbleParameters.getParameter< real_t >("initial_height");
         const real_t Length         = TaylorBubbleParameters.getParameter< real_t >("length");
         init_Taylor_bubble_cylindric(blocks, phase_field_ID, BubbleRadius, InitialHeight, Length,
                                      real_c(interface_thickness));
         break;
      }
      default:
         WALBERLA_ABORT("Scenario is not defined.")
      }
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the phase-field done")

      /////////////////
      // ADD SWEEPS //
      ///////////////

      pystencils::initialize_phase_field_distributions init_h(allen_cahn_PDFs_ID, phase_field_ID, velocity_field_ID,
                                                              interface_thickness);
      pystencils::initialize_velocity_based_distributions init_g(density_field_ID, hydrodynamic_PDFs_ID,
                                                                 velocity_field_ID);

      pystencils::phase_field_LB_step phase_field_LB_step(allen_cahn_PDFs_ID, phase_field_ID, velocity_field_ID,
                                                          mobility, interface_thickness);
      pystencils::hydro_LB_step hydro_LB_step(
         hydrodynamic_PDFs_ID, phase_field_ID, velocity_field_ID, gravitational_acceleration, interface_thickness,
         density_liquid, density_gas, surface_tension, relaxation_time_liquid, relaxation_time_gas);

      ////////////////////////
      // ADD COMMUNICATION //
      //////////////////////
      auto UniformBufferedSchemeVelocityDistributions =
         make_shared< blockforest::communication::UniformBufferedScheme< Stencil_hydro_T > >(blocks);
      auto generatedPackInfo_velocity_based_distributions =
         make_shared< lbm::PackInfo_velocity_based_distributions >(hydrodynamic_PDFs_ID);
      UniformBufferedSchemeVelocityDistributions->addPackInfo(generatedPackInfo_velocity_based_distributions);
      auto Comm_velocity_based_distributions =
         std::function< void() >([&]() { UniformBufferedSchemeVelocityDistributions->communicate(); });
      auto Comm_velocity_based_distributions_start =
         std::function< void() >([&]() { UniformBufferedSchemeVelocityDistributions->startCommunication(); });
      auto Comm_velocity_based_distributions_wait =
         std::function< void() >([&]() { UniformBufferedSchemeVelocityDistributions->wait(); });

      auto UniformBufferedSchemePhaseField =
         make_shared< blockforest::communication::UniformBufferedScheme< Stencil_hydro_T > >(blocks);
      auto generatedPackInfo_phase_field = make_shared< pystencils::PackInfo_phase_field >(phase_field_ID);
      UniformBufferedSchemePhaseField->addPackInfo(generatedPackInfo_phase_field);
      auto Comm_phase_field = std::function< void() >([&]() { UniformBufferedSchemePhaseField->communicate(); });
      auto Comm_phase_field_start =
         std::function< void() >([&]() { UniformBufferedSchemePhaseField->startCommunication(); });
      auto Comm_phase_field_wait = std::function< void() >([&]() { UniformBufferedSchemePhaseField->wait(); });

      auto UniformBufferedSchemePhaseFieldDistributions =
         make_shared< blockforest::communication::UniformBufferedScheme< Stencil_hydro_T > >(blocks);
      auto generatedPackInfo_phase_field_distributions =
         make_shared< lbm::PackInfo_phase_field_distributions >(allen_cahn_PDFs_ID);
      UniformBufferedSchemePhaseFieldDistributions->addPackInfo(generatedPackInfo_phase_field_distributions);
      auto Comm_phase_field_distributions =
         std::function< void() >([&]() { UniformBufferedSchemePhaseFieldDistributions->communicate(); });
      auto Comm_phase_field_distributions_start =
         std::function< void() >([&]() { UniformBufferedSchemePhaseFieldDistributions->startCommunication(); });
      auto Comm_phase_field_distributions_wait =
         std::function< void() >([&]() { UniformBufferedSchemePhaseFieldDistributions->wait(); });

      ////////////////////////
      // BOUNDARY HANDLING //
      //////////////////////

      const FlagUID fluidFlagUID("Fluid");
      const FlagUID wallFlagUID("NoSlip");

      auto boundariesConfig = config->getBlock("Boundaries");
      if (boundariesConfig)
      {
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldID, boundariesConfig);
         if (tube)
         {
            // initialize cylindrical domain walls
            const Vector3< real_t > domainCylinderBottomEnd =
               Vector3< real_t >(real_c(domainSize[0]) * real_c(0.5), real_c(0), real_c(domainSize[2]) * real_c(0.5));
            const Vector3< real_t > domainCylinderTopEnd = Vector3< real_t >(
               real_c(domainSize[0]) * real_c(0.5), real_c(domainSize[1]), real_c(domainSize[2]) * real_c(0.5));
            const real_t radius = real_c(domainSize[0]) * real_c(0.5);
            const geometry::Cylinder cylinderTube(domainCylinderBottomEnd, domainCylinderTopEnd, radius);

            for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
            {
               FlagField_T* flagField = blockIt->template getData< FlagField_T >(flagFieldID);
               auto wallFlag          = flagField->getOrRegisterFlag(wallFlagUID);
               WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flagField, {
                  Cell globalCell = Cell(x, y, z);
                  blocks->transformBlockLocalToGlobalCell(globalCell, *blockIt);
                  const Vector3< real_t > globalPoint = blocks->getCellCenter(globalCell);

                  if (!geometry::contains(cylinderTube, globalPoint)) { addFlag(flagField->get(x, y, z), wallFlag); }
               }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
            }
         }
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID);
      }

      lbm::phase_field_LB_NoSlip phase_field_LB_NoSlip(blocks, allen_cahn_PDFs_ID);
      lbm::hydro_LB_NoSlip hydro_LB_NoSlip(blocks, hydrodynamic_PDFs_ID);
      pystencils::ContactAngle contact_angle(blocks, phase_field_ID, interface_thickness);

      phase_field_LB_NoSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldID, FlagUID("NoSlip"), fluidFlagUID);
      hydro_LB_NoSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldID, FlagUID("NoSlip"), fluidFlagUID);
      contact_angle.fillFromFlagField< FlagField_T >(blocks, flagFieldID, wallFlagUID, fluidFlagUID);

      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the normals-field done")

      ////////////////
      // TIME LOOP //
      //////////////

      SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

      timeloop.add() << BeforeFunction(Comm_velocity_based_distributions_start, "Start Hydro PDFs Communication")
                     << Sweep(phase_field_LB_NoSlip, "NoSlip Phase");
      timeloop.add() << Sweep(phase_field_LB_step, "Phase LB Step")
                     << AfterFunction(Comm_velocity_based_distributions_wait, "Wait Hydro PDFs Communication");
      timeloop.add() << Sweep(contact_angle, "Contact Angle")
                     << AfterFunction(Comm_phase_field, "Communication Phase Field");

      timeloop.add() << BeforeFunction(Comm_phase_field_distributions_start, "Start Phase PDFs Communication")
                     << Sweep(hydro_LB_NoSlip, "NoSlip Hydro");
      timeloop.add() << Sweep(hydro_LB_step, "Hydro LB Step")
                     << AfterFunction(Comm_phase_field_distributions_wait, "Wait Phase PDFs Communication");

      if (scenario == 4)
      {
         python_coupling::PythonCallback smear_interface("interface_diffusion");
         if (smear_interface.isCallable())
         {
            smear_interface.data().exposePtr("blocks", blocks);
            smear_interface();
         }
      }

      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the PDFs")
      for (auto& block : *blocks)
      {
         init_h(&block);
         init_g(&block);
      }
      Comm_phase_field_distributions();
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the PDFs done")
      uint_t dbWriteFrequency = parameters.getParameter< uint_t >("dbWriteFrequency", 10000000);

      timeloop.addFuncAfterTimeStep(
         [&]() {
            if (timeloop.getCurrentTimeStep() % dbWriteFrequency == 0)
            {
               python_coupling::PythonCallback callback("at_end_of_time_step");
               if (callback.isCallable())
               {
                  callback.data().exposeValue("blocks", blocks);
                  callback.data().exposeValue("timeStep", timeloop.getCurrentTimeStep());
                  callback.data().exposeValue("stencil_phase", StencilNamePhase);
                  callback.data().exposeValue("stencil_hydro", StencilNameHydro);
                  callback();
               }
            }
         },
         "Python callback");

      int meshWriteFrequency = parameters.getParameter< int >("meshWriteFrequency", 0);
      int counter            = 0;
      int targetRank         = 0;
      if (meshWriteFrequency > 0)
      {
         timeloop.addFuncAfterTimeStep(
            [&]() {
               if (timeloop.getCurrentTimeStep() % uint_t(meshWriteFrequency) == 0)
               {
                  auto mesh = postprocessing::realFieldToSurfaceMesh< PhaseField_T >(blocks, phase_field_ID, 0.5, 0,
                                                                                     true, targetRank, MPI_COMM_WORLD);
                  WALBERLA_EXCLUSIVE_WORLD_SECTION(targetRank)
                  {
                     std::string path = "";
                     std::ostringstream out;
                     out << std::internal << std::setfill('0') << std::setw(6) << counter;
                     geometry::writeMesh(path + out.str() + ".obj", *mesh);
                     counter++;
                  }
               }
            },
            "Mesh writer");
      }

      // remaining time logger
      timeloop.addFuncAfterTimeStep(
         timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
         "remaining time logger");

      uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput   = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                           "simulation_step", false, true, true, false, 0);
         auto phaseWriter = make_shared< field::VTKWriter< PhaseField_T > >(phase_field_ID, "PhaseField");
         vtkOutput->addCellDataWriter(phaseWriter);

         auto velWriter = make_shared< field::VTKWriter< VelocityField_T > >(velocity_field_ID, "Velocity");
         vtkOutput->addCellDataWriter(velWriter);

         timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      lbm::PerformanceEvaluation< FlagField_T > performance(blocks, flagFieldID, fluidFlagUID);
      WcTimingPool timeloopTiming;
      WcTimer simTimer;

      WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
      WALBERLA_MPI_WORLD_BARRIER()
      simTimer.start();
      timeloop.run(timeloopTiming);
      WALBERLA_MPI_WORLD_BARRIER()
      simTimer.end();

      auto time = real_c(simTimer.max());
      WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
      performance.logResultOnRoot(timesteps, time);

      const auto reducedTimeloopTiming = timeloopTiming.getReduced();
      WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)
   }
   return EXIT_SUCCESS;
}
