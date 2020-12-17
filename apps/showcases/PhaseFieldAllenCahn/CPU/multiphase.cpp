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
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"
#include "python_coupling/export/FieldExports.h"

#include "timeloop/SweepTimeloop.h"

#include "CalculateNormals.h"
#include "InitializerFunctions.h"
#include "PythonExports.h"
#include "contact.h"

//////////////////////////////
// INCLUDE GENERATED FILES //
////////////////////////////

#include "GenDefines.h"
#include "hydro_LB_NoSlip.h"
#include "hydro_LB_step.h"
#include "initialize_phase_field_distributions.h"
#include "initialize_velocity_based_distributions.h"
#include "phase_field_LB_NoSlip.h"
#include "phase_field_LB_step.h"
#include "stream_hydro.h"

////////////
// USING //
//////////

using namespace walberla;

using PdfField_phase_T = GhostLayerField< real_t, Stencil_phase_T::Size >;
using PdfField_hydro_T = GhostLayerField< real_t, Stencil_hydro_T::Size >;
using VelocityField_T  = GhostLayerField< real_t, Stencil_hydro_T::Dimension >;
using NormalsField_T   = GhostLayerField< int8_t, Stencil_hydro_T::Dimension >;
using PhaseField_T     = GhostLayerField< real_t, 1 >;
using flag_t           = walberla::uint8_t;
using FlagField_T      = FlagField< flag_t >;

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

      auto domainSetup                = config->getOneBlock("DomainSetup");
      Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");

      ////////////////////////////////////////
      // ADD GENERAL SIMULATION PARAMETERS //
      //////////////////////////////////////

      auto parameters                    = config->getOneBlock("Parameters");
      const std::string timeStepStrategy = parameters.getParameter< std::string >("timeStepStrategy", "normal");
      const uint_t timesteps             = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const real_t remainingTimeLoggerFrequency =
         parameters.getParameter< real_t >("remainingTimeLoggerFrequency", 3.0);
      const uint_t scenario  = parameters.getParameter< uint_t >("scenario", uint_c(1));
      const real_t alpha     = parameters.getParameter< real_t >("contactAngle", real_c(90));
      const real_t alpha_rad = alpha * (math::pi / 180);
      Vector3< int > overlappingWidth =
         parameters.getParameter< Vector3< int > >("overlappingWidth", Vector3< int >(1, 1, 1));

      /////////////////////////
      // ADD DATA TO BLOCKS //
      ///////////////////////

      BlockDataID lb_phase_field =
         field::addToStorage< PdfField_phase_T >(blocks, "lb phase field", real_t(0), field::fzyx);
      BlockDataID lb_velocity_field =
         field::addToStorage< PdfField_hydro_T >(blocks, "lb velocity field", real_t(0), field::fzyx);
      BlockDataID vel_field   = field::addToStorage< VelocityField_T >(blocks, "vel", real_t(0), field::fzyx);
      BlockDataID phase_field = field::addToStorage< PhaseField_T >(blocks, "phase", real_t(0), field::fzyx);

      BlockDataID normals     = field::addToStorage< NormalsField_T >(blocks, "normals", int8_t(0), field::fzyx);
      BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the phase-field")
      if (scenario == 1)
      {
         auto bubbleParameters                  = config->getOneBlock("Bubble");
         const Vector3< real_t > bubbleMidPoint = bubbleParameters.getParameter< Vector3< real_t > >("bubbleMidPoint");
         const real_t bubbleRadius              = bubbleParameters.getParameter< real_t >("bubbleRadius", 20.0);
         const bool bubble                      = bubbleParameters.getParameter< bool >("bubble", true);
         initPhaseField_sphere(blocks, phase_field, bubbleRadius, bubbleMidPoint, bubble);
      }
      else if (scenario == 2)
      {
         initPhaseField_RTI(blocks, phase_field);
      }
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the phase-field done")

      /////////////////
      // ADD SWEEPS //
      ///////////////

      auto physical_parameters     = config->getOneBlock("PhysicalParameters");
      const real_t density_liquid  = physical_parameters.getParameter< real_t >("density_liquid", real_c(1.0));
      const real_t density_gas     = physical_parameters.getParameter< real_t >("density_gas");
      const real_t surface_tension = physical_parameters.getParameter< real_t >("surface_tension");
      const real_t mobility        = physical_parameters.getParameter< real_t >("mobility");
      const real_t gravitational_acceleration =
         physical_parameters.getParameter< real_t >("gravitational_acceleration");
      const real_t relaxation_time_liquid = physical_parameters.getParameter< real_t >("relaxation_time_liquid");
      const real_t relaxation_time_gas    = physical_parameters.getParameter< real_t >("relaxation_time_gas");

      pystencils::initialize_phase_field_distributions init_h(lb_phase_field, phase_field, vel_field);
      pystencils::initialize_velocity_based_distributions init_g(lb_velocity_field, vel_field);

      pystencils::phase_field_LB_step phase_field_LB_step(
         lb_phase_field, phase_field, vel_field, mobility,
         Cell(overlappingWidth[0], overlappingWidth[1], overlappingWidth[2]));
      pystencils::hydro_LB_step hydro_LB_step(lb_velocity_field, phase_field, vel_field, gravitational_acceleration,
                                              density_liquid, density_gas, surface_tension, relaxation_time_liquid,
                                              relaxation_time_gas,
                                              Cell(overlappingWidth[0], overlappingWidth[1], overlappingWidth[2]));
      pystencils::stream_hydro stream_hydro(lb_velocity_field,
                                            Cell(overlappingWidth[0], overlappingWidth[1], overlappingWidth[2]));

      ////////////////////////
      // ADD COMMUNICATION //
      //////////////////////

      blockforest::communication::UniformBufferedScheme< Stencil_hydro_T > Comm_phase_field(blocks);
      Comm_phase_field.addPackInfo(make_shared< field::communication::PackInfo< PhaseField_T > >(phase_field));

      blockforest::communication::UniformBufferedScheme< Stencil_hydro_T > Comm_velocity_based_distributions(blocks);
      Comm_velocity_based_distributions.addPackInfo(
         make_shared< field::communication::PackInfo< PdfField_hydro_T > >(lb_velocity_field));

      blockforest::communication::UniformBufferedScheme< Stencil_hydro_T > Comm_phase_field_distributions(blocks);
      Comm_phase_field_distributions.addPackInfo(
         make_shared< field::communication::PackInfo< PdfField_phase_T > >(lb_phase_field));

      ////////////////////////
      // BOUNDARY HANDLING //
      //////////////////////

      const FlagUID fluidFlagUID("Fluid");
      const FlagUID wallFlagUID("NoSlip");

      auto boundariesConfig = config->getBlock("Boundaries");
      if (boundariesConfig)
      {
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldID, boundariesConfig);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID);
      }

      calculate_normals(blocks, normals, flagFieldID, fluidFlagUID, wallFlagUID);
      lbm::phase_field_LB_NoSlip phase_field_LB_NoSlip(blocks, lb_phase_field);
      lbm::hydro_LB_NoSlip hydro_LB_NoSlip(blocks, lb_velocity_field);
      lbm::contact contact_angle(blocks, phase_field, alpha_rad);

      phase_field_LB_NoSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldID, FlagUID("NoSlip"), fluidFlagUID);
      hydro_LB_NoSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldID, FlagUID("NoSlip"), fluidFlagUID);
      contact_angle.fillFromNormalField< NormalsField_T >(blocks, normals);

      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the normals-field done")

      ////////////////
      // TIME LOOP //
      //////////////

      auto timeLoop       = make_shared< SweepTimeloop >(blocks->getBlockStorage(), timesteps);
      auto normalTimeStep = [&]() {
         for (auto& block : *blocks)
         {
            phase_field_LB_NoSlip(&block);
            Comm_phase_field_distributions();

            phase_field_LB_step(&block);

            Comm_phase_field();
            contact_angle(&block);

            hydro_LB_step(&block);
            hydro_LB_NoSlip(&block);
            Comm_velocity_based_distributions();

            stream_hydro(&block);
         }
      };
      auto simpleOverlapTimeStep = [&]() {
         for (auto& block : *blocks)
            phase_field_LB_NoSlip(&block);

         Comm_phase_field_distributions.startCommunication();
         for (auto& block : *blocks)
            phase_field_LB_step.inner(&block);
         Comm_phase_field_distributions.wait();
         for (auto& block : *blocks)
            phase_field_LB_step.outer(&block);

         Comm_phase_field.startCommunication();
         for (auto& block : *blocks)
            hydro_LB_step.inner(&block);
         Comm_phase_field.wait();
         for (auto& block : *blocks)
            hydro_LB_step.outer(&block);

         for (auto& block : *blocks)
            hydro_LB_NoSlip(&block);

         Comm_velocity_based_distributions.startCommunication();
         for (auto& block : *blocks)
            stream_hydro.inner(&block);
         Comm_velocity_based_distributions.wait();
         for (auto& block : *blocks)
            stream_hydro.outer(&block);
      };
      std::function< void() > timeStep;
      if (timeStepStrategy == "overlap")
      {
         timeStep = std::function< void() >(simpleOverlapTimeStep);
         WALBERLA_LOG_INFO_ON_ROOT("overlapping timestep")
      }
      else
      {
         timeStep = std::function< void() >(normalTimeStep);
         WALBERLA_LOG_INFO_ON_ROOT("normal timestep with no overlapping")
      }

      timeLoop->add() << BeforeFunction(timeStep) << Sweep([](IBlock*) {}, "time step");

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
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation of the PDFs done")
      uint_t dbWriteFrequency = parameters.getParameter< uint_t >("dbWriteFrequency", 10000000);

      timeLoop->addFuncAfterTimeStep(
         [&]() {
            if (timeLoop->getCurrentTimeStep() % dbWriteFrequency == 0)
            {
               python_coupling::PythonCallback callback("at_end_of_time_step");
               if (callback.isCallable())
               {
                  callback.data().exposeValue("blocks", blocks);
                  callback.data().exposeValue( "timeStep", timeLoop->getCurrentTimeStep());
                  callback();
               }
            }
         },
         "Python callback");

      // remaining time logger
      timeLoop->addFuncAfterTimeStep(
         timing::RemainingTimeLogger(timeLoop->getNrOfTimeSteps(), remainingTimeLoggerFrequency),
         "remaining time logger");

      uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput   = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
         auto phaseWriter = make_shared< field::VTKWriter< PhaseField_T > >(phase_field, "PhaseField");
         vtkOutput->addCellDataWriter(phaseWriter);

         // auto normlasWriter = make_shared<field::VTKWriter<NormalsField_T>>(normals, "Normals");
         // vtkOutput->addCellDataWriter(normlasWriter);

         // auto velWriter = make_shared<field::VTKWriter<VelocityField_T>>(vel_field, "Velocity");
         // vtkOutput->addCellDataWriter(velWriter);

         timeLoop->addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")
      WcTimer simTimer;
      simTimer.start();
      timeLoop->run();
      simTimer.end();
      WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
      auto time            = simTimer.last();
      auto nrOfCells       = real_c(cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2]);
      auto mlupsPerProcess = nrOfCells * real_c(timesteps) / time * 1e-6;
      WALBERLA_LOG_RESULT_ON_ROOT("MLUPS per process: " << mlupsPerProcess)
      WALBERLA_LOG_RESULT_ON_ROOT("Time per time step: " << time / real_c(timesteps) << " s")
   }
   return EXIT_SUCCESS;
}
