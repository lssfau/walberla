//
// Created by dy94rovu on 6/24/24.
//
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
//! \file MaterialTransport.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Ravi Ayyala Somayajula <ravi.k.ayyala@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/grid_generator/SCIterator.h"
#include "core/logging/all.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/vtk/all.h"

#include "geometry/InitBoundaryHandling.h"

#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/communication/UniformGPUScheme.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/vtk/all.h"

#include "lbm_mesapd_coupling/DataTypesCodegen.h"
#include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMSweepCollection.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"

#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "sqlite/SQLite.h"

#include "vtk/all.h"

#include "ConcentrationMacroGetter.h"
#include "FluidMacroGetter.h"
#include "GeneralInfoHeader.h"
#include "../../utilities/InitializerFunctions.h"
#include "PackInfoConcentration.h"
#include "PackInfoFluid.h"

namespace MaterialTransport{
///////////
// USING //
///////////

using namespace walberla;
using namespace lbm_mesapd_coupling::psm::gpu;
typedef pystencils::PackInfoFluid PackInfoFluid_T;
typedef pystencils::PackInfoConcentration PackInfoConcentration_T;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

///////////
// FLAGS //
///////////

// Fluid Flags
const FlagUID Fluid_Flag("Fluid");
const FlagUID Density_Fluid_Flag("Density_Fluid");
const FlagUID NoSlip_Fluid_Flag("NoSlip_Fluid");
const FlagUID FreeSlip_Fluid_Flag("FreeSlip_Fluid");
const FlagUID Inflow_Fluid_Flag("Inflow_Fluid");
const FlagUID Inflow_Fluid_Flag_Poisuelle("Inflow_Fluid_Poisuelle");
const FlagUID Neumann_Fluid_Flag("Neumann_Fluid");

// Concentration Flags
const FlagUID Concentration_Flag("Concentration");
const FlagUID Density_Concentration_Flag_west("Density_Concentration_west");
const FlagUID Density_Concentration_Flag_south("Density_Concentration_south");
const FlagUID Density_Concentration_Flag_north("Density_Concentration_north");
const FlagUID NoSlip_Concentration_Flag("NoSlip_Concentration");
const FlagUID Inflow_Concentration_Flag("Inflow_Concentration");
//const FlagUID Dirichlet_Concentration_Flag("Dirichlet_Concentration");
const FlagUID Neumann_Concentration_Flag("Neumann_Concentration");


//////////
// MAIN //
//////////


int main(int argc, char** argv) {
   Environment env(argc, argv);
   auto cfgFile = env.config();
   if (!cfgFile) { WALBERLA_ABORT("Usage: " << argv[0] << " path-to-configuration-file \n"); }

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   gpu::selectDeviceBasedOnMpiRank();
#endif

   WALBERLA_LOG_INFO_ON_ROOT("waLBerla revision: " << std::string(WALBERLA_GIT_SHA1).substr(0, 8));
   WALBERLA_LOG_INFO_ON_ROOT("compiler flags: " << std::string(WALBERLA_COMPILER_FLAGS));
   WALBERLA_LOG_INFO_ON_ROOT("build machine: " << std::string(WALBERLA_BUILD_MACHINE));
   WALBERLA_LOG_INFO_ON_ROOT(*cfgFile);

   // Read config file
   Config::BlockHandle numericalSetup = cfgFile->getBlock("NumericalSetup");
   const std::string simulationName   = numericalSetup.getParameter< std::string >("simName");
   const uint_t numXBlocks            = numericalSetup.getParameter< uint_t >("numXBlocks");
   const uint_t numYBlocks            = numericalSetup.getParameter< uint_t >("numYBlocks");
   const uint_t numZBlocks            = numericalSetup.getParameter< uint_t >("numZBlocks");
   WALBERLA_CHECK_EQUAL(numXBlocks * numYBlocks * numZBlocks, uint_t(MPIManager::instance()->numProcesses()),
                        "When using GPUs, the number of blocks ("
                           << numXBlocks * numYBlocks * numZBlocks << ") has to match the number of MPI processes ("
                           << uint_t(MPIManager::instance()->numProcesses()) << ")");
   const bool periodicInX                 = numericalSetup.getParameter< bool >("periodicInX");
   const bool periodicInY                 = numericalSetup.getParameter< bool >("periodicInY");
   const bool periodicInZ                 = numericalSetup.getParameter< bool >("periodicInZ");
   const uint_t numXCellsPerBlock         = numericalSetup.getParameter< uint_t >("numXCellsPerBlock");
   const uint_t numYCellsPerBlock         = numericalSetup.getParameter< uint_t >("numYCellsPerBlock");
   const uint_t numZCellsPerBlock         = numericalSetup.getParameter< uint_t >("numZCellsPerBlock");
   const bool sendDirectlyFromGPU         = numericalSetup.getParameter< bool >("sendDirectlyFromGPU");
   const bool useCommunicationHiding      = numericalSetup.getParameter< bool >("useCommunicationHiding");
   const Vector3< uint_t > frameWidth     = numericalSetup.getParameter< Vector3< uint_t > >("frameWidth");
   const uint_t timeSteps                 = numericalSetup.getParameter< uint_t >("timeSteps");

   const real_t uInflow_LBM        = numericalSetup.getParameter< real_t >("uInflowL");
   const real_t uInflow_SI        = numericalSetup.getParameter< real_t >("uInflowSI");
   const real_t relaxationRateFluid = numericalSetup.getParameter< real_t >("relaxationRateFLuid");
   const real_t relaxationRateConcentration = numericalSetup.getParameter< real_t >("relaxationRateConcentration");
   const real_t dx_SI = numericalSetup.getParameter< real_t >("dx");

   const Vector3< real_t > Uinitialize (uInflow_LBM,0,0);
   if ((periodicInY && numYBlocks == 1) || (periodicInZ && numZBlocks == 1))
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Using only 1 block in periodic dimensions can lead to unexpected behavior.")
   }
   const real_t viscosity = lbm::collision_model::viscosityFromOmega(relaxationRateFluid);
   const real_t diffusivity = lbm::collision_model::viscosityFromOmega(relaxationRateConcentration);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(viscosity)

   Config::BlockHandle outputSetup      = cfgFile->getBlock("Output");
   const uint_t vtkSpacing              = outputSetup.getParameter< uint_t >("vtkSpacing");
   const std::string vtkFolder          = outputSetup.getParameter< std::string >("vtkFolder");
   const uint_t performanceLogFrequency = outputSetup.getParameter< uint_t >("performanceLogFrequency");


   // Calculation and Conversions
   Vector3< uint_t > domainSize;
   domainSize[0] = numXBlocks * numXCellsPerBlock;
   domainSize[1] = numYBlocks * numYCellsPerBlock;
   domainSize[2] = numZBlocks * numZCellsPerBlock;

   const real_t dt_SI                       = uInflow_LBM / uInflow_SI * dx_SI;
   const real_t dx                          = real_t(1);  // LBM
   const real_t dt                          = real_t(1);  // LBM
   const real_t tinitial                    = real_t(0);
   real_t sigma_0                           = real_t(10*dx);
   real_t sigma_D                           = real_t(std::sqrt(2*diffusivity*tinitial));
   const Vector3< real_t > x_0(real_t(200),real_t(200),real_t(1));

   WALBERLA_LOG_INFO_ON_ROOT("Simulation setup:");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - dt_SI = " << dt_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - dx_SI = " << dx_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: kinematic. visc (LBM) = " << viscosity << ", relaxation rate (omega) = " << relaxationRateFluid);
   WALBERLA_LOG_INFO_ON_ROOT(" - concentration: diffusivity (LBM) = " << diffusivity << ", relaxation rate (omega_c)= " << relaxationRateConcentration);
   WALBERLA_LOG_INFO_ON_ROOT(" - Fluid Lattice Velocity = " << uInflow_LBM);
   WALBERLA_LOG_INFO_ON_ROOT(" - Relaxation time (Tau) Fluid = " << real_c(1/real_t(relaxationRateFluid)));
   WALBERLA_LOG_INFO_ON_ROOT(" - Relaxation time (Tau) Concentration = " << real_c(1/real_t(relaxationRateConcentration)));
   WALBERLA_LOG_INFO_ON_ROOT("Peclet number is " << (uInflow_LBM*domainSize[1])/diffusivity);


   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid(
      numXBlocks, numYBlocks, numZBlocks, numXCellsPerBlock, numYCellsPerBlock, numZCellsPerBlock, real_t(1), uint_t(0),
      false, false, periodicInX, periodicInY, periodicInZ, // periodicity
      false);

   auto simulationDomain = blocks->getDomain();

   // Setting initial PDFs to nan helps to detect bugs in the initialization/BC handling
   // Depending on WALBERLA_BUILD_WITH_GPU_SUPPORT, pdfFieldCPUGPUID is either a CPU or a CPU field
   BlockDataID velFieldFluidID;
   BlockDataID densityConcentrationFieldID;

   #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      // Fluid PDFs on GPU
      BlockDataID pdfFieldFluidID =
         field::addToStorage< PdfField_fluid_T >(blocks, "pdf fluid field (fzyx)", real_c(std::nan("")), field::fzyx);
      BlockDataID pdfFieldFluidCPUGPUID = gpu::addGPUFieldToStorage< PdfField_fluid_T >(blocks, pdfFieldFluidID, "pdf fluid field GPU");

      // Concentration PDFs on GPU
      BlockDataID pdfFieldConcentrationID =
         field::addToStorage< PdfField_concentration_T >(blocks, "pdf concentration field (fzyx)", real_c(std::nan("")), field::fzyx);
      BlockDataID pdfFieldConcentrationCPUGPUID = gpu::addGPUFieldToStorage< PdfField_concentration_T >(blocks, pdfFieldConcentrationID, "pdf concentration field GPU");

      // Fluid velocity field on GPU
      velFieldFluidID  = field::addToStorage< VelocityField_fluid_T >(blocks, "velocity fluid field", real_t(0), field::fzyx);
      BlockDataID velFieldFluidCPUGPUID = gpu::addGPUFieldToStorage< VelocityField_fluid_T >(blocks, velFieldFluidID, "velocity fluid field GPU");

      // Concentration Density on GPU
      densityConcentrationFieldID = field::addToStorage< DensityField_concentration_T >(blocks, "density concentration field", real_t(0), field::fzyx);
      BlockDataID densityConcentrationFieldCPUGPUID = gpu::addGPUFieldToStorage< DensityField_concentration_T >(blocks, densityConcentrationFieldID, "density concentration field GPU");

   #else

      // Fluid PDFs on CPU
      BlockDataID pdfFieldFluidCPUGPUID =
         field::addToStorage< PdfField_fluid_T >(blocks, "pdf fluid field CPU", real_c(std::nan("")), field::fzyx);

      BlockDataID velFieldFluidCPUGPUID = field::addToStorage< VelocityField_fluid_T >(blocks, "velocity fluid field CPU", real_t(0), field::fzyx);
      // Concentration PDFs on CPU
      BlockDataID pdfFieldConcentrationCPUGPUID =
         field::addToStorage< PdfField_concentration_T >(blocks, "pdf concentration field CPU", real_c(std::nan("")), field::fzyx);

      BlockDataID densityConcentrationFieldCPUGPUID = field::addToStorage< DensityField_concentration_T >(blocks, "density concentration field", real_t(0), field::fzyx);

   #endif
      BlockDataID densityFluidFieldID = field::addToStorage< DensityField_fluid_T >(blocks, "density fluid field", real_t(0), field::fzyx);
      //BlockDataID densityConcentrationFieldID = field::addToStorage< DensityField_concentration_T >(blocks, "density concentration field", real_t(0), field::fzyx);
      BlockDataID flagFieldFluidID = field::addFlagFieldToStorage< FlagField_T >(blocks, "fluid flag field");
      BlockDataID flagFieldConcentrationID = field::addFlagFieldToStorage< FlagField_T >(blocks, "concentration flag field");


      // Assemble boundary block string
      std::string boundariesBlockString = " BoundariesFluid";

      if (!periodicInX && simulationName == "2d_poisuelle" )
      {
         boundariesBlockString += "{"
                                  "Border { direction W;    walldistance -1;  flag Inflow_Fluid_Poisuelle; }"
                                  "Border { direction E;    walldistance -1;  flag Neumann_Fluid; }";
      }

      if(!periodicInX && simulationName == "2d_plate"){
         boundariesBlockString += "{"
                                  "Border { direction W;    walldistance -1;  flag Inflow_Fluid; }"
                                  "Border { direction E;    walldistance -1;  flag Density_Fluid; }";
      }

      if (!periodicInY && simulationName == "2d_poisuelle")
      {
         boundariesBlockString += "Border { direction S;    walldistance -1;  flag Neumann_Fluid; }"
                                  "Border { direction N;    walldistance -1;  flag NoSlip_Fluid; }";
      }

      if (!periodicInY && simulationName == "2d_plate")
      {
         boundariesBlockString += "Border { direction S;    walldistance -1;  flag Neumann_Fluid; }"
                                  "Border { direction N;    walldistance -1;  flag Neumann_Fluid; }";
      }

      if (!periodicInZ)
      {
         boundariesBlockString += "Border { direction T;    walldistance -1;  flag NoSlip_Fluid; }"
                                  "Border { direction B;    walldistance -1;  flag NoSlip_Fluid; }";
      }

      boundariesBlockString += "}";

      boundariesBlockString += "\n BoundariesConcentration";

      if(!periodicInX)
      {
         boundariesBlockString += "{"
                                  "Border { direction W;    walldistance -1;  flag Density_Concentration_west; }"
                                  "Border { direction E;    walldistance -1;  flag Neumann_Concentration; }";
      }
      if (!periodicInY && simulationName == "2d_plate")
      {
         boundariesBlockString += "Border { direction S;    walldistance -1;  flag Density_Concentration_south; }"
                                  "Border { direction N;    walldistance -1;  flag Neumann_Concentration; }";
      }
      if (!periodicInY && simulationName == "2d_poisuelle")
      {
         boundariesBlockString += "Border { direction S;    walldistance -1;  flag Neumann_Concentration; }"
                                  "Border { direction N;    walldistance -1;  flag Density_Concentration_north; }";
      }

      if (!periodicInZ)
      {
         boundariesBlockString += "Border { direction T;    walldistance -1;  flag Neumann_Concentration; }"
                                  "Border { direction B;    walldistance -1;  flag Neumann_Concentration; }";
      }
      boundariesBlockString += "}";
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream boundariesFile("boundaries.prm");
         boundariesFile << boundariesBlockString;
         boundariesFile.close();
      }
      WALBERLA_MPI_BARRIER()

      auto boundariesCfgFile = Config();
      boundariesCfgFile.readParameterFile("boundaries.prm");
      auto boundariesConfigFluid = boundariesCfgFile.getBlock("BoundariesFluid");
      auto boundariesConfigConcentration = boundariesCfgFile.getBlock("BoundariesConcentration");


      // map boundaries into the fluid field simulation
      lbm::BC_Fluid_Density density_fluid_bc(blocks, pdfFieldFluidCPUGPUID, real_t(1.0));
      lbm::BC_Fluid_NoSlip noSlip_fluid_bc(blocks, pdfFieldFluidCPUGPUID);
      lbm::BC_Fluid_FreeSlip freeSlip_fluid_bc(blocks, pdfFieldFluidCPUGPUID);
      lbm::BC_Fluid_UBB ubb_fluid_bc(blocks, pdfFieldFluidCPUGPUID, uInflow_LBM, real_t(0));

      // velocity call back for poisuelle flow:

      auto VelocityCallback = [](const Cell& pos, const shared_ptr< StructuredBlockForest >& SbF, IBlock& block) {


            Cell globalCell;
            const uint_t FieldGhostLayers = 1;
            cell_idx_t ghost = cell_idx_c(FieldGhostLayers);
            CellInterval domainBB = SbF->getDomainCellBB();
            domainBB.xMin() -= ghost;
            domainBB.xMax() += ghost;

            // WEST - Inflow
            CellInterval west(domainBB.xMin(), domainBB.yMin(),
                              domainBB.zMin(), domainBB.xMin(),
                              domainBB.yMax(), domainBB.zMin());

            //const auto cellAABB = SbF->getBlockLocalCellAABB(block, pos);
            //auto cellCenter = cellAABB.center();

            auto h_y          = real_c(west.ySize());
            SbF->transformBlockLocalToGlobalCell(globalCell, block, pos);
            //WALBERLA_LOG_INFO_ON_ROOT("h_y is " << h_y << " " << "y coord is " << globalCell[1]);
            auto y1 = real_c(globalCell[1]);

            real_t u = 0.05*(1- std::pow((y1/h_y),2));

            Vector3< real_t > result(u, 0.0, 0.0);
            return result;

      };

         std::function< Vector3< real_t >(const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
         velocity_initialisation = std::bind(VelocityCallback,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);


      lbm::BC_PoisuelleUBB ubb_fluid_poisuelle (blocks,pdfFieldFluidCPUGPUID,velocity_initialisation);
      lbm::BC_Fluid_Neumann neumann_fluid_bc(blocks, pdfFieldFluidCPUGPUID);

      lbm::BC_Concentration_Density density_concentration_bc_west(blocks, pdfFieldConcentrationCPUGPUID, real_t(0));
      lbm::BC_Concentration_Density density_concentration_bc_south(blocks, pdfFieldConcentrationCPUGPUID, real_t(1));
      lbm::BC_Concentration_Density density_concentration_bc_north(blocks, pdfFieldConcentrationCPUGPUID, real_t(1));
      lbm::BC_Concentration_Neumann neumann_concentration_bc(blocks, pdfFieldConcentrationCPUGPUID);

      if (!periodicInX && !periodicInY && !periodicInZ &&
          simulationName == "2d_plate") // executes for the 2d plate validation
      {
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldFluidID, boundariesConfigFluid);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldFluidID, Fluid_Flag);
         density_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Density_Fluid_Flag, Fluid_Flag);
         noSlip_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, NoSlip_Fluid_Flag, Fluid_Flag);
         ubb_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Inflow_Fluid_Flag, Fluid_Flag);
         neumann_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Neumann_Fluid_Flag, Fluid_Flag);
         // map boundaries into the concentration field simulation
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldConcentrationID,
                                                       boundariesConfigConcentration);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldConcentrationID, Concentration_Flag);

         density_concentration_bc_south.fillFromFlagField< FlagField_T >(
            blocks, flagFieldConcentrationID, Density_Concentration_Flag_south, Concentration_Flag);

         density_concentration_bc_west.fillFromFlagField< FlagField_T >(
            blocks, flagFieldConcentrationID, Density_Concentration_Flag_west, Concentration_Flag);

         neumann_concentration_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldConcentrationID,
                                                                   Neumann_Concentration_Flag, Concentration_Flag);
      }
      if (periodicInX && periodicInY && periodicInZ && simulationName == "2d_gaussian")
      { // executes for the 2d gaussian curve validation

         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldFluidID, Fluid_Flag);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldConcentrationID, Concentration_Flag);
      }

      if (!periodicInX && !periodicInY && !periodicInZ &&
          simulationName == "2d_poisuelle") // executes for the 2d plate validation
      {
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldFluidID, boundariesConfigFluid);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldFluidID, Fluid_Flag);
         density_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Density_Fluid_Flag, Fluid_Flag);
         noSlip_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, NoSlip_Fluid_Flag, Fluid_Flag);
         ubb_fluid_poisuelle.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Inflow_Fluid_Flag_Poisuelle, Fluid_Flag);
         neumann_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Neumann_Fluid_Flag, Fluid_Flag);
         // map boundaries into the concentration field simulation
         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldConcentrationID,
                                                       boundariesConfigConcentration);
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldConcentrationID, Concentration_Flag);

         density_concentration_bc_north.fillFromFlagField< FlagField_T >(
            blocks, flagFieldConcentrationID, Density_Concentration_Flag_north, Concentration_Flag);

         density_concentration_bc_west.fillFromFlagField< FlagField_T >(
            blocks, flagFieldConcentrationID, Density_Concentration_Flag_west, Concentration_Flag);

         neumann_concentration_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldConcentrationID,
                                                                   Neumann_Concentration_Flag, Concentration_Flag);
      }


      ///////////////
      // TIME LOOP //
      ///////////////
      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      initConcentrationField(blocks,densityConcentrationFieldID,simulationDomain,domainSize);
      initFluidField(blocks, velFieldFluidID , Uinitialize);

      #else

      if (simulationName == "2d_gaussian")
      {
         initConcentrationFieldGaussian(blocks, densityConcentrationFieldCPUGPUID, simulationDomain, domainSize,
                                        sigma_0, sigma_D, Uinitialize, x_0);
      }
      if(simulationName == "2d_poisuelle"){
         WALBERLA_LOG_INFO_ON_ROOT("Initializing Poisuelle Fluid Domain");
         initFluidFieldPoiseuille(blocks,velFieldFluidCPUGPUID,Uinitialize,domainSize);
      }
      if (simulationName == "2d_plate")
      {
         WALBERLA_LOG_INFO_ON_ROOT("Initializing 2d plate Fluid Domain");
         initFluidField(blocks, velFieldFluidCPUGPUID, Uinitialize);
      }
#endif
      pystencils::InitializeFluidDomain initializeFluidDomain(pdfFieldFluidCPUGPUID, velFieldFluidCPUGPUID, real_t(1));
      pystencils::InitializeConcentrationDomain initializeConcentrationDomain(
         densityConcentrationFieldCPUGPUID, pdfFieldConcentrationCPUGPUID, velFieldFluidCPUGPUID);

      for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
      {
         initializeFluidDomain(&(*blockIt));

         initializeConcentrationDomain(&(*blockIt));
      }

      ///////////////////////
      // ADD COMMUNICATION //
      //////////////////////

      // Setup of the fluid LBM communication for synchronizing the fluid pdf field between neighboring blocks
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
gpu::communication::UniformGPUScheme< Stencil_Fluid_T > com_fluid(blocks, sendDirectlyFromGPU, false);
#else
walberla::blockforest::communication::UniformBufferedScheme< Stencil_Fluid_T > com_fluid(blocks);
#endif
com_fluid.addPackInfo(make_shared< PackInfoFluid_T >(pdfFieldFluidCPUGPUID));
auto communication_fluid = std::function< void() >([&]() { com_fluid.communicate(); });

// Setup of the concentration LBM communication for synchronizing the concentration pdf field between neighboring blocks
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   gpu::communication::UniformGPUScheme< Stencil_Concentration_T > com_concentration(blocks, sendDirectlyFromGPU, false);
#else
   walberla::blockforest::communication::UniformBufferedScheme< Stencil_Concentration_T > com_concentration(blocks);
#endif
   com_concentration.addPackInfo(make_shared< PackInfoConcentration_T >(pdfFieldConcentrationCPUGPUID));
   auto communication_concentration = std::function< void() >([&]() { com_concentration.communicate(); });

   // time loop objects for communication with and without hiding

   SweepTimeloop commTimeloop(blocks->getBlockStorage(), timeSteps);
   SweepTimeloop timeloop(blocks->getBlockStorage(), timeSteps);

   timeloop.addFuncBeforeTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");

   // objects to get the macroscopic quantities

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   pystencils:: FluidMacroGetter getterSweep_fluid(densityFluidFieldID,pdfFieldFluidID,velFieldFluidID,real_t(0),real_t(0),real_t(0));
#else
   pystencils::FluidMacroGetter getterSweep_fluid(densityFluidFieldID,pdfFieldFluidCPUGPUID,velFieldFluidCPUGPUID);
#endif

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   pystencils:: ConcentrationMacroGetter getterSweep_concentration(densityConcentrationFieldID,pdfFieldConcentrationID,velFieldFluidID,real_t(0),real_t(0),real_t(0));
#else
   pystencils::ConcentrationMacroGetter getterSweep_concentration(densityConcentrationFieldCPUGPUID,pdfFieldConcentrationCPUGPUID);
#endif


   // VTK output -> comeback later to this part
   if (vtkSpacing != uint_t(0))
   {

      // Fields
      auto vtkOutput_Fluid = vtk::createVTKOutput_BlockData(blocks, "vtk files fluid", vtkSpacing, 0, false, vtkFolder);

      vtkOutput_Fluid->addBeforeFunction(communication_fluid);


      vtkOutput_Fluid->addBeforeFunction([&]() {
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         gpu::fieldCpy< PdfField_fluid_T, gpu::GPUField< real_t > >(blocks, pdfFieldFluidID, pdfFieldFluidCPUGPUID);
         gpu::fieldCpy< PdfField_concentration_T, gpu::GPUField< real_t > >(blocks, pdfFieldConcentrationID, pdfFieldConcentrationCPUGPUID);
         gpu::fieldCpy< VelocityField_fluid_T, gpu::GPUField< real_t > >(blocks, velFieldFluidID, velFieldFluidCPUGPUID);
         gpu::fieldCpy< DensityField_concentration_T, gpu::GPUField< real_t > >(blocks, densityConcentrationFieldID, densityConcentrationFieldCPUGPUID);

#endif
         for (auto& block : *blocks)
         {
            getterSweep_fluid(&block);
         }
      });

      auto vtkOutput_Concentration = vtk::createVTKOutput_BlockData(blocks, "vtk files concentration", vtkSpacing, 0, false, vtkFolder);

      vtkOutput_Concentration->addBeforeFunction(communication_concentration);

      vtkOutput_Concentration->addBeforeFunction([&](){
         for (auto& block : *blocks)
         {
            getterSweep_concentration(&block);

         }
      });

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      vtkOutput_Fluid->addCellDataWriter(make_shared< field::VTKWriter< VelocityField_fluid_T > >(velFieldFluidID, "Fluid Velocity"));
#else
      vtkOutput_Fluid->addCellDataWriter(make_shared< field::VTKWriter< VelocityField_fluid_T > >(velFieldFluidCPUGPUID, "Fluid Velocity"));
#endif
      vtkOutput_Fluid->addCellDataWriter(make_shared< field::VTKWriter< DensityField_fluid_T > >(densityFluidFieldID, "Fluid Density"));
      vtkOutput_Fluid->addCellDataWriter(make_shared< field::VTKWriter< FlagField_T > >(flagFieldFluidID, "FluidFlagField"));

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      vtkOutput_Concentration->addCellDataWriter(make_shared< field::VTKWriter< DensityField_concentration_T > >(densityConcentrationFieldID, "Concentration"));
#else
      vtkOutput_Concentration->addCellDataWriter(make_shared< field::VTKWriter< DensityField_concentration_T > >(densityConcentrationFieldCPUGPUID, "Concentration"));
#endif

      vtkOutput_Concentration->addCellDataWriter(make_shared< field::VTKWriter< FlagField_T > >(flagFieldConcentrationID, "ConcentrationFlagField"));

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput_Fluid), "VTK output Fluid");
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput_Concentration), "VTK output Concentration");
   }

   if (vtkSpacing != uint_t(0)) { vtk::writeDomainDecomposition(blocks, "domain_decomposition", vtkFolder); }


   // Add performance logging
   lbm::PerformanceLogger< FlagField_T > performanceLoggerFluid(blocks, flagFieldFluidID, Fluid_Flag, performanceLogFrequency);
   lbm::PerformanceLogger< FlagField_T > performanceLoggerConcentration(blocks, flagFieldConcentrationID, Concentration_Flag, performanceLogFrequency);
   if (performanceLogFrequency > 0)
   {
      timeloop.addFuncAfterTimeStep(performanceLoggerFluid, "Evaluate performance logging of fluid");
   }

   // Add LBM (fluid and concentration) communication function and boundary handling sweep

   pystencils::LBMFluidSweep lbmFluidSweep(pdfFieldFluidCPUGPUID, velFieldFluidCPUGPUID,relaxationRateFluid);
   pystencils::LBMConcentrationSweep lbmConcentrationSweep(pdfFieldConcentrationCPUGPUID, velFieldFluidCPUGPUID,
                                                           relaxationRateConcentration);

   if (simulationName == "2d_plate")
   {
      timeloop.add() << BeforeFunction(communication_fluid, "LBM fluid Communication")
                     << Sweep(deviceSyncWrapper(density_fluid_bc.getSweep()), "Boundary Handling (Fluid Density)");
      timeloop.add() << Sweep(deviceSyncWrapper(neumann_fluid_bc.getSweep()), "Boundary Handling (Fluid Neumann)");
      timeloop.add() << Sweep(deviceSyncWrapper(ubb_fluid_bc.getSweep()), "Boundary Handling (Fluid UBB)");

      timeloop.add() << BeforeFunction(communication_concentration, "LBM concentration Communication")
                     << Sweep(deviceSyncWrapper(neumann_concentration_bc.getSweep()),
                              "Boundary Handling (Concentration Neumann)");

      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc_south.getSweep()),
                              "Boundary Handling (Concentration Density south)");

      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc_west.getSweep()),
                              "Boundary Handling (Concentration Density west)");
   }
   auto tracker = make_shared< lbm::TimestepTracker >(0);
   if (simulationName == "2d_poisuelle")
   {
      timeloop.add() << BeforeFunction(communication_fluid, "LBM fluid Communication")
                     << Sweep(deviceSyncWrapper(density_fluid_bc.getSweep()), "Boundary Handling (Fluid Density)");
      timeloop.add() << Sweep(deviceSyncWrapper(neumann_fluid_bc.getSweep()), "Boundary Handling (Fluid Neumann)");
      timeloop.add() << Sweep(deviceSyncWrapper(noSlip_fluid_bc.getSweep()), "Boundary Handling (no slip fluid)");
      timeloop.add() << Sweep(ubb_fluid_poisuelle.getSweep(tracker), "Boundary Handling (ubb poisuelle)");

      timeloop.add() << BeforeFunction(communication_concentration, "LBM concentration Communication")
                     << Sweep(deviceSyncWrapper(neumann_concentration_bc.getSweep()),
                              "Boundary Handling (Concentration Neumann)");

      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc_north.getSweep()),
                              "Boundary Handling (Concentration Density north)");

      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc_west.getSweep()),
                              "Boundary Handling (Concentration Density west)");
   }

   if (periodicInX && periodicInY && periodicInZ &&
       simulationName == "2d_gaussian") // executes for the 2d gaussian curve validation
   {
      timeloop.add() << BeforeFunction(communication_fluid, "LBM fluid Communication")
                     << Sweep(deviceSyncWrapper(lbmFluidSweep), "LBM Fluid sweep");
      timeloop.add() << BeforeFunction(communication_concentration, "LBM concentration Communication")
                     << Sweep(deviceSyncWrapper(lbmConcentrationSweep), "LBM Concentration sweep");
   }

   if (!periodicInX && !periodicInY && !periodicInZ &&
       simulationName == "2d_plate") // executes for the 2d plate validation
   {
      timeloop.add() << Sweep(deviceSyncWrapper(lbmFluidSweep), "LBM Fluid sweep");
      timeloop.add() << Sweep(deviceSyncWrapper(lbmConcentrationSweep), "LBM Concentration sweep");
   }
   if (!periodicInX && !periodicInY && !periodicInZ &&
       simulationName == "2d_poisuelle") // executes for the 2d plate validation
   {
      timeloop.add() << Sweep(deviceSyncWrapper(lbmFluidSweep), "LBM Fluid sweep");
      timeloop.add() << Sweep(deviceSyncWrapper(lbmConcentrationSweep), "LBM Concentration sweep");
   }

   WcTimingPool timeloopTiming;
   // TODO: maybe add warmup phase
   for (uint_t timeStep = 0; timeStep < timeSteps; ++timeStep)
   {
      if (useCommunicationHiding) { commTimeloop.singleStep(timeloopTiming); }
      timeloop.singleStep(timeloopTiming);
   }
   timeloopTiming.logResultOnRoot();
   auto timeloopTimingReduced = timeloopTiming.getReduced();











      return EXIT_SUCCESS;
}
} // namespace MaterialTransport



int main(int argc, char** argv) { MaterialTransport::main(argc, argv); }