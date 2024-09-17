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

#include "../utilities/InitializerFunctions.h"
#include "ConcentrationMacroGetter.h"
#include "FluidMacroGetter.h"
#include "GeneralInfoHeader.h"
#include "PackInfoConcentration.h"
#include "PackInfoFluid.h"

namespace MaterialTransport
{
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
const FlagUID Inflow_Fluid_Flag("Inflow_Fluid");
const FlagUID Neumann_Fluid_Flag("Neumann_Fluid");

// Concentration Flags
const FlagUID Concentration_Flag("Concentration");
const FlagUID Density_Concentration_Flag_west("Density_Concentration_west");
const FlagUID Density_Concentration_Flag_south("Density_Concentration_south");
const FlagUID NoSlip_Concentration_Flag("NoSlip_Concentration");
const FlagUID Inflow_Concentration_Flag("Inflow_Concentration");
// const FlagUID Dirichlet_Concentration_Flag("Dirichlet_Concentration");
const FlagUID Neumann_Concentration_Flag("Neumann_Concentration");

//////////
// MAIN //
//////////

int main(int argc, char** argv)
{
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

   Config::BlockHandle NumericalSetup = cfgFile->getBlock("NumericalSetup");
   const real_t xSize_SI              = NumericalSetup.getParameter< real_t >("xSize");
   const real_t ySize_SI              = NumericalSetup.getParameter< real_t >("ySize");
   const real_t zSize_SI              = NumericalSetup.getParameter< real_t >("zSize");
   const uint_t numXBlocks            = NumericalSetup.getParameter< uint_t >("numXBlocks");
   const uint_t numYBlocks            = NumericalSetup.getParameter< uint_t >("numYBlocks");
   const uint_t numZBlocks            = NumericalSetup.getParameter< uint_t >("numZBlocks");
   WALBERLA_CHECK_EQUAL(numXBlocks * numYBlocks * numZBlocks, uint_t(MPIManager::instance()->numProcesses()),
                        "When using GPUs, the number of blocks ("
                           << numXBlocks * numYBlocks * numZBlocks << ") has to match the number of MPI processes ("
                           << uint_t(MPIManager::instance()->numProcesses()) << ")");
   const real_t tSI                   = NumericalSetup.getParameter< real_t >("tSI");
   const bool periodicInX             = NumericalSetup.getParameter< bool >("periodicInX");
   const bool periodicInY             = NumericalSetup.getParameter< bool >("periodicInY");
   const bool periodicInZ             = NumericalSetup.getParameter< bool >("periodicInZ");
   const bool sendDirectlyFromGPU     = NumericalSetup.getParameter< bool >("sendDirectlyFromGPU");
   const bool useCommunicationHiding  = NumericalSetup.getParameter< bool >("useCommunicationHiding");
   const Vector3< uint_t > frameWidth = NumericalSetup.getParameter< Vector3< uint_t > >("frameWidth");

   const real_t uInflow_SI = NumericalSetup.getParameter< real_t >("uInflowL");
   const real_t kappa_SI    = NumericalSetup.getParameter< real_t >("kappa_SI");
   const real_t kappa_L     = NumericalSetup.getParameter< real_t >("kappa_L");
   const real_t dx_SI       = NumericalSetup.getParameter< real_t >("dx_SI");
   //const real_t uInflow_SI  = NumericalSetup.getParameter< real_t >("uInflowSI");
   const real_t dt_SI              = real_c((dx_SI * dx_SI * kappa_L) / (kappa_SI));
   const real_t uInflow_LBM  = (dt_SI/dx_SI)*uInflow_SI;

   const Vector3< real_t > Uinitialize(uInflow_LBM, uInflow_LBM, uInflow_LBM);
   //const Vector3< real_t > Uinitialize(uInflow_LBM,0,0);
   if ((periodicInY && numYBlocks == 1) || (periodicInZ && numZBlocks == 1))
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Using only 1 block in periodic dimensions can lead to unexpected behavior.")
   }

   Config::BlockHandle outputSetup      = cfgFile->getBlock("Output");
   const uint_t vtkSpacing              = outputSetup.getParameter< uint_t >("vtkSpacing");
   const std::string vtkFolder          = outputSetup.getParameter< std::string >("vtkFolder");
   const uint_t performanceLogFrequency = outputSetup.getParameter< uint_t >("performanceLogFrequency");


   const real_t omegaConcentration = lbm::collision_model::omegaFromViscosity(kappa_L);

   Vector3< uint_t > domainSize(uint_c((xSize_SI / dx_SI)), uint_c((ySize_SI / dx_SI)), uint_c((zSize_SI/dx_SI)));
   const uint_t timeSteps   = uint_c(std::ceil(tSI / dt_SI));
   const uint_t tinitial    = uint_c(0);
   const real_t sigma_0     = real_t(1);
   const real_t diffusivity = kappa_L;
   const real_t sigma_D     = real_t(std::sqrt(2 * diffusivity * tinitial));

   const Vector3< real_t > x_0(real_t(domainSize[0] / 2), real_t(domainSize[1] / 2), real_t(domainSize[2] / 2));
   //const Vector3< real_t > x_0(0.5,0.5,0.5);
   Vector3< uint_t > cellsPerBlockPerDirection(domainSize[0] / numXBlocks, domainSize[1] / numYBlocks,
                                               domainSize[2] / numZBlocks);
   const real_t courant_number = (uInflow_SI*dt_SI)/dx_SI;
   const real_t peclet_number = (xSize_SI*uInflow_SI)/kappa_SI;
   WALBERLA_LOG_INFO_ON_ROOT("Simulation setup:");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - dt_SI = " << dt_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - dx_SI = " << dx_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - concentration: diffusivity (LBM) = " << diffusivity);
   WALBERLA_LOG_INFO_ON_ROOT(" - Fluid Lattice Velocity = " << uInflow_LBM);
   WALBERLA_LOG_INFO_ON_ROOT(" - Number of time steps = " << timeSteps);
   WALBERLA_LOG_INFO_ON_ROOT(" - omega (concentration) = " << omegaConcentration);
   WALBERLA_LOG_INFO_ON_ROOT(" - relaxation time (tau) = " << 1/(omegaConcentration));
   WALBERLA_LOG_INFO_ON_ROOT(" - Courant number = " << courant_number);
   WALBERLA_LOG_INFO_ON_ROOT(" - Peclet number = " <<  peclet_number);


   WALBERLA_CHECK_EQUAL(domainSize[0], cellsPerBlockPerDirection[0] * numXBlocks,
                        "number of cells in x of " << domainSize[0]
                                                   << " is not divisible by given number of blocks in x direction");
   WALBERLA_CHECK_EQUAL(domainSize[1], cellsPerBlockPerDirection[1] * numYBlocks,
                        "number of cells in y of " << domainSize[1]
                                                   << " is not divisible by given number of blocks in y direction");
   WALBERLA_CHECK_EQUAL(domainSize[2], cellsPerBlockPerDirection[2] * numZBlocks,
                        "number of cells in z of " << domainSize[2]
                                                   << " is not divisible by given number of blocks in z direction");

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   shared_ptr< StructuredBlockForest > blocks =
      blockforest::createUniformBlockGrid(numXBlocks, numYBlocks, numZBlocks, cellsPerBlockPerDirection[0],
                                          cellsPerBlockPerDirection[1], cellsPerBlockPerDirection[2], real_t(1),
                                          uint_t(0), false, false, periodicInX, periodicInY, periodicInZ, // periodicity
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
   BlockDataID pdfFieldFluidCPUGPUID =
      gpu::addGPUFieldToStorage< PdfField_fluid_T >(blocks, pdfFieldFluidID, "pdf fluid field GPU");

   // Concentration PDFs on GPU
   BlockDataID pdfFieldConcentrationID = field::addToStorage< PdfField_concentration_T >(
      blocks, "pdf concentration field (fzyx)", real_c(std::nan("")), field::fzyx);
   BlockDataID pdfFieldConcentrationCPUGPUID = gpu::addGPUFieldToStorage< PdfField_concentration_T >(
      blocks, pdfFieldConcentrationID, "pdf concentration field GPU");

   // Fluid velocity field on GPU
   velFieldFluidID =
      field::addToStorage< VelocityField_fluid_T >(blocks, "velocity fluid field", real_t(0), field::fzyx);
   BlockDataID velFieldFluidCPUGPUID =
      gpu::addGPUFieldToStorage< VelocityField_fluid_T >(blocks, velFieldFluidID, "velocity fluid field GPU");

   // Concentration Density on GPU
   densityConcentrationFieldID = field::addToStorage< DensityField_concentration_T >(
      blocks, "density concentration field", real_t(0), field::fzyx);
   BlockDataID densityConcentrationFieldCPUGPUID = gpu::addGPUFieldToStorage< DensityField_concentration_T >(
      blocks, densityConcentrationFieldID, "density concentration field GPU");

#else

   // Fluid PDFs on CPU
   BlockDataID pdfFieldFluidCPUGPUID =
      field::addToStorage< PdfField_fluid_T >(blocks, "pdf fluid field CPU", real_c(std::nan("")), field::fzyx);

   BlockDataID velFieldFluidCPUGPUID =
      field::addToStorage< VelocityField_fluid_T >(blocks, "velocity fluid field CPU", real_t(0), field::fzyx);

   // Concentration PDFs on CPU
   BlockDataID pdfFieldConcentrationCPUGPUID = field::addToStorage< PdfField_concentration_T >(
      blocks, "pdf concentration field CPU", real_c(std::nan("")), field::fzyx);

   BlockDataID densityConcentrationFieldCPUGPUID = field::addToStorage< DensityField_concentration_T >(
      blocks, "density concentration field", real_t(0), field::fzyx);

#endif
   BlockDataID densityFluidFieldID =
      field::addToStorage< DensityField_fluid_T >(blocks, "density fluid field", real_t(0), field::fzyx);
   densityConcentrationFieldID = field::addToStorage< DensityField_concentration_T >(
      blocks, "density concentration field", real_t(0), field::fzyx);
   BlockDataID flagFieldFluidID = field::addFlagFieldToStorage< FlagField_T >(blocks, "fluid flag field");
   BlockDataID flagFieldConcentrationID =
      field::addFlagFieldToStorage< FlagField_T >(blocks, "concentration flag field");
   BlockDataID analyticalConcentrationFieldID = field::addToStorage< DensityField_concentration_T >(
      blocks, "density concentration field analytical", real_t(0), field::fzyx);
   BlockDataID ErrorFieldID = field::addToStorage< DensityField_concentration_T >(
      blocks, "density concentration field error", real_t(0), field::fzyx);

   // Assemble boundary block string
   std::string boundariesBlockString = " BoundariesFluid";

   if (!periodicInX)
   {
      boundariesBlockString += "{"
                               "Border { direction W;    walldistance -1;  flag Inflow_Fluid; }"
                               "Border { direction E;    walldistance -1;  flag Density_Fluid; }";
   }

   if (!periodicInY)
   {
      boundariesBlockString += "Border { direction S;    walldistance -1;  flag Neumann_Fluid; }"
                               "Border { direction N;    walldistance -1;  flag Neumann_Fluid; }";
   }

   if (!periodicInZ)
   {
      boundariesBlockString += "Border { direction T;    walldistance -1;  flag Neumann_Fluid; }"
                               "Border { direction B;    walldistance -1;  flag Neumann_Fluid; }";
   }

   boundariesBlockString += "}";

   boundariesBlockString += "\n BoundariesConcentration";

   if (!periodicInX)
   {
      boundariesBlockString += "{"
                               "Border { direction W;    walldistance -1;  flag Density_Concentration_west; }"
                               "Border { direction E;    walldistance -1;  flag Neumann_Concentration; }";
   }
   if (!periodicInY)
   {
      boundariesBlockString += "Border { direction S;    walldistance -1;  flag Density_Concentration_south; }"
                               "Border { direction N;    walldistance -1;  flag Neumann_Concentration; }";
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
   auto boundariesConfigFluid         = boundariesCfgFile.getBlock("BoundariesFluid");
   auto boundariesConfigConcentration = boundariesCfgFile.getBlock("BoundariesConcentration");

   // map boundaries into the fluid field simulation
   lbm::BC_Fluid_Density density_fluid_bc(blocks, pdfFieldFluidCPUGPUID, real_t(1.0));
   lbm::BC_Fluid_NoSlip noSlip_fluid_bc(blocks, pdfFieldFluidCPUGPUID);
   lbm::BC_Fluid_UBB ubb_fluid_bc(blocks, pdfFieldFluidCPUGPUID, uInflow_LBM, real_t(0), real_t(0));
   lbm::BC_Fluid_Neumann neumann_fluid_bc(blocks, pdfFieldFluidCPUGPUID);
   lbm::BC_Concentration_Density density_concentration_bc_west(blocks, pdfFieldConcentrationCPUGPUID, real_t(0));
   lbm::BC_Concentration_Density density_concentration_bc_south(blocks, pdfFieldConcentrationCPUGPUID, real_t(1));
   lbm::BC_Concentration_Neumann neumann_concentration_bc(blocks, pdfFieldConcentrationCPUGPUID);

   if (periodicInX && periodicInY && periodicInZ)
   { // executes for the 2d gaussian curve validation

      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldFluidID, Fluid_Flag);
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldConcentrationID, Concentration_Flag);
   }

   ///////////////
   // TIME LOOP //
   ///////////////
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   initConcentrationField(blocks, densityConcentrationFieldID, simulationDomain, domainSize);
   initFluidField(blocks, velFieldFluidID, Uinitialize);

#else


   /*initConcentrationFieldGaussian(blocks, densityConcentrationFieldCPUGPUID, simulationDomain, domainSize, sigma_0,
                                  sigma_D, Uinitialize, x_0);*/
   // initConcentrationField(blocks,densityConcentrationFieldCPUGPUID,simulationDomain,domainSize);
   initConcentrationFieldSinusoidal(blocks,
                                    densityConcentrationFieldCPUGPUID,simulationDomain,
                                    domainSize,sigma_0,sigma_D,Uinitialize,x_0,dx_SI,dt_SI);

   /*initConcentrationFieldPacket(blocks,
                                densityConcentrationFieldCPUGPUID,simulationDomain,
                                domainSize,sigma_0,sigma_D,
                                Uinitialize,x_0,dx_SI,dt_SI,diffusivity);*/


   initFluidField(blocks, velFieldFluidCPUGPUID, Uinitialize);

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
   gpu::communication::UniformGPUScheme< Stencil_Concentration_T > com_concentration(blocks, sendDirectlyFromGPU,
                                                                                     false);
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
   pystencils::FluidMacroGetter getterSweep_fluid(densityFluidFieldID, pdfFieldFluidID, velFieldFluidID, real_t(0),
                                                  real_t(0), real_t(0));
#else
   pystencils::FluidMacroGetter getterSweep_fluid(densityFluidFieldID, pdfFieldFluidCPUGPUID, velFieldFluidCPUGPUID);
#endif

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   pystencils::ConcentrationMacroGetter getterSweep_concentration(densityConcentrationFieldID, pdfFieldConcentrationID,
                                                                  velFieldFluidID, real_t(0), real_t(0), real_t(0));
#else
   pystencils::ConcentrationMacroGetter getterSweep_concentration(densityConcentrationFieldCPUGPUID,
                                                                  pdfFieldConcentrationCPUGPUID);
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
         gpu::fieldCpy< PdfField_concentration_T, gpu::GPUField< real_t > >(blocks, pdfFieldConcentrationID,
                                                                            pdfFieldConcentrationCPUGPUID);
         gpu::fieldCpy< VelocityField_fluid_T, gpu::GPUField< real_t > >(blocks, velFieldFluidID,
                                                                         velFieldFluidCPUGPUID);
         gpu::fieldCpy< DensityField_concentration_T, gpu::GPUField< real_t > >(blocks, densityConcentrationFieldID,
                                                                                densityConcentrationFieldCPUGPUID);

#endif
         for (auto& block : *blocks)
         {
            getterSweep_fluid(&block);
         }
      });

      auto vtkOutput_Concentration =
         vtk::createVTKOutput_BlockData(blocks, "vtk files concentration", vtkSpacing, 0, false, vtkFolder);

      vtkOutput_Concentration->addBeforeFunction(communication_concentration);

      vtkOutput_Concentration->addBeforeFunction([&]() {
         for (auto& block : *blocks)
         {
            getterSweep_concentration(&block);
         }
      });

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      vtkOutput_Fluid->addCellDataWriter(
         make_shared< field::VTKWriter< VelocityField_fluid_T > >(velFieldFluidID, "Fluid Velocity"));
#else
      vtkOutput_Fluid->addCellDataWriter(
         make_shared< field::VTKWriter< VelocityField_fluid_T > >(velFieldFluidCPUGPUID, "Fluid Velocity"));
#endif
      vtkOutput_Fluid->addCellDataWriter(
         make_shared< field::VTKWriter< DensityField_fluid_T > >(densityFluidFieldID, "Fluid Density"));
      vtkOutput_Fluid->addCellDataWriter(
         make_shared< field::VTKWriter< FlagField_T > >(flagFieldFluidID, "FluidFlagField"));

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      vtkOutput_Concentration->addCellDataWriter(
         make_shared< field::VTKWriter< DensityField_concentration_T > >(densityConcentrationFieldID, "Concentration"));
#else
      vtkOutput_Concentration->addCellDataWriter(make_shared< field::VTKWriter< DensityField_concentration_T > >(
         densityConcentrationFieldCPUGPUID, "Concentration"));
#endif

      vtkOutput_Concentration->addCellDataWriter(
         make_shared< field::VTKWriter< FlagField_T > >(flagFieldConcentrationID, "ConcentrationFlagField"));
      vtkOutput_Concentration->addCellDataWriter(make_shared< field::VTKWriter< DensityField_concentration_T > >(
         analyticalConcentrationFieldID, "Analytical Concentration"));
      vtkOutput_Concentration->addCellDataWriter(
         make_shared< field::VTKWriter< DensityField_concentration_T > >(ErrorFieldID, "Error Concentration"));

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput_Fluid), "VTK output Fluid");
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput_Concentration), "VTK output Concentration");
   }

   if (vtkSpacing != uint_t(0)) { vtk::writeDomainDecomposition(blocks, "domain_decomposition", vtkFolder); }

   // Add performance logging
   lbm::PerformanceLogger< FlagField_T > performanceLoggerFluid(blocks, flagFieldFluidID, Fluid_Flag,
                                                                performanceLogFrequency);
   lbm::PerformanceLogger< FlagField_T > performanceLoggerConcentration(blocks, flagFieldConcentrationID,
                                                                        Concentration_Flag, performanceLogFrequency);
   if (performanceLogFrequency > 0)
   {
      timeloop.addFuncAfterTimeStep(performanceLoggerFluid, "Evaluate performance logging of fluid");
   }

   // Add LBM (fluid and concentration) communication function and boundary handling sweep

   // pystencils::LBMFluidSweep lbmFluidSweep(pdfFieldFluidCPUGPUID, relaxationRateFluid);
   pystencils::LBMConcentrationSweep lbmConcentrationSweep(pdfFieldConcentrationCPUGPUID, velFieldFluidCPUGPUID,
                                                           omegaConcentration);

   if (periodicInX && periodicInY && periodicInZ) // executes for the 2d gaussian curve validation
   {
      /*timeloop.add() << BeforeFunction(communication_fluid, "LBM fluid Communication")
                     << Sweep(deviceSyncWrapper(lbmFluidSweep), "LBM Fluid sweep");*/
      timeloop.add() << BeforeFunction(communication_concentration, "LBM concentration Communication")
                     << Sweep(deviceSyncWrapper(lbmConcentrationSweep), "LBM Concentration sweep");
   }

   auto calculate_macroscopic_quantities = [&]() {
      for (auto& block : *blocks)
      {
         getterSweep_concentration(&block);
      }
   };
   timeloop.addFuncAfterTimeStep(calculate_macroscopic_quantities, "calculate macroscopic quantities");

   WcTimingPool timeloopTiming;
   // TODO: maybe add warmup phase

   /*analyticalSolGaussian(blocks, analyticalConcentrationFieldID, simulationDomain, domainSize, sigma_0, diffusivity,
                         Uinitialize, x_0, 0,0);*/
   analyticalSolSinusoidal(blocks,
                           analyticalConcentrationFieldID,simulationDomain,
                           domainSize,sigma_0,diffusivity,
                           Uinitialize,x_0,0,dx_SI,dt_SI);
   /*analyticalSolPacket(blocks,
                       analyticalConcentrationFieldID,simulationDomain,
                       domainSize,sigma_0,diffusivity,
                       Uinitialize,x_0,0,dx_SI,dt_SI);*/
   real_t cummulative_norm = 0;

   for (uint_t timeStep = 1; timeStep <= timeSteps; ++timeStep)
   {
      if (useCommunicationHiding) { commTimeloop.singleStep(timeloopTiming); }
      timeloop.singleStep(timeloopTiming);

      /*analyticalSolGaussian(blocks, analyticalConcentrationFieldID, simulationDomain, domainSize, sigma_0, diffusivity,
                            Uinitialize, x_0, timeStep,0);*/
      analyticalSolSinusoidal(blocks,
                                   analyticalConcentrationFieldID,simulationDomain,
                                   domainSize,sigma_0,diffusivity,
                                   Uinitialize,x_0,timeStep,dx_SI,dt_SI);

      /*analyticalSolPacket(blocks,
                          analyticalConcentrationFieldID,simulationDomain,
                          domainSize,sigma_0,diffusivity,
                          Uinitialize,x_0,timeStep,dx_SI,dt_SI);*/

      std::vector< real_t > norms = computeErrorL2(blocks, densityConcentrationFieldCPUGPUID,
                                                   analyticalConcentrationFieldID, ErrorFieldID, simulationDomain);
      WALBERLA_LOG_INFO_ON_ROOT("Infinity norm is " << norms[0] << "\n L1 norm is " << norms[1] << "\n L2 norm is "
                                                    << norms[2]);
      cummulative_norm += norms[2];
   }
   cummulative_norm  = cummulative_norm/timeSteps;
   WALBERLA_LOG_INFO_ON_ROOT("cummulative L2 norm is "<< cummulative_norm);
   timeloopTiming.logResultOnRoot();
   auto timeloopTimingReduced = timeloopTiming.getReduced();

   return EXIT_SUCCESS;
}
} // namespace MaterialTransport

int main(int argc, char** argv) { MaterialTransport::main(argc, argv); }