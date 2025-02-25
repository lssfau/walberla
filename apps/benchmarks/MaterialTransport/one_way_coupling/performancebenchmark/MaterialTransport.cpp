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
const FlagUID Inflow_Fluid_Flag("Inflow_Fluid");

// Concentration Flags
const FlagUID Concentration_Flag("Concentration");
const FlagUID Density_Concentration_Flag("Density_Concentration");
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
   const uint_t numXBlocks            = numericalSetup.getParameter< uint_t >("numXBlocks");
   const uint_t numYBlocks            = numericalSetup.getParameter< uint_t >("numYBlocks");
   const uint_t numZBlocks            = numericalSetup.getParameter< uint_t >("numZBlocks");
   WALBERLA_CHECK_EQUAL(numXBlocks * numYBlocks * numZBlocks, uint_t(MPIManager::instance()->numProcesses()),
                        "When using GPUs, the number of blocks ("
                           << numXBlocks * numYBlocks * numZBlocks << ") has to match the number of MPI processes ("
                           << uint_t(MPIManager::instance()->numProcesses()) << ")");
   const bool periodicInY                 = numericalSetup.getParameter< bool >("periodicInY");
   const bool periodicInZ                 = numericalSetup.getParameter< bool >("periodicInZ");
   const uint_t numXCellsPerBlock         = numericalSetup.getParameter< uint_t >("numXCellsPerBlock");
   const uint_t numYCellsPerBlock         = numericalSetup.getParameter< uint_t >("numYCellsPerBlock");
   const uint_t numZCellsPerBlock         = numericalSetup.getParameter< uint_t >("numZCellsPerBlock");
   const bool sendDirectlyFromGPU         = numericalSetup.getParameter< bool >("sendDirectlyFromGPU");
   const bool useCommunicationHiding      = numericalSetup.getParameter< bool >("useCommunicationHiding");
   const Vector3< uint_t > frameWidth     = numericalSetup.getParameter< Vector3< uint_t > >("frameWidth");
   const uint_t timeSteps                 = numericalSetup.getParameter< uint_t >("timeSteps");

   const real_t uInflow        = numericalSetup.getParameter< real_t >("uInflow");
   const real_t relaxationRate = numericalSetup.getParameter< real_t >("relaxationRate");
   const Vector3< real_t > Uinitialize (uInflow,0,0);
   if ((periodicInY && numYBlocks == 1) || (periodicInZ && numZBlocks == 1))
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Using only 1 block in periodic dimensions can lead to unexpected behavior.")
   }
   const real_t viscosity = lbm::collision_model::viscosityFromOmega(relaxationRate);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(viscosity)

   Config::BlockHandle outputSetup      = cfgFile->getBlock("Output");
   const uint_t vtkSpacing              = outputSetup.getParameter< uint_t >("vtkSpacing");
   const std::string vtkFolder          = outputSetup.getParameter< std::string >("vtkFolder");
   const uint_t performanceLogFrequency = outputSetup.getParameter< uint_t >("performanceLogFrequency");

   Vector3< uint_t > domainSize;
   domainSize[0] = numXBlocks * numXCellsPerBlock;
   domainSize[1] = numYBlocks * numYCellsPerBlock;
   domainSize[2] = numZBlocks * numZCellsPerBlock;

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const bool periodicInX                     = false;
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
      velFieldFluidID  = field::addToStorage< VelocityField_fluid_T >(blocks, "velocity fluid field", real_t(0), field::fzyx);
      // Concentration PDFs on CPU
      BlockDataID pdfFieldConcentrationCPUGPUID =
         field::addToStorage< PdfField_concentration_T >(blocks, "pdf concentration field CPU", real_c(std::nan("")), field::fzyx);

      BlockDataID densityConcentrationFieldCPUGPUID = field::addToStorage< DensityField_concentration_T >(blocks, "density concentration field", real_t(0), field::fzyx);

   #endif
      BlockDataID densityFluidFieldID = field::addToStorage< DensityField_fluid_T >(blocks, "density fluid field", real_t(0), field::fzyx);
      densityConcentrationFieldID = field::addToStorage< DensityField_concentration_T >(blocks, "density concentration field", real_t(0), field::fzyx);
      BlockDataID flagFieldFluidID = field::addFlagFieldToStorage< FlagField_T >(blocks, "fluid flag field");
      BlockDataID flagFieldConcentrationID = field::addFlagFieldToStorage< FlagField_T >(blocks, "concentration flag field");


      // Assemble boundary block string
      std::string boundariesBlockString = " BoundariesFluid"
                                          "{"
                                          "Border { direction W;    walldistance -1;  flag Inflow_Fluid; }"
                                          "Border { direction E;    walldistance -1;  flag Density_Fluid; }";

      if (!periodicInY)
      {
         boundariesBlockString += "Border { direction S;    walldistance -1;  flag NoSlip_Fluid; }"
                                  "Border { direction N;    walldistance -1;  flag NoSlip_Fluid; }";
      }

      if (!periodicInZ)
      {
         boundariesBlockString += "Border { direction T;    walldistance -1;  flag NoSlip_Fluid; }"
                                  "Border { direction B;    walldistance -1;  flag NoSlip_Fluid; }";
      }

      boundariesBlockString += "}";

      boundariesBlockString += "\n BoundariesConcentration";

      boundariesBlockString += "{"
                                "Border { direction W;    walldistance -1;  flag Neumann_Concentration; }"
                                "Border { direction E;    walldistance -1;  flag Neumann_Concentration; }";

      if (!periodicInY)
      {
         boundariesBlockString += "Border { direction S;    walldistance -1;  flag Neumann_Concentration; }"
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
      auto boundariesConfigFluid = boundariesCfgFile.getBlock("BoundariesFluid");
      auto boundariesConfigConcentration = boundariesCfgFile.getBlock("BoundariesConcentration");


      // map boundaries into the fluid field simulation
      geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldFluidID, boundariesConfigFluid);
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldFluidID, Fluid_Flag);
      lbm::BC_Fluid_Density density_fluid_bc(blocks, pdfFieldFluidCPUGPUID, real_t(1.0));
      density_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Density_Fluid_Flag, Fluid_Flag);
      lbm::BC_Fluid_NoSlip noSlip_fluid_bc(blocks, pdfFieldFluidCPUGPUID);
      noSlip_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, NoSlip_Fluid_Flag, Fluid_Flag);
      lbm::BC_Fluid_UBB ubb_fluid_bc(blocks, pdfFieldFluidCPUGPUID, uInflow, real_t(0), real_t(0));
      ubb_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Inflow_Fluid_Flag, Fluid_Flag);

      // map boundaries into the concentration field simulation
      geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldConcentrationID, boundariesConfigConcentration);
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldConcentrationID, Concentration_Flag);
      lbm::BC_Concentration_Density density_concentration_bc(blocks, pdfFieldConcentrationCPUGPUID, real_t(1.0));
      density_concentration_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldConcentrationID, Density_Concentration_Flag, Concentration_Flag);
      lbm::BC_Concentration_Neumann neumann_concentration_bc(blocks,pdfFieldConcentrationCPUGPUID);
      neumann_concentration_bc.fillFromFlagField<FlagField_T>(blocks,flagFieldConcentrationID,Neumann_Concentration_Flag,Concentration_Flag);

      ///////////////
      // TIME LOOP //
      ///////////////
      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      initConcentrationField(blocks,densityConcentrationFieldID,simulationDomain,domainSize);
      initFluidField(blocks, velFieldFluidID , Uinitialize,domainSize);
      gpu::fieldCpy<gpu::GPUField< real_t >, DensityField_concentration_T>(blocks, densityConcentrationFieldCPUGPUID,densityConcentrationFieldID);
      WALBERLA_LOG_INFO_ON_ROOT("code reached here on gpu");
      gpu::fieldCpy< gpu::GPUField< real_t >, VelocityField_fluid_T >(blocks, velFieldFluidCPUGPUID, velFieldFluidID);
      pystencils::InitializeFluidDomain initializeFluidDomain(pdfFieldFluidCPUGPUID,velFieldFluidCPUGPUID,real_t(0),real_t(0),real_t(0),real_t(1));
      pystencils::InitializeConcentrationDomain initializeConcentrationDomain(densityConcentrationFieldCPUGPUID ,pdfFieldConcentrationCPUGPUID,velFieldFluidCPUGPUID,real_t(0),real_t(0),real_t(0));
      #else
      initConcentrationField(blocks,densityConcentrationFieldCPUGPUID,simulationDomain,domainSize);
      initFluidField(blocks, velFieldFluidCPUGPUID, Uinitialize,domainSize);
      pystencils::InitializeFluidDomain initializeFluidDomain(pdfFieldFluidCPUGPUID,velFieldFluidCPUGPUID,real_t(0),real_t(0),real_t(0),real_t(1));
      pystencils::InitializeConcentrationDomain initializeConcentrationDomain(densityConcentrationFieldCPUGPUID ,pdfFieldConcentrationCPUGPUID,velFieldFluidCPUGPUID,real_t(0),real_t(0),real_t(0));

      #endif

      //pystencils::InitializeFluidDomain initializeFluidDomain(pdfFieldFluidCPUGPUID,velFieldFluidCPUGPUID,real_t(0),real_t(0),real_t(0),real_t(1));
      //pystencils::InitializeConcentrationDomain initializeConcentrationDomain(densityConcentrationFieldCPUGPUID ,pdfFieldConcentrationCPUGPUID,velFieldFluidCPUGPUID,real_t(0),real_t(0),real_t(0));

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
   pystencils::FluidMacroGetter getterSweep_fluid(densityFluidFieldID,pdfFieldFluidCPUGPUID,velFieldFluidID,real_t(0),real_t(0),real_t(0));
#endif

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   pystencils:: ConcentrationMacroGetter getterSweep_concentration(densityConcentrationFieldID,pdfFieldConcentrationID);
#else
   pystencils::ConcentrationMacroGetter getterSweep_concentration(densityConcentrationFieldID,pdfFieldConcentrationCPUGPUID);
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
      vtkOutput_Concentration->addCellDataWriter(make_shared< field::VTKWriter< DensityField_concentration_T > >(densityConcentrationFieldID, "Concentration"));
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
      //timeloop.addFuncAfterTimeStep(performanceLoggerConcentration, "Evaluate performance logging of Concentration");
   }

   // Add LBM (fluid and concentration) communication function and boundary handling sweep
   if (useCommunicationHiding)
   {
      timeloop.add() << Sweep(deviceSyncWrapper(density_fluid_bc.getSweep()), "Boundary Handling (Fluid Density)");
      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc.getSweep()), "Boundary Handling (Concentration Density)");
   }
   else
   {
      timeloop.add() << BeforeFunction(communication_fluid, "LBM fluid Communication")
                     << Sweep(deviceSyncWrapper(density_fluid_bc.getSweep()), "Boundary Handling (Fluid Density)");
      timeloop.add() << BeforeFunction(communication_concentration, "LBM concentration Communication")
                     << Sweep(deviceSyncWrapper(neumann_concentration_bc.getSweep()), "Boundary Handling (Concentration Neumann)");
   }
   timeloop.add() << Sweep(deviceSyncWrapper(ubb_fluid_bc.getSweep()), "Boundary Handling (Fluid UBB)");
   //timeloop.add() << Sweep(deviceSyncWrapper(neumann_concentration_bc.getSweep()), "Boundary Handling (Concentration Neumann)");
   //timeloop.add() << Sweep(deviceSyncWrapper(ubb_concentration_bc.getSweep()), "Boundary Handling (Concentration UBB)");
   if (!periodicInY || !periodicInZ)
   {
      timeloop.add() << Sweep(deviceSyncWrapper(noSlip_fluid_bc.getSweep()), "Boundary Handling (Fluid NoSlip)");
      //timeloop.add() << Sweep(deviceSyncWrapper(noSlip_concentration_bc.getSweep()), "Boundary Handling (Concentration NoSlip)");
   }

   pystencils::LBMFluidSweep lbmFluidSweep(pdfFieldFluidCPUGPUID, velFieldFluidCPUGPUID,real_t(0.0), real_t(0.0), real_t(0.0), relaxationRate);
   pystencils::LBMConcentrationSweep lbmConcentrationSweep(pdfFieldConcentrationCPUGPUID,velFieldFluidCPUGPUID,real_t(0.0), real_t(0.0), real_t(0.0), relaxationRate);
   pystencils::LBMFluidSplitSweep lbmFluidSplitSweep(pdfFieldFluidCPUGPUID,velFieldFluidCPUGPUID,real_t(0.0), real_t(0.0), real_t(0.0), relaxationRate,frameWidth);
   pystencils::LBMConcentrationSplitSweep lbmConcentrationSplitSweep(pdfFieldConcentrationCPUGPUID,velFieldFluidCPUGPUID,real_t(0.0), real_t(0.0), real_t(0.0), relaxationRate,frameWidth);

   if (useCommunicationHiding)
   {
      commTimeloop.add() << BeforeFunction([&]() { com_fluid.startCommunication(); }, "LBM fluid Communication (start)")
                         << Sweep(deviceSyncWrapper(lbmFluidSplitSweep.getInnerSweep()), "LBM fluid inner sweep")
                         << AfterFunction([&]() { com_fluid.wait(); }, "LBM fluid Communication (wait)");
      timeloop.add() << Sweep(deviceSyncWrapper(lbmFluidSplitSweep.getOuterSweep()), "LBM fluid outer sweep");


      commTimeloop.add() << BeforeFunction([&]() { com_concentration.startCommunication(); }, "LBM concentration Communication (start)")
                         << Sweep(deviceSyncWrapper(lbmConcentrationSplitSweep.getInnerSweep()), "LBM concentration inner sweep")
                         << AfterFunction([&]() { com_concentration.wait(); }, "LBM concentration Communication (wait)");
      timeloop.add() << Sweep(deviceSyncWrapper(lbmConcentrationSplitSweep.getOuterSweep()), "LBM concentration outer sweep");

   }
   else
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