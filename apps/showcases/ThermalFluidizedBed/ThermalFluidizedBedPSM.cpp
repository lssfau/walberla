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
//! \file ThermalFluidizedBedPSM.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Ravi Ayyala Somayajula <ravi.k.ayyala@fau.de>
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Modification of showcases/FluidizedBed/FluidizedBedPSM.cpp for incorporating temperature
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/Debug.h"
#include "core/grid_generator/SCIterator.h"
#include "core/logging/all.h"
#include "core/math/all.h"
#include "core/mpi/Broadcast.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/waLBerlaBuildInfo.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/vtk/all.h"

#include "geometry/InitBoundaryHandling.h"

#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/communication/UniformGPUScheme.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/vtk/all.h"


#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/LinkedCells.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/AssocToBlock.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/InsertParticleIntoLinkedCells.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"


// codegen folder file includes
#include "lbm_mesapd_coupling/DataTypesCodegen.h"
#include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMSweepCollection.h"


// code generated files
#include "ConcentrationMacroGetter.h"
#include "FluidMacroGetter.h"
#include "GeneralInfoHeader.h"
#include "PSMFluidSweep.h"
#include "PackInfoConcentration.h"
#include "PackInfoFluid.h"

namespace thermalfluidized_bed
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

const FlagUID Fluid_Flag("fluid");
const FlagUID NoSlip_Flag("no slip");
const FlagUID Inflow_Flag("inflow");
const FlagUID Outflow_Flag("outflow");


///////////////////
// FLAGS Codegen //
//////////////////

// Fluid Flags
const FlagUID Density_Fluid_Flag("Density_Fluid");
const FlagUID NoSlip_Fluid_Flag("NoSlip_Fluid");
const FlagUID Inflow_Fluid_Flag("Inflow_Fluid");
const FlagUID OutFlow_Fluid_Flag("OutFlow_Fluid");

// Concentration Flags
const FlagUID Concentration_Flag("Concentration");
const FlagUID Density_Concentration_Flag_west("Density_Concentration_west");
const FlagUID Density_Concentration_Flag_east("Density_Concentration_east");
const FlagUID NoSlip_Concentration_Flag("NoSlip_Concentration");
const FlagUID Inflow_Concentration_Flag("Inflow_Concentration");
const FlagUID Neumann_Concentration_Flag("Neumann_Concentration");


//*******************************************************************************************************************

void createPlane(const shared_ptr< mesa_pd::data::ParticleStorage >& ps,
                 const shared_ptr< mesa_pd::data::ShapeStorage >& ss, Vector3< real_t > position,
                 Vector3< real_t > normal)
{
   mesa_pd::data::Particle&& p0 = *ps->create(true);
   p0.setPosition(position);
   p0.setInteractionRadius(std::numeric_limits< real_t >::infinity());
   p0.setShapeID(ss->create< mesa_pd::data::HalfSpace >(normal));
   p0.setOwner(mpi::MPIManager::instance()->rank());
   p0.setType(0);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
}

void createPlaneSetup(const shared_ptr< mesa_pd::data::ParticleStorage >& ps,
                      const shared_ptr< mesa_pd::data::ShapeStorage >& ss, const math::AABB& simulationDomain,
                      bool periodicInX, bool periodicInY, real_t offsetAtInflow, real_t offsetAtOutflow)
{
   createPlane(ps, ss, simulationDomain.minCorner() + Vector3< real_t >(0, 0, offsetAtInflow),
               Vector3< real_t >(0, 0, 1));
   createPlane(ps, ss, simulationDomain.maxCorner() + Vector3< real_t >(0, 0, offsetAtOutflow),
               Vector3< real_t >(0, 0, -1));

   if (!periodicInX)
   {
      createPlane(ps, ss, simulationDomain.minCorner(), Vector3< real_t >(1, 0, 0));
      createPlane(ps, ss, simulationDomain.maxCorner(), Vector3< real_t >(-1, 0, 0));
   }

   if (!periodicInY)
   {
      createPlane(ps, ss, simulationDomain.minCorner(), Vector3< real_t >(0, 1, 0));
      createPlane(ps, ss, simulationDomain.maxCorner(), Vector3< real_t >(0, -1, 0));
   }
}

struct ParticleInfo
{
   real_t averageVelocity = 0_r;
   real_t maximumVelocity = 0_r;
   uint_t numParticles    = 0;
   real_t maximumHeight   = 0_r;
   real_t particleVolume  = 0_r;
   real_t heightOfMass    = 0_r;

   void allReduce()
   {
      walberla::mpi::allReduceInplace(numParticles, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(averageVelocity, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(maximumVelocity, walberla::mpi::MAX);
      walberla::mpi::allReduceInplace(maximumHeight, walberla::mpi::MAX);
      walberla::mpi::allReduceInplace(particleVolume, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(heightOfMass, walberla::mpi::SUM);

      averageVelocity /= real_c(numParticles);
      heightOfMass /= particleVolume;
   }
};

std::ostream& operator<<(std::ostream& os, ParticleInfo const& m)
{
   return os << "Particle Info: uAvg = " << m.averageVelocity << ", uMax = " << m.maximumVelocity
             << ", numParticles = " << m.numParticles << ", zMax = " << m.maximumHeight << ", Vp = " << m.particleVolume
             << ", zMass = " << m.heightOfMass;
}

template< typename Accessor_T >
ParticleInfo evaluateParticleInfo(const Accessor_T& ac)
{
   static_assert(std::is_base_of< mesa_pd::data::IAccessor, Accessor_T >::value, "Provide a valid accessor");

   ParticleInfo info;
   for (uint_t i = 0; i < ac.size(); ++i)
   {
      if (isSet(ac.getFlags(i), mesa_pd::data::particle_flags::GHOST)) continue;
      if (isSet(ac.getFlags(i), mesa_pd::data::particle_flags::GLOBAL)) continue;

      ++info.numParticles;
      real_t velMagnitude   = ac.getLinearVelocity(i).length();
      real_t particleVolume = ac.getShape(i)->getVolume();
      real_t height         = ac.getPosition(i)[2];
      info.averageVelocity += velMagnitude;
      info.maximumVelocity = std::max(info.maximumVelocity, velMagnitude);
      info.maximumHeight   = std::max(info.maximumHeight, height);
      info.particleVolume += particleVolume;
      info.heightOfMass += particleVolume * height;
   }

   info.allReduce();

   return info;
}

struct FluidInfo
{
   uint_t numFluidCells   = 0;
   real_t averageVelocity = 0_r;
   real_t maximumVelocity = 0_r;
   real_t averageDensity  = 0_r;
   real_t maximumDensity  = 0_r;

   void allReduce()
   {
      walberla::mpi::allReduceInplace(numFluidCells, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(averageVelocity, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(maximumVelocity, walberla::mpi::MAX);
      ;
      walberla::mpi::allReduceInplace(averageDensity, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(maximumDensity, walberla::mpi::MAX);

      averageVelocity /= real_c(numFluidCells);
      averageDensity /= real_c(numFluidCells);
   }
};

std::ostream& operator<<(std::ostream& os, FluidInfo const& m)
{
   return os << "Fluid Info: numFluidCells = " << m.numFluidCells << ", uAvg = " << m.averageVelocity
             << ", uMax = " << m.maximumVelocity << ", densityAvg = " << m.averageDensity
             << ", densityMax = " << m.maximumDensity;
}

FluidInfo evaluateFluidInfo(const shared_ptr< StructuredBlockStorage >& blocks, const BlockDataID& densityFieldID,
                            const BlockDataID& velocityFieldID)
{
   FluidInfo info;

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto densityField  = blockIt->getData< DensityField_fluid_T >(densityFieldID);
      auto velocityField = blockIt->getData< VelocityField_fluid_T >(velocityFieldID);

      WALBERLA_FOR_ALL_CELLS_XYZ(
         densityField, ++info.numFluidCells; Vector3< real_t > velocity(
            velocityField->get(x, y, z, 0), velocityField->get(x, y, z, 1), velocityField->get(x, y, z, 2));
         real_t density = densityField->get(x, y, z); real_t velMagnitude = velocity.length();
         info.averageVelocity += velMagnitude; info.maximumVelocity = std::max(info.maximumVelocity, velMagnitude);
         info.averageDensity += density; info.maximumDensity        = std::max(info.maximumDensity, density);)
   }
   info.allReduce();
   return info;
}

//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Basic simulation of a fluidization setup
 *
 * Initially, the mono-sized sphere are created on a structured grid inside the domain.
 * The domain is either periodic or bounded by walls in the horizontal directions (x and y).
 * In z-direction, a constant inflow from below is provided
 * and a pressure boundary condition is set at the top, resembling an outflow boundary.
 *
 * The simulation is run for the given number of seconds (runtime).
 *
 * All parameters should be set via the input file.
 *
 * For the overall algorithm and the different model parameters, see
 * Rettinger, Rüde - An efficient four-way coupled lattice Boltzmann - discrete element method for
 * fully resolved simulations of particle-laden flows (2020, preprint: https://arxiv.org/abs/2003.01490)
 *
 */
//*******************************************************************************************************************
int main(int argc, char** argv)
{
   Environment env(argc, argv);

   auto cfgFile = env.config();
   if (!cfgFile) { WALBERLA_ABORT("Usage: " << argv[0] << " path-to-configuration-file \n"); }

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   walberla::gpu::selectDeviceBasedOnMpiRank();
#endif

   WALBERLA_LOG_INFO_ON_ROOT("waLBerla revision: " << std::string(WALBERLA_GIT_SHA1).substr(0, 8));
   WALBERLA_LOG_INFO_ON_ROOT("compiler flags: " << std::string(WALBERLA_COMPILER_FLAGS));
   WALBERLA_LOG_INFO_ON_ROOT("build machine: " << std::string(WALBERLA_BUILD_MACHINE));
   WALBERLA_LOG_INFO_ON_ROOT(*cfgFile);

   // read all parameters from the config file

   ///////////////////////////////////////////////////////////////
   // parameters reading and printing for thermal fluidized bed //
   ///////////////////////////////////////////////////////////////

   Config::BlockHandle physicalSetup         = cfgFile->getBlock("PhysicalSetup");
   const real_t xSize_SI                     = physicalSetup.getParameter< real_t >("xSize");
   const real_t ySize_SI                     = physicalSetup.getParameter< real_t >("ySize");
   const real_t zSize_SI                     = physicalSetup.getParameter< real_t >("zSize");
   const bool periodicInX                    = physicalSetup.getParameter< bool >("periodicInX");
   const bool periodicInY                    = physicalSetup.getParameter< bool >("periodicInY");
   const real_t runtime_SI                   = physicalSetup.getParameter< real_t >("runtime");
   const real_t uInflow_SI                   = physicalSetup.getParameter< real_t >("uInflow");
   const real_t gravitationalAcceleration_SI = physicalSetup.getParameter< real_t >("gravitationalAcceleration");
   const real_t kinematicViscosityFluid_SI   = physicalSetup.getParameter< real_t >("kinematicViscosityFluid");
   const real_t densityFluid_SI              = physicalSetup.getParameter< real_t >("densityFluid");
   const real_t particleDiameter_SI          = physicalSetup.getParameter< real_t >("particleDiameter");
   const real_t densityParticle_SI           = physicalSetup.getParameter< real_t >("densityParticle");
   const real_t dynamicFrictionCoefficient   = physicalSetup.getParameter< real_t >("dynamicFrictionCoefficient");
   const real_t coefficientOfRestitution     = physicalSetup.getParameter< real_t >("coefficientOfRestitution");
   const real_t collisionTimeFactor          = physicalSetup.getParameter< real_t >("collisionTimeFactor");
   const real_t particleGenerationSpacing_SI = physicalSetup.getParameter< real_t >("particleGenerationSpacing");

   Config::BlockHandle numericalSetup     = cfgFile->getBlock("NumericalSetup");
   const real_t dx_SI                     = numericalSetup.getParameter< real_t >("dx");
   const real_t uInflow                   = numericalSetup.getParameter< real_t >("uInflow");
   const uint_t numXBlocks                = numericalSetup.getParameter< uint_t >("numXBlocks");
   const uint_t numYBlocks                = numericalSetup.getParameter< uint_t >("numYBlocks");
   const uint_t numZBlocks                = numericalSetup.getParameter< uint_t >("numZBlocks");
   const bool useLubricationForces        = numericalSetup.getParameter< bool >("useLubricationForces");
   const uint_t numberOfParticleSubCycles = numericalSetup.getParameter< uint_t >("numberOfParticleSubCycles");
   const bool sendDirectlyFromGPU     = numericalSetup.getParameter< bool >("sendDirectlyFromGPU");
   const bool useCommunicationHiding  = numericalSetup.getParameter< bool >("useCommunicationHiding");
   const Vector3< uint_t > frameWidth = numericalSetup.getParameter< Vector3< uint_t > >("frameWidth");
   const real_t Thot           = numericalSetup.getParameter< real_t >("Thot");
   const real_t Tcold          = numericalSetup.getParameter< real_t >("Tcold");
   const real_t Pr             = numericalSetup.getParameter< real_t >("PrandtlNumber");
   const real_t Ra             = numericalSetup.getParameter< real_t >("RayleighNumber");
   const real_t Ma             = numericalSetup.getParameter< real_t >("Mach");
   const real_t relaxationRate = numericalSetup.getParameter< real_t >("relaxationRate");
   const uint_t timeSteps      = numericalSetup.getParameter< uint_t >("timeSteps");
   const Vector3< uint_t > particleSubBlockSize =
      numericalSetup.getParameter< Vector3< uint_t > >("particleSubBlockSize");
   const real_t linkedCellWidthRation = numericalSetup.getParameter< real_t >("linkedCellWidthRation");
   const bool particleBarriers        = numericalSetup.getParameter< bool >("particleBarriers");

   Config::BlockHandle outputSetup      = cfgFile->getBlock("Output");
   const real_t infoSpacing_SI          = outputSetup.getParameter< real_t >("infoSpacing");
   const real_t vtkSpacingParticles_SI  = outputSetup.getParameter< real_t >("vtkSpacingParticles");
   const real_t vtkSpacingFluid_SI      = outputSetup.getParameter< real_t >("vtkSpacingFluid");
   const std::string vtkFolder          = outputSetup.getParameter< std::string >("vtkFolder");
   const uint_t performanceLogFrequency = outputSetup.getParameter< uint_t >("performanceLogFrequency");

   // convert SI units to simulation (LBM) units and check setup

   Vector3< uint_t > domainSize(uint_c((xSize_SI / dx_SI)), uint_c((ySize_SI / dx_SI)),
                                uint_c((zSize_SI / dx_SI)));
   WALBERLA_LOG_INFO_ON_ROOT("domain size in x is " << domainSize[0]);
   WALBERLA_CHECK_FLOAT_EQUAL(real_t(domainSize[0]) * dx_SI, xSize_SI, "domain size in x is not divisible by given dx");
   WALBERLA_CHECK_FLOAT_EQUAL(real_t(domainSize[1]) * dx_SI, ySize_SI, "domain size in y is not divisible by given dx");
   WALBERLA_CHECK_FLOAT_EQUAL(real_t(domainSize[2]) * dx_SI, zSize_SI, "domain size in z is not divisible by given dx");

   Vector3< uint_t > cellsPerBlockPerDirection(domainSize[0] / numXBlocks, domainSize[1] / numYBlocks,
                                               domainSize[2] / numZBlocks);

   WALBERLA_CHECK_EQUAL(domainSize[0], cellsPerBlockPerDirection[0] * numXBlocks,
                        "number of cells in x of " << domainSize[0]
                                                   << " is not divisible by given number of blocks in x direction");
   WALBERLA_CHECK_EQUAL(domainSize[1], cellsPerBlockPerDirection[1] * numYBlocks,
                        "number of cells in y of " << domainSize[1]
                                                   << " is not divisible by given number of blocks in y direction");
   WALBERLA_CHECK_EQUAL(domainSize[2], cellsPerBlockPerDirection[2] * numZBlocks,
                        "number of cells in z of " << domainSize[2]
                                                   << " is not divisible by given number of blocks in z direction");

   WALBERLA_CHECK_GREATER(
      particleDiameter_SI / dx_SI, 5_r,
      "Your numerical resolution is below 5 cells per diameter and thus too small for such simulations!");

   const real_t densityRatio           = densityParticle_SI / densityFluid_SI;
   const real_t ReynoldsNumberParticle = uInflow_SI * particleDiameter_SI / kinematicViscosityFluid_SI;
   const real_t GalileiNumber = std::sqrt((densityRatio - 1_r) * particleDiameter_SI * gravitationalAcceleration_SI) *
                                particleDiameter_SI / kinematicViscosityFluid_SI;

   // in simulation units: dt = 1, dx = 1, densityFluid = 1

   const real_t dt_SI                     = uInflow / uInflow_SI * dx_SI;
   const real_t diameter                  = particleDiameter_SI / dx_SI;
   const real_t particleGenerationSpacing = particleGenerationSpacing_SI / dx_SI;
   const real_t viscosity                 = kinematicViscosityFluid_SI * dt_SI / (dx_SI * dx_SI);
   const real_t omega                     = lbm::collision_model::omegaFromViscosity(viscosity);
   const real_t gravitationalAcceleration = gravitationalAcceleration_SI * dt_SI * dt_SI / dx_SI;
   const real_t particleVolume            = math::pi / 6_r * diameter * diameter * diameter;

   const real_t densityFluid    = real_t(1);
   const real_t densityParticle = densityRatio;
   const real_t dx              = real_t(1);

   const uint_t numTimeSteps        = uint_c(std::ceil(runtime_SI / dt_SI));
   const uint_t infoSpacing         = uint_c(std::ceil(infoSpacing_SI / dt_SI));
   const uint_t vtkSpacingParticles = uint_c(std::ceil(vtkSpacingParticles_SI / dt_SI));
   const uint_t vtkSpacingFluid     = uint_c(std::ceil(vtkSpacingFluid_SI / dt_SI));

   const Vector3< real_t > inflowVec(0_r, 0_r, uInflow);

   const real_t poissonsRatio         = real_t(0.22);
   const real_t kappa                 = real_t(2) * (real_t(1) - poissonsRatio) / (real_t(2) - poissonsRatio);
   const real_t particleCollisionTime = collisionTimeFactor * diameter;

   WALBERLA_LOG_INFO_ON_ROOT("Simulation setup:");
   WALBERLA_LOG_INFO_ON_ROOT(" - particles: diameter = " << diameter << ", densityRatio = " << densityRatio);
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: kin. visc = " << viscosity << ", relaxation rate = " << omega);
   WALBERLA_LOG_INFO_ON_ROOT(" - grav. acceleration = " << gravitationalAcceleration);
   WALBERLA_LOG_INFO_ON_ROOT(" - Galileo number = " << GalileiNumber);
   WALBERLA_LOG_INFO_ON_ROOT(" - particle Reynolds number = " << ReynoldsNumberParticle);
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - cells per blocks per direction = " << cellsPerBlockPerDirection);
   WALBERLA_LOG_INFO_ON_ROOT(" - dx = " << dx_SI << " m");
   WALBERLA_LOG_INFO_ON_ROOT(" - dt = " << dt_SI << " s");
   WALBERLA_LOG_INFO_ON_ROOT(" - total time steps = " << numTimeSteps);
   WALBERLA_LOG_INFO_ON_ROOT(" - particle generation spacing = " << particleGenerationSpacing);
   WALBERLA_LOG_INFO_ON_ROOT(" - info spacing = " << infoSpacing);
   WALBERLA_LOG_INFO_ON_ROOT(" - vtk spacing particles = " << vtkSpacingParticles
                                                           << ", fluid slice = " << vtkSpacingFluid);

   /*WALBERLA_CHECK_EQUAL(numXBlocks * numYBlocks * numZBlocks, uint_t(MPIManager::instance()->numProcesses()),
                        "When using GPUs, the number of blocks ("
                           << numXBlocks * numYBlocks * numZBlocks << ") has to match the number of MPI processes ("
                           << uint_t(MPIManager::instance()->numProcesses()) << ")");*/


   // Necessary Calculations for temperature field

   const real_t length_conversion = (numZBlocks * cellsPerBlockPerDirection[2]);
   Vector3< uint_t > domainSizeLB(uint_c(length_conversion), uint_c(length_conversion), uint_c(length_conversion));

   const real_t uCharacteristicLB = Ma * (std::sqrt(real_c(1 / 3.0)));
   const Vector3< real_t > Uinitialize(uCharacteristicLB, 0, 0);

   const real_t delta_T              = Thot - Tcold;
   const real_t T0                   = (Thot + Tcold) / (2 * delta_T);
   const real_t kinematicViscosityLB = (Ma / std::sqrt(3)) * std::sqrt(Pr / Ra) * length_conversion;
   // std::min((Ma / std::sqrt(3)) * std::sqrt(Pr / Ra) * length_conversion, std::sqrt(3.6 / Ra));
   const real_t thermalDiffusivityLB = kinematicViscosityLB / Pr;

   const real_t time_conversion = (length_conversion * length_conversion) / (thermalDiffusivityLB);
   const real_t gravityLB       = gravitationalAcceleration_SI * time_conversion * time_conversion / length_conversion;
   const real_t alphaLB         = (uCharacteristicLB * uCharacteristicLB) / (gravityLB * delta_T * length_conversion);

   const real_t rho_0 = 1.0;
   const real_t omega_f = lbm::collision_model::omegaFromViscosity(kinematicViscosityLB);
   const real_t omega_t = lbm::collision_model::omegaFromViscosity(thermalDiffusivityLB);
   const real_t qk = 1/(4*thermalDiffusivityLB + 0.5);
   const real_t qe = 1.0;
   //  calculations for verification and correctness

   const real_t RayleighNumber =
      (Pr * alphaLB * gravityLB * delta_T * length_conversion * length_conversion * length_conversion) /
      (kinematicViscosityLB * kinematicViscosityLB);
   const real_t uchar                 = std::sqrt(Ra / Pr) * (kinematicViscosityLB / length_conversion);
   const real_t machnumber            = uchar * std::sqrt(3);
   const real_t thermal_diffusivity_2 = (Ma / std::sqrt(3)) / (std::sqrt(Ra * Pr)) * length_conversion;
   const real_t ratio                 = length_conversion / thermalDiffusivityLB;

   // Print temperature field related quantities

   WALBERLA_LOG_INFO_ON_ROOT("length conversion is " << length_conversion << " " << "time conversion is "
                                                     << time_conversion);
   WALBERLA_LOG_INFO_ON_ROOT("kinematic viscosity is " << kinematicViscosityLB);
   WALBERLA_LOG_INFO_ON_ROOT("Omega fluid is " << omega_f << "omega temperature is " << omega_t);
   WALBERLA_LOG_INFO_ON_ROOT("Thermal Diffusivity LB is " << thermalDiffusivityLB << " " << thermal_diffusivity_2);
   WALBERLA_LOG_INFO_ON_ROOT("length to diffusivity ratio is " << length_conversion / thermalDiffusivityLB);
   WALBERLA_LOG_INFO_ON_ROOT("Rayleigh number is " << RayleighNumber);
   WALBERLA_LOG_INFO_ON_ROOT("Characteristic velocity is " << uchar << " " << uCharacteristicLB);
   WALBERLA_LOG_INFO_ON_ROOT("Mach number is " << machnumber);
   WALBERLA_LOG_INFO_ON_ROOT("Domain size in lattice is " << domainSizeLB[0]);
   WALBERLA_LOG_INFO_ON_ROOT("alpha LB is " << alphaLB);
   WALBERLA_LOG_INFO_ON_ROOT("gravity LB is " << gravityLB);
   WALBERLA_LOG_INFO_ON_ROOT("Temperature difference delta_t is " << delta_T);
   //WALBERLA_LOG_INFO_ON_ROOT("Number of time steps is " << timeSteps);


   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const bool periodicInZ                     = false;
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid(
      numXBlocks, numYBlocks, numZBlocks, cellsPerBlockPerDirection[0], cellsPerBlockPerDirection[1],
      cellsPerBlockPerDirection[2], dx, 0, false, false, periodicInX, periodicInY, periodicInZ, // periodicity
      false);

   auto simulationDomain = blocks->getDomain();

   //////////////////
   // RPD COUPLING //
   //////////////////

   auto rpdDomain = std::make_shared< mesa_pd::domain::BlockForestDomain >(blocks->getBlockForestPointer());

   // init data structures
   auto ps                  = walberla::make_shared< mesa_pd::data::ParticleStorage >(1);
   auto ss                  = walberla::make_shared< mesa_pd::data::ShapeStorage >();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor            = walberla::make_shared< ParticleAccessor_T >(ps, ss);

   // prevent particles from interfering with inflow and outflow by putting the bounding planes slightly in front
   const real_t planeOffsetFromInflow  = dx;
   const real_t planeOffsetFromOutflow = dx;
   createPlaneSetup(ps, ss, simulationDomain, periodicInX, periodicInY, planeOffsetFromInflow, planeOffsetFromOutflow);

   auto sphereShape = ss->create< mesa_pd::data::Sphere >(diameter * real_t(0.5));
   ss->shapes[sphereShape]->updateMassAndInertia(densityParticle);

   // create spheres
   auto generationDomain = simulationDomain.getExtended(-particleGenerationSpacing * 0.5_r);
   for (auto pt : grid_generator::SCGrid(generationDomain, generationDomain.center(), particleGenerationSpacing))
   {
      if (rpdDomain->isContainedInProcessSubdomain(uint_c(mpi::MPIManager::instance()->rank()), pt))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pt);
         p.setInteractionRadius(diameter * real_t(0.5));
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         p.setType(1);
         p.setLinearVelocity(0.1_r * Vector3< real_t >(math::realRandom(
                                        -uInflow, uInflow))); // set small initial velocity to break symmetries
      }
   }



   //////////////////////////////////////////////////////////
   // ADD DATA TO BLOCKS and Boundary Handling for Codegen //
   /////////////////////////////////////////////////////////

   // Setting initial PDFs to nan helps to detect bugs in the initialization/BC handling
   // Depending on WALBERLA_BUILD_WITH_GPU_SUPPORT, pdfFieldCPUGPUID is either a CPU or a CPU field
   BlockDataID velFieldFluidID;
   BlockDataID densityConcentrationFieldID;

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   // Fluid PDFs on GPU
   BlockDataID pdfFieldFluidID =
      field::addToStorage< PdfField_fluid_T >(blocks, "pdf fluid field (fzyx) GPU", real_c(std::nan("")), field::fzyx);
   BlockDataID pdfFieldFluidCPUGPUID =
      gpu::addGPUFieldToStorage< PdfField_fluid_T >(blocks, pdfFieldFluidID, "pdf fluid field GPU");

   // Concentration PDFs on GPU
   BlockDataID pdfFieldConcentrationID = field::addToStorage< PdfField_concentration_T >(
      blocks, "pdf concentration field (fzyx) GPU", real_c(std::nan("")), field::fzyx);
   BlockDataID pdfFieldConcentrationCPUGPUID = gpu::addGPUFieldToStorage< PdfField_concentration_T >(
      blocks, pdfFieldConcentrationID, "pdf concentration field GPU");

   // Fluid velocity field on GPU
   velFieldFluidID =
      field::addToStorage< VelocityField_fluid_T >(blocks, "velocity fluid field GPU", real_t(0), field::fzyx);
   BlockDataID velFieldFluidCPUGPUID =
      gpu::addGPUFieldToStorage< VelocityField_fluid_T >(blocks, velFieldFluidID, "velocity fluid field GPU");

   // Concentration Density on GPU
   densityConcentrationFieldID = field::addToStorage< DensityField_concentration_T >(
      blocks, "density concentration field GPU", real_t(0), field::fzyx);
   BlockDataID densityConcentrationFieldCPUGPUID = gpu::addGPUFieldToStorage< DensityField_concentration_T >(
      blocks, densityConcentrationFieldID, "density concentration field GPU");

   // fraction field on GPU
   BlockDataID BFieldID = field::addToStorage< GhostLayerField< real_t, 1 > >(blocks, "B field GPU", 0, field::fzyx, 1);
#else

   // Fluid PDFs on CPU

   BlockDataID pdfFieldFluidCPUGPUID =
      field::addToStorage< PdfField_fluid_T >(blocks, "pdf fluid field CPU", real_c(std::nan("")), field::fzyx);

   BlockDataID velFieldFluidCPUGPUID =
      field::addToStorage< VelocityField_fluid_T >(blocks, "velocity fluid field CPU", real_t(0), field::fzyx);
   velFieldFluidID =
      field::addToStorage< VelocityField_fluid_T >(blocks, "velocity fluid field CPU", real_t(0), field::fzyx);
   // Concentration PDFs on CPU
   BlockDataID pdfFieldConcentrationCPUGPUID = field::addToStorage< PdfField_concentration_T >(
      blocks, "pdf concentration field CPU", real_c(std::nan("")), field::fzyx);

   BlockDataID densityConcentrationFieldCPUGPUID = field::addToStorage< DensityField_concentration_T >(
      blocks, "density concentration field CPU", real_t(0), field::fzyx);

#endif
   BlockDataID densityFluidFieldID =
      field::addToStorage< DensityField_fluid_T >(blocks, "density fluid field CPU", real_t(0), field::fzyx);
   //densityConcentrationFieldID = field::addToStorage< DensityField_concentration_T >(
     // blocks, "density concentration field CPU", real_t(0), field::fzyx);
   BlockDataID flagFieldFluidID = field::addFlagFieldToStorage< FlagField_T >(blocks, "fluid flag field CPU");
   BlockDataID flagFieldConcentrationID =
      field::addFlagFieldToStorage< FlagField_T >(blocks, "concentration flag field CPU");

   // Assemble boundary block string
   std::string boundariesBlockString = " BoundariesFluid";
      boundariesBlockString+=    "{";

   if (!periodicInX)
   {
      boundariesBlockString += "Border { direction W;    walldistance -1;  flag NoSlip_Fluid; }"
                               "Border { direction E;    walldistance -1;  flag NoSlip_Fluid; }";
   }

   if (!periodicInY)
   {
      boundariesBlockString += "Border { direction S;    walldistance -1;  flag NoSlip_Fluid; }"
                               "Border { direction N;    walldistance -1;  flag NoSlip_Fluid; }";
   }

   if (!periodicInZ)
   {
      boundariesBlockString += "Border { direction T;    walldistance -1;  flag Density_Fluid; }"
                               "Border { direction B;    walldistance -1;  flag Inflow_Fluid; }";
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
      boundariesBlockString += "Border { direction T;    walldistance -1;  flag Density_Concentration_west; }"
                               "Border { direction B;    walldistance -1;  flag Density_Concentration_west; }";
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
   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldFluidID, boundariesConfigFluid);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldFluidID, Fluid_Flag);
   lbm::BC_Fluid_Density density_fluid_bc(blocks,pdfFieldFluidCPUGPUID,real_t(1));
   density_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Density_Fluid_Flag, Fluid_Flag);
   lbm::BC_Fluid_NoSlip noSlip_fluid_bc(blocks, pdfFieldFluidCPUGPUID);
   noSlip_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, NoSlip_Fluid_Flag, Fluid_Flag);
   lbm::BC_Fluid_UBB ubb_fluid_bc(blocks, pdfFieldFluidCPUGPUID, real_t(0), real_t(0), uInflow);
   ubb_fluid_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldFluidID, Inflow_Fluid_Flag, Fluid_Flag);

   lbm::BC_Fluid_Outflow outflow_fluid_bc(blocks,pdfFieldFluidCPUGPUID,pdfFieldFluidID);
   outflow_fluid_bc.fillFromFlagField<FlagField_T>(blocks, flagFieldFluidID, OutFlow_Fluid_Flag, Fluid_Flag);

   // map boundaries into the concentration field simulation
   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldConcentrationID, boundariesConfigConcentration);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldConcentrationID, Concentration_Flag);
   lbm::BC_Concentration_Density density_concentration_bc_west(blocks, pdfFieldConcentrationCPUGPUID, Thot);
   density_concentration_bc_west.fillFromFlagField< FlagField_T >(blocks, flagFieldConcentrationID,
                                                                  Density_Concentration_Flag_west, Concentration_Flag);
   lbm::BC_Concentration_Density density_concentration_bc_east(blocks, pdfFieldConcentrationCPUGPUID, Tcold);
   density_concentration_bc_east.fillFromFlagField< FlagField_T >(blocks, flagFieldConcentrationID,
                                                                  Density_Concentration_Flag_east, Concentration_Flag);
   lbm::BC_Concentration_Neumann neumann_concentration_bc(blocks, pdfFieldConcentrationCPUGPUID);
   neumann_concentration_bc.fillFromFlagField< FlagField_T >(blocks, flagFieldConcentrationID,
                                                             Neumann_Concentration_Flag, Concentration_Flag);   // Assemble boundary block string

   // set up RPD functionality
   std::function< void(void) > syncCall = [&ps, &rpdDomain]() {
      // keep overlap for lubrication
      const real_t overlap = real_t(1.5);
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   };

   syncCall();
   real_t timeStepSizeRPD = real_t(1) / real_t(numberOfParticleSubCycles);
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce(timeStepSizeRPD);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(timeStepSizeRPD);
   mesa_pd::kernel::LinearSpringDashpot collisionResponse(2);
   collisionResponse.setFrictionCoefficientDynamic(0, 1, dynamicFrictionCoefficient);
   collisionResponse.setFrictionCoefficientDynamic(1, 1, dynamicFrictionCoefficient);
   real_t massSphere       = densityParticle * particleVolume;
   real_t meffSpherePlane  = massSphere;
   real_t meffSphereSphere = massSphere * massSphere / (real_t(2) * massSphere);
   collisionResponse.setStiffnessAndDamping(0, 1, coefficientOfRestitution, particleCollisionTime, kappa,
                                            meffSpherePlane);
   collisionResponse.setStiffnessAndDamping(1, 1, coefficientOfRestitution, particleCollisionTime, kappa,
                                            meffSphereSphere);
   mesa_pd::kernel::AssocToBlock assoc(blocks->getBlockForestPointer());
   mesa_pd::mpi::ReduceProperty reduceProperty;
   mesa_pd::mpi::ReduceContactHistory reduceAndSwapContactHistory;

   // set up coupling functionality
   Vector3< real_t > gravitationalForce(real_t(0), real_t(0),
                                        -(densityParticle - densityFluid) * gravitationalAcceleration * particleVolume);
   lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel(
      viscosity, [](real_t r) { return (real_t(0.001 + real_t(0.00007) * r)) * r; });

   ///////////////
   // TIME LOOP //
   ///////////////

   /////////////////////////////////////////////////////////
   // map particles and initialize PSM domain for codegen //
   ////////////////////////////////////////////////////////

   // Map particles into the fluid domain
   ParticleAndVolumeFractionSoA_T< Weighting > particleAndVolumeFractionSoA(blocks, omega);
   PSMSweepCollection psmSweepCollection(blocks, accessor, lbm_mesapd_coupling::RegularParticlesSelector(),
                                         particleAndVolumeFractionSoA, particleSubBlockSize);

   // Initialize PDFs

   pystencils::InitializeFluidDomain pdfSetterFluid(particleAndVolumeFractionSoA.BsFieldID,particleAndVolumeFractionSoA.BFieldID,densityConcentrationFieldCPUGPUID,particleAndVolumeFractionSoA.particleVelocitiesFieldID,pdfFieldFluidCPUGPUID,T0,alphaLB,gravityLB,real_t(1),real_t(0),real_t(1));
   pystencils::InitializeConcentrationDomain pdfSetterConcentration(
      densityConcentrationFieldCPUGPUID, pdfFieldConcentrationCPUGPUID, velFieldFluidCPUGPUID);


   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      psmSweepCollection.particleMappingSweep(&(*blockIt));
      psmSweepCollection.setParticleVelocitiesSweep(&(*blockIt));
      pdfSetterFluid(&(*blockIt));
      pdfSetterConcentration(&(*blockIt));

   }

   ////////////////////////////////////////
   // Communication setting for codegen //
   ///////////////////////////////////////

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

   // create the timeloop
   SweepTimeloop timeloop(blocks->getBlockStorage(), numTimeSteps);
   timeloop.addFuncBeforeTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");
   SweepTimeloop commTimeloop(blocks->getBlockStorage(), numTimeSteps);
   // objects to get the macroscopic quantities

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   pystencils::FluidMacroGetter getterSweep_fluid(BFieldID,densityConcentrationFieldID, densityFluidFieldID,
                                                  pdfFieldFluidID, velFieldFluidID, T0, alphaLB, gravityLB,
                                                  rho_0);
#else
   pystencils::FluidMacroGetter getterSweep_fluid(particleAndVolumeFractionSoA.BFieldID,densityConcentrationFieldID, densityFluidFieldID,
                                                  pdfFieldFluidCPUGPUID, velFieldFluidID, T0, alphaLB, gravityLB,
                                                  rho_0);

#endif

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   pystencils::ConcentrationMacroGetter getterSweep_concentration(densityConcentrationFieldID, pdfFieldConcentrationID);
#else
   pystencils::ConcentrationMacroGetter getterSweep_concentration(densityConcentrationFieldID,
                                                                  pdfFieldConcentrationCPUGPUID);
#endif

   // vtk output
   if (vtkSpacingParticles != uint_t(0))
   {
      // particles
      auto particleVtkOutput = make_shared< mesa_pd::vtk::ParticleVtkOutput >(ps);
      particleVtkOutput->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      particleVtkOutput->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      particleVtkOutput->addOutput< mesa_pd::data::SelectParticleInteractionRadius >("radius");
      // limit output to process-local spheres
      particleVtkOutput->setParticleSelector([sphereShape](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         return pIt->getShapeID() == sphereShape &&
                !(mesa_pd::data::particle_flags::isSet(pIt->getFlags(), mesa_pd::data::particle_flags::GHOST));
      });
      auto particleVtkWriter =
         vtk::createVTKOutput_PointData(particleVtkOutput, "particles", vtkSpacingParticles, vtkFolder);
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(particleVtkWriter), "VTK (sphere data)");
   }

   if (vtkSpacingFluid != uint_t(0))
   {
      // vtk files for fluid field
      auto vtkOutput_Fluid = vtk::createVTKOutput_BlockData(blocks, "fluid", vtkSpacingFluid, 0, false, vtkFolder);

      vtkOutput_Fluid->addBeforeFunction(communication_fluid);

      AABB fluidSliceAABB(real_t(0), real_c(domainSize[1]) * real_t(0.5) - real_t(1), real_t(0), real_c(domainSize[0]),
                     real_c(domainSize[1]) * real_t(0.5) + real_t(1), real_c(domainSize[2]));

      vtkOutput_Fluid->addBeforeFunction([&]() {
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         gpu::fieldCpy< PdfField_fluid_T, gpu::GPUField< real_t > >(blocks, pdfFieldFluidID, pdfFieldFluidCPUGPUID);
         gpu::fieldCpy< PdfField_concentration_T, gpu::GPUField< real_t > >(blocks, pdfFieldConcentrationID,
                                                                            pdfFieldConcentrationCPUGPUID);
         gpu::fieldCpy< VelocityField_fluid_T, gpu::GPUField< real_t > >(blocks, velFieldFluidID,
                                                                         velFieldFluidCPUGPUID);
         gpu::fieldCpy< DensityField_concentration_T, gpu::GPUField< real_t > >(blocks, densityConcentrationFieldID,
                                                                                densityConcentrationFieldCPUGPUID);
         gpu::fieldCpy< GhostLayerField< real_t, 1 >, BFieldGPU_T >(blocks, BFieldID,
                                                                    particleAndVolumeFractionSoA.BFieldID);
#endif
         for (auto& block : *blocks)
         {
            getterSweep_fluid(&block);
         }
      });

       //////////////////////////////
      // Write a fluid field slice //
      //////////////////////////////

      vtk::AABBCellFilter aabbSliceFilterFluid(fluidSliceAABB);

      field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldFluidID);
      fluidFilter.addFlag(Fluid_Flag);

      vtk::ChainedFilter combinedSliceFilter;
      combinedSliceFilter.addFilter(fluidFilter);
      combinedSliceFilter.addFilter(aabbSliceFilterFluid);
      vtkOutput_Fluid->addCellInclusionFilter(combinedSliceFilter);
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


// VTK files for the temperature field

      auto vtkOutput_Concentration =
         vtk::createVTKOutput_BlockData(blocks, "vtk files concentration", vtkSpacingFluid, 0, false, vtkFolder);

      vtkOutput_Concentration->addBeforeFunction(communication_concentration);

      vtkOutput_Concentration->addBeforeFunction([&]() {
         for (auto& block : *blocks)
         {
            getterSweep_concentration(&block);
         }
      });

      ///////////////////////////////////////
      // Write a concentration field slice //
      //////////////////////////////////////

      AABB concentrationSliceAABB(real_t(0), real_c(domainSize[1]) * real_t(0.5) - real_t(1), real_t(0), real_c(domainSize[0]),
                     real_c(domainSize[1]) * real_t(0.5) + real_t(1), real_c(domainSize[2]));

      vtk::AABBCellFilter aabbSliceFilterConcentration(concentrationSliceAABB);

      field::FlagFieldCellFilter< FlagField_T > ConcentrationFilter(flagFieldConcentrationID);
      ConcentrationFilter.addFlag(Concentration_Flag);


      combinedSliceFilter.addFilter(ConcentrationFilter);
      combinedSliceFilter.addFilter(aabbSliceFilterConcentration);
      vtkOutput_Concentration->addCellInclusionFilter(combinedSliceFilter);

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      vtkOutput_Concentration->addCellDataWriter(
         make_shared< field::VTKWriter< DensityField_concentration_T > >(densityConcentrationFieldID, "Concentration"));
#else
      vtkOutput_Concentration->addCellDataWriter(
         make_shared< field::VTKWriter< DensityField_concentration_T > >(densityConcentrationFieldCPUGPUID, "Concentration"));
#endif

      vtkOutput_Concentration->addCellDataWriter(
         make_shared< field::VTKWriter< FlagField_T > >(flagFieldConcentrationID, "ConcentrationFlagField"));

      // fraction mapping field, only a slice
      auto fractionFieldVTK =
         vtk::createVTKOutput_BlockData(blocks, "fraction_field", vtkSpacingFluid, 0, false, vtkFolder);

      fractionFieldVTK->addCellInclusionFilter(combinedSliceFilter);

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      fractionFieldVTK->addCellDataWriter(
         make_shared< field::VTKWriter< BField_T > >(BFieldID, "OverlapFraction"));
#else
      fractionFieldVTK->addCellDataWriter(
         make_shared< field::VTKWriter< BField_T > >(particleAndVolumeFractionSoA.BFieldID, "OverlapFraction"));
#endif


      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(fractionFieldVTK), "VTK (fraction field data");
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput_Fluid), "VTK output Fluid");
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput_Concentration), "VTK output Concentration");
   }

   if (vtkSpacingFluid != uint_t(0) || vtkSpacingParticles != uint_t(0))
   {
      vtk::writeDomainDecomposition(blocks, "domain_decomposition", vtkFolder);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////////
   // add LBM communication, boundary handling and the LBM sweeps to the time loop  for codegen //
   //////////////////////////////////////////////////////////////////////////////////////////////
   pystencils::PSMFluidSweep psmFluidSweep(
      particleAndVolumeFractionSoA.BsFieldID, particleAndVolumeFractionSoA.BFieldID, densityConcentrationFieldCPUGPUID,
      particleAndVolumeFractionSoA.particleForcesFieldID, particleAndVolumeFractionSoA.particleVelocitiesFieldID,
      pdfFieldFluidCPUGPUID, velFieldFluidCPUGPUID, T0, alphaLB, gravityLB, omega_f, rho_0);

   pystencils::PSMFluidSweepSplit psmFluidSplitSweep(
      particleAndVolumeFractionSoA.BsFieldID, particleAndVolumeFractionSoA.BFieldID, densityConcentrationFieldCPUGPUID,particleAndVolumeFractionSoA.particleForcesFieldID,
      particleAndVolumeFractionSoA.particleVelocitiesFieldID,
      pdfFieldFluidCPUGPUID,velFieldFluidCPUGPUID,T0,alphaLB,gravityLB,omega_f,rho_0,frameWidth);

   pystencils::LBMConcentrationSweep lbmConcentrationSweep(densityConcentrationFieldCPUGPUID, pdfFieldConcentrationCPUGPUID,
                                                           velFieldFluidCPUGPUID,qe,qk);

   pystencils::LBMConcentrationSplitSweep lbmConcentrationSplitSweep(densityConcentrationFieldCPUGPUID, pdfFieldConcentrationCPUGPUID,
                                                                     velFieldFluidCPUGPUID,qe,qk,frameWidth);

   if(useCommunicationHiding){
      timeloop.add() << Sweep(deviceSyncWrapper(density_fluid_bc.getSweep()), "Boundary Handling (Outflow Fluid)");
      timeloop.add() << Sweep(deviceSyncWrapper(ubb_fluid_bc.getSweep()),"Boundary Handling (UBB inflow fluid)");
      timeloop.add() << Sweep(deviceSyncWrapper(neumann_concentration_bc.getSweep()),
               "Boundary Handling (Concentration Neumann)");
      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc_west.getSweep()),
                              "Boundary Handling (Concentration Density west)");
      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc_east.getSweep()),
                              "Boundary Handling (Concentration Density east)");

      commTimeloop.add() << BeforeFunction([&]() { com_fluid.startCommunication(); })
                         << Sweep(deviceSyncWrapper(psmSweepCollection.particleMappingSweep), "Particle mapping");
      commTimeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.setParticleVelocitiesSweep),
                                  "Set particle velocities");
      commTimeloop.add() << Sweep(deviceSyncWrapper(psmFluidSplitSweep.getInnerSweep()), "PSM inner sweep")
                         << AfterFunction([&]() { com_fluid.wait(); }, "LBM Communication (wait)");
      timeloop.add() << Sweep(deviceSyncWrapper(psmFluidSplitSweep.getOuterSweep()), "PSM outer sweep");

      commTimeloop.add() << BeforeFunction([&]() { com_concentration.startCommunication(); }, "LBM concentration Communication (start)")
                         << Sweep(deviceSyncWrapper(lbmConcentrationSplitSweep.getInnerSweep()), "LBM concentration inner sweep")
                         << AfterFunction([&]() { com_concentration.wait(); }, "LBM concentration Communication (wait)");
      timeloop.add() << Sweep(deviceSyncWrapper(lbmConcentrationSplitSweep.getOuterSweep()), "LBM concentration outer sweep");

      // after both the sweeps, reduce the particle forces.
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.reduceParticleForcesSweep),
                              "Reduce particle forces");
   }

   else{
      timeloop.add() << BeforeFunction(communication_fluid, "LBM fluid Communication")
                     << Sweep(deviceSyncWrapper(density_fluid_bc.getSweep()), "Boundary Handling (Outflow Fluid)");
      timeloop.add() << Sweep(deviceSyncWrapper(ubb_fluid_bc.getSweep()),"Boundary Handling (UBB inflow fluid)");

      timeloop.add() << BeforeFunction(communication_concentration, "LBM concentration Communication")
                     << Sweep(deviceSyncWrapper(neumann_concentration_bc.getSweep()),
                              "Boundary Handling (Concentration Neumann)");
      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc_west.getSweep()),
                              "Boundary Handling (Concentration Density west)");
      timeloop.add() << Sweep(deviceSyncWrapper(density_concentration_bc_east.getSweep()),
                              "Boundary Handling (Concentration Density east)");

      timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.particleMappingSweep), "Particle mapping");
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.setParticleVelocitiesSweep),
                              "Set particle velocities");
      timeloop.add() << Sweep(deviceSyncWrapper(psmFluidSweep), "PSM Fluid sweep");
      timeloop.add() << Sweep(deviceSyncWrapper(lbmConcentrationSweep), "LBM Concentration sweep");

      // after both the sweeps, reduce the particle forces.
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.reduceParticleForcesSweep),
                              "Reduce particle forces");
   }


   // Add performance logging
   lbm::PerformanceLogger< FlagField_T > performanceLogger(blocks, flagFieldFluidID, Fluid_Flag, performanceLogFrequency);
   if (performanceLogFrequency > 0)
   {
      timeloop.addFuncAfterTimeStep(performanceLogger, "Evaluate performance logging");
   }
   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;
   const bool useOpenMP = true;

   real_t linkedCellWidth = linkedCellWidthRation * diameter;
   mesa_pd::data::LinkedCells linkedCells(rpdDomain->getUnionOfLocalAABBs().getExtended(linkedCellWidth),
                                          linkedCellWidth);
   mesa_pd::kernel::InsertParticleIntoLinkedCells ipilc;

   // time loop
   for (uint_t timeStep = 0; timeStep < numTimeSteps; ++timeStep)
   {
      // perform a single simulation step -> this contains LBM and setting of the hydrodynamic interactions
      if(useCommunicationHiding) { commTimeloop.singleStep(timeloopTiming); }
      timeloop.singleStep(timeloopTiming);

      if (particleBarriers) WALBERLA_MPI_BARRIER();
      timeloopTiming["RPD forEachParticle assoc"].start();
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, assoc, *accessor);
      if (particleBarriers) WALBERLA_MPI_BARRIER();
      timeloopTiming["RPD forEachParticle assoc"].end();
      timeloopTiming["RPD reduceProperty HydrodynamicForceTorqueNotification"].start();
      reduceProperty.operator()< mesa_pd::HydrodynamicForceTorqueNotification >(*ps);
      if (particleBarriers) WALBERLA_MPI_BARRIER();
      timeloopTiming["RPD reduceProperty HydrodynamicForceTorqueNotification"].end();

      if (timeStep == 0)
      {
         lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel
            initializeHydrodynamicForceTorqueForAveragingKernel;
         timeloopTiming["RPD forEachParticle initializeHydrodynamicForceTorqueForAveragingKernel"].start();
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
                             initializeHydrodynamicForceTorqueForAveragingKernel, *accessor);
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD forEachParticle initializeHydrodynamicForceTorqueForAveragingKernel"].end();
      }
      timeloopTiming["RPD forEachParticle averageHydrodynamicForceTorque"].start();
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque,
                          *accessor);
      if (particleBarriers) WALBERLA_MPI_BARRIER();
      timeloopTiming["RPD forEachParticle averageHydrodynamicForceTorque"].end();

      for (auto subCycle = uint_t(0); subCycle < numberOfParticleSubCycles; ++subCycle)
      {
         timeloopTiming["RPD forEachParticle vvIntegratorPreForce"].start();
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD forEachParticle vvIntegratorPreForce"].end();
         timeloopTiming["RPD syncCall"].start();
         syncCall();
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD syncCall"].end();

         timeloopTiming["RPD linkedCells.clear"].start();
         linkedCells.clear();
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD linkedCells.clear"].end();
         timeloopTiming["RPD forEachParticle ipilc"].start();
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, ipilc, *accessor, linkedCells);
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD forEachParticle ipilc"].end();

         if (useLubricationForces)
         {
            // lubrication correction
            timeloopTiming["RPD forEachParticlePairHalf lubricationCorrectionKernel"].start();
            linkedCells.forEachParticlePairHalf(
               useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
               [&lubricationCorrectionKernel, &rpdDomain](const size_t idx1, const size_t idx2, auto& ac) {
                  mesa_pd::collision_detection::AnalyticContactDetection acd;
                  acd.getContactThreshold() = lubricationCorrectionKernel.getNormalCutOffDistance();
                  mesa_pd::kernel::DoubleCast double_cast;
                  mesa_pd::mpi::ContactFilter contact_filter;
                  if (double_cast(idx1, idx2, ac, acd, ac))
                  {
                     if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                     {
                        double_cast(acd.getIdx1(), acd.getIdx2(), ac, lubricationCorrectionKernel, ac,
                                    acd.getContactNormal(), acd.getPenetrationDepth());
                     }
                  }
               },
               *accessor);
            if (particleBarriers) WALBERLA_MPI_BARRIER();
            timeloopTiming["RPD forEachParticlePairHalf lubricationCorrectionKernel"].end();
         }

         // collision response
         timeloopTiming["RPD forEachParticlePairHalf collisionResponse"].start();
         linkedCells.forEachParticlePairHalf(
            useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
            [&collisionResponse, &rpdDomain, timeStepSizeRPD](const size_t idx1, const size_t idx2, auto& ac) {
               mesa_pd::collision_detection::AnalyticContactDetection acd;
               mesa_pd::kernel::DoubleCast double_cast;
               mesa_pd::mpi::ContactFilter contact_filter;
               if (double_cast(idx1, idx2, ac, acd, ac))
               {
                  if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                  {
                     collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(),
                                       acd.getPenetrationDepth(), timeStepSizeRPD);
                  }
               }
            },
            *accessor);
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD forEachParticlePairHalf collisionResponse"].end();

         timeloopTiming["RPD reduceProperty reduceAndSwapContactHistory"].start();
         reduceAndSwapContactHistory(*ps);
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD reduceProperty reduceAndSwapContactHistory"].end();

         // add hydrodynamic force
         lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
         timeloopTiming["RPD forEachParticle addHydrodynamicInteraction + addGravitationalForce"].start();
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction,
                             *accessor);

         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor);
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD forEachParticle addHydrodynamicInteraction + addGravitationalForce"].end();

         timeloopTiming["RPD reduceProperty ForceTorqueNotification"].start();
         reduceProperty.operator()< mesa_pd::ForceTorqueNotification >(*ps);
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD reduceProperty ForceTorqueNotification"].end();

         timeloopTiming["RPD forEachParticle vvIntegratorPostForce"].start();
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
         if (particleBarriers) WALBERLA_MPI_BARRIER();
         timeloopTiming["RPD forEachParticle vvIntegratorPostForce"].end();
      }

      timeloopTiming["RPD syncCall"].start();
      syncCall();
      if (particleBarriers) WALBERLA_MPI_BARRIER();
      timeloopTiming["RPD syncCall"].end();

      timeloopTiming["RPD forEachParticle resetHydrodynamicForceTorque"].start();
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor);
      if (particleBarriers) WALBERLA_MPI_BARRIER();
      timeloopTiming["RPD forEachParticle resetHydrodynamicForceTorque"].end();

      if (timeStep % infoSpacing == 0)
      {
         timeloopTiming["Evaluate infos"].start();

         auto particleInfo = evaluateParticleInfo(*accessor);
         WALBERLA_LOG_INFO_ON_ROOT(particleInfo);


         ////////////////////////////////////////////
         // evaluation of fluid info using codegen //
         ///////////////////////////////////////////
         //auto fluidInfo = evaluateFluidInfoCodegen(blocks,velFieldFluidCPUGPUID,densityFluidFieldID);
         gpu::fieldCpy< PdfField_fluid_T, gpu::GPUField< real_t > >(blocks, pdfFieldFluidID, pdfFieldFluidCPUGPUID);
         gpu::fieldCpy< PdfField_concentration_T, gpu::GPUField< real_t > >(blocks, pdfFieldConcentrationID,
                                                                            pdfFieldConcentrationCPUGPUID);
         gpu::fieldCpy< VelocityField_fluid_T, gpu::GPUField< real_t > >(blocks, velFieldFluidID,
                                                                         velFieldFluidCPUGPUID);
         gpu::fieldCpy< DensityField_concentration_T, gpu::GPUField< real_t > >(blocks, densityConcentrationFieldID,
                                                                                densityConcentrationFieldCPUGPUID);
         gpu::fieldCpy< GhostLayerField< real_t, 1 >, BFieldGPU_T >(blocks, BFieldID,
                                                                    particleAndVolumeFractionSoA.BFieldID);
         for (auto& block : *blocks)
         {
            getterSweep_fluid(&block);
         }
         auto fluidInfo = evaluateFluidInfo(blocks, densityFluidFieldID, velFieldFluidID);
         WALBERLA_LOG_INFO_ON_ROOT(fluidInfo);

         timeloopTiming["Evaluate infos"].end();
      }
   }

   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace thermalfluidized_bed

int main(int argc, char** argv) { thermalfluidized_bed::main(argc, argv); }
