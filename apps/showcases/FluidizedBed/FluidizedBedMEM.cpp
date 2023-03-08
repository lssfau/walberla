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
//! \file FluidizedBedMEM.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
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
#include "field/communication/PackInfo.h"
#include "field/vtk/all.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/all.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/DoubleCast.h"
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

namespace fluidized_bed
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

using LatticeModel_T = lbm::D3Q19< lbm::collision_model::SRT >;

using Stencil_T  = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField< LatticeModel_T >;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

using ScalarField_T = GhostLayerField< real_t, 1 >;

const uint_t FieldGhostLayers = 1;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag("fluid");
const FlagUID NoSlip_Flag("no slip");
const FlagUID MO_Flag("moving obstacle");
const FlagUID FormerMO_Flag("former moving obstacle");
const FlagUID Inflow_Flag("inflow");
const FlagUID Outflow_Flag("outflow");

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
template< typename ParticleAccessor_T >
class MyBoundaryHandling
{
 public:
   using NoSlip_T  = lbm::NoSlip< LatticeModel_T, flag_t >;
   using MO_T      = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using Inflow_T  = lbm::SimpleUBB< LatticeModel_T, flag_t >;
   using Outflow_T = lbm::SimplePressure< LatticeModel_T, flag_t >;
   using Type      = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, MO_T, Inflow_T, Outflow_T >;

   MyBoundaryHandling(const BlockDataID& flagFieldID, const BlockDataID& pdfFieldID, const BlockDataID& particleFieldID,
                      const shared_ptr< ParticleAccessor_T >& ac, Vector3< real_t > inflowVelocity, bool periodicInX,
                      bool periodicInY)
      : flagFieldID_(flagFieldID), pdfFieldID_(pdfFieldID), particleFieldID_(particleFieldID), ac_(ac),
        inflowVelocity_(inflowVelocity), periodicInX_(periodicInX), periodicInY_(periodicInY)
   {}

   Type* operator()(IBlock* const block, const StructuredBlockStorage* const storage) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block);
      WALBERLA_ASSERT_NOT_NULLPTR(storage);

      auto* flagField     = block->getData< FlagField_T >(flagFieldID_);
      auto* pdfField      = block->getData< PdfField_T >(pdfFieldID_);
      auto* particleField = block->getData< lbm_mesapd_coupling::ParticleField_T >(particleFieldID_);

      const auto fluid =
         flagField->flagExists(Fluid_Flag) ? flagField->getFlag(Fluid_Flag) : flagField->registerFlag(Fluid_Flag);

      Type* handling =
         new Type("moving obstacle boundary handling", flagField, fluid, NoSlip_T("NoSlip", NoSlip_Flag, pdfField),
                  MO_T("MO", MO_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block),
                  Inflow_T("Inflow", Inflow_Flag, pdfField, inflowVelocity_),
                  Outflow_T("Outflow", Outflow_Flag, pdfField, real_t(1)));

      const auto inflow  = flagField->getFlag(Inflow_Flag);
      const auto outflow = flagField->getFlag(Outflow_Flag);
      const auto noslip  = flagField->getFlag(NoSlip_Flag);

      CellInterval domainBB = storage->getDomainCellBB();

      domainBB.xMin() -= cell_idx_c(FieldGhostLayers);
      domainBB.xMax() += cell_idx_c(FieldGhostLayers);

      domainBB.zMin() -= cell_idx_c(FieldGhostLayers);
      domainBB.zMax() += cell_idx_c(FieldGhostLayers);

      domainBB.yMin() -= cell_idx_c(FieldGhostLayers);
      domainBB.yMax() += cell_idx_c(FieldGhostLayers);

      // inflow at bottom
      CellInterval bottom(domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(),
                          domainBB.zMin());
      storage->transformGlobalToBlockLocalCellInterval(bottom, *block);
      handling->forceBoundary(inflow, bottom);

      // outflow at top
      CellInterval top(domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(),
                       domainBB.zMax());
      storage->transformGlobalToBlockLocalCellInterval(top, *block);
      handling->forceBoundary(outflow, top);

      if (!periodicInX_)
      {
         // left
         CellInterval left(domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(),
                           domainBB.zMax());
         storage->transformGlobalToBlockLocalCellInterval(left, *block);
         handling->forceBoundary(noslip, left);

         // right
         CellInterval right(domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(),
                            domainBB.zMax());
         storage->transformGlobalToBlockLocalCellInterval(right, *block);
         handling->forceBoundary(noslip, right);
      }

      if (!periodicInY_)
      {
         // front
         CellInterval front(domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(),
                            domainBB.zMax());
         storage->transformGlobalToBlockLocalCellInterval(front, *block);
         handling->forceBoundary(noslip, front);

         // back
         CellInterval back(domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(),
                           domainBB.zMax());
         storage->transformGlobalToBlockLocalCellInterval(back, *block);
         handling->forceBoundary(noslip, back);
      }

      handling->fillWithDomain(FieldGhostLayers);

      return handling;
   }

 private:
   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr< ParticleAccessor_T > ac_;

   Vector3< real_t > inflowVelocity_;
   bool periodicInX_;
   bool periodicInY_;
};
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

template< typename BoundaryHandling_T >
FluidInfo evaluateFluidInfo(const shared_ptr< StructuredBlockStorage >& blocks, const BlockDataID& pdfFieldID,
                            const BlockDataID& boundaryHandlingID)
{
   FluidInfo info;

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto pdfField         = blockIt->getData< PdfField_T >(pdfFieldID);
      auto boundaryHandling = blockIt->getData< BoundaryHandling_T >(boundaryHandlingID);

      WALBERLA_FOR_ALL_CELLS_XYZ(
         pdfField, if (!boundaryHandling->isDomain(x, y, z)) continue; ++info.numFluidCells;
         Vector3< real_t > velocity(0_r); real_t density = pdfField->getDensityAndVelocity(velocity, x, y, z);
         real_t velMagnitude                             = velocity.length(); info.averageVelocity += velMagnitude;
         info.maximumVelocity = std::max(info.maximumVelocity, velMagnitude); info.averageDensity += density;
         info.maximumDensity  = std::max(info.maximumDensity, density);)
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

   WALBERLA_LOG_INFO_ON_ROOT("waLBerla revision: " << std::string(WALBERLA_GIT_SHA1).substr(0, 8));
   WALBERLA_LOG_INFO_ON_ROOT(*cfgFile);

   // read all parameters from the config file

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
   const real_t particleGenerationSpacing_SI = physicalSetup.getParameter< real_t >("particleGenerationSpacing");

   Config::BlockHandle numericalSetup     = cfgFile->getBlock("NumericalSetup");
   const real_t dx_SI                     = numericalSetup.getParameter< real_t >("dx");
   const real_t uInflow                   = numericalSetup.getParameter< real_t >("uInflow");
   const uint_t numXBlocks                = numericalSetup.getParameter< uint_t >("numXBlocks");
   const uint_t numYBlocks                = numericalSetup.getParameter< uint_t >("numYBlocks");
   const uint_t numZBlocks                = numericalSetup.getParameter< uint_t >("numZBlocks");
   const bool useLubricationForces        = numericalSetup.getParameter< bool >("useLubricationForces");
   const uint_t numberOfParticleSubCycles = numericalSetup.getParameter< uint_t >("numberOfParticleSubCycles");

   Config::BlockHandle outputSetup     = cfgFile->getBlock("Output");
   const real_t infoSpacing_SI         = outputSetup.getParameter< real_t >("infoSpacing");
   const real_t vtkSpacingParticles_SI = outputSetup.getParameter< real_t >("vtkSpacingParticles");
   const real_t vtkSpacingFluid_SI     = outputSetup.getParameter< real_t >("vtkSpacingFluid");
   const std::string vtkFolder         = outputSetup.getParameter< std::string >("vtkFolder");

   // convert SI units to simulation (LBM) units and check setup

   Vector3< uint_t > domainSize(uint_c(std::ceil(xSize_SI / dx_SI)), uint_c(std::ceil(ySize_SI / dx_SI)),
                                uint_c(std::ceil(zSize_SI / dx_SI)));
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
   const real_t particleCollisionTime = 4_r * diameter;

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
         p.setType(0);
         p.setLinearVelocity(0.1_r * Vector3< real_t >(math::realRandom(
                                        -uInflow, uInflow))); // set small initial velocity to break symmetries
      }
   }

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   LatticeModel_T latticeModel = LatticeModel_T(omega);

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >(blocks, "pdf field", latticeModel, inflowVec,
                                                                        densityFluid, uint_t(1), field::zyxf);

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   // add particle field
   BlockDataID particleFieldID = field::addToStorage< lbm_mesapd_coupling::ParticleField_T >(
      blocks, "particle field", accessor->getInvalidUid(), field::zyxf, FieldGhostLayers);

   // add boundary handling
   using BoundaryHandling_T       = MyBoundaryHandling< ParticleAccessor_T >::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
      MyBoundaryHandling< ParticleAccessor_T >(flagFieldID, pdfFieldID, particleFieldID, accessor, inflowVec,
                                               periodicInX, periodicInY),
      "boundary handling");

   // set up RPD functionality
   std::function< void(void) > syncCall = [&ps, &rpdDomain]() {
      const real_t overlap = real_t(1.5);
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   };

   syncCall();

   real_t timeStepSizeRPD = real_t(1) / real_t(numberOfParticleSubCycles);
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce(timeStepSizeRPD);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(timeStepSizeRPD);
   mesa_pd::kernel::LinearSpringDashpot collisionResponse(1);
   collisionResponse.setFrictionCoefficientDynamic(0, 0, dynamicFrictionCoefficient);
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
   lbm_mesapd_coupling::ParticleMappingKernel< BoundaryHandling_T > particleMappingKernel(blocks, boundaryHandlingID);
   lbm_mesapd_coupling::MovingParticleMappingKernel< BoundaryHandling_T > movingParticleMappingKernel(
      blocks, boundaryHandlingID, particleFieldID);

   ///////////////
   // TIME LOOP //
   ///////////////

   // map particles into the LBM simulation
   // note: planes are not mapped and are thus only visible to the particles, not to the fluid
   // instead, the respective boundary conditions for the fluid are explicitly set, see the boundary handling
   ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, movingParticleMappingKernel,
                       *accessor, MO_Flag);

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   blockforest::communication::UniformBufferedScheme< Stencil_T > optimizedPDFCommunicationScheme(blocks);
   optimizedPDFCommunicationScheme.addPackInfo(
      make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >(pdfFieldID)); // optimized sync

   blockforest::communication::UniformBufferedScheme< Stencil_T > fullPDFCommunicationScheme(blocks);
   fullPDFCommunicationScheme.addPackInfo(
      make_shared< field::communication::PackInfo< PdfField_T > >(pdfFieldID)); // full sync

   // create the timeloop

   SweepTimeloop timeloop(blocks->getBlockStorage(), numTimeSteps);

   timeloop.addFuncBeforeTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");

   // vtk output
   if (vtkSpacingParticles != uint_t(0))
   {
      // sphere
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
      // velocity field, only a slice
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData(blocks, "fluid", vtkSpacingFluid, 0, false, vtkFolder);

      pdfFieldVTK->addBeforeFunction(fullPDFCommunicationScheme);

      AABB sliceAABB(real_t(0), real_c(domainSize[1]) * real_t(0.5) - real_t(1), real_t(0), real_c(domainSize[0]),
                     real_c(domainSize[1]) * real_t(0.5) + real_t(1), real_c(domainSize[2]));
      vtk::AABBCellFilter aabbSliceFilter(sliceAABB);

      field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldID);
      fluidFilter.addFlag(Fluid_Flag);

      vtk::ChainedFilter combinedSliceFilter;
      combinedSliceFilter.addFilter(fluidFilter);
      combinedSliceFilter.addFilter(aabbSliceFilter);

      pdfFieldVTK->addCellInclusionFilter(combinedSliceFilter);

      pdfFieldVTK->addCellDataWriter(
         make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "Velocity"));
      pdfFieldVTK->addCellDataWriter(
         make_shared< lbm::DensityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "Density"));

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(pdfFieldVTK), "VTK (fluid field data)");
   }

   if (vtkSpacingFluid != uint_t(0) || vtkSpacingParticles != uint_t(0))
   {
      vtk::writeDomainDecomposition(blocks, "domain_decomposition", vtkFolder);
   }

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip
   // treatment)
   timeloop.add() << BeforeFunction(optimizedPDFCommunicationScheme, "LBM Communication")
                  << Sweep(BoundaryHandling_T::getBlockSweep(boundaryHandlingID), "Boundary Handling");

   // stream + collide LBM step
   auto lbmSweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >(pdfFieldID, flagFieldID, Fluid_Flag);
   timeloop.add() << Sweep(makeSharedSweep(lbmSweep), "LBM stream / collide");

   // this is carried out after the particle integration, it corrects the flag field and restores missing PDF
   // information then, the checkpointing file can be written, as otherwise some cells are invalid and can not be
   // recovered
   SweepTimeloop timeloopAfterParticles(blocks->getBlockStorage(), numTimeSteps);

   // sweep for updating the particle mapping into the LBM simulation
   bool strictlyConserveMomentum = false;
   timeloopAfterParticles.add() << Sweep(
      lbm_mesapd_coupling::makeMovingParticleMapping< PdfField_T, BoundaryHandling_T >(
         blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag,
         lbm_mesapd_coupling::RegularParticlesSelector(), strictlyConserveMomentum),
      "Particle Mapping");

   // sweep for restoring PDFs in cells previously occupied by particles
   bool reconstruction_recomputeTargetDensity = false;
   bool reconstruction_useCentralDifferences  = true;
   auto gradReconstructor = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor< BoundaryHandling_T >(
      blocks, boundaryHandlingID, omega, reconstruction_recomputeTargetDensity, reconstruction_useCentralDifferences);
   timeloopAfterParticles.add()
      << BeforeFunction(fullPDFCommunicationScheme, "PDF Communication")
      << Sweep(makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager< PdfField_T, BoundaryHandling_T >(
                  blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag,
                  gradReconstructor, strictlyConserveMomentum)),
               "PDF Restore");

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;
   const bool useOpenMP = false;

   // time loop
   for (uint_t timeStep = 0; timeStep < numTimeSteps; ++timeStep)
   {
      // perform a single simulation step -> this contains LBM and setting of the hydrodynamic interactions
      timeloop.singleStep(timeloopTiming);

      reduceProperty.operator()< mesa_pd::HydrodynamicForceTorqueNotification >(*ps);

      if (timeStep == 0)
      {
         lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel
            initializeHydrodynamicForceTorqueForAveragingKernel;
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
                             initializeHydrodynamicForceTorqueForAveragingKernel, *accessor);
      }
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque,
                          *accessor);

      for (auto subCycle = uint_t(0); subCycle < numberOfParticleSubCycles; ++subCycle)
      {
         timeloopTiming["RPD"].start();

         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
         syncCall();

         if (useLubricationForces)
         {
            // lubrication correction
            ps->forEachParticlePairHalf(
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
         }

         // collision response
         ps->forEachParticlePairHalf(
            useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
            [&collisionResponse, &rpdDomain, timeStepSizeRPD, coefficientOfRestitution, particleCollisionTime,
             kappa](const size_t idx1, const size_t idx2, auto& ac) {
               mesa_pd::collision_detection::AnalyticContactDetection acd;
               mesa_pd::kernel::DoubleCast double_cast;
               mesa_pd::mpi::ContactFilter contact_filter;
               if (double_cast(idx1, idx2, ac, acd, ac))
               {
                  if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                  {
                     auto meff = real_t(1) / (ac.getInvMass(idx1) + ac.getInvMass(idx2));
                     collisionResponse.setStiffnessAndDamping(0, 0, coefficientOfRestitution, particleCollisionTime,
                                                              kappa, meff);
                     collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(),
                                       acd.getPenetrationDepth(), timeStepSizeRPD);
                  }
               }
            },
            *accessor);

         reduceAndSwapContactHistory(*ps);

         // add hydrodynamic force
         lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction,
                             *accessor);

         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor);

         reduceProperty.operator()< mesa_pd::ForceTorqueNotification >(*ps);

         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
         syncCall();

         timeloopTiming["RPD"].end();
      }

      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor);

      // update particle mapping
      timeloopAfterParticles.singleStep(timeloopTiming);

      if (timeStep % infoSpacing == 0)
      {
         timeloopTiming["Evaluate infos"].start();

         auto particleInfo = evaluateParticleInfo(*accessor);
         WALBERLA_LOG_INFO_ON_ROOT(particleInfo);

         auto fluidInfo = evaluateFluidInfo< BoundaryHandling_T >(blocks, pdfFieldID, boundaryHandlingID);
         WALBERLA_LOG_INFO_ON_ROOT(fluidInfo);

         timeloopTiming["Evaluate infos"].end();
      }
   }

   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace fluidized_bed

int main(int argc, char** argv) { fluidized_bed::main(argc, argv); }
