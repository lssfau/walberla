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
//! \file AntiDunes.cpp
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Jonas Plewinski <jonas.plewinski@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
// This showcase simulates antidunes, i.e., particulate dunes that travel in opposite stream-wise direction.
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/StabilityChecker.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/BlockStateDetectorSweep.h"
#include "lbm/free_surface/SurfaceMeshWriter.h"
#include "lbm/free_surface/TotalMassComputer.h"
#include "lbm/free_surface/VtkWriter.h"
#include "lbm/free_surface/dynamics/CellConversionSweep.h"
#include "lbm/free_surface/dynamics/ConversionFlagsResetSweep.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionSweep.h"
#include "lbm/free_surface/dynamics/PdfRefillingSweep.h"
#include "lbm/free_surface/dynamics/StreamReconstructAdvectSweep.h"
#include "lbm/free_surface/surface_geometry/CurvatureSweep.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"
#include "lbm/free_surface/surface_geometry/SmoothingSweep.h"
#include "lbm/free_surface/surface_geometry/Utility.h"

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/ParticleAccessorWithBaseShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ClearGhostOwnerSync.h"
#include "mesa_pd/mpi/ClearNextNeighborSync.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/SyncGhostOwners.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "vtk/VTKOutput.h"

#include <core/waLBerlaBuildInfo.h>

#include "AntidunesBoundaryHandling.h"
#include "AntidunesLatticeModel.h"
#include "PIDController.h"
#include "Utility.h"

namespace walberla
{
namespace antidunes
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using LatticeModel_T        = lbm::AntidunesLatticeModel;
using LatticeModelStencil_T = LatticeModel_T::Stencil;
using PdfField_T            = lbm::PdfField< LatticeModel_T >;
using PdfCommunication_T    = blockforest::SimpleCommunication< LatticeModelStencil_T >;

// the geometry computations in SurfaceGeometryHandler require meaningful values in the ghost layers in corner
// directions (flag field and fill level field); this holds, even if the lattice model uses a D3Q19 stencil
using CommunicationStencil_T =
   typename std::conditional< LatticeModel_T::Stencil::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;
using Communication_T = blockforest::SimpleCommunication< CommunicationStencil_T >;

using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithBaseShape;

using flag_t      = uint32_t;
using FlagField_T = FlagField< flag_t >;
using AntidunesBoundaryHandling_T =
   free_surface::AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >;

using StateSweep = walberla::free_surface::BlockStateDetectorSweep< FlagField_T >;

const FlagUID FormerMO_Flag("former moving obstacle");

// empty sweep required for using selectors (e.g. StateSweep::fullFreeSurface)
struct emptySweep
{
   void operator()(IBlock*) {}
};

// data handling for loading a field of type ScalarField_T from file
template< typename ScalarField_T >
class ScalarFieldHandling : public field::BlockDataHandling< ScalarField_T >
{
 public:
   ScalarFieldHandling(const weak_ptr< StructuredBlockStorage >& blocks, uint_t numberGhostLayer)
      : blocks_(blocks), numberGhostLayer_(numberGhostLayer)
   {}

 protected:
   ScalarField_T* allocate(IBlock* const block) override { return allocateDispatch(block); }

   ScalarField_T* reallocate(IBlock* const block) override { return allocateDispatch(block); }

 private:
   weak_ptr< StructuredBlockStorage > blocks_;
   uint_t numberGhostLayer_;

   ScalarField_T* allocateDispatch(IBlock* const block)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block);

      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blocks);

      return new ScalarField_T(blocks->getNumberOfXCells(*block), blocks->getNumberOfYCells(*block),
                               blocks->getNumberOfZCells(*block), numberGhostLayer_, real_c(0), field::fzyx);
   }
}; // class ScalarFieldHandling

// function describing the global initialization profile
inline real_t initializationProfile(real_t x, real_t amplitude, real_t offset, real_t wavelength)
{
   return amplitude * std::cos(x / wavelength * real_c(2) * math::pi + math::pi) + offset;
}

real_t getHydrostaticDensity(real_t height, real_t referenceHeight, real_t gravitationalAcceleration)
{
   return real_c(1) + real_c(3) * gravitationalAcceleration * (height - referenceHeight);
}

void initializePoiseuilleProfile(StructuredBlockForest& forest, const BlockDataID& pdfFieldID,
                                 const ConstBlockDataID& fillFieldID, const real_t& averageBedHeight,
                                 const real_t& averageFluidHeight, const real_t& forcingX, const real_t& viscosity,
                                 real_t amplitude, real_t wavelength)
{
   WALBERLA_LOG_INFO_ON_ROOT("Initializing Poiseuille velocity profile");

   const real_t rho = real_t(1);

   for (auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt)
   {
      auto pdfField                        = blockIt->getData< PdfField_T >(pdfFieldID);
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(fillFieldID);

      WALBERLA_FOR_ALL_CELLS_XYZ(
         pdfField, const Vector3< real_t > coord = forest.getBlockLocalCellCenter(*blockIt, Cell(x, y, z));

         Vector3< real_t > velocity(real_t(0));

         auto localBedHeight = initializationProfile(coord[0], amplitude, averageBedHeight, wavelength);
         auto heightAboveBed = coord[2] - localBedHeight;

         const real_t fillLevel = fillField->get(x, y, z);

         if (heightAboveBed >= 0_r && fillLevel > 0_r) {
            velocity[0] =
               forcingX / (real_t(2) * viscosity) * heightAboveBed * (2_r * averageFluidHeight - heightAboveBed);
         } pdfField->setToEquilibrium(x, y, z, velocity, rho);)
   }
}

/***********************************************************************************************************************
 * Initialize the hydrostatic pressure in the direction in which a force is acting in ALL cells (regardless of a cell's
 * flag). The velocity remains unchanged.
 *
 * The force vector must have only one component, i.e., the direction of the force can only be in x-, y- or z-axis.
 * The variable fluidHeight determines the height at which the density is equal to reference density (=1).
 **********************************************************************************************************************/
template< typename PdfField_T >
void initHydrostaticPressure(const std::weak_ptr< StructuredBlockForest >& blockForestPtr,
                             const BlockDataID& pdfFieldID,
                             std::function< real_t(const Vector3< real_t >&) > hydrostaticDensityFct)
{
   WALBERLA_LOG_INFO_ON_ROOT("Initializing hydrostatic pressure");

   const auto blockForest = blockForestPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      PdfField_T* const pdfField = blockIt->getData< PdfField_T >(pdfFieldID);

      CellInterval local = pdfField->xyzSizeWithGhostLayer(); // block-, i.e., process-local cell interval

      for (auto cellIt = local.begin(); cellIt != local.end(); ++cellIt)
      {
         // initialize the (hydrostatic) pressure, i.e., LBM density
         // Bernoulli: p = p0 + density * gravity * height
         // => LBM (density=1): rho = rho0 + gravity * height = 1 + 1/cs^2 * g * h = 1 + 3 * g * h
         // shift global cell by 0.5 since density is set for cell center

         Vector3< real_t > cellCenter = blockForest->getBlockLocalCellCenter(*blockIt, *cellIt);
         const real_t rho             = hydrostaticDensityFct(cellCenter);

         const Vector3< real_t > velocity = pdfField->getVelocity(*cellIt);

         pdfField->setDensityAndVelocity(*cellIt, velocity, rho);
      }
   }
}

template< typename FreeSurfaceBoundaryHandling_T, typename PdfField_T, typename FlagField_T >
class MeanVelocityComputer
{
 public:
   MeanVelocityComputer(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                        const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                        const ConstBlockDataID& pdfFieldID, const std::shared_ptr< Vector3< real_t > >& meanVelocity,
                        real_t averagingFactor)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), pdfFieldID_(pdfFieldID),
        meanVelocity_(meanVelocity), averagingFactor_(averagingFactor)
   {}

   void operator()()
   {
      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      getMeanVelocity(blockForest, freeSurfaceBoundaryHandling);
   }

   void getMeanVelocity(const std::shared_ptr< const StructuredBlockForest >& blockForest,
                        const std::shared_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling)
   {
      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      // use separate variables for the velocity in each direction; syntax meanVelocity[0] does not work in OMP-macro
      real_t meanVelocityX = real_c(0);
      real_t meanVelocityY = real_c(0);
      real_t meanVelocityZ = real_c(0);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         const FlagField_T* const flagField = blockIt->template getData< const FlagField_T >(flagFieldID);
         const PdfField_T* const pdfField   = blockIt->template getData< const PdfField_T >(pdfFieldID_);

         WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, pdfFieldIt, pdfField,
                                    omp parallel for schedule(static) reduction(+:meanVelocityX)
                                       reduction(+:meanVelocityY) reduction(+:meanVelocityZ),
                                    {
            if (flagInfo.isLiquid(flagFieldIt) || flagInfo.isInterface(flagFieldIt))
            {
               const Vector3< real_t > velocity = pdfField->getVelocity(pdfFieldIt.cell());

               meanVelocityX += velocity[0];
               meanVelocityY += velocity[1];
               meanVelocityZ += velocity[2];
            }
                                    }) // WALBERLA_FOR_ALL_CELLS_OMP
      }

      Vector3< real_t > meanVelocity(meanVelocityX, meanVelocityY, meanVelocityZ);
      mpi::allReduceInplace< real_t >(meanVelocity, mpi::SUM);

      meanVelocity *= averagingFactor_;
      *meanVelocity_ = meanVelocity;
   };

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;

   const ConstBlockDataID pdfFieldID_;

   std::shared_ptr< Vector3< real_t > > meanVelocity_;
   real_t averagingFactor_;
}; // class MeanVelocityComputer

class ForcingAdjuster
{
 public:
   ForcingAdjuster(const shared_ptr< StructuredBlockStorage >& blocks, const BlockDataID& pdfFieldID,
                   real_t targetVelocity, real_t externalForcing, real_t proportionalGain, real_t derivativeGain,
                   real_t integralGain, real_t maxRamp, real_t minActuatingVariable, real_t maxActuatingVariable)
      : blocks_(blocks), pdfFieldID_(pdfFieldID), currentExternalForcing_(externalForcing),
        pid_(targetVelocity, externalForcing, proportionalGain, derivativeGain, integralGain, maxRamp,
             minActuatingVariable, maxActuatingVariable)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Creating PID controller with pg = " << pid_.getProportionalGain()
                                                                     << ", dg = " << pid_.getDerivateGain()
                                                                     << ", ig = " << pid_.getIntegralGain());
   }

   void operator()(const real_t currentMeanVelocity)
   {
      // compute new forcing value on root (since flow rate only known on root)
      WALBERLA_ROOT_SECTION()
      {
         real_t newExternalForcing = pid_.update(currentMeanVelocity);
         currentExternalForcing_   = newExternalForcing;
      }

      // send updated external forcing to all other processes
      mpi::broadcastObject(currentExternalForcing_);

      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         // get the field data out of the block
         auto pdf                     = block->getData< PdfField_T >(pdfFieldID_);
         pdf->latticeModel().force_0_ = currentExternalForcing_;
      }
   }

   real_t getExternalForcing() { return currentExternalForcing_; }
   void storePIDSnapshot(std::string filename)
   {
      WALBERLA_ROOT_SECTION() { pid_.writeStateToFile(filename); }
   }
   void loadPIDSnapshot(std::string filename) { pid_.readStateFromFile(filename); }

 private:
   shared_ptr< StructuredBlockStorage > blocks_;
   BlockDataID pdfFieldID_;

   real_t currentExternalForcing_;
   PIDController pid_;
}; // ForcingAdjuster

int main(int argc, char** argv)
{
   Environment walberlaEnv(argc, argv);

   WALBERLA_LOG_INFO_ON_ROOT("waLBerla Revision: " << std::string(WALBERLA_GIT_SHA1).substr(0, 8))

   if (argc < 2) { WALBERLA_ABORT("Please specify a parameter file as input argument.") }

   WALBERLA_LOG_DEVEL_ON_ROOT("Using generated lattice model.");
   auto configPtr = walberlaEnv.config();

   // print content of parameter file
   WALBERLA_LOG_INFO_ON_ROOT(*configPtr);

   WALBERLA_ROOT_SECTION()
   {
      std::ofstream file;
      file.open("parameterConfiguration.cfg");
      file << *configPtr;
      file.close();
   }

   // get block forest parameters from parameter file
   auto blockForestParameters            = configPtr->getOneBlock("BlockForestParameters");
   const Vector3< uint_t > cellsPerBlock = blockForestParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");
   const Vector3< bool > periodicity     = blockForestParameters.getParameter< Vector3< bool > >("periodicity");
   const bool loadSnapshot               = blockForestParameters.getParameter< bool >("loadSnapshot");
   const bool storeSnapshot              = blockForestParameters.getParameter< bool >("storeSnapshot");
   const uint_t snapshotFrequency        = blockForestParameters.getParameter< uint_t >("snapshotFrequency");
   const std::string snapshotBaseFolder  = blockForestParameters.getParameter< std::string >("snapshotBaseFolder");

   // get domain parameters from parameter file
   auto domainParameters              = configPtr->getOneBlock("DomainParameters");
   const Vector3< uint_t > domainSize = domainParameters.getParameter< Vector3< uint_t > >("domainSize");
   const uint_t wavePeriods           = domainParameters.getParameter< uint_t >("wavePeriods");
   const real_t liquidHeightFactor    = domainParameters.getParameter< real_t >("liquidHeightFactor");
   const real_t floorHeightFactor     = domainParameters.getParameter< real_t >("floorHeightFactor");
   const real_t initialAmplitude      = domainParameters.getParameter< real_t >("initialAmplitude");

   // compute number of blocks as defined by domainSize and cellsPerBlock
   Vector3< uint_t > numBlocks;
   for (uint i = 0; i <= 2; ++i)
   {
      numBlocks[i] = domainSize[i] / cellsPerBlock[i];
      WALBERLA_CHECK_EQUAL(numBlocks[i] * cellsPerBlock[i], domainSize[i],
                           "Domain size in direction " << i << " is not a multiple of cells per block.")
   }

   // get number of (MPI) processes
   const uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
   WALBERLA_CHECK_EQUAL(numProcesses, numBlocks[0] * numBlocks[1] * numBlocks[2],
                        "The number of MPI processes is different from the number of blocks as defined by "
                        "\"domainSize/cellsPerBlock\".")

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numProcesses);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numBlocks);

   // get PID controller parameters
   auto PIDParameters                       = configPtr->getOneBlock("PIDParameters");
   const real_t targetMeanVelocityMagnitude = PIDParameters.getParameter< real_t >("targetMeanVelocityMagnitude");
   const real_t proportionalGain            = PIDParameters.getParameter< real_t >("proportionalGain");
   const real_t derivativeGain              = PIDParameters.getParameter< real_t >("derivativeGain");
   const real_t integralGain                = PIDParameters.getParameter< real_t >("integralGain");
   const real_t maxRamp                     = PIDParameters.getParameter< real_t >("maxRamp");
   const real_t minActuatingVariable        = PIDParameters.getParameter< real_t >("minActuatingVariable");
   const real_t maxActuatingVariable        = PIDParameters.getParameter< real_t >("maxActuatingVariable");

   // read particle infos
   const auto particleParameters               = configPtr->getOneBlock("ParticleParameters");
   const std::string particleInFileName        = particleParameters.getParameter< std::string >("inFileName");
   const uint_t bedCopiesInX                   = particleParameters.getParameter< uint_t >("bedCopiesInX");
   const uint_t bedCopiesInY                   = particleParameters.getParameter< uint_t >("bedCopiesInY");
   const real_t particleDensityRatio           = particleParameters.getParameter< real_t >("densityRatio");
   const real_t particleFixingHeightFactor     = particleParameters.getParameter< real_t >("fixingHeightFactor");
   const real_t particleFrictionCoefficient    = particleParameters.getParameter< real_t >("frictionCoefficient");
   const real_t particleRestitutionCoefficient = particleParameters.getParameter< real_t >("restitutionCoefficient");
   const uint_t particleNumSubCycles           = particleParameters.getParameter< uint_t >("numSubCycles");
   const bool useLubricationCorrection         = particleParameters.getParameter< bool >("useLubricationCorrection");
   const bool useNoSlipParticles               = particleParameters.getParameter< bool >("useNoSlipParticles");
   const real_t particlePoissonsRatio          = 0.22_r;
   const real_t particleKappa = real_t(2) * (real_t(1) - particlePoissonsRatio) / (real_t(2) - particlePoissonsRatio);
   real_t particleCollisionTimeNonDim = 4_r;
   bool useOpenMP                     = false;
   const uint_t vtkSpacingParticles =
      configPtr->getOneBlock("VTK").getOneBlock("fluid_field").getParameter< uint_t >("writeFrequency");
   const std::string vtkFolder =
      configPtr->getOneBlock("VTK").getOneBlock("fluid_field").getParameter< std::string >("baseFolder");

   // get physics parameters from parameter file
   auto physicsParameters   = configPtr->getOneBlock("PhysicsParameters");
   const bool enableWetting = physicsParameters.getParameter< bool >("enableWetting");
   const uint_t timesteps   = physicsParameters.getParameter< uint_t >("timesteps");
   const real_t Re          = physicsParameters.getParameter< real_t >("Re");
   const real_t Fr          = physicsParameters.getParameter< real_t >("Fr");
   const real_t We          = physicsParameters.getParameter< real_t >("We");

   // get avgDiameter and scaling factor
   real_t avgParticleDiameter   = 0_r;
   real_t particleScalingFactor = 0_r;
   getAvgDiameterScalingFactor(particleInFileName, domainSize, bedCopiesInX, bedCopiesInY, avgParticleDiameter,
                               particleScalingFactor);
   const real_t liquidHeight         = avgParticleDiameter * liquidHeightFactor;
   const real_t floorHeight          = avgParticleDiameter * floorHeightFactor;
   const real_t particleFixingHeight = particleFixingHeightFactor * avgParticleDiameter;

   WALBERLA_CHECK_FLOAT_UNEQUAL(liquidHeight, 0.0)
   WALBERLA_CHECK_FLOAT_UNEQUAL(floorHeight, 0.0)

   const real_t absoluteLiquidHeight = liquidHeight + floorHeight;

   const real_t viscosity      = targetMeanVelocityMagnitude * liquidHeight / Re;
   const real_t relaxationRate = real_t(1.0) / (real_t(3) * viscosity + real_t(0.5));

   const real_t gravity = (targetMeanVelocityMagnitude / Fr) * (targetMeanVelocityMagnitude / Fr) / liquidHeight;
   const real_t fx      = real_t(3) * targetMeanVelocityMagnitude * viscosity / (liquidHeight * liquidHeight);
   Vector3< real_t > force(fx, real_t(0.0), -gravity);

   const real_t surfaceTension =
      real_t(1.0) * targetMeanVelocityMagnitude * targetMeanVelocityMagnitude * liquidHeight / We;

   // compute SI dx and dt
   const real_t viscosity_SI = real_t(1.0016e-6); // kinemtic viscosity of water at 20°C at 1 bar
   const real_t dx_SI        = 1_r / particleScalingFactor;
   const real_t dt_SI        = viscosity / viscosity_SI * dx_SI * dx_SI;

   WALBERLA_LOG_INFO_ON_ROOT("\nPhysical parameters:")
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(liquidHeight);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(floorHeight);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(dx_SI);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(dt_SI);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(force);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(relaxationRate);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(viscosity);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(absoluteLiquidHeight);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(avgParticleDiameter);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleScalingFactor);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleFixingHeight);
   WALBERLA_LOG_INFO_ON_ROOT("\nFree surface physical parameters")
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(surfaceTension);

   if ((periodicity[0] && numBlocks[0] < uint_c(3)) || (periodicity[1] && numBlocks[1] < uint_c(3)) ||
       (periodicity[2] && numBlocks[2] < uint_c(3)))
   {
      WALBERLA_ABORT("When using particles, use at least three blocks per periodic direction.");
   }

   // read model parameters from parameter file
   const auto modelParameters               = configPtr->getOneBlock("ModelParameters");
   const std::string pdfReconstructionModel = modelParameters.getParameter< std::string >("pdfReconstructionModel");
   const std::string pdfRefillingModel      = modelParameters.getParameter< std::string >("pdfRefillingModel");
   const std::string excessMassDistributionModel =
      modelParameters.getParameter< std::string >("excessMassDistributionModel");
   const std::string curvatureModel          = modelParameters.getParameter< std::string >("curvatureModel");
   const bool useSimpleMassExchange          = modelParameters.getParameter< bool >("useSimpleMassExchange");
   const real_t cellConversionThreshold      = modelParameters.getParameter< real_t >("cellConversionThreshold");
   const real_t cellConversionForceThreshold = modelParameters.getParameter< real_t >("cellConversionForceThreshold");

   // read evaluation parameters from parameter file
   const auto evaluationParameters      = configPtr->getOneBlock("EvaluationParameters");
   const uint_t performanceLogFrequency = evaluationParameters.getParameter< uint_t >("performanceLogFrequency");
   const uint_t evaluationFrequency     = evaluationParameters.getParameter< uint_t >("evaluationFrequency");
   const std::string baseFolderName     = evaluationParameters.getParameter< std::string >("baseFolderName");

   uint_t beginTimeStep = 0;
   const std::string checkpointConfigFile("antidunesCheckpointConfig.file");
   if (loadSnapshot)
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ifstream file;
         file.open(snapshotBaseFolder + "/" + checkpointConfigFile);
         if (file.fail()) WALBERLA_ABORT("Error: " << checkpointConfigFile << " could not be opened!");
         file >> beginTimeStep;
         file >> force[0];
         file.close();
      }
      mpi::broadcastObject(beginTimeStep);
      mpi::broadcastObject(force);

      WALBERLA_LOG_INFO_ON_ROOT("Successfully read config parameters from checkpoint config file:")
      WALBERLA_LOG_INFO_ON_ROOT(" - beginTimeStep = " << beginTimeStep)
      WALBERLA_LOG_INFO_ON_ROOT(" - force = < " << force[0] << ", " << force[1] << ", " << force[2] << " >")
   }

   if (loadSnapshot)
   {
      // modify config file to start VTK output from "loadFromTimestep" rather than from 0
      std::vector< config::Config::Block* > configVTKBlock;
      configPtr->getWritableGlobalBlock().getWritableBlocks("VTK", configVTKBlock, 1, 1);
      std::vector< config::Config::Block* > configVTKFluidFieldBlock;
      configVTKBlock[0]->getWritableBlocks("fluid_field", configVTKFluidFieldBlock, 1, 1);
      configVTKFluidFieldBlock[0]->setOrAddParameter("initialExecutionCount", std::to_string(beginTimeStep));
   }

   WALBERLA_ROOT_SECTION()
   {
      // create base directories if they do not yet exist
      filesystem::path tpath(baseFolderName);
      if (!filesystem::exists(tpath)) filesystem::create_directory(tpath);

      filesystem::path snapshotPath(snapshotBaseFolder);
      if (!filesystem::exists(snapshotPath)) filesystem::create_directory(snapshotPath);
   }

   std::shared_ptr< StructuredBlockForest > blockForest(nullptr);
   const std::string blockForestFile("blockForest.file");

   if (loadSnapshot)
   {
      // load block forest from file
      MPIManager::instance()->useWorldComm();

      blockForest = make_shared< StructuredBlockForest >(
         make_shared< BlockForest >(uint_c(MPIManager::instance()->rank()),
                                    (std::string(snapshotBaseFolder + "/" + blockForestFile)).c_str(), true, false),
         cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
      blockForest->createCellBoundingBoxes();
   }
   else
   {
      // create uniform block forest
      blockForest = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                                        cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                                        real_c(1.0),                                          // dx
                                                        true, // one block per process
                                                        periodicity[0], periodicity[1], periodicity[2]); // periodicity
   }

   // save block forest to file (but do not overwrite existing file if snapshot is loaded)
   if (storeSnapshot && !loadSnapshot)
   {
      blockForest->getBlockForest().saveToFile(snapshotBaseFolder + "/" + blockForestFile);
   }

   const auto vtkParameters      = configPtr->getOneBlock("VTK");
   const auto vtkFluidParameters = vtkParameters.getOneBlock("fluid_field");

   LatticeModel_T latticeModel = LatticeModel_T(force[0], force[1], force[2], relaxationRate);

   BlockDataID pdfFieldID;
   const std::string pdfFieldFile("pdfField.file");
   BlockDataID fillFieldID;
   const std::string fillFieldFile("fillField.file");

   if (loadSnapshot)
   {
      // load PDF field from file
      shared_ptr< lbm::internal::PdfFieldHandling< LatticeModel_T > > pdfFieldDataHandling =
         make_shared< lbm::internal::PdfFieldHandling< LatticeModel_T > >(
            blockForest, latticeModel, false, Vector3< real_t >(real_c(0)), real_c(1), uint_c(1), field::fzyx);
      pdfFieldID = (blockForest->getBlockStorage())
                      .loadBlockData(snapshotBaseFolder + "/" + pdfFieldFile, pdfFieldDataHandling, "PDF field");

      // load fill level field from file
      std::shared_ptr< ScalarFieldHandling< ScalarField_T > > fillFieldDataHandling =
         std::make_shared< ScalarFieldHandling< ScalarField_T > >(blockForest, uint_c(2));
      fillFieldID =
         (blockForest->getBlockStorage())
            .loadBlockData(snapshotBaseFolder + "/" + fillFieldFile, fillFieldDataHandling, "Fill level field");
   }
   else
   {
      // add PDF field
      pdfFieldID =
         lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel,
                                   Vector3< real_t >(targetMeanVelocityMagnitude, 0, 0), real_t(1.0), field::fzyx);

      // add fill level field (initialized with 0, i.e., gas everywhere)
      fillFieldID =
         field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(0.0), field::fzyx, uint_c(2));
   }

   // MesaPD data structures
   auto particleStorage  = std::make_shared< mesa_pd::data::ParticleStorage >(1);
   auto particleAccessor = std::make_shared< mesa_pd::data::ParticleAccessorWithBaseShape >(particleStorage);
   auto mesapdDomain     = std::make_shared< mesa_pd::domain::BlockForestDomain >(blockForest->getBlockForestPointer());

   BlockDataID particleStorageID;
   const std::string particleStorageFile("particleStorageFile.file");
   if (loadSnapshot)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Initializing particles from checkpointing file!");
      particleStorageID = blockForest->loadBlockData(snapshotBaseFolder + "/" + particleStorageFile,
                                                     mesa_pd::domain::createBlockForestDataHandling(particleStorage),
                                                     "Particle Storage");
      mesa_pd::mpi::ClearNextNeighborSync CNNS;
      CNNS(*particleAccessor);

      mesa_pd::mpi::ClearGhostOwnerSync CGOS;
      CGOS(*particleAccessor);
   }
   else
   {
      particleStorageID =
         blockForest->addBlockData(mesa_pd::domain::createBlockForestDataHandling(particleStorage), "Particle Storage");
   }

   BlockDataID particleFieldID = field::addToStorage< lbm_mesapd_coupling::ParticleField_T >(
      blockForest, "Particle field", particleAccessor->getInvalidUid(), field::fzyx, uint_c(2));

   auto densityReferenceHeight = absoluteLiquidHeight;
   auto hydrostaticDensityFct  = [force, densityReferenceHeight](const Vector3< real_t >& position) {
      uint_t forceComponent = 2; // gravity is here strictly only acting in z direction!
      return getHydrostaticDensity(position[forceComponent], densityReferenceHeight, force[forceComponent]);
   };

   // add boundary handling
   const std::shared_ptr< AntidunesBoundaryHandling_T > antidunesBoundaryHandling =
      std::make_shared< AntidunesBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID, particleFieldID,
                                                      particleAccessor, hydrostaticDensityFct);
   const BlockDataID flagFieldID                                    = antidunesBoundaryHandling->getFlagFieldID();
   const typename AntidunesBoundaryHandling_T::FlagInfo_T& flagInfo = antidunesBoundaryHandling->getFlagInfo();

   real_t sinusAmplitude  = real_c(0.5) * initialAmplitude;
   real_t sinusOffset     = floorHeight;
   real_t sinusWavelength = real_c(domainSize[0]) / real_c(wavePeriods);

   if (!loadSnapshot)
   {
      // samples used in the Monte-Carlo like estimation of the fill level
      const uint_t fillLevelInitSamples = uint_c(100); // actually there will be 101 since 0 is also included

      const uint_t numTotalPoints = (fillLevelInitSamples + uint_c(1)) * (fillLevelInitSamples + uint_c(1));
      const real_t stepsize       = real_c(1) / real_c(fillLevelInitSamples);

      // initialize free-surface sine profile
      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(fillFieldID);

         WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
            // cell in block-local coordinates
            const Cell localCell = fillFieldIt.cell();

            // get cell in global coordinates
            Cell globalCell = localCell;
            blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

            // Monte-Carlo like estimation of the fill level:
            // create uniformly-distributed sample points in each cell and count the number of points below the sine
            // profile; this fraction of points is used as the fill level to initialize the profile
            uint_t numPointsBelow = uint_c(0);

            for (uint_t xSample = uint_c(0); xSample <= fillLevelInitSamples; ++xSample)
            {
               // Pascal et al. (2021) defined the amplitude to span from minimum peak to maximum peak; in
               // initializationProfile(), the amplitude is defined to range from the average to the maximum peak
               const real_t functionValue =
                  initializationProfile(real_c(globalCell[0]) + real_c(xSample) * stepsize, sinusAmplitude,
                                        absoluteLiquidHeight + real_c(0.5), sinusWavelength);

               for (uint_t zSample = uint_c(0); zSample <= fillLevelInitSamples; ++zSample)
               {
                  const real_t zPoint = real_c(globalCell[2]) + real_c(zSample) * stepsize;
                  // with operator <, a fill level of 1 can not be reached when the line is equal to the cell's top
                  // border; with operator <=, a fill level of 0 can not be reached when the line is equal to the cell's
                  // bottom border
                  if (zPoint < functionValue) { ++numPointsBelow; }
               }
            }

            // fill level is fraction of points below sine profile
            fillField->get(localCell) = real_c(numPointsBelow) / real_c(numTotalPoints);
         }) // WALBERLA_FOR_ALL_CELLS
      }

      initializePoiseuilleProfile(*blockForest, pdfFieldID, fillFieldID, floorHeight, liquidHeight + real_c(0.5),
                                  force[0], viscosity, sinusAmplitude, sinusWavelength);
   }

   // initialize domain boundary conditions from config file
   const auto boundaryParameters = configPtr->getOneBlock("BoundaryParameters");
   antidunesBoundaryHandling->initFromConfig(boundaryParameters);

   std::function< void(void) > syncCall;
   auto simulationDomainAABB = blockForest->getDomain();

   lbm_mesapd_coupling::ParticleMappingKernel< AntidunesBoundaryHandling_T::BoundaryHandling_T > particleMappingKernel(
      blockForest, antidunesBoundaryHandling->getHandlingID());
   lbm_mesapd_coupling::MovingParticleMappingKernel< AntidunesBoundaryHandling_T::BoundaryHandling_T >
      movingParticleMappingKernel(blockForest, antidunesBoundaryHandling->getHandlingID(), particleFieldID);

   uint_t numParticles = 0;
   // initialize bottom solid sine profile
   if (!loadSnapshot)
   {
      auto createParticleFct = [sinusAmplitude, sinusOffset, sinusWavelength](Vector3< real_t > pos) {
         return pos[2] < initializationProfile(pos[0], sinusAmplitude, sinusOffset, sinusWavelength);
      };

      real_t maxParticleHeight = 0_r;
      initSpheresFromFile(particleInFileName, *particleStorage, *mesapdDomain, particleDensityRatio, domainSize,
                          createParticleFct, simulationDomainAABB, bedCopiesInX, bedCopiesInY, numParticles,
                          maxParticleHeight, particleScalingFactor);
      WALBERLA_LOG_INFO_ON_ROOT("Max particle height " << maxParticleHeight);
      if ((sinusOffset + sinusAmplitude) > maxParticleHeight)
         WALBERLA_ABORT("Created particle bed is below desired sinus shape!");
      if (2_r * sinusAmplitude > (maxParticleHeight - particleFixingHeight))
         WALBERLA_ABORT("Created mobile particle bed is not high enough for desired sinus shape!");
      if (useNoSlipParticles && (particleFixingHeight < maxParticleHeight))
         WALBERLA_ABORT("You are using no-slip BCs on particles (which does not set hydrodynamic forces) but do not "
                        "fix all particles - this leads to wrong behavior and is not permitted!")

      // fix lower particles
      particleStorage->forEachParticle(
         useOpenMP, mesa_pd::kernel::SelectAll(), *particleAccessor,
         [particleFixingHeight](const size_t idx, auto& ac) {
            if (ac.getPosition(idx)[2] < particleFixingHeight)
               mesa_pd::data::particle_flags::set(ac.getFlagsRef(idx), mesa_pd::data::particle_flags::FIXED);
         },
         *particleAccessor);
   }
   else
   {
      real_t avgParticleDiameterTest = 0_r;
      particleStorage->forEachParticle(
         false, mesa_pd::kernel::SelectLocal(), *particleAccessor,
         [&numParticles, &avgParticleDiameterTest](const size_t idx, auto& ac) {
            auto sp = static_cast< mesa_pd::data::Sphere* >(ac.getBaseShape(idx).get());
            ++numParticles;
            avgParticleDiameterTest += 2_r * sp->getRadius();
         },
         *particleAccessor);
      mpi::allReduceInplace(numParticles, mpi::SUM);
      mpi::allReduceInplace(avgParticleDiameterTest, mpi::SUM);
      avgParticleDiameterTest /= real_c(numParticles);
      WALBERLA_LOG_INFO_ON_ROOT("Read particles from check pointing file with avg diameter of "
                                << avgParticleDiameterTest)
      if (std::abs(avgParticleDiameterTest - avgParticleDiameter) / avgParticleDiameterTest > 0.05)
      {
         WALBERLA_ABORT("Particle diameters not correct.")
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT("Created " << numParticles << " particles");

   // create planes
   createPlane(*particleStorage, simulationDomainAABB.minCorner(), Vector3< real_t >(0, 0, 1));
   createPlane(*particleStorage, simulationDomainAABB.maxCorner(), Vector3< real_t >(0, 0, -1));

   const real_t blockSyncExtension    = real_t(2.5);
   real_t maxPossibleParticleDiameter = avgParticleDiameter * 1.1_r;
   if (maxPossibleParticleDiameter < 2_r * real_c(cellsPerBlock.min()) - blockSyncExtension)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Using next neighbor sync for particles");
      syncCall = [particleStorage, mesapdDomain, blockSyncExtension]() {
         mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
         syncNextNeighborFunc(*particleStorage, *mesapdDomain, blockSyncExtension);
      };
      syncCall();
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT("Using ghost owner sync for particles")
      syncCall = [particleStorage, mesapdDomain, blockSyncExtension]() {
         mesa_pd::mpi::SyncGhostOwners syncGhostOwnersFunc;
         syncGhostOwnersFunc(*particleStorage, *mesapdDomain, blockSyncExtension);
      };
      for (uint_t i = 0; i < uint_c(std::ceil(maxPossibleParticleDiameter / real_c(cellsPerBlock.min()))); ++i)
         syncCall();
   }

   if (useNoSlipParticles)
   {
      particleStorage->forEachParticle(useOpenMP, SphereSelector(), *particleAccessor, particleMappingKernel,
                                       *particleAccessor, AntidunesBoundaryHandling_T::noSlipFlagID);
   }
   else
   {
      particleStorage->forEachParticle(useOpenMP, SphereSelector(), *particleAccessor, movingParticleMappingKernel,
                                       *particleAccessor, AntidunesBoundaryHandling_T::movingObstacleFlagID);
   }

   // IMPORTANT REMARK: this must be only called after every solid flag has been set; otherwise, the boundary handling
   // might not detect solid flags correctly
   antidunesBoundaryHandling->initFlagsFromFillLevel();

   // communication after initialization
   Communication_T communication(blockForest, flagFieldID, fillFieldID);
   communication();

   PdfCommunication_T pdfCommunication(blockForest, pdfFieldID);
   pdfCommunication();

   // add bubble model
   std::shared_ptr< walberla::free_surface::bubble_model::BubbleModelBase > bubbleModel =
      std::make_shared< walberla::free_surface::bubble_model::BubbleModelConstantPressure >(real_c(1));

   // initialize hydrostatic pressure
   if (!loadSnapshot) { initHydrostaticPressure< PdfField_T >(blockForest, pdfFieldID, hydrostaticDensityFct); }

   // set density in non-liquid or non-interface cells to 1 (after initializing with hydrostatic pressure)
   // setDensityInNonFluidCellsToOne< FlagField_T, PdfField_T >(blockForest, flagInfo, flagFieldID, pdfFieldID);

   // create timeloop
   SweepTimeloop timeloop(blockForest, timesteps);
   timeloop.setCurrentTimeStep(beginTimeStep);

   timeloop.addFuncBeforeTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");

   // Laplace pressure = 2 * surface tension * curvature; curvature computation is not necessary with no surface
   // tension
   bool computeCurvature = false;
   if (!realIsEqual(surfaceTension, real_c(0), real_c(1e-14))) { computeCurvature = true; }

   auto blockStateUpdate = StateSweep(blockForest, flagInfo, flagFieldID);

   // add surface geometry handler
   BlockDataID curvatureFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Curvature field", real_c(0), field::fzyx, uint_c(1));
   BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID obstacleNormalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Obstacle normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   // add field for smoothed fill levels
   BlockDataID smoothFillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Smooth fill level field", real_c(0), field::fzyx, uint_c(1));

   // smooth fill level field for decreasing error in finite difference normal and curvature computation (see
   // dissertation of S. Bogner, 2017 (section 4.4.2.1))
   walberla::free_surface::SmoothingSweep< CommunicationStencil_T, FlagField_T, ScalarField_T, VectorField_T >
      smoothingSweep(smoothFillFieldID, fillFieldID, flagFieldID,
                     walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, flagInfo.getObstacleIDSet(),
                     enableWetting);
   // IMPORTANT REMARK: SmoothingSweep must be executed on all blocks, because the algorithm works on all liquid,
   // interface and gas cells. This is necessary since the normals are not only computed in interface cells, but also
   // in the neighborhood of interface cells. Therefore, meaningful values for the fill levels of the second
   // neighbors of interface cells are also required in NormalSweep.
   timeloop.add() << Sweep(smoothingSweep, "Sweep: fill level smoothing")
                  << AfterFunction(Communication_T(blockForest, smoothFillFieldID),
                                   "Communication: after smoothing sweep");

   // compute interface normals (using smoothed fill level field)
   walberla::free_surface::NormalSweep< CommunicationStencil_T, FlagField_T, ScalarField_T, VectorField_T > normalSweep(
      normalFieldID, smoothFillFieldID, flagFieldID, walberla::free_surface::flagIDs::interfaceFlagID,
      walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, flagInfo.getObstacleIDSet(), true, false, true,
      false);
   timeloop.add() << Sweep(normalSweep, "Sweep: normal computation", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: normal")
                  << AfterFunction(Communication_T(blockForest, normalFieldID), "Communication: after normal sweep");

   if (computeCurvature)
   {
      // compute interface curvature using finite differences according to Brackbill et al.
      walberla::free_surface::CurvatureSweepFiniteDifferences< CommunicationStencil_T, FlagField_T, ScalarField_T,
                                                               VectorField_T >
         curvSweep(curvatureFieldID, normalFieldID, obstacleNormalFieldID, flagFieldID,
                   walberla::free_surface::flagIDs::interfaceFlagID,
                   walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, flagInfo.getObstacleIDSet(), false,
                   real_c(0));
      timeloop.add() << Sweep(curvSweep, "Sweep: curvature computation (finite difference method)",
                              StateSweep::fullFreeSurface)
                     << Sweep(emptySweep(), "Empty sweep: curvature")
                     << AfterFunction(Communication_T(blockForest, curvatureFieldID),
                                      "Communication: after curvature sweep");
   }

   // add surface dynamics handler

   // add standard waLBerla boundary handling
   timeloop.add() << Sweep(antidunesBoundaryHandling->getBoundarySweep(), "Sweep: boundary handling",
                           Set< SUID >::emptySet(), StateSweep::onlyGasAndBoundary)
                  << Sweep(emptySweep(), "Empty sweep: boundary handling", StateSweep::onlyGasAndBoundary);

   // sweep for
   // - reconstruction of PDFs in interface cells
   // - streaming of PDFs in interface cells (and liquid cells on the same block)
   // - advection of mass
   // - update bubble volumes
   // - marking interface cells for conversion
   const walberla::free_surface::StreamReconstructAdvectSweep<
      LatticeModel_T, typename AntidunesBoundaryHandling_T::BoundaryHandling_T, FlagField_T,
      typename AntidunesBoundaryHandling_T::FlagInfo_T, ScalarField_T, VectorField_T, true >
      streamReconstructAdvectSweep(surfaceTension, antidunesBoundaryHandling->getHandlingID(), fillFieldID, flagFieldID,
                                   pdfFieldID, normalFieldID, curvatureFieldID, flagInfo, bubbleModel.get(),
                                   pdfReconstructionModel, useSimpleMassExchange, cellConversionThreshold,
                                   cellConversionForceThreshold);
   // sweep acts only on blocks with at least one interface cell (due to StateSweep::fullFreeSurface)
   timeloop.add() << Sweep(streamReconstructAdvectSweep, "Sweep: StreamReconstructAdvect", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: StreamReconstructAdvect")
                  // do not communicate PDFs here:
                  // - stream on blocks with "StateSweep::fullFreeSurface" was performed here using post-collision PDFs
                  // - stream on other blocks is performed below and should also use post-collision PDFs
                  // => if PDFs were communicated here, the ghost layer of other blocks would have post-stream PDFs
                  << AfterFunction(Communication_T(blockForest, fillFieldID, flagFieldID),
                                   "Communication: after StreamReconstructAdvect sweep")
                  << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(blockForest, fillFieldID),
                                   "Second ghost layer update: after StreamReconstructAdvect sweep (fill level field)")
                  << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(blockForest, flagFieldID),
                                   "Second ghost layer update: after StreamReconstructAdvect sweep (flag field)");

   auto lbmSweepGenerated = typename LatticeModel_T::Sweep(pdfFieldID);

   // temporary class for being able to call the LBM collision with operator()
   class CollideSweep
   {
    public:
      CollideSweep(const typename LatticeModel_T::Sweep& sweep) : sweep_(sweep){};

      void operator()(IBlock* const block, const uint_t numberOfGhostLayersToInclude = uint_t(0))
      {
         sweep_.collide(block, numberOfGhostLayersToInclude);
      }

    private:
      typename LatticeModel_T::Sweep sweep_;
   };

   timeloop.add() << Sweep(CollideSweep(lbmSweepGenerated), "Sweep: collision (generated)", StateSweep::fullFreeSurface)
                  << Sweep(lbmSweepGenerated, "Sweep: streamCollide (generated)", StateSweep::onlyLBM)
                  << Sweep(emptySweep(), "Empty sweep: streamCollide (generated)")
                  << AfterFunction(PdfCommunication_T(blockForest, pdfFieldID),
                                   "Communication: after streamCollide (generated)");

   // convert cells
   // - according to the flags from StreamReconstructAdvectSweep (interface -> gas/liquid)
   // - to ensure a closed layer of interface cells (gas/liquid -> interface)
   // - detect and register bubble merges/splits (bubble volumes are already updated in StreamReconstructAdvectSweep)
   // - convert cells and initialize PDFs near inflow boundaries
   const walberla::free_surface::CellConversionSweep< LatticeModel_T, AntidunesBoundaryHandling_T::BoundaryHandling_T,
                                                      ScalarField_T >
      cellConvSweep(antidunesBoundaryHandling->getHandlingID(), pdfFieldID, flagInfo, bubbleModel.get());
   timeloop.add() << Sweep(cellConvSweep, "Sweep: cell conversion", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: cell conversion")
                  //<< AfterFunction(PdfCommunication_T(blockForest, pdfFieldID),
                  //
                  //                 "Communication: after cell conversion sweep (PDF field)")
                  // communicate the flag field also in corner directions
                  << AfterFunction(Communication_T(blockForest, flagFieldID),
                                   "Communication: after cell conversion sweep (flag field)")
                  << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(blockForest, flagFieldID),
                                   "Second ghost layer update: after cell conversion sweep (flag field)");

   // reinitialize PDFs, i.e., refill cells that were converted from gas to interface
   // - when the flag "convertedFromGasToInterface" has been set (by CellConversionSweep)
   // - according to the method specified with pdfRefillingModel_
   const walberla::free_surface::EquilibriumRefillingSweep< LatticeModel_T, FlagField_T > equilibriumRefillingSweep(
      pdfFieldID, flagFieldID, flagInfo, true);
   timeloop.add() << Sweep(equilibriumRefillingSweep, "Sweep: EquilibriumRefilling", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: EquilibriumRefilling")
                  << AfterFunction(PdfCommunication_T(blockForest, pdfFieldID),
                                   "Communication: after EquilibriumRefilling sweep");

   // distribute excess mass:
   // - excess mass: mass that is free after conversion from interface to gas/liquid cells
   // - update the bubble model
   // IMPORTANT REMARK: this sweep computes the mass via the density, i.e., the PDF field must be up-to-date and the
   // PdfRefillingSweep must have been performed
   const walberla::free_surface::ExcessMassDistributionSweepInterfaceEvenly< LatticeModel_T, FlagField_T, ScalarField_T,
                                                                             VectorField_T >
      distributeMassSweep(excessMassDistributionModel, fillFieldID, flagFieldID, pdfFieldID, flagInfo);
   timeloop.add() << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: distribute excess mass")
                  << AfterFunction(Communication_T(blockForest, fillFieldID),
                                   "Communication: after excess mass distribution sweep")
                  << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(blockForest, fillFieldID),
                                   "Second ghost layer update: after excess mass distribution sweep (fill level field)")
                  // update bubble model, i.e., perform registered bubble merges/splits; bubble merges/splits are
                  // already detected and registered by CellConversionSweep
                  << AfterFunction(
                        std::bind(&walberla::free_surface::bubble_model::BubbleModelBase::update, bubbleModel),
                        "Sweep: bubble model update");

   // reset all flags that signal cell conversions (except "keepInterfaceForWettingFlag")
   walberla::free_surface::ConversionFlagsResetSweep< FlagField_T > resetConversionFlagsSweep(flagFieldID, flagInfo);
   timeloop.add() << Sweep(resetConversionFlagsSweep, "Sweep: conversion flag reset", StateSweep::fullFreeSurface)
                  << Sweep(emptySweep(), "Empty sweep: conversion flag reset")
                  << AfterFunction(Communication_T(blockForest, flagFieldID),
                                   "Communication: after excess mass distribution sweep")
                  << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(blockForest, flagFieldID),
                                   "Second ghost layer update: after excess mass distribution sweep (flag field)");

   // update block states
   timeloop.add() << Sweep(blockStateUpdate, "Sweep: block state update");

   // add VTK output
   walberla::free_surface::addVTKOutput< LatticeModel_T, AntidunesBoundaryHandling_T, PdfField_T, FlagField_T,
                                         ScalarField_T, VectorField_T >(
      blockForest, timeloop, configPtr, flagInfo, pdfFieldID, flagFieldID, fillFieldID, BlockDataID(), curvatureFieldID,
      normalFieldID, obstacleNormalFieldID);

   // add triangle mesh output of free surface
   walberla::free_surface::SurfaceMeshWriter< ScalarField_T, FlagField_T > surfaceMeshWriter(
      blockForest, fillFieldID, flagFieldID, walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, real_c(0),
      configPtr);
   surfaceMeshWriter(); // write initial mesh
   timeloop.addFuncAfterTimeStep(surfaceMeshWriter, "Writer: surface mesh");

   if (vtkSpacingParticles != uint_t(0))
   {
      // particle field
      auto particleFieldVTK =
         vtk::createVTKOutput_BlockData(blockForest, "particle_field", vtkSpacingParticles, 0, false, vtkFolder);
      auto cellBB_filterParameters             = vtkFluidParameters.getOneBlock("CellBB_filter");
      const Vector3< uint_t > cellBB_filterMin = cellBB_filterParameters.getParameter< Vector3< uint_t > >("min");
      const Vector3< uint_t > cellBB_filterMax = cellBB_filterParameters.getParameter< Vector3< uint_t > >("max");
      AABB sliceAABB(real_c(cellBB_filterMin[0]), real_c(cellBB_filterMin[1]), real_c(cellBB_filterMin[2]),
                     real_c(cellBB_filterMax[0] + uint_t(1)), real_c(cellBB_filterMax[1] + uint_t(1)),
                     real_c(cellBB_filterMax[2] + uint_t(1)));

      particleFieldVTK->addCellInclusionFilter(vtk::AABBCellFilter(sliceAABB));
      particleFieldVTK->addCellDataWriter(
         make_shared< field::VTKWriter< GhostLayerField< walberla::id_t, 1 > > >(particleFieldID, "particleField"));
      particleFieldVTK->setSamplingResolution(vtkFluidParameters.getParameter< real_t >("samplingResolution"));
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(particleFieldVTK), "VTK (particle field data");
   }

   if (vtkSpacingParticles != uint_t(0))
   {
      // sphere
      auto particleVtkOutput = make_shared< mesa_pd::vtk::ParticleVtkOutput >(particleStorage);
      particleVtkOutput->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      particleVtkOutput->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      particleVtkOutput->addOutput< mesa_pd::data::SelectParticleInteractionRadius >("radius");
      // limit output to process-local spheres
      particleVtkOutput->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return (pIt->getBaseShape()->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE) &&
                !isSet(pIt->getFlags(), GHOST);
      });
      auto particleVtkWriter =
         vtk::createVTKOutput_PointData(particleVtkOutput, "particles", vtkSpacingParticles, vtkFolder,
                                        std::string("simulation_step"), false, true, true, true, beginTimeStep);
      timeloop.addFuncAfterTimeStep(vtk::writeFiles(particleVtkWriter), "VTK (sphere data)");
   }

   // add logging for computational performance
   const lbm::PerformanceLogger< FlagField_T > performanceLogger(
      blockForest, flagFieldID, walberla::free_surface::flagIDs::liquidInterfaceFlagIDs, performanceLogFrequency);
   timeloop.addFuncAfterTimeStep(performanceLogger, "Evaluator: performance logging");

   // LBM stability check
   timeloop.addFuncAfterTimeStep(makeSharedFunctor(field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                    walberlaEnv.config(), blockForest, pdfFieldID, flagFieldID,
                                    walberla::free_surface::flagIDs::liquidInterfaceFlagIDs)),
                                 "LBM stability check");

   // add sweep for evaluating the fluid's mean velocity
   const std::shared_ptr< Vector3< real_t > > meanVelocity = std::make_shared< Vector3< real_t > >(real_c(0));
   const real_t velocityAveragingFactor = 1_r / (liquidHeight * real_c(domainSize[0]) * real_c(domainSize[1]));
   MeanVelocityComputer< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T > meanVelocityComputer(
      blockForest, antidunesBoundaryHandling, pdfFieldID, meanVelocity, velocityAveragingFactor);

   // PID Controller
   shared_ptr< ForcingAdjuster > forcingAdjuster =
      make_shared< ForcingAdjuster >(blockForest, pdfFieldID, targetMeanVelocityMagnitude, force[0], proportionalGain,
                                     derivativeGain, integralGain, maxRamp, minActuatingVariable, maxActuatingVariable);

   if (loadSnapshot) { forcingAdjuster->loadPIDSnapshot(snapshotBaseFolder + "/" + "pidState.file"); }

   WcTimingPool timingPool;

   // this is carried out after the particle integration, it corrects the flag field and restores missing PDF
   // information then, the checkpointing file can be written, as otherwise some cells are invalid and can not be
   // recovered
   SweepTimeloop timeloopAfterParticles(blockForest, timesteps);
   timeloopAfterParticles.setCurrentTimeStep(beginTimeStep);

   // sweep for updating the particle mapping into the LBM simulation
   bool strictlyConserveMomentum = false;
   timeloopAfterParticles.add() << Sweep(
      lbm_mesapd_coupling::makeMovingParticleMapping< PdfField_T, AntidunesBoundaryHandling_T::BoundaryHandling_T >(
         blockForest, pdfFieldID, antidunesBoundaryHandling->getHandlingID(), particleFieldID, particleAccessor,
         AntidunesBoundaryHandling_T::movingObstacleFlagID, FormerMO_Flag,
         lbm_mesapd_coupling::RegularParticlesSelector(), strictlyConserveMomentum),
      "Particle Mapping");

   // sweep for restoring PDFs in cells previously occupied by particles
   bool reconstruction_recomputeTargetDensity = false;
   bool reconstruction_useCentralDifferences  = true;
   auto gradReconstructor =
      lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor< AntidunesBoundaryHandling_T::BoundaryHandling_T >(
         blockForest, antidunesBoundaryHandling->getHandlingID(), relaxationRate, reconstruction_recomputeTargetDensity,
         reconstruction_useCentralDifferences);

   timeloopAfterParticles.add()
      << Sweep(makeSharedSweep(
                  lbm_mesapd_coupling::makePdfReconstructionManager< PdfField_T,
                                                                     AntidunesBoundaryHandling_T::BoundaryHandling_T >(
                     blockForest, pdfFieldID, antidunesBoundaryHandling->getHandlingID(), particleFieldID,
                     particleAccessor, FormerMO_Flag, walberla::free_surface::flagIDs::liquidFlagID, gradReconstructor,
                     strictlyConserveMomentum)),
               "PDF Restore")
      << AfterFunction(Communication_T(blockForest, flagFieldID, particleFieldID),
                       "Communication: after PDF reconstruction sweep") // unsure if necessary but added for consistency
      << AfterFunction(pdfCommunication, "PDF Communication");

   real_t timeStepSizeParticles = real_t(1) / real_t(particleNumSubCycles);
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce(timeStepSizeParticles);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(timeStepSizeParticles);
   mesa_pd::kernel::LinearSpringDashpot collisionResponse(1);
   collisionResponse.setFrictionCoefficientDynamic(0, 0, particleFrictionCoefficient);
   mesa_pd::mpi::ReduceProperty reduceProperty;
   mesa_pd::mpi::ReduceContactHistory reduceAndSwapContactHistory;
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;
   real_t particleCollisionTime = particleCollisionTimeNonDim * avgParticleDiameter;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel(
      viscosity, [](real_t r) { return (real_t(0.001 + real_t(0.00007) * r)) * r; });

   WALBERLA_LOG_INFO_ON_ROOT("Will use particle time step size of "
                             << timeStepSizeParticles << " and collision time of " << particleCollisionTime);

   AverageDataSliceEvaluator< PdfField_T, AntidunesBoundaryHandling_T, FlagField_T, ScalarField_T >
      averageDataSliceEvaluator(blockForest, flagFieldID, fillFieldID, pdfFieldID);

   std::shared_ptr< real_t > totalFluidMass = std::make_shared< real_t >(real_c(0));
   walberla::free_surface::TotalMassComputer< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T >
      totalFluidMassEvaluator(blockForest, antidunesBoundaryHandling, pdfFieldID, fillFieldID, evaluationFrequency,
                              totalFluidMass);

   BedloadTransportEvaluator< ParticleAccessor_T > bedloadTransportEvaluator(
      particleAccessor, 1_r / real_c(domainSize[0] * domainSize[1]), numParticles);
   auto bedLoadTransportFileName = baseFolderName + "/bedload.txt";
   WALBERLA_LOG_INFO_ON_ROOT("Writing bedload info to file " << bedLoadTransportFileName);

   auto fluidInfoFileName = baseFolderName + "/fluidInfo.txt";
   WALBERLA_LOG_INFO_ON_ROOT("Writing fluid info to file " << fluidInfoFileName);

   // write info file
   WALBERLA_ROOT_SECTION()
   {
      std::ofstream evalInfoFile(baseFolderName + "/info.txt");
      evalInfoFile << evaluationFrequency << "\n";
      evalInfoFile << gravity << "\n";
      evalInfoFile << viscosity << "\n";
      evalInfoFile << particleDensityRatio << "\n";
      evalInfoFile << avgParticleDiameter << "\n";
      evalInfoFile << domainSize[0] << "\n";
      evalInfoFile << domainSize[1] << "\n";
      evalInfoFile << domainSize[2] << "\n";
      evalInfoFile << numParticles << "\n";
      evalInfoFile << dx_SI << "\n";
      evalInfoFile << dt_SI << "\n";
      evalInfoFile.close();
   }

   Vector3< real_t > totalHydrodynamicForceOnParticles(0_r); // only root will have valid values

   for (uint_t t = beginTimeStep; t != timesteps; ++t)
   {
      timeloop.singleStep(timingPool, true);

      timingPool["Mesa_pd"].start();

      reduceProperty.operator()< mesa_pd::HydrodynamicForceTorqueNotification >(*particleStorage);

      if (t == 0)
      {
         lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel
            initializeHydrodynamicForceTorqueForAveragingKernel;
         particleStorage->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *particleAccessor,
                                          initializeHydrodynamicForceTorqueForAveragingKernel, *particleAccessor);
      }
      particleStorage->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *particleAccessor,
                                       averageHydrodynamicForceTorque, *particleAccessor);

      for (auto subCycle = uint_t(0); subCycle < particleNumSubCycles; ++subCycle)
      {
         particleStorage->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *particleAccessor,
                                          vvIntegratorPreForce, *particleAccessor);
         syncCall();

         if (useLubricationCorrection)
         {
            // lubrication correction
            particleStorage->forEachParticlePairHalf(
               useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *particleAccessor,
               [&lubricationCorrectionKernel, &mesapdDomain](const size_t idx1, const size_t idx2, auto& ac) {
                  mesa_pd::collision_detection::AnalyticContactDetection acd;
                  acd.getContactThreshold() = lubricationCorrectionKernel.getNormalCutOffDistance();
                  mesa_pd::kernel::DoubleCast double_cast;
                  mesa_pd::mpi::ContactFilter contact_filter;
                  if (double_cast(idx1, idx2, ac, acd, ac))
                  {
                     if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *mesapdDomain))
                     {
                        double_cast(acd.getIdx1(), acd.getIdx2(), ac, lubricationCorrectionKernel, ac,
                                    acd.getContactNormal(), acd.getPenetrationDepth());
                     }
                  }
               },
               *particleAccessor);
         }

         // collision response
         particleStorage->forEachParticlePairHalf(
            useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *particleAccessor,
            [&collisionResponse, &mesapdDomain, timeStepSizeParticles, particleRestitutionCoefficient,
             particleCollisionTime, particleKappa](const size_t idx1, const size_t idx2, auto& ac) {
               mesa_pd::collision_detection::AnalyticContactDetection acd;
               mesa_pd::kernel::DoubleCast double_cast;
               mesa_pd::mpi::ContactFilter contact_filter;
               if (double_cast(idx1, idx2, ac, acd, ac))
               {
                  if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *mesapdDomain))
                  {
                     auto meff = real_t(1) / (ac.getInvMass(idx1) + ac.getInvMass(idx2));
                     collisionResponse.setStiffnessAndDamping(0, 0, particleRestitutionCoefficient,
                                                              particleCollisionTime, particleKappa, meff);
                     collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(),
                                       acd.getPenetrationDepth(), timeStepSizeParticles);
                  }
               }
            },
            *particleAccessor);

         reduceAndSwapContactHistory(*particleStorage);

         // add hydrodynamic force
         lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
         particleStorage->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *particleAccessor,
                                          addHydrodynamicInteraction, *particleAccessor);

         // add external forces
         particleStorage->forEachParticle(
            useOpenMP, mesa_pd::kernel::SelectLocal(), *particleAccessor,
            [particleDensityRatio, force](const size_t idx, auto& ac) {
               mesa_pd::addForceAtomic(
                  idx, ac,
                  ac.getVolume(idx) *
                     Vector3< real_t >(force[0], force[1], (particleDensityRatio - real_t(1)) * force[2]));
            },
            *particleAccessor);

         reduceProperty.operator()< mesa_pd::ForceTorqueNotification >(*particleStorage);

         particleStorage->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *particleAccessor,
                                          vvIntegratorPostForce, *particleAccessor);
         syncCall();
      }

      // has to be evaluated here before the force info is erased from particles
      if (t % evaluationFrequency == uint_c(0))
         totalHydrodynamicForceOnParticles = getTotalHydrodynamicForceOnParticles(particleAccessor);

      particleStorage->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *particleAccessor,
                                       resetHydrodynamicForceTorque, *particleAccessor);
      timingPool["Mesa_pd"].end();

      // update particle mapping
      timeloopAfterParticles.singleStep(timingPool, true);

      timingPool["Evaluation"].start();

      if (t % evaluationFrequency == uint_c(0))
      {
         averageDataSliceEvaluator();
         totalFluidMassEvaluator.computeMass(blockForest, antidunesBoundaryHandling);
         bedloadTransportEvaluator();
         meanVelocityComputer();

         WALBERLA_ROOT_SECTION()
         {
            write2DVectorToFile(averageDataSliceEvaluator.getSolidVolumeFractionVector(),
                                averageDataSliceEvaluator.getXLen(), averageDataSliceEvaluator.getZLen(),
                                baseFolderName + "/svfSlice_" + std::to_string(t) + ".txt");
            write2DVectorToFile(averageDataSliceEvaluator.getFillLevelVector(), averageDataSliceEvaluator.getXLen(),
                                averageDataSliceEvaluator.getZLen(),
                                baseFolderName + "/fillSlice_" + std::to_string(t) + ".txt");
            write2DVectorToFile(averageDataSliceEvaluator.getVelocityXVector(), averageDataSliceEvaluator.getXLen(),
                                averageDataSliceEvaluator.getZLen(),
                                baseFolderName + "/velXSlice_" + std::to_string(t) + ".txt");

            std::ofstream bedloadFile(bedLoadTransportFileName, std::ofstream::app);
            bedloadFile << t << " " << bedloadTransportEvaluator.getTransportRate() << " "
                        << bedloadTransportEvaluator.getAverageVelocity() << " " << totalHydrodynamicForceOnParticles[0]
                        << " " << totalHydrodynamicForceOnParticles[1] << " " << totalHydrodynamicForceOnParticles[2]
                        << "\n";
            bedloadFile.close();

            WALBERLA_LOG_DEVEL("____________________________________________________________________");
            WALBERLA_LOG_DEVEL("time step = " << t);
            const real_t froudeNumber = (*meanVelocity)[0] / real_c(std::sqrt(liquidHeight * std::abs(force[2])));

            const real_t reynoldsNumber = (*meanVelocity)[0] * liquidHeight / viscosity;

            const real_t weberNumber =
               real_t(1.0) * (*meanVelocity)[0] * (*meanVelocity)[0] * liquidHeight / surfaceTension;

            WALBERLA_LOG_DEVEL(" - Total fluid mass = " << std::setprecision(16) << (*totalFluidMass));
            auto maxFluidZPos = averageDataSliceEvaluator.getMaxFluidZPos();
            WALBERLA_LOG_DEVEL(" - Max fluid z-position = " << maxFluidZPos);
            WALBERLA_LOG_DEVEL(" - Froude number = " << froudeNumber);
            WALBERLA_LOG_DEVEL(" - Reynolds number = " << reynoldsNumber);
            WALBERLA_LOG_DEVEL(" - We = " << weberNumber);

            WALBERLA_LOG_DEVEL(" - meanVelocity = " << *meanVelocity);

            std::ofstream fluidInfoFile(fluidInfoFileName, std::ofstream::app);
            fluidInfoFile << t << " " << force[0] << " " << (*meanVelocity)[0] << " " << maxFluidZPos << " "
                          << std::setprecision(16) << (*totalFluidMass) << "\n";
            fluidInfoFile.close();

            if (std::isnan(reynoldsNumber)) WALBERLA_ABORT("reynoldsNumber is inf!")
         }

         WALBERLA_LOG_DEVEL_ON_ROOT(" -> CurrentExternalForce in x-direction before update = " << force[0]);
         (*forcingAdjuster)(meanVelocity->length());
         force[0] = forcingAdjuster->getExternalForcing();
         WALBERLA_LOG_DEVEL_ON_ROOT(" -> CurrentExternalForce in x-direction after update  = " << force[0]);
      }
      timingPool["Evaluation"].end();

      if (storeSnapshot)
      {
         if (t % snapshotFrequency == uint_c(0) && t > uint_c(0))
         {
            WALBERLA_LOG_INFO_ON_ROOT("Writing checkpointing file in time step " << t)

            blockForest->saveBlockData(snapshotBaseFolder + "/tmp_" + pdfFieldFile, pdfFieldID);
            blockForest->saveBlockData(snapshotBaseFolder + "/tmp_" + fillFieldFile, fillFieldID);
            blockForest->saveBlockData(snapshotBaseFolder + "/tmp_" + particleStorageFile, particleStorageID);

            WALBERLA_ROOT_SECTION()
            {
               std::string tmpCheckpointConfigFile = snapshotBaseFolder + "/tmp_" + checkpointConfigFile;
               std::ofstream file;
               file.open(tmpCheckpointConfigFile.c_str());

               file << std::setprecision(16);
               file << t + 1 << "\n";
               file << force[0] << "\n";
               file.close();
            }

            forcingAdjuster->storePIDSnapshot(snapshotBaseFolder + "/" + "pidState.file");

            WALBERLA_MPI_BARRIER();

            // rename checkpoint files to "real" ones
            // otherwise, the checkpointed state might be incomplete if the simulation stops due to over time during
            // checkpointing
            WALBERLA_ROOT_SECTION()
            {
               renameFile(snapshotBaseFolder + "/tmp_" + pdfFieldFile, snapshotBaseFolder + "/" + pdfFieldFile);
               renameFile(snapshotBaseFolder + "/tmp_" + fillFieldFile, snapshotBaseFolder + "/" + fillFieldFile);
               renameFile(snapshotBaseFolder + "/tmp_" + particleStorageFile,
                          snapshotBaseFolder + "/" + particleStorageFile);
               renameFile(snapshotBaseFolder + "/tmp_" + checkpointConfigFile,
                          snapshotBaseFolder + "/" + checkpointConfigFile);
            }
         }
      }

      if (t % performanceLogFrequency == uint_c(0) && t > uint_c(0)) { timingPool.logResultOnRoot(); }
   }

   return EXIT_SUCCESS;
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv) { return walberla::antidunes::main(argc, argv); }
