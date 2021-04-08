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
//! \file SphereWallCollision.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/waLBerlaBuildInfo.h"
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/logging/all.h"
#include "core/mpi/Broadcast.h"
#include "core/mpi/Reduce.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "field/interpolators/FieldInterpolatorCreators.h"
#include "field/interpolators/NearestNeighborFieldInterpolator.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/ExtrapolationDirectionFinder.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/OmegaBulkAdaption.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/kernel/NonLinearSpringDashpot.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ClearNextNeighborSync.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <functional>

#ifdef WALBERLA_BUILD_WITH_CODEGEN
#include "GeneratedLBM.h"
#endif

namespace sphere_wall_collision
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

#ifdef WALBERLA_BUILD_WITH_CODEGEN
using LatticeModel_T = lbm::GeneratedLBM;
#else
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::D3Q19MRT>;
#endif

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

using ScalarField_T = GhostLayerField< real_t, 1>;

const uint_t FieldGhostLayers = 1;

#define USE_TRTLIKE
//#define USE_D3Q27TRTLIKE
//#define USE_CUMULANT

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID NoSlip_Flag( "no slip" );
const FlagUID MO_Flag( "moving obstacle" );
const FlagUID FormerMO_Flag( "former moving obstacle" );
const FlagUID SimplePressure_Flag("simple pressure");

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
template <typename ParticleAccessor_T>
class MyBoundaryHandling
{
public:

   using NoSlip_T = lbm::NoSlip< LatticeModel_T, flag_t >;
   using MO_T = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using SimplePressure_T = lbm::SimplePressure< LatticeModel_T, flag_t >;
   using Type = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, MO_T, SimplePressure_T >;

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID,
                       const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor_T>& ac,
                       bool applyOutflowBCAtTop) :
         flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), particleFieldID_( particleFieldID ), ac_( ac ), applyOutflowBCAtTop_(applyOutflowBCAtTop) {}

   Type * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( storage );

      auto * flagField     = block->getData< FlagField_T >( flagFieldID_ );
      auto *  pdfField     = block->getData< PdfField_T > ( pdfFieldID_ );
      auto * particleField = block->getData< lbm_mesapd_coupling::ParticleField_T > ( particleFieldID_ );

      const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

      Type * handling = new Type( "moving obstacle boundary handling", flagField, fluid,
                                  NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                  MO_T( "MO", MO_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ),
                                  SimplePressure_T( "SimplePressure", SimplePressure_Flag, pdfField, real_t(1) ) );

      if(applyOutflowBCAtTop_)
      {
         const auto simplePressure = flagField->getFlag( SimplePressure_Flag );

         CellInterval domainBB = storage->getDomainCellBB();
         domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
         domainBB.xMax() += cell_idx_c( FieldGhostLayers );

         domainBB.yMin() -= cell_idx_c( FieldGhostLayers );
         domainBB.yMax() += cell_idx_c( FieldGhostLayers );

         domainBB.zMin() -= cell_idx_c( FieldGhostLayers );
         domainBB.zMax() += cell_idx_c( FieldGhostLayers );

         // TOP
         CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
         storage->transformGlobalToBlockLocalCellInterval( top, *block );
         handling->forceBoundary( simplePressure, top );
      }

      handling->fillWithDomain( FieldGhostLayers );

      return handling;
   }

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;

   bool applyOutflowBCAtTop_;
};
//*******************************************************************************************************************

//*******************************************************************************************************************
/*!\brief Evaluating the position and velocity of the sphere
 *
 */
//*******************************************************************************************************************
template< typename ParticleAccessor_T>
class SpherePropertyLogger
{
public:
   SpherePropertyLogger( const shared_ptr< ParticleAccessor_T > & ac, walberla::id_t sphereUid,
                         const std::string & fileName, bool fileIO,
                         real_t diameter, real_t uIn, uint_t numberOfSubSteps) :
      ac_( ac ), sphereUid_( sphereUid ), fileName_( fileName ), fileIO_(fileIO),
      diameter_( diameter ), uIn_( uIn ), numberOfSubSteps_(numberOfSubSteps),
      gap_( real_t(0) ), settlingVelocity_( real_t(0) )
   {
      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str() );
            file << "#\t t/tref\t gapSize\t gapSize/D\t velZ\t velZ/uIn\n";
            file.close();
         }
      }
   }

   void operator()(const uint_t timeStep, const uint_t subStep )
   {
      real_t curTimestep = real_c(timeStep) + real_c(subStep) / real_c(numberOfSubSteps_);

      Vector3<real_t> pos(real_t(0));
      Vector3<real_t> transVel(real_t(0));

      size_t idx = ac_->uidToIdx(sphereUid_);
      if( idx != std::numeric_limits<size_t>::max())
      {
         if(!mesa_pd::data::particle_flags::isSet( ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST))
         {
            pos = ac_->getPosition(idx);
            transVel = ac_->getLinearVelocity(idx);
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( pos[0], mpi::SUM );
         mpi::allReduceInplace( pos[1], mpi::SUM );
         mpi::allReduceInplace( pos[2], mpi::SUM );
         mpi::allReduceInplace( transVel[2], mpi::SUM );
      }

      position_ = pos;
      gap_ = pos[2] - real_t(0.5) * diameter_;
      settlingVelocity_ = transVel[2];

      WALBERLA_ROOT_SECTION()
      {
         // coarsen the output to the result file
         if(subStep == 0)
         {
            gapLogging_.push_back(gap_);
            velocityLogging_.push_back(settlingVelocity_);
         }
      }

      if( fileIO_ )
         writeToFile( curTimestep, gap_, settlingVelocity_);
   }

   real_t getGapSize() const
   {
      return gap_;
   }

   real_t getSettlingVelocity() const
   {
      return settlingVelocity_;
   }

   Vector3<real_t> getPosition() const
   {
      return position_;
   }

   void writeResult(const std::string & filename, uint_t tImpact)
   {
      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_CHECK_EQUAL(gapLogging_.size(), velocityLogging_.size());
         std::ofstream file;
         file.open( filename.c_str());

         file << "#\t t\t t/tref\t gapSize\t gapSize/D\t velZ\t velZ/uIn\n";
         real_t tref = diameter_ / uIn_;

         for( uint_t i = uint_t(0); i < gapLogging_.size(); ++i)
         {
            int timestep = int(i) - int(tImpact);
            file << timestep << "\t" << real_c(timestep) / tref<< "\t"
                 << "\t" << gapLogging_[i] << "\t" << gapLogging_[i] / diameter_
                 << "\t" << velocityLogging_[i] << "\t" << velocityLogging_[i] / uIn_
                 << "\n";
         }

         file.close();
      }
   }


private:
   void writeToFile( real_t timestep, real_t gap, real_t velocity )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName_.c_str(), std::ofstream::app );

         real_t tref = diameter_ / uIn_;

         file << timestep << "\t" << timestep / tref<< "\t"
              << "\t" << gap << "\t" << gap / diameter_
              << "\t" << velocity << "\t" << velocity / uIn_
              << "\n";
         file.close();
      }
   }


   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereUid_;
   std::string fileName_;
   bool fileIO_;
   real_t diameter_, uIn_;
   uint_t numberOfSubSteps_;

   real_t gap_, settlingVelocity_;
   Vector3<real_t> position_;

   std::vector<real_t> gapLogging_, velocityLogging_;

};


void createPlaneSetup(const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss, const math::AABB & simulationDomain, bool applyOutflowBCAtTop )
{
   // create bounding planes
   mesa_pd::data::Particle&& p0 = *ps->create(true);
   p0.setPosition(simulationDomain.minCorner());
   p0.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p0.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,0,1) ));
   p0.setOwner(mpi::MPIManager::instance()->rank());
   p0.setType(0);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   if(!applyOutflowBCAtTop)
   {
      //only create top plane when no outflow BC should be set there
      mesa_pd::data::Particle&& p1 = *ps->create(true);
      p1.setPosition(simulationDomain.maxCorner());
      p1.setInteractionRadius(std::numeric_limits<real_t>::infinity());
      p1.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,0,-1) ));
      p1.setOwner(mpi::MPIManager::instance()->rank());
      p1.setType(0);
      mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
      mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   }

}

struct NonCollidingSphereSelector
{
   explicit NonCollidingSphereSelector( real_t radius ) : radius_(radius) {}

   template< typename ParticleAccessor_T >
   bool inline operator()(const size_t particleIdx, const ParticleAccessor_T & ac) const
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");
      return !mesa_pd::data::particle_flags::isSet( ac.getFlags(particleIdx), mesa_pd::data::particle_flags::FIXED) &&
             !mesa_pd::data::particle_flags::isSet( ac.getFlags(particleIdx), mesa_pd::data::particle_flags::GLOBAL) &&
             ac.getPosition(particleIdx)[2] > radius_;
   }
private:
   real_t radius_;
};


class ForcesOnSphereLogger
{
public:
   ForcesOnSphereLogger( const std::string & fileName, bool fileIO, real_t weightForce, uint_t numberOfSubSteps ) :
         fileName_( fileName ), fileIO_( fileIO ), weightForce_( weightForce ), numberOfSubSteps_( numberOfSubSteps )
   {
      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str() );
            file << "#t\t BuoyAndGravF\t HydInteractionF\t LubricationF\t CollisionForce\n";
            file.close();
         }
      }
   }

   void operator()(const uint_t timeStep, const uint_t subStep, real_t hydForce, real_t lubForce, real_t collisionForce )
   {
      real_t curTimestep = real_c(timeStep) + real_c(subStep) / real_c(numberOfSubSteps_);

      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str(), std::ofstream::app );

            file << std::setprecision(10) << curTimestep << "\t" << weightForce_ << "\t" << hydForce << "\t" << lubForce << "\t" << collisionForce << "\n";
            file.close();
         }
      }

   }


private:

   std::string fileName_;
   bool fileIO_;

   real_t weightForce_;
   uint_t numberOfSubSteps_;
};

template<typename ParticleAccessor_T>
real_t getForce(walberla::id_t uid, ParticleAccessor_T & ac)
{
   auto idx = ac.uidToIdx(uid);
   real_t force = real_t(0);
   if(idx != ac.getInvalidIdx())
   {
      force = ac.getForce(idx)[2];
   }
   WALBERLA_MPI_SECTION()
   {
      mpi::allReduceInplace( force, mpi::SUM );
   }
   return force;
}

template< typename DensityInterpolator_T>
void logDensityToFile(const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID densityInterpolatorID, real_t surfaceDistance,
                      Vector3<real_t> spherePosition, real_t diameter, const std::string & fileName, uint_t timestep, real_t averageDensityInSystem )
{
   real_t evalRadius = real_t(0.5) * diameter + surfaceDistance;
   std::vector< Vector3<real_t> > evaluationPositions;
   evaluationPositions.push_back(spherePosition + Vector3<real_t>(real_t(0), real_t(0), evalRadius));
   evaluationPositions.push_back(spherePosition + Vector3<real_t>(evalRadius / real_c(sqrt(real_t(2))), real_t(0), real_c(evalRadius / sqrt(real_t(2)))));
   evaluationPositions.push_back(spherePosition + Vector3<real_t>(evalRadius, real_t(0), real_t(0)));
   evaluationPositions.push_back(spherePosition + Vector3<real_t>(evalRadius / real_c(sqrt(real_t(2))), real_t(0), -evalRadius / real_c(sqrt(real_t(2)))));
   evaluationPositions.push_back(spherePosition + Vector3<real_t>(real_t(0), real_t(0), -evalRadius));

   std::vector<real_t> densityAtPos(evaluationPositions.size(), real_t(0));

   for( auto & block: *blocks)
   {
      auto densityInterpolator = block.getData<DensityInterpolator_T>(densityInterpolatorID);
      for(auto i = uint_t(0); i < evaluationPositions.size(); ++i)
      {
         auto pos = evaluationPositions[i];
         if(block.getAABB().contains(pos))
         {
            densityInterpolator->get(pos, &densityAtPos[i]);
         }
      }
   }

   for(auto i = uint_t(0); i < densityAtPos.size(); ++i)
   {
      // reduce to root
      mpi::reduceInplace(densityAtPos[i], mpi::SUM);
   }

   WALBERLA_ROOT_SECTION()
   {
      std::ofstream file;
      file.open( fileName.c_str(), std::ofstream::app );

      file << std::setprecision(8) << timestep << "\t" << averageDensityInSystem;
      for(auto density: densityAtPos)
      {
         file << std::setprecision(8) << "\t" << density;
      }
      file << "\n";
      file.close();
   }
}

template< typename BoundaryHandling_T>
real_t getAverageDensityInSystem(const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID pdfFieldID, BlockDataID boundaryHandlingID)
{
   real_t totalMass = real_t(0);
   uint_t count = uint_t(0);
   for( auto & block: *blocks)
   {
      auto pdfField = block.getData<PdfField_T>(pdfFieldID);
      auto boundaryHandling = block.getData<BoundaryHandling_T>(boundaryHandlingID);
      WALBERLA_FOR_ALL_CELLS_XYZ(pdfField,
                                 if(boundaryHandling->isDomain(x,y,z))
                                 {
                                    totalMass += pdfField->getDensity(x,y,z);
                                    ++count;
                                 }
      );
   }

   // reduce to root
   mpi::reduceInplace(totalMass, mpi::SUM);
   mpi::reduceInplace(count, mpi::SUM);

   return totalMass/real_c(count);
}


//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief PHYSICAL test case of a sphere settling and colliding with a wall submerged in a viscous fluid.
 *
 * The trajectory of the bouncing sphere is logged and compared to reference experiments.
 *
 * for experiments see: Gondret, Lance, Petit - "Bouncing motion of spherical particles in fluids" 2002
 * for simulations see e.g.: Biegert, Vowinckel, Meiburg - "A collision model for grain-resolving simulations of flows over
 * dense, mobile, polydisperse granular sediment beds" 2017
 */
//*******************************************************************************************************************


int main( int argc, char **argv )
{
   Environment env( argc, argv );

   if( !env.config() )
   {
      WALBERLA_ABORT("Usage: " << argv[0] << " path-to-configuration-file \n");
   }

   Config::BlockHandle simulationInput = env.config()->getBlock( "BouncingSphere" );

   // setup description from separate block
   const std::string setup = simulationInput.getParameter<std::string>("setup");
   Config::BlockHandle setupDescription = env.config()->getBlock( setup );

   const real_t Re = setupDescription.getParameter<real_t>("Re");
   const real_t densityFluid_SI = setupDescription.getParameter<real_t>("densityFluid_SI");
   const real_t dynamicViscosityFluid_SI = setupDescription.getParameter<real_t>("dynamicViscosityFluid_SI");
   const real_t densitySphere_SI = setupDescription.getParameter<real_t>("densitySphere_SI");
   const real_t diameter_SI = setupDescription.getParameter<real_t>("diameter_SI");
   const real_t gravitationalAcceleration_SI = setupDescription.getParameter<real_t>("gravitationalAcceleration_SI");
   const real_t restitutionCoeff = setupDescription.getParameter<real_t>("restitutionCoeff");
   const uint_t numberOfBouncesUntilTermination = setupDescription.getParameter<uint_t>("numberOfBouncesUntilTermination");
   const Vector3<real_t> domainSizeNonDim = setupDescription.getParameter<Vector3<real_t> >("domainSize");
   const Vector3<uint_t> numberOfBlocksPerDirection = setupDescription.getParameter<Vector3<uint_t> >("numberOfBlocksPerDirection");

   // numerical parameters

   // control
   const bool randomizeInitialSpherePosition = simulationInput.getParameter<bool>("randomizeInitialSpherePosition");
   const bool initializeSphereVelocity = simulationInput.getParameter<bool>("initializeSphereVelocity");
   const bool applyOutflowBCAtTop = simulationInput.getParameter<bool>("applyOutflowBCAtTop");
   bool artificiallyAccelerateSphere = simulationInput.getParameter<bool>("artificiallyAccelerateSphere");

   // LBM
   const real_t uIn = simulationInput.getParameter<real_t>("uIn");
   const real_t diameter = simulationInput.getParameter<real_t>("diameter");
   const real_t magicNumber = simulationInput.getParameter<real_t>("magicNumber");
   const real_t bulkViscRateFactor = simulationInput.getParameter<real_t>("bulkViscRateFactor");
   const bool useOmegaBulkAdaption = simulationInput.getParameter<bool>("useOmegaBulkAdaption");

   const uint_t numRPDSubCycles = simulationInput.getParameter<uint_t>("numRPDSubCycles");
   const bool useLubricationCorrection = simulationInput.getParameter<bool>("useLubricationCorrection");
   const real_t lubricationCutOffDistance = simulationInput.getParameter<real_t>("lubricationCutOffDistance");
   const real_t lubricationMinimalGapSizeNonDim = simulationInput.getParameter<real_t>("lubricationMinimalGapSizeNonDim");
   const bool useVelocityVerlet = simulationInput.getParameter<bool>("useVelocityVerlet");

   // Collision Response
   const real_t collisionTime = simulationInput.getParameter<real_t>("collisionTime");
   const real_t StCrit = simulationInput.getParameter<real_t>("StCrit");
   const bool useACTM = simulationInput.getParameter<bool>("useACTM");

   const bool averageForceTorqueOverTwoTimeSteps = simulationInput.getParameter<bool>("averageForceTorqueOverTwoTimeSteps");
   const bool conserveMomentum = simulationInput.getParameter<bool>("conserveMomentum");
   const bool disableFluidForceDuringContact = simulationInput.getParameter<bool>("disableFluidForceDuringContact");
   const std::string reconstructorType = simulationInput.getParameter<std::string>("reconstructorType");

   const bool fileIO = simulationInput.getParameter<bool>("fileIO");
   const uint_t vtkIOFreq = simulationInput.getParameter<uint_t>("vtkIOFreq");
   const std::string baseFolder = simulationInput.getParameter<std::string>("baseFolder");
   bool vtkOutputOnCollision = simulationInput.getParameter<bool>("vtkOutputOnCollision");

   bool writeCheckPointFile = simulationInput.getParameter<bool>("writeCheckPointFile", false);
   bool readFromCheckPointFile = simulationInput.getParameter<bool>("readFromCheckPointFile", false);

   WALBERLA_ROOT_SECTION()
   {
      if( fileIO )
      {
         // create base directory if it does not yet exist
         filesystem::path tpath( baseFolder );
         if( !filesystem::exists( tpath ) )
            filesystem::create_directory( tpath );
      }
   }


   //////////////////////////////////////
   // SIMULATION PROPERTIES in SI units//
   //////////////////////////////////////

   const real_t densityRatio = densitySphere_SI / densityFluid_SI;
   const real_t kinematicViscosityFluid_SI = dynamicViscosityFluid_SI / densityFluid_SI;
   // here, we do the analogy with experiments via the Galileo number which is a-priori well-defined in contrast to the Stokes number which is a result of the simulation
   // -> no, this is not a good idea as the sphere in the experiments has in some cases not yet reached the terminal settling velocity
   // thus, it is better to pose the reported resulting Re as an input parameter
   const real_t uIn_SI = Re * kinematicViscosityFluid_SI / diameter_SI;
   const real_t ug_SI = std::sqrt((densityRatio-real_t(1)) * gravitationalAcceleration_SI * diameter_SI);
   const real_t Ga = ug_SI * diameter_SI / kinematicViscosityFluid_SI;

   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////

   const real_t dx_SI = diameter_SI / diameter;
   const real_t sphereVolume = math::pi / real_t(6) * diameter * diameter * diameter;
   const real_t dt_SI = uIn / uIn_SI * dx_SI;
   const real_t viscosity =  kinematicViscosityFluid_SI * dt_SI / ( dx_SI * dx_SI );
   const real_t gravitationalAcceleration = gravitationalAcceleration_SI * dt_SI * dt_SI / dx_SI;
   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);

   const real_t Re_p = diameter * uIn / viscosity;
   const real_t St = densityRatio * uIn * diameter / ( real_t(9) * viscosity );

   const real_t omegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, bulkViscRateFactor);

   const real_t densityFluid = real_t(1);
   const real_t densitySphere = densityRatio;

   const real_t dx = real_t(1);

   const real_t responseTime = densityRatio * diameter * diameter / ( real_t(18) * viscosity );
   const real_t accelerationFactor = real_t(1) / (real_t(0.1) * responseTime);
   const real_t tref = diameter / uIn;

   const real_t particleMass = densityRatio * sphereVolume;
   const real_t Mij = particleMass; // * particleMass / ( real_t(2) * particleMass ); // Mij = M for sphere-wall collision

   const real_t uInCrit = real_t(9) * StCrit * viscosity / ( densityRatio * diameter);

   Vector3<uint_t> domainSize( uint_c(domainSizeNonDim[0] * diameter ), uint_c(domainSizeNonDim[1] * diameter ), uint_c(domainSizeNonDim[2] * diameter ) );

   real_t initialSpherePosition = real_c(domainSize[2]) - real_t(1.5) * diameter;

   if( randomizeInitialSpherePosition )
   {
      uint_t seed1 = uint_c(std::chrono::system_clock::now().time_since_epoch().count());
      mpi::broadcastObject(seed1); // root process chooses seed and broadcasts it
      std::mt19937 g1 (static_cast<unsigned int>(seed1));

      real_t initialPositionOffset = real_t(0.5)*math::realRandom<real_t>(real_t(-1), real_t(1), g1);
      initialSpherePosition += initialPositionOffset;
   }


   WALBERLA_LOG_INFO_ON_ROOT("Setup (in SI units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter_SI << ", densityRatio = " << densityRatio);
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid_SI << ", dyn. visc = " << dynamicViscosityFluid_SI << ", kin. visc = " << kinematicViscosityFluid_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - ug = " << ug_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - Galileo number = " << Ga );
   WALBERLA_LOG_INFO_ON_ROOT(" - target Reynolds number = " << Re );
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = <" << real_c(domainSizeNonDim[0]) * diameter_SI << "," << real_c(domainSizeNonDim[1]) * diameter_SI  << ","<< real_c(domainSizeNonDim[2]) * diameter_SI  << ">");
   WALBERLA_LOG_INFO_ON_ROOT(" - dx = " << dx_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - dt = " << dt_SI);

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   if(applyOutflowBCAtTop) WALBERLA_LOG_INFO_ON_ROOT(" - outflow BC at top");
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter << ", density = " << densitySphere );
   WALBERLA_LOG_INFO_ON_ROOT(" - initial sphere position = " << initialSpherePosition );
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid << ", relaxation time (tau) = " << real_t(1) / omega << ", omega = " << omega << ", kin. visc = " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - magic number = " << magicNumber);
   WALBERLA_LOG_INFO_ON_ROOT(" - omega bulk = " << omegaBulk << ", bulk visc = " << lbm_mesapd_coupling::bulkViscosityFromOmegaBulk(omegaBulk) << " ( bulk visc rate factor (conste) = " << bulkViscRateFactor << ")");
   if(useOmegaBulkAdaption) WALBERLA_LOG_INFO_ON_ROOT(" - using omega bulk adaption in vicinity of sphere");
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration );
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = " << uIn );
   WALBERLA_LOG_INFO_ON_ROOT(" - target Reynolds number = " << Re_p );
   WALBERLA_LOG_INFO_ON_ROOT(" - target Stokes number = " << St );
   WALBERLA_LOG_INFO_ON_ROOT(" - tref = " << tref );
   WALBERLA_LOG_INFO_ON_ROOT(" - Stokes response time = " << responseTime );
   if( artificiallyAccelerateSphere )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - artificially accelerating sphere with factor " << accelerationFactor)
   }
   WALBERLA_LOG_INFO_ON_ROOT(" - integrator = " << (useVelocityVerlet ? "Velocity Verlet" : "Explicit Euler") );
   if( vtkIOFreq > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files to folder \"" << baseFolder << "\" with frequency " << vtkIOFreq);
   }
   if(vtkOutputOnCollision)
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files during collision events");
   }

   WALBERLA_LOG_INFO_ON_ROOT("Collision Response properties:" );
   WALBERLA_LOG_INFO_ON_ROOT(" - collision time = " << collisionTime );
   if(useACTM)
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using nonlinear collision model with ACTM" );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using linear collision model with fixed parameters" );
   }
   WALBERLA_LOG_INFO_ON_ROOT(" - coeff of restitution = " << restitutionCoeff );
   WALBERLA_LOG_INFO_ON_ROOT(" - number of RPD sub cycles = " << numRPDSubCycles );
   WALBERLA_LOG_INFO_ON_ROOT(" - lubrication correction = " << (useLubricationCorrection ? "yes" : "no") );
   if( useLubricationCorrection )
   {
      WALBERLA_LOG_INFO_ON_ROOT("  - lubrication correction cut off = " << lubricationCutOffDistance );
      WALBERLA_LOG_INFO_ON_ROOT("  - lubrication correction minimal gap size non dim = " << lubricationMinimalGapSizeNonDim  << " ( = " << lubricationMinimalGapSizeNonDim * diameter * real_t(0.5) << " cells )");
   }
   WALBERLA_LOG_INFO_ON_ROOT("Coupling properties:" );
   WALBERLA_LOG_INFO_ON_ROOT(" - disable hydrodynamic forces during contact = " << (disableFluidForceDuringContact ? "yes" : "no") );
   if(disableFluidForceDuringContact) WALBERLA_LOG_INFO_ON_ROOT("  - StCrit = " << StCrit << ", uInCrit = " << uInCrit );
   WALBERLA_LOG_INFO_ON_ROOT(" - average forces over two LBM time steps = " << (averageForceTorqueOverTwoTimeSteps ? "yes" : "no") );
   WALBERLA_LOG_INFO_ON_ROOT(" - conserve momentum = " << (conserveMomentum ? "yes" : "no") );
   WALBERLA_LOG_INFO_ON_ROOT(" - reconstructor = " << reconstructorType );

   std::string checkPointFileName = "checkPointingFile_" + setup
                                    + "_uIn" + std::to_string(uIn)
                                    + "_d" + std::to_string(diameter)
                                    + "_bulkViscFactor" + std::to_string(bulkViscRateFactor)
                                    + "_mn" + std::to_string(magicNumber);
   if(useOmegaBulkAdaption) checkPointFileName +="_omegaBulkAdaption";
   if(applyOutflowBCAtTop) checkPointFileName +="_outflowBC";
#if defined(USE_TRTLIKE)
   checkPointFileName += "_trtlike";
#elif defined(USE_CUMULANT)
   checkPointFileName += "_cumulant";
#elif defined(USE_D3Q27TRTLIKE)
   checkPointFileName += "_d3q27trtlike";
#endif

   WALBERLA_LOG_INFO_ON_ROOT("Checkpointing:");
   WALBERLA_LOG_INFO_ON_ROOT(" - read from checkpoint file = " << (readFromCheckPointFile ? "yes" : "no"));
   WALBERLA_LOG_INFO_ON_ROOT(" - write checkpoint file = " << (writeCheckPointFile ? "yes" : "no"));
   if( readFromCheckPointFile || writeCheckPointFile) WALBERLA_LOG_INFO_ON_ROOT(" - checkPointingFileName = " << checkPointFileName );

   if(readFromCheckPointFile && writeCheckPointFile )
   {
      // decide which option to choose
      if(filesystem::exists( checkPointFileName+"_lbm.txt"))
      {
         WALBERLA_LOG_INFO_ON_ROOT("Checkpoint file already exists! Will skip writing check point file and start from this check point!");
         writeCheckPointFile = false;
      } else
      {
         WALBERLA_LOG_INFO_ON_ROOT("Checkpoint file does not exists yet! Will skip reading check point file and just regularly start the simulation from the beginning!");
         readFromCheckPointFile = false;
      }
   }

   if(readFromCheckPointFile && artificiallyAccelerateSphere)
   {
      // accelerating after reading from a checkpoint file does not make sense, as current actual runtime is not known
      WALBERLA_LOG_INFO_ON_ROOT("NOTE: switching off artificial acceleration of sphere due to reading from check point file");
      artificiallyAccelerateSphere = false;
   }


   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////


   Vector3<uint_t> cellsPerBlockPerDirection( domainSize[0] / numberOfBlocksPerDirection[0],
                                              domainSize[1] / numberOfBlocksPerDirection[1],
                                              domainSize[2] / numberOfBlocksPerDirection[2] );
   for( uint_t i = 0; i < 3; ++i ) {
      WALBERLA_CHECK_EQUAL(cellsPerBlockPerDirection[i] * numberOfBlocksPerDirection[i], domainSize[i],
                           "Unmatching domain decomposition in direction " << i << "!");
   }

   shared_ptr< StructuredBlockForest > blocks;
   if( readFromCheckPointFile )
   {
      WALBERLA_LOG_INFO_ON_ROOT("Reading block forest from file!");
      blocks = blockforest::createUniformBlockGrid(checkPointFileName+"_forest.txt", cellsPerBlockPerDirection[0], cellsPerBlockPerDirection[1], cellsPerBlockPerDirection[2], false);
   }
   else
   {
      blocks = blockforest::createUniformBlockGrid( numberOfBlocksPerDirection[0], numberOfBlocksPerDirection[1], numberOfBlocksPerDirection[2],
                                                      cellsPerBlockPerDirection[0], cellsPerBlockPerDirection[1], cellsPerBlockPerDirection[2], dx,
                                                      0, false, false,
                                                      true, true, false, //periodicity
                                                      false );
   
      if(writeCheckPointFile)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Writing block forest to file!");
         blocks->getBlockForest().saveToFile(checkPointFileName+"_forest.txt");
      }
   }
   
   WALBERLA_LOG_INFO_ON_ROOT("Domain partitioning:");
   WALBERLA_LOG_INFO_ON_ROOT(" - blocks per direction = " << numberOfBlocksPerDirection );
   WALBERLA_LOG_INFO_ON_ROOT(" - cells per block = " << cellsPerBlockPerDirection );

   //write domain decomposition to file
   if( vtkIOFreq > 0 )
   {
      vtk::writeDomainDecomposition( blocks, "initial_domain_decomposition", baseFolder );
   }


   //////////////////
   // RPD COUPLING //
   //////////////////

   auto rpdDomain = std::make_shared<mesa_pd::domain::BlockForestDomain>(blocks->getBlockForestPointer());

   //init data structures
   auto ps = walberla::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = walberla::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);
   BlockDataID particleStorageID;
   if(readFromCheckPointFile)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing particles from checkpointing file!" );
      particleStorageID = blocks->loadBlockData( checkPointFileName + "_mesa.txt", mesa_pd::domain::createBlockForestDataHandling(ps), "Particle Storage" );
      mesa_pd::mpi::ClearNextNeighborSync CNNS;
      CNNS(*accessor);
   } else {
      particleStorageID = blocks->addBlockData(mesa_pd::domain::createBlockForestDataHandling(ps), "Particle Storage");
   }

   // bounding planes
   createPlaneSetup(ps,ss,blocks->getDomain(), applyOutflowBCAtTop);

   // create sphere and store Uid
   Vector3<real_t> initialPosition( real_t(0.5) * real_c(domainSize[0]), real_t(0.5) * real_c(domainSize[1]), initialSpherePosition );
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( diameter * real_t(0.5) );
   ss->shapes[sphereShape]->updateMassAndInertia(densitySphere);

   walberla::id_t sphereUid = 0;
   if( readFromCheckPointFile )
   {
      for( auto pIt = ps->begin(); pIt != ps->end(); ++pIt )
      {
         // find sphere in loaded data structure and store uid for later reference
         if( pIt->getShapeID() == sphereShape )
         {
            sphereUid = pIt->getUid();
         }
      }
   } else
   {
      // create sphere
      if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), initialPosition ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(initialPosition);
         p.setInteractionRadius(diameter * real_t(0.5));
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         p.setType(0);
         if( initializeSphereVelocity ) p.setLinearVelocity(Vector3<real_t>(real_t(0), real_t(0), -real_t(0.1) * uIn));
         sphereUid = p.getUid();
      }
   }
   mpi::allReduceInplace(sphereUid, mpi::SUM);

   if(sphereUid == 0) WALBERLA_ABORT("No sphere present - aborting!"); //something went wrong in the checkpointing probably

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add omega bulk field
   BlockDataID omegaBulkFieldID = field::addToStorage<ScalarField_T>( blocks, "omega bulk field", omegaBulk, field::fzyx);

   // create the lattice model
   real_t lambda_e = lbm::collision_model::TRT::lambda_e( omega );
   real_t lambda_d = lbm::collision_model::TRT::lambda_d( omega, magicNumber );
#ifdef WALBERLA_BUILD_WITH_CODEGEN
#ifdef USE_CUMULANT
   WALBERLA_LOG_INFO_ON_ROOT("Using generated cumulant lattice model!");
   LatticeModel_T latticeModel = LatticeModel_T(omega);
#elif defined(USE_TRTLIKE) || defined(USE_D3Q27TRTLIKE)
   WALBERLA_LOG_INFO_ON_ROOT("Using generated TRT-like lattice model!");
   LatticeModel_T latticeModel = LatticeModel_T(omegaBulkFieldID, lambda_d, lambda_e);
#endif
#else
   WALBERLA_LOG_INFO_ON_ROOT("Using waLBerla built-in MRT lattice model and ignoring omega bulk field since not supported!");
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::D3Q19MRT( omegaBulk, omegaBulk, lambda_d, lambda_e, lambda_e, lambda_d ));
#endif

   // add PDF field
   BlockDataID pdfFieldID;
   if( readFromCheckPointFile )
   {
      // add PDF field
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing PDF Field from checkpointing file!" );
      shared_ptr< lbm::internal::PdfFieldHandling< LatticeModel_T > > dataHandling =
            make_shared< lbm::internal::PdfFieldHandling< LatticeModel_T > >( blocks, latticeModel, false,
                                                                              Vector3< real_t >( real_t(0) ), real_t(1),
                                                                              uint_t(1), field::fzyx );

      pdfFieldID = blocks->loadBlockData( checkPointFileName+"_lbm.txt", dataHandling, "pdf field" );

   } else {
      // add PDF field
      pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field", latticeModel,
                                                                Vector3< real_t >( real_t(0) ), real_t(1),
                                                                uint_t(1), field::fzyx );
   }


   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::fzyx, FieldGhostLayers );

   // add boundary handling
   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor, applyOutflowBCAtTop), "boundary handling" );

   // interpolation functionality
   using DensityAdaptor_T = typename lbm::Adaptor< LatticeModel_T >::Density;
   BlockDataID densityAdaptorID = field::addFieldAdaptor<  DensityAdaptor_T >( blocks, pdfFieldID, "density adaptor" );

   using DensityInterpolator_T = typename field::NearestNeighborFieldInterpolator<DensityAdaptor_T, FlagField_T>;
   BlockDataID densityInterpolatorID = field::addFieldInterpolator< DensityInterpolator_T, FlagField_T >( blocks, densityAdaptorID, flagFieldID, Fluid_Flag );

   // set up RPD functionality
   std::function<void(void)> syncCall = [&ps,&rpdDomain](){
      const real_t overlap = real_t( 1.5 );
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   };

   syncCall();

   real_t timeStepSizeRPD = real_t(1)/real_t(numRPDSubCycles);
   mesa_pd::kernel::ExplicitEuler explEulerIntegrator(timeStepSizeRPD);
   mesa_pd::kernel::VelocityVerletPreForceUpdate  vvIntegratorPreForce(timeStepSizeRPD);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(timeStepSizeRPD);

   // linear model
   mesa_pd::kernel::LinearSpringDashpot linearCollisionResponse(1);
   linearCollisionResponse.setStiffnessAndDamping(0,0,restitutionCoeff,collisionTime,real_t(0),Mij); // no response in tangential direction

   // nonlinear model for ACTM
   mesa_pd::kernel::NonLinearSpringDashpot nonLinearCollisionResponse(1, collisionTime);
   nonLinearCollisionResponse.setLnCORsqr(0,0,std::log(restitutionCoeff) * std::log(restitutionCoeff));
   nonLinearCollisionResponse.setMeff(0,0,Mij);

   mesa_pd::mpi::ReduceProperty reduceProperty;

   mesa_pd::mpi::ReduceContactHistory reduceAndSwapContactHistory;

   // set up coupling functionality
   Vector3<real_t> gravitationalForce( real_t(0), real_t(0), -(densitySphere - densityFluid) * gravitationalAcceleration * sphereVolume );
   lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel(viscosity, [lubricationMinimalGapSizeNonDim](real_t r){ return lubricationMinimalGapSizeNonDim * r;}, lubricationCutOffDistance, real_t(0), real_t(0) ); // no tangential components needed
   lbm_mesapd_coupling::ParticleMappingKernel<BoundaryHandling_T> particleMappingKernel(blocks, boundaryHandlingID);
   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);

   ///////////////
   // TIME LOOP //
   ///////////////

   // map planes into the LBM simulation -> act as no-slip boundaries
   ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), *accessor, particleMappingKernel, *accessor, NoSlip_Flag);

   // map particles into the LBM simulation
   ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   blockforest::communication::UniformBufferedScheme< Stencil_T > optimizedPDFCommunicationScheme( blocks );
   optimizedPDFCommunicationScheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) ); // optimized sync

   blockforest::communication::UniformBufferedScheme< Stencil_T > fullPDFCommunicationScheme( blocks );
   fullPDFCommunicationScheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) ); // full sync

   // create the timeloop
   const uint_t timesteps = uint_c( real_t(3) * initialPosition[2] / uIn);

   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   timeloop.addFuncBeforeTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

   // vtk output
   if( vtkIOFreq != uint_t(0) )
   {
      // sphere
      auto particleVtkOutput = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
      particleVtkOutput->addOutput<mesa_pd::data::SelectParticleOwner>("owner");
      particleVtkOutput->addOutput<mesa_pd::data::SelectParticleUid>("uid");
      particleVtkOutput->addOutput<mesa_pd::data::SelectParticleLinearVelocity>("velocity");
      particleVtkOutput->setParticleSelector( [sphereShape](const mesa_pd::data::ParticleStorage::iterator& pIt) {return pIt->getShapeID() == sphereShape;} ); //limit output to sphere
      auto particleVtkWriter = vtk::createVTKOutput_PointData(particleVtkOutput, "Particles", vtkIOFreq, baseFolder, "simulation_step");
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( particleVtkWriter ), "VTK (sphere data)" );

      // flag field (written only once in the first time step, ghost layers are also written)

      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", timesteps, FieldGhostLayers, false, baseFolder );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );


      // pdf field, as slice
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", vtkIOFreq, 0, false, baseFolder );

      pdfFieldVTK->addBeforeFunction( fullPDFCommunicationScheme );

      AABB sliceAABB( real_t(0), real_c(domainSize[1])*real_t(0.5)-real_t(1), real_t(0),
                      real_c(domainSize[0]), real_c(domainSize[1])*real_t(0.5)+real_t(1), real_c(domainSize[2]) );
      vtk::AABBCellFilter aabbSliceFilter( sliceAABB );

      field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
      fluidFilter.addFlag( Fluid_Flag );

      vtk::ChainedFilter combinedSliceFilter;
      combinedSliceFilter.addFilter( fluidFilter );
      combinedSliceFilter.addFilter( aabbSliceFilter );

      pdfFieldVTK->addCellInclusionFilter( combinedSliceFilter );

      pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( pdfFieldVTK ), "VTK (fluid field data)" );

      // omega bulk field
      //timeloop.addFuncBeforeTimeStep( field::createVTKOutput<ScalarField_T, float>( omegaBulkFieldID, *blocks, "omega_bulk_field", vtkIOFreq, uint_t(0), false, baseFolder ), "VTK (omega bulk field)" );
   }



   // update bulk omega in all cells to adapt to changed particle position
   if(useOmegaBulkAdaption)
   {
      using OmegaBulkAdapter_T = lbm_mesapd_coupling::OmegaBulkAdapter<ParticleAccessor_T, lbm_mesapd_coupling::RegularParticlesSelector>;
      real_t defaultOmegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, real_t(1));
      real_t adaptionLayerSize = real_t(2);
      shared_ptr<OmegaBulkAdapter_T> omegaBulkAdapter = make_shared<OmegaBulkAdapter_T>(blocks, omegaBulkFieldID, accessor, defaultOmegaBulk, omegaBulk, adaptionLayerSize, lbm_mesapd_coupling::RegularParticlesSelector());

      // initial adaption
      for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
         (*omegaBulkAdapter)(blockIt.get());
      }

      timeloop.add() << Sweep( makeSharedSweep(omegaBulkAdapter), "Omega Bulk Adapter");
   }


   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   timeloop.add() << BeforeFunction( optimizedPDFCommunicationScheme, "LBM Communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   // stream + collide LBM step
#ifdef WALBERLA_BUILD_WITH_CODEGEN
   auto lbmSweep = LatticeModel_T::Sweep( pdfFieldID );
   timeloop.add() << Sweep( lbmSweep, "LB sweep" );
#else
   auto lbmSweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   timeloop.add() << Sweep( makeSharedSweep( lbmSweep ), "cell-wise LB sweep" );
#endif

   // this is carried out after the particle integration, it corrects the flag field and restores missing PDF information
   // then, the checkpointing file can be written, as otherwise some cells are invalid and can not be recovered
   SweepTimeloop timeloopAfterParticles( blocks->getBlockStorage(), timesteps );

   // sweep for updating the particle mapping into the LBM simulation
   timeloopAfterParticles.add() << Sweep( lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID,boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag, lbm_mesapd_coupling::RegularParticlesSelector(), conserveMomentum), "Particle Mapping" );

   // sweep for restoring PDFs in cells previously occupied by particles
   if( reconstructorType == "EAN" )
   {
      auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(blocks);
      auto equilibriumAndNonEquilibriumSphereNormalReconstructor = lbm_mesapd_coupling::makeEquilibriumAndNonEquilibriumReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, uint_t(3), true);
      timeloopAfterParticles.add() << BeforeFunction( fullPDFCommunicationScheme, "LBM Communication" )
                                   << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, equilibriumAndNonEquilibriumSphereNormalReconstructor, conserveMomentum)), "PDF Restore" );
   }else if( reconstructorType == "Grad" )
   {
      bool recomputeTargetDensity = false;
      auto gradReconstructor = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, omega, recomputeTargetDensity,true);

      timeloopAfterParticles.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                                   << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, gradReconstructor, conserveMomentum) ), "PDF Restore" );
   }else if( reconstructorType == "Eq" )
   {
      timeloopAfterParticles.add() << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, conserveMomentum)), "PDF Restore" );
   }else if( reconstructorType == "Ext")
   {
      auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(blocks);
      auto extrapolationReconstructor = lbm_mesapd_coupling::makeExtrapolationReconstructor<BoundaryHandling_T, lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder, false>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, uint_t(3), true);
      timeloopAfterParticles.add() << BeforeFunction( fullPDFCommunicationScheme, "LBM Communication" )
                                   << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, extrapolationReconstructor, conserveMomentum)), "PDF Restore" );
   } else
   {
      WALBERLA_ABORT("Reconstructor Type " << reconstructorType << " not implemented!");
   }


   // evaluation functionality
   std::string loggingFileName( baseFolder + "/LoggingSphereWallCollision.txt");
   if( fileIO  )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing logging output to file \"" << loggingFileName << "\"");
   }
   SpherePropertyLogger<ParticleAccessor_T> logger( accessor, sphereUid, loggingFileName, fileIO, diameter, uIn, numRPDSubCycles );

   std::string forceLoggingFileName( baseFolder + "/ForceLoggingSphereWallCollision.txt");
   if( fileIO  )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing force logging output to file \"" << forceLoggingFileName << "\"");
   }
   ForcesOnSphereLogger sphereForceLogger(forceLoggingFileName, fileIO, gravitationalForce[2], numRPDSubCycles);

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   // evaluation quantities
   uint_t numBounces = uint_t(0);
   uint_t tImpact = uint_t(0);
   std::vector<uint_t> impactTimes;

   real_t curVel(real_t(0));
   real_t oldVel(real_t(0));
   real_t maxSettlingVel(real_t(0));
   std::vector<real_t> maxSettlingVelBetweenBounces;

   real_t minGapSize(real_t(0));
   real_t maxGapSize(real_t(0));
   std::vector<real_t> maxGapSizeBetweenBounces;

   real_t actualSt(real_t(0));
   real_t actualRe(real_t(0));

   WALBERLA_LOG_INFO_ON_ROOT("Running for maximum of " << timesteps << " timesteps!");
   
   const bool useOpenMP = false;

   // special vtk output
   auto pdfFieldVTKCollision = vtk::createVTKOutput_BlockData( blocks, "collision_fluid_field", 1, 1, false, baseFolder );

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) );
   pdfFieldVTKCollision->addBeforeFunction( pdfGhostLayerSync );

   AABB boxAABB( initialPosition[0] - diameter, initialPosition[1] - diameter, -real_t(1),
                 initialPosition[0] + diameter, initialPosition[1] + diameter, real_t(2) * diameter);
   vtk::AABBCellFilter aabbBoxFilter( boxAABB );

   field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
   fluidFilter.addFlag( Fluid_Flag );

   vtk::ChainedFilter combinedSliceFilter;
   combinedSliceFilter.addFilter( fluidFilter );
   combinedSliceFilter.addFilter( aabbBoxFilter );

   pdfFieldVTKCollision->addCellInclusionFilter( combinedSliceFilter );

   pdfFieldVTKCollision->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
   pdfFieldVTKCollision->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );


   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {
      // perform a single simulation step -> this contains LBM and setting of the hydrodynamic interactions
      timeloop.singleStep( timeloopTiming );
         
      reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);

      if( averageForceTorqueOverTwoTimeSteps )
      {
         if( i == 0 )
         {
            lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel initializeHydrodynamicForceTorqueForAveragingKernel;
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, initializeHydrodynamicForceTorqueForAveragingKernel, *accessor );
         }
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque, *accessor );
      }

      real_t hydForce(real_t(0));
      real_t lubForce(real_t(0));
      real_t collisionForce(real_t(0));

      for(auto subCycle = uint_t(0); subCycle < numRPDSubCycles; ++subCycle )
      {

         timeloopTiming["RPD"].start();

         if( useVelocityVerlet )
         {

            Vector3<real_t> oldForce;

            if(artificiallyAccelerateSphere)
            {
               // since the pre-force step of VV updates the position based on velocity and old force, we set oldForce to zero here for this step while artificially accelerating
               // to not perturb the step after which artificial acceleration is switched off (which requires valid oldForce values then) we store the old force and then re-apply it
               auto idx = accessor->uidToIdx(sphereUid);
               if(idx != accessor->getInvalidIdx())
               {
                  oldForce = accessor->getOldForce(idx);
                  accessor->setOldForce(idx, Vector3<real_t>(real_t(0)));
               }
            }

            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
            syncCall();

            if(artificiallyAccelerateSphere)
            {
               // re-apply old force
               auto idx = accessor->uidToIdx(sphereUid);
               if(idx != accessor->getInvalidIdx())
               {
                  accessor->setOldForce(idx, oldForce);
               }
            }
         }

         // lubrication correction
         ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                     [&lubricationCorrectionKernel,&rpdDomain](const size_t idx1, const size_t idx2, auto& ac)
                                     {
                                        mesa_pd::collision_detection::AnalyticContactDetection acd;
                                        acd.getContactThreshold() = lubricationCorrectionKernel.getNormalCutOffDistance();
                                        mesa_pd::kernel::DoubleCast double_cast;
                                        mesa_pd::mpi::ContactFilter contact_filter;
                                        if (double_cast(idx1, idx2, ac, acd, ac ))
                                        {
                                           if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                           {
                                              double_cast(acd.getIdx1(), acd.getIdx2(), ac, lubricationCorrectionKernel, ac, acd.getContactNormal(), acd.getPenetrationDepth());
                                           }
                                        }
                                     },
                                     *accessor );

         lubForce = getForce(sphereUid,*accessor);

         // one could add linked cells here

         // collision response
         ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                     [&linearCollisionResponse, &nonLinearCollisionResponse, &rpdDomain, timeStepSizeRPD, useACTM]
                                     (const size_t idx1, const size_t idx2, auto& ac)
                                     {
                                        mesa_pd::collision_detection::AnalyticContactDetection acd;
                                        mesa_pd::kernel::DoubleCast double_cast;
                                        mesa_pd::mpi::ContactFilter contact_filter;
                                        if (double_cast(idx1, idx2, ac, acd, ac ))
                                        {
                                           if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                           {
                                              if(useACTM) nonLinearCollisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth(), timeStepSizeRPD);
                                              else        linearCollisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth(), timeStepSizeRPD);
                                           }
                                        }
                                     },
                                     *accessor );

         collisionForce = getForce(sphereUid,*accessor) - lubForce;

         reduceAndSwapContactHistory(*ps);

         // add hydrodynamic force
         if( disableFluidForceDuringContact )
         {
            lbm_mesapd_coupling::StokesNumberBasedSphereSelector sphereSelector(StCrit, densityFluid, densitySphere, viscosity);
            lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
            ps->forEachParticle(useOpenMP, sphereSelector, *accessor, addHydrodynamicInteraction, *accessor );
         }
         else
         {
            lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction, *accessor );
         }

         hydForce = getForce(sphereUid,*accessor) - lubForce - collisionForce;

         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor );

         reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);

         // integration
         if( useVelocityVerlet ) ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
         else ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, explEulerIntegrator, *accessor);

         syncCall();

         if( artificiallyAccelerateSphere)
         {
            // overwrite velocity of particle with prescribed one
            lbm_mesapd_coupling::RegularParticlesSelector sphereSelector;
            real_t newSphereVel = uIn * (std::exp(-accelerationFactor * real_t(i) ) - real_t(1));
            ps->forEachParticle(useOpenMP, sphereSelector, *accessor, [newSphereVel](const size_t idx, ParticleAccessor_T& ac){ ac.setLinearVelocity(idx, Vector3<real_t>(real_t(0), real_t(0), newSphereVel));}, *accessor);
         }

         timeloopTiming["RPD"].end();

         // logging
         timeloopTiming["Logging"].start();
         logger(i, subCycle);
         sphereForceLogger(i, subCycle, hydForce, lubForce, collisionForce);
         timeloopTiming["Logging"].end();

      }

      ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );


      // update particle mapping
      timeloopAfterParticles.singleStep(timeloopTiming);

      // check for termination
      oldVel = curVel;
      curVel = logger.getSettlingVelocity();

      maxSettlingVel = std::min(maxSettlingVel, curVel);
      minGapSize = std::min(minGapSize, logger.getGapSize());

      if(numBounces >= uint_t(1))
      {
         maxGapSize = std::max(maxGapSize, logger.getGapSize());
      }

      // detect a bounce
      if( oldVel < real_t(0) && curVel > real_t(0) && logger.getGapSize() < real_t(1))
      {

         ++numBounces;
         WALBERLA_LOG_INFO_ON_ROOT("Detected bounce number " << numBounces);
         WALBERLA_LOG_INFO_ON_ROOT("Max settling velocity " << maxSettlingVel << " -> St = "  << densityRatio * std::abs(maxSettlingVel) * diameter / ( real_t(9) * viscosity ) );

         if(numBounces == uint_t(1))
         {
            actualSt = densityRatio * std::abs(maxSettlingVel) * diameter / ( real_t(9) * viscosity );
            actualRe = std::abs(maxSettlingVel) * diameter / viscosity;
         }

         // reset and store quantities
         maxSettlingVelBetweenBounces.push_back(maxSettlingVel);
         maxSettlingVel = real_t(0);

         if(numBounces > uint_t(1))
         {
            maxGapSizeBetweenBounces.push_back(maxGapSize);
            maxGapSize = real_t(0);
         }

         if( numBounces >= numberOfBouncesUntilTermination )
         {
            break;
         }
      }

      // impact times are measured when the contact between sphere and wall is broken up again
      if( tImpact == uint_t(0) && numBounces == uint_t(1) && logger.getGapSize() > real_t(0) )
      {
         tImpact = i;
         WALBERLA_LOG_INFO_ON_ROOT("Detected impact time at time step " << tImpact);

         // switch off special vtk output after first bounce
         if(vtkOutputOnCollision) vtkOutputOnCollision = false;
      }

      if( numBounces > impactTimes.size() && logger.getGapSize() > real_t(0))
      {
         impactTimes.push_back(i);
      }

      // evaluate density around sphere
      if(fileIO)
      {
         std::string densityLoggingFileName(baseFolder + "/DensityLogging.txt");
         if(i == uint_t(0))
         {
            WALBERLA_ROOT_SECTION()
            {
               std::ofstream file;
               file.open( densityLoggingFileName.c_str());
               file << "# t \t avgDensity \t densitiesAtDifferentPositions \n";
               file.close();
            }
         }

         real_t surfaceDistance = real_t(2);
         real_t avgDensity = getAverageDensityInSystem<BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID);
         logDensityToFile<DensityInterpolator_T>(blocks, densityInterpolatorID, surfaceDistance, logger.getPosition(), diameter, densityLoggingFileName, i, avgDensity);
      }


      // check if sphere is close to bottom plane
      if( logger.getGapSize() < real_t(1) * diameter ) {

         // write a single checkpointing file before any collision relevant parameters take effect
         if( writeCheckPointFile )
         {
            WALBERLA_LOG_INFO_ON_ROOT("Writing checkpointing file in time step " << timeloop.getCurrentTimeStep());

            blocks->saveBlockData( checkPointFileName + "_lbm.txt", pdfFieldID );
            blocks->saveBlockData( checkPointFileName + "_mesa.txt", particleStorageID );

            writeCheckPointFile = false;
         }

         // switch off acceleration after writing of check point file
         artificiallyAccelerateSphere = false;

      }

      if(vtkOutputOnCollision && logger.getGapSize() < real_t(0))
      {
         pdfFieldVTKCollision->write();
      }

   }

   WALBERLA_LOG_INFO_ON_ROOT("Detected " << numBounces << " bounces, terminating simulation.");
   WALBERLA_LOG_INFO_ON_ROOT("Maximum settling velocities: " );
   for( auto vel : maxSettlingVelBetweenBounces)
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - vel = " << vel << " -> St = "  << densityRatio * std::abs(vel) * diameter / ( real_t(9) * viscosity ) << ", Re = " << std::abs(vel) * diameter / viscosity );
   }
   WALBERLA_LOG_INFO_ON_ROOT("Maximum gap sizes between bounces:")
   for( auto gapSize : maxGapSizeBetweenBounces)
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - gap size = " << gapSize << ", gap size / diameter = "  << gapSize / diameter );
   }


   std::string resultFile(baseFolder + "/ResultBouncingSphere.txt");
   WALBERLA_LOG_INFO_ON_ROOT("Writing logging file " << resultFile);
   logger.writeResult(resultFile, tImpact);

   std::string summaryFile(baseFolder + "/Summary.txt");
   WALBERLA_LOG_INFO_ON_ROOT("Writing summary file " << summaryFile);
   WALBERLA_ROOT_SECTION()
   {
      std::ofstream file;
      file.open( summaryFile.c_str());

      file << "waLBerla Revision = " << std::string(WALBERLA_GIT_SHA1).substr(0,8) << "\n";
      file << "\nInput parameters:\n";
      file << "Ga = " << Ga << "\n";
      file << "uIn = " << uIn << "\n";
      file << "LBM parameters:\n";
      file << " - magic number = " << magicNumber << "\n";
      file << " - bulk viscosity rate factor = " << bulkViscRateFactor << "\n";
      file << " - use omega bulk adaption = " << useOmegaBulkAdaption << " (layer size = 2)\n";
      file << "Collision parameters:\n";
      file << " - subCycles = " << numRPDSubCycles << "\n";
      file << " - collision time (Tc) = " << collisionTime << "\n";
      file << " - use ACTM = " << useACTM << "\n";
      file << "use lubrication correction = " << useLubricationCorrection << "\n";
      file << " - minimum gap size non dim = " << lubricationMinimalGapSizeNonDim << "\n";
      file << " - cut off distance = " << lubricationCutOffDistance << "\n";
      file << "switch off fluid interaction force during contact = " << disableFluidForceDuringContact << "\n";
      file << " - St_Crit = " << StCrit << "\n";
      file << "reconstructor type " << reconstructorType << "\n";
      file << "apply outflow BC at top = " << applyOutflowBCAtTop << "\n";
      file << "started from checkpoint file = " << readFromCheckPointFile << " ( " << checkPointFileName << ")\n";

      file << "\nOutput quantities:\n";
      file << "actual St = " << actualSt << "\n";
      file << "actual Re = " << actualRe << "\n";
      file << "impact times:\n";
      uint_t bounce = uint_t(0);
      for( auto impact : impactTimes)
      {
         file << " - " << bounce++ << ": impact time = " << impact << "\n";
      }
      file << "settling velocities:\n";
      bounce = uint_t(0);
      for( auto vel : maxSettlingVelBetweenBounces)
      {
        file << " - " << bounce++ << ": vel = " << vel << " -> St = "  << densityRatio * std::abs(vel) * diameter / ( real_t(9) * viscosity ) << ", Re = " << std::abs(vel) * diameter / viscosity << "\n";
      }
      file << "maximal overlap = " << std::abs(minGapSize) << " (" << std::abs(minGapSize)/diameter * real_t(100) << "%)\n";
      file << "maximum gap sizes:\n";
      bounce = uint_t(1);
      for( auto gapSize : maxGapSizeBetweenBounces)
      {
         file << " - " << bounce++ << ": gap size = " << gapSize << ", gap size / diameter = "  << gapSize / diameter << "\n";
      }

      file.close();
   }




   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace sphere_wall_collision

int main( int argc, char **argv ){
   sphere_wall_collision::main(argc, argv);
}
