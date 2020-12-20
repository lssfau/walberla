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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/domain/BlockForestDataHandling.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/domain/InfoCollection.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <blockforest/Initialization.h>
#include <blockforest/loadbalancing/DynamicCurve.h>
#include <blockforest/loadbalancing/InfoCollection.h>
#include <core/debug/TestSubsystem.h>
#include <core/logging/Logging.h>
#include <pe/amr/level_determination/MinMaxLevelDetermination.h>
#include <pe/amr/weight_assignment/WeightAssignmentFunctor.h>


namespace walberla {

using namespace walberla::mesa_pd;

///switches between level 0 and 1 in for every call
class ReGrid_Switcher
{
public:
   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > &, const BlockForest & forest );
};

void ReGrid_Switcher::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                  std::vector< const Block * > &, const BlockForest & /*forest*/ )
{
   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      it->second = it->first->getLevel() == 0 ? 1 : 0;
      WALBERLA_LOG_DETAIL("setting block to level: " << it->second);
   }
}

void createSphere(data::ParticleStorage& ps, domain::IDomain& domain, const Vec3& pos)
{
   auto owned = domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), pos );
   if (owned)
   {
      data::Particle&& p          = *ps.create();
      p.getPositionRef()          = pos;
      p.getInteractionRadiusRef() = real_t(1);
      p.getOwnerRef()             = mpi::MPIManager::instance()->rank();
   }
}

/// This tests starts with one block located on one process.
/// It then generates 1 particle in every quadrant.
/// Load balancing is setup such that it will subdivide the block into 8 subblocks.
/// Correct transmission of all particles is checked
/// Blocks are coarsened again and particles are checked.
int main( bool simple )
{
   //   logging::Logging::instance()->setStreamLogLevel(logging::Logging::DETAIL);

   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(0,0,0,10,10,10), // simulation domain
                                                 Vector3<uint_t>(1,1,1), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );

   forest->recalculateBlockLevelsInRefresh( true );
   forest->alwaysRebalanceInRefresh( false );
   forest->reevaluateMinTargetLevelsAfterForcedRefinement( false );
   forest->allowRefreshChangingDepth( true );

   forest->allowMultipleRefreshCycles( false );
   forest->checkForEarlyOutInRefresh( true );
   forest->checkForLateOutInRefresh( true );

   auto infoCollection = make_shared<blockforest::InfoCollection>();

   if (simple)
   {
      // use ReGrid_Switcher
      forest->setRefreshMinTargetLevelDeterminationFunction( ReGrid_Switcher() );
      forest->setRefreshPhantomBlockMigrationPreparationFunction(
               blockforest::DynamicCurveBalance< blockforest::NoPhantomData >( false, true, false ) );
   } else
   {
      forest->setRefreshMinTargetLevelDeterminationFunction( pe::amr::MinMaxLevelDetermination(infoCollection, 2, 5) );

      forest->setRefreshPhantomBlockDataAssignmentFunction( pe::amr::WeightAssignmentFunctor( infoCollection ) );
      forest->setRefreshPhantomBlockDataPackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      forest->setRefreshPhantomBlockDataUnpackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

      forest->setRefreshPhantomBlockMigrationPreparationFunction(
               blockforest::DynamicCurveBalance< pe::amr::WeightAssignmentFunctor::PhantomBlockWeight >( false, true, false ) );
   }

   auto ps = std::make_shared<data::ParticleStorage> (100);
   data::ParticleAccessor ac(ps);
   forest->addBlockData(domain::createBlockForestDataHandling(ps), "Storage");

   domain::BlockForestDomain domain(forest);
   // creating particles
   createSphere(*ps, domain, Vec3(real_t(2.5), real_t(2.5), real_t(2.5)));
   createSphere(*ps, domain, Vec3(real_t(2.5), real_t(2.5), real_t(7.5)));
   createSphere(*ps, domain, Vec3(real_t(2.5), real_t(7.5), real_t(2.5)));
   createSphere(*ps, domain, Vec3(real_t(2.5), real_t(7.5), real_t(7.5)));
   createSphere(*ps, domain, Vec3(real_t(7.5), real_t(2.5), real_t(2.5)));
   createSphere(*ps, domain, Vec3(real_t(7.5), real_t(2.5), real_t(7.5)));
   createSphere(*ps, domain, Vec3(real_t(7.5), real_t(7.5), real_t(2.5)));
   createSphere(*ps, domain, Vec3(real_t(7.5), real_t(7.5), real_t(7.5)));

   //***** TESTING *****
   WALBERLA_CHECK_EQUAL( forest->size(), mpi::MPIManager::instance()->rank() == 0 ? 1 : 0 );
   WALBERLA_CHECK_EQUAL( ps->size(), mpi::MPIManager::instance()->rank() == 0 ? 8 : 0 );

   domain::createWithNeighborhood(ac, *forest, *infoCollection);
   forest->refresh();
   WALBERLA_CHECK_EQUAL( forest->size(), 1 );
   WALBERLA_CHECK_EQUAL( ps->size(), 1 );
   WALBERLA_CHECK_FLOAT_EQUAL(ps->getPosition(0), forest->begin()->getAABB().center());

   if (!simple)
   {
      forest->setRefreshMinTargetLevelDeterminationFunction( pe::amr::MinMaxLevelDetermination(infoCollection, 2, 9) );
   }
   domain::createWithNeighborhood(ac, *forest, *infoCollection);
   forest->refresh();
   WALBERLA_CHECK_EQUAL( forest->size(), mpi::MPIManager::instance()->rank() == 0 ? 1 : 0 );
   WALBERLA_CHECK_EQUAL( ps->size(), mpi::MPIManager::instance()->rank() == 0 ? 8 : 0 );

   return EXIT_SUCCESS;
}

} // namespace walberla

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   WALBERLA_UNUSED(argc);
   WALBERLA_UNUSED(argv);

   walberla::main( true );
   walberla::main( false );
   return EXIT_SUCCESS;
}
