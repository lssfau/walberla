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
//! \file DynamicRefinement.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"
#include "pe/synchronization/ClearSynchronization.h"
#include "pe/utility/GetBody.h"
#include "pe/utility/DestroyBody.h"

#include "blockforest/Initialization.h"
#include <blockforest/loadbalancing/DynamicCurve.h>
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "core/debug/TestSubsystem.h"

namespace walberla {
using namespace walberla::pe;

using BodyTuple = std::tuple<Sphere> ;

class ReGrid
{
public:

   ReGrid( const BlockDataID storageID, const size_t minParticles, const size_t maxParticles) :
      storageID_( storageID ), minParticles_(minParticles), maxParticles_(maxParticles)
   {}

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > &, const BlockForest & forest );

private:
   const BlockDataID storageID_;
   const size_t      minParticles_;
   const size_t      maxParticles_;
};

void ReGrid::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                         std::vector< const Block * > &, const BlockForest & /*forest*/ )
{
   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      const auto numberOfParticles = (*(it->first->getData< Storage >( storageID_ )))[0].size();
      //WALBERLA_LOG_DEVEL("storage size: " << localBodyStorage.size());

      it->second = it->first->getLevel(); //keep everything as it is
      if (numberOfParticles < minParticles_)
      {
         WALBERLA_LOG_DEVEL(it->first->getLevel() << " -> " << it->first->getLevel() - uint_t(1) << " (" << numberOfParticles << ")" );
         if (it->first->getLevel() > 0)
            it->second = it->first->getLevel() - uint_t(1);
      } else if (numberOfParticles > maxParticles_)
      {
         it->second = it->first->getLevel() + uint_t(1);
         WALBERLA_LOG_DEVEL(it->first->getLevel() << " -> " << it->first->getLevel() + uint_t(1) << " (" << numberOfParticles << ")" );
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            math::AABB(0,0,0,20,20,20),
            uint_c( 1), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            false,                              // max blocks per process
            false, false, false,                // full periodicity
            false);

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "HCCD");
   forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   auto & blockforest = forest->getBlockForest();
   blockforest.recalculateBlockLevelsInRefresh( true );
   blockforest.alwaysRebalanceInRefresh( false );
   blockforest.reevaluateMinTargetLevelsAfterForcedRefinement( false );
   blockforest.allowRefreshChangingDepth( true );

   blockforest.allowMultipleRefreshCycles( false );
   blockforest.checkForEarlyOutInRefresh( true );
   blockforest.checkForLateOutInRefresh( true );

   ReGrid regrid( storageID, 20, 20 );

   blockforest.setRefreshMinTargetLevelDeterminationFunction( regrid );

   blockforest.setRefreshPhantomBlockMigrationPreparationFunction(
            blockforest::DynamicCurveBalance< blockforest::NoPhantomData >( true, true ) );

   real_t spacing(2.5);
   for (auto blkIt = forest->begin(); blkIt != forest->end(); ++blkIt)
   {
      IBlock & currentBlock = *blkIt;
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB(), Vector3<real_t>(spacing) * real_t(0.5), spacing); it != grid_generator::SCIterator(); ++it)
      {
         createSphere( *globalBodyStorage, forest->getBlockStorage(), storageID, 0, *it, 1 );
      }
   }
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);

   clearSynchronization( forest->getBlockForest(), storageID );
   forest->refresh();
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);

   WALBERLA_ASSERT_EQUAL( forest->size(), 8 );
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
//      IBlock & currentBlock = *blockIt;
//      Storage * storage = currentBlock.getData< Storage >( storageID );
//      BodyStorage& localStorage = (*storage)[0];
//      BodyStorage& shadowStorage = (*storage)[1];

      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt)
      {
         WALBERLA_ASSERT( blockIt->getAABB().contains(bodyIt->getPosition()) );

//         WALBERLA_LOG_DEVEL( blockIt->getAABB() );
//         WALBERLA_LOG_DEVEL(*bodyIt );
      }
   }

   WALBERLA_LOG_DEVEL("========================================================");

   clearSynchronization( forest->getBlockForest(), storageID );
   forest->refresh();
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);

   WALBERLA_ASSERT_EQUAL( forest->size(), 64 );
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
//      IBlock & currentBlock = *blockIt;
//      Storage * storage = currentBlock.getData< Storage >( storageID );
//      BodyStorage& localStorage = (*storage)[0];
//      BodyStorage& shadowStorage = (*storage)[1];

      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt)
      {
         WALBERLA_ASSERT( blockIt->getAABB().contains(bodyIt->getPosition()) );

//         WALBERLA_LOG_DEVEL( blockIt->getAABB() );
//         WALBERLA_LOG_DEVEL(*bodyIt );
      }
   }

   WALBERLA_LOG_DEVEL("========================================================");

   clearSynchronization( forest->getBlockForest(), storageID );
   forest->refresh();
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);

   WALBERLA_ASSERT_EQUAL( forest->size(), 8 );
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
//      IBlock & currentBlock = *blockIt;
//      Storage * storage = currentBlock.getData< Storage >( storageID );
//      BodyStorage& localStorage = (*storage)[0];
//      BodyStorage& shadowStorage = (*storage)[1];

      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt)
      {
         WALBERLA_ASSERT( blockIt->getAABB().contains(bodyIt->getPosition()) );

//         WALBERLA_LOG_DEVEL( blockIt->getAABB() );
//         WALBERLA_LOG_DEVEL(*bodyIt );
      }
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
