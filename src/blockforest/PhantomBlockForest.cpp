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
//! \file PhantomBlockForest.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockForest.h"
#include "BlockReconstruction.h"
#include "PhantomBlockForest.h"

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/mpi/BufferSizeTrait.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"

#include <memory>
#include <set>



namespace walberla {
namespace blockforest {



PhantomBlockForest::PhantomBlockForest( BlockForest & blockforest ) :
   blockforest_( blockforest ), depth_( uint_t(0) ), neighborhood_( blockforest.getNeighborhood() )
{}



void PhantomBlockForest::initialize( const BlockStateDeterminationFunction & function, const bool recalculateDepth )
{
   BlockReconstruction::AABBReconstruction aabbReconstruction( blockforest_.getDomain(),
                                                               blockforest_.getXSize(), blockforest_.getYSize(), blockforest_.getZSize(),
                                                               blockforest_.getTreeIdDigits() );

   std::map< uint_t, std::map< BlockID, std::pair< uint_t, Set<SUID> > > > blockNeighborhood;
   std::set< mpi::MPIRank > ranksToRecvFrom;

   const uint_t process = blockforest_.getProcess();

   const auto & blocks = blockforest_.getBlockMap();

   for( auto it = blocks.begin(); it != blocks.end(); ++it )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( it->second );
      auto & block = it->second;

      const auto & id = block->getId();
      const auto & state = block->getState();
      const auto & aabb = block->getAABB();

      const uint_t level = block->getLevel();
      const uint_t targetLevel = block->getTargetLevel();

      block->clearTargetProcess();

      if( targetLevel == level )
      {
         auto cstate = function ? function( std::vector< std::pair< BlockID, Set<SUID> > >( 1, std::make_pair( id, state ) ), id ) : state;

         auto phantom = std::make_shared< PhantomBlock >( *this, id, cstate, aabb, targetLevel, level,
                                                                      std::vector< uint_t >( uint_t(1), process ), process );
         blocks_[ phantom->getId() ] = phantom;

         block->addTargetProcess( process );

         blockNeighborhood[ process ][ id ] = std::make_pair( process, cstate );
      }
      else if( targetLevel > level )
      {
         WALBERLA_ASSERT_EQUAL( targetLevel, (level + uint_t(1)) );

         const real_t xMid = ( aabb.xMin() + aabb.xMax() ) / real_c(2);
         const real_t yMid = ( aabb.yMin() + aabb.yMax() ) / real_c(2);
         const real_t zMid = ( aabb.zMin() + aabb.zMax() ) / real_c(2);

         for( uint_t c = 0; c != 8; ++c )
         {
            BlockID cid( id, c );

            auto cstate = function ? function( std::vector< std::pair< BlockID, Set<SUID> > >( 1, std::make_pair( id, state ) ), cid ) : state;

            AABB caabb( ( ( c & 1 ) ? xMid : aabb.xMin() ),
                        ( ( c & 2 ) ? yMid : aabb.yMin() ),
                        ( ( c & 4 ) ? zMid : aabb.zMin() ),
                        ( ( c & 1 ) ? aabb.xMax() : xMid ),
                        ( ( c & 2 ) ? aabb.yMax() : yMid ),
                        ( ( c & 4 ) ? aabb.zMax() : zMid ) );

            auto phantom = std::make_shared< PhantomBlock >( *this, cid, cstate, caabb, targetLevel, level,
                                                                         std::vector< uint_t >( uint_t(1), process ), process );
            blocks_[ phantom->getId() ] = phantom;

            block->addTargetProcess( process );

            blockNeighborhood[ process ][ cid ] = std::make_pair( process, cstate );
         }
         WALBERLA_ASSERT_EQUAL( block->getTargetProcess().size(), uint_t(8) );
      }
      else
      {
         WALBERLA_ASSERT_EQUAL( targetLevel, (level - uint_t(1)) );

         BlockID fid = id.getFatherId();
         const uint_t branchId = id.getBranchId();

         if( branchId == uint_t(0) )
         {
            std::vector< std::pair< BlockID, Set<SUID> > > sourceStates( 8, std::make_pair( id, state ) );

            AABB faabb;
            aabbReconstruction( faabb, fid );

            std::vector< uint_t > sourceProcesses( uint_t(8), process );
            const auto & neighborhood = block->getNeighborhood();
            for( auto neighbor = neighborhood.begin(); neighbor != neighborhood.end(); ++neighbor )
            {
               auto & nid = neighbor->getId();
               if( blockforest_.getLevelFromBlockId( nid ) == level && nid.getFatherId() == fid )
               {
                  sourceProcesses[ neighbor->getId().getBranchId() ] = neighbor->getProcess();
                  sourceStates[ neighbor->getId().getBranchId() ] = std::make_pair( nid, neighbor->getState() );
               }
            }

            auto cstate = function ? function( sourceStates, fid ) : ( sourceStates[0].second + sourceStates[1].second + sourceStates[2].second +
                                                                       sourceStates[3].second + sourceStates[4].second + sourceStates[5].second +
                                                                       sourceStates[6].second + sourceStates[7].second );

            auto phantom = std::make_shared< PhantomBlock >( *this, fid, cstate, faabb, targetLevel, level, sourceProcesses, process );
            blocks_[ phantom->getId() ] = phantom;

            block->addTargetProcess( process );

            blockNeighborhood[ process ][ fid ] = std::make_pair( process, cstate );
         }
         else
         {
            const auto & neighborhood = block->getNeighborhood();
            for( auto neighbor = neighborhood.begin(); neighbor != neighborhood.end(); ++neighbor )
            {
               const auto & nid = neighbor->getId();
               if( blockforest_.getLevelFromBlockId( nid ) == level && nid.getFatherId() == fid && nid.getBranchId() == uint_t(0) )
               {
                  WALBERLA_ASSERT_EQUAL( block->getTargetProcess().size(), uint_t(0) );
                  block->addTargetProcess( neighbor->getProcess() );
#ifdef NDEBUG
                  break;
#endif
               }
            }
         }
         WALBERLA_ASSERT_EQUAL( block->getTargetProcess().size(), uint_t(1) );
      }

      const auto & neighborhood = block->getNeighborhood();
      for( auto neighbor = neighborhood.begin(); neighbor != neighborhood.end(); ++neighbor )
         if( neighbor->getProcess() != process )
            ranksToRecvFrom.insert( numeric_cast< mpi::MPIRank >( neighbor->getProcess() ) );

      if( recalculateDepth )
         depth_ = std::max( depth_, targetLevel );
   }

   if( recalculateDepth )
      mpi::allReduceInplace( depth_, mpi::MAX );
   else
      depth_ = blockforest_.getDepth();

#ifdef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID
   WALBERLA_CHECK_LESS_EQUAL( blockforest_.getTreeIdDigits() + depth_ * uint_t(3), UINT_BITS );
#endif

   // phantom block neighborhood construction

   // TODO: use OpenMP buffer system?
   mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 1941 ); // phantomblockforest = 112 104 97 110 116 111 109 98 108 111 99 107 102 111 114 101 115 116
   bufferSystem.setReceiverInfo( ranksToRecvFrom, true ); // ATTENTION: true = the size of a message from A to B varies

   for( int i = 0; i != 2; ++i )
   {
      for( auto rank = ranksToRecvFrom.begin(); rank != ranksToRecvFrom.end(); ++rank )
      {
         WALBERLA_ASSERT_UNEQUAL( uint_c(*rank), process );
         bufferSystem.sendBuffer( *rank ) << blockNeighborhood[ process ];
      }

      bufferSystem.sendAll();

      for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
      {
         WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process );
         WALBERLA_ASSERT( ( i == 0 && blockNeighborhood.find( uint_c( recvIt.rank() ) ) == blockNeighborhood.end() ) ||
                          ( i == 1 && blockNeighborhood.find( uint_c( recvIt.rank() ) ) != blockNeighborhood.end() ) );
         recvIt.buffer() >> blockNeighborhood[ uint_c( recvIt.rank() ) ];
      }

      std::map< BlockID, std::pair< uint_t, Set<SUID> > > & localMap = blockNeighborhood[ process ];
      for( auto it = blockNeighborhood.begin(); it != blockNeighborhood.end(); ++it )
      {
         if( it->first != process )
         {
            // localMap.insert( it->second.begin(), it->second.end() ); // using 'insert' doesn't allow to assert consistency ... -> manual for loop
            for( auto blockProcessPair = it->second.begin(); blockProcessPair != it->second.end(); ++blockProcessPair )
            {
#ifndef NDEBUG
               if( localMap.find( blockProcessPair->first ) != localMap.end() )
                  WALBERLA_ASSERT_EQUAL( localMap[ blockProcessPair->first ], blockProcessPair->second );
#endif
               localMap[ blockProcessPair->first ] = blockProcessPair->second;
            }
            it->second.clear();
         }
      }
   }

   std::vector< BlockReconstruction::NeighborhoodReconstructionBlock > neighbors;

   std::map< BlockID, std::pair< uint_t, Set<SUID> > > & localMap = blockNeighborhood[ process ];
   for( auto it = localMap.begin(); it != localMap.end(); ++it )
      neighbors.emplace_back( it->first, it->second.first, it->second.second, aabbReconstruction );

   BlockReconstruction::NeighborhoodReconstruction< PhantomBlock > neighborhoodReconstruction( blockforest_.getDomain(),
                                                                                               blockforest_.isXPeriodic(),
                                                                                               blockforest_.isYPeriodic(),
                                                                                               blockforest_.isZPeriodic() );

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it ) // TODO: can be done in parallel with OpenMP
      neighborhoodReconstruction( it->second.get(), neighbors );

   updateNeighborhood();
}



void PhantomBlockForest::assignBlockData( const PhantomBlockDataAssignmentFunction & function )
{
   if( function )
   {
      std::vector< std::pair< const PhantomBlock *, walberla::any > > blockData;
   
      for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      {
         auto & block = it->second;
         WALBERLA_ASSERT_NOT_NULLPTR( block.get() );
         blockData.emplace_back( block.get(), walberla::any() );
      }
      
      function( blockData, *this );
      WALBERLA_CHECK_EQUAL( blockData.size(), blocks_.size() );
      
      auto bit = blocks_.begin();
      for( auto it = blockData.begin(); it != blockData.end(); ++it )
      {
         WALBERLA_CHECK_NOT_NULLPTR( it->first );
         WALBERLA_CHECK( bit != blocks_.end() );
         WALBERLA_CHECK( bit->second.get() == it->first );

         bit->second->addData( it->second );

         ++bit;
      }
   }
}



bool PhantomBlockForest::calculateMigrationInformation( const MigrationPreparationFunction & function, const uint_t iteration )
{
   processesToRecvFrom_.clear();

   if( function )
   {
      std::vector< std::pair< const PhantomBlock *, uint_t > > targetProcess;
   
      for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      {
         auto & block = it->second;
         WALBERLA_ASSERT_NOT_NULLPTR( block.get() );
         targetProcess.emplace_back( block.get(), block->getTargetProcess() );
      }
      
      bool runAgain = function( targetProcess, processesToRecvFrom_, *this, iteration );

      WALBERLA_CHECK_EQUAL( targetProcess.size(), blocks_.size() );
      for( auto it = processesToRecvFrom_.begin(); it != processesToRecvFrom_.end(); ++it )
         WALBERLA_CHECK_LESS( *it, uint_c( MPIManager::instance()->numProcesses() ) );
      
      auto bit = blocks_.begin();
      for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it )
      {
         WALBERLA_CHECK_NOT_NULLPTR( it->first );
         WALBERLA_CHECK( bit != blocks_.end() );
         WALBERLA_CHECK( bit->second.get() == it->first );
         
         WALBERLA_CHECK_LESS( it->second, uint_c( MPIManager::instance()->numProcesses() ) );

         bit->second->setTargetProcess( it->second );

         ++bit;
      }
      
      return runAgain;
   }

   return false;
}



void PhantomBlockForest::migrate( const PhantomBlockDataPackFunction & packBlockData,
                                  const PhantomBlockDataUnpackFunction & unpackBlockData )
{
   WALBERLA_CHECK( ( packBlockData && unpackBlockData ) || ( !packBlockData && !unpackBlockData ) );

   // exchange target process information with neighbors so that phantom block neighbor data can be updated

   prepareMigration();

   // migrate phantom blocks

   const uint_t process = blockforest_.getProcess();

   std::map< uint_t, std::vector< PhantomBlock * > > blocksToSend;

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      auto & block = it->second;
      WALBERLA_ASSERT_NOT_NULLPTR( block.get() );
      if( block->getTargetProcess() != process )
         blocksToSend[ block->getTargetProcess() ].push_back( block.get() );
   }

   std::set< mpi::MPIRank > ranksToRecvFrom;
   for( auto p = processesToRecvFrom_.begin(); p != processesToRecvFrom_.end(); ++p )
      if( *p != process )
         ranksToRecvFrom.insert( static_cast< mpi::MPIRank >(*p) );

   // TODO: use OpenMP buffer system?
   mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 1944 ); // phantomblockforest = 112 104 97 110 116 111 109 98 108 111 99 107 102 111 114 101 115 116 + 3
   bufferSystem.setReceiverInfo( ranksToRecvFrom, true ); // true = the size of a message from A to B varies

   WALBERLA_LOG_PROGRESS( "PhantomBlockForest migration: pack phantom blocks" );

   for( auto p = blocksToSend.begin(); p != blocksToSend.end(); ++p )
   {
      WALBERLA_ASSERT_UNEQUAL( p->first, process );

      auto & buffer = bufferSystem.sendBuffer( p->first );
      auto & blocks = p->second;

      // TODO: Is the data amount reduction really necessary? Is it even a good idea?!

      buffer << uint_c( blocks.size() );
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
      {
         auto * pBlock = *block;
         WALBERLA_ASSERT_NOT_NULLPTR( pBlock );

         buffer << pBlock->getId() << pBlock->getState();

         if( pBlock->getSourceLevel() == pBlock->getLevel() )
         {
            WALBERLA_ASSERT_EQUAL( pBlock->getSourceProcess().size(), uint_t(1) );
            buffer << int8_t(0) << uint32_c( pBlock->getSourceProcess()[0] );
         }
         else if( pBlock->getSourceLevel() > pBlock->getLevel() )
         {
            WALBERLA_ASSERT_EQUAL( pBlock->getSourceLevel(), pBlock->getLevel() + uint_t(1) );
            WALBERLA_ASSERT_EQUAL( pBlock->getSourceProcess().size(), uint_t(8) );
            buffer << int8_t(1) << uint32_c( pBlock->getSourceProcess()[0] ) << uint32_c( pBlock->getSourceProcess()[1] )
                                << uint32_c( pBlock->getSourceProcess()[2] ) << uint32_c( pBlock->getSourceProcess()[3] )
                                << uint32_c( pBlock->getSourceProcess()[4] ) << uint32_c( pBlock->getSourceProcess()[5] )
                                << uint32_c( pBlock->getSourceProcess()[6] ) << uint32_c( pBlock->getSourceProcess()[7] );
         }
         else
         {
            WALBERLA_ASSERT_EQUAL( pBlock->getSourceLevel(), pBlock->getLevel() - uint_t(1) );
            WALBERLA_ASSERT_EQUAL( pBlock->getSourceProcess().size(), uint_t(1) );
            buffer << int8_t(-1) << uint32_c( pBlock->getSourceProcess()[0] );
         }

         buffer << uint8_c( pBlock->getNeighborhoodSize() );
         for( uint_t n = uint_t(0); n != pBlock->getNeighborhoodSize(); ++n )
            buffer << pBlock->getNeighborId(n) << uint32_c( pBlock->getNeighborProcess(n) ) << pBlock->getNeighborState(n);

         if( packBlockData )
         {
            WALBERLA_ASSERT( static_cast<bool>(unpackBlockData) );
            packBlockData( buffer, *pBlock );
         }

         blocks_.erase( pBlock->getId() );
      }
   }

   WALBERLA_LOG_PROGRESS( "PhantomBlockForest migration: send phantom blocks" );

   bufferSystem.sendAll();

   WALBERLA_LOG_PROGRESS( "PhantomBlockForest migration: receive and unpack phantom blocks" );

   BlockReconstruction::NeighborhoodSectionReconstruction< PhantomBlock > neighborhoodSectionReconstruction(
            blockforest_.getDomain(), blockforest_.getXSize(), blockforest_.getYSize(), blockforest_.getZSize(),
            blockforest_.isXPeriodic(), blockforest_.isYPeriodic(), blockforest_.isZPeriodic(), blockforest_.getTreeIdDigits() );

   for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
   {
      WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process );
      auto & buffer = recvIt.buffer();

      uint_t nBlocks( uint_t(0) );
      buffer >> nBlocks;
      for( uint_t i = 0; i != nBlocks; ++ i )
      {
         BlockID id;
         Set<SUID> state;
         buffer >> id >> state;
         math::AABB aabb;
         const uint_t level = blockforest_.getAABBFromBlockId( aabb, id );
         uint_t sourceLevel = level;

         int8_t offset( uint8_t(0) );
         buffer >> offset;
         std::vector< uint_t > sp;
         if( offset == int8_t(0) || offset == int8_t(-1) )
         {
            uint32_t s( uint32_t(0) );
            buffer >> s;
            sp.push_back( s );
            if( offset == int8_t(-1) )
            {
               WALBERLA_ASSERT_GREATER( sourceLevel, uint_t(0) );
               --sourceLevel;
            }
         }
         else
         {
            WALBERLA_ASSERT_EQUAL( offset, int8_t(1) );

            uint32_t s0( uint32_t(0) );
            uint32_t s1( uint32_t(0) );
            uint32_t s2( uint32_t(0) );
            uint32_t s3( uint32_t(0) );
            uint32_t s4( uint32_t(0) );
            uint32_t s5( uint32_t(0) );
            uint32_t s6( uint32_t(0) );
            uint32_t s7( uint32_t(0) );
            buffer >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7;
            sp.push_back( s0 ); sp.push_back( s1 ); sp.push_back( s2 ); sp.push_back( s3 );
            sp.push_back( s4 ); sp.push_back( s5 ); sp.push_back( s6 ); sp.push_back( s7 );
            ++sourceLevel;
         }

         auto phantom = std::make_shared< PhantomBlock >( *this, id, state, aabb, level, sourceLevel, sp, process );
         blocks_[ id ] = phantom;

         phantom->clearNeighborhood();
         uint8_t neighbors( uint8_t(0) );
         buffer >> neighbors;
         for( uint8_t n = uint8_t(0); n != neighbors; ++n )
         {
            BlockID nId;
            uint32_t nProcess;
            Set<SUID> nState;
            buffer >> nId >> nProcess >> nState;
            phantom->addNeighbor( nId, nProcess, nState );
         }
         neighborhoodSectionReconstruction( phantom.get() );

         if( unpackBlockData )
         {
            WALBERLA_ASSERT( static_cast<bool>(packBlockData) );
            walberla::any data;
            unpackBlockData( buffer, *phantom, data );
            phantom->addData( data );
         }
      }
   }

   // update process neighborhood

   WALBERLA_LOG_PROGRESS( "PhantomBlockForest migration: update process neighborhood" );

   updateNeighborhood();
   processesToRecvFrom_.clear();
}



void PhantomBlockForest::updateNeighborhood()
{
   std::set< uint_t > neighbors;

   if( blockforest_.insertBuffersIntoProcessNetwork() )
   {
      std::map< mpi::MPIRank, mpi::MPISize > ranksToRecvFrom;
      for( auto n = neighborhood_.begin(); n != neighborhood_.end(); ++n )
         ranksToRecvFrom[ static_cast< mpi::MPIRank >(*n) ] = mpi::BufferSizeTrait<int>::size;

      mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 1942 ); // phantomblockforest = 112 104 97 110 116 111 109 98 108 111 99 107 102 111 114 101 115 116 + 1
      bufferSystem.setReceiverInfo( ranksToRecvFrom );

      const int empty = blocks_.empty() ? 1 : 0;
      for( auto rank = ranksToRecvFrom.begin(); rank != ranksToRecvFrom.end(); ++rank )
      {
         WALBERLA_ASSERT_UNEQUAL( uint_c(rank->first), blockforest_.getProcess() );
         bufferSystem.sendBuffer( rank->first ) << empty;
      }

      bufferSystem.sendAll();

      for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
      {
         uint_t np = uint_c( recvIt.rank() );
         int nEmpty( 0 );
         WALBERLA_ASSERT_UNEQUAL( np, blockforest_.getProcess() );
         recvIt.buffer() >> nEmpty;
         if( nEmpty == 1 )
            neighbors.insert( np );
      }

      if( blocks_.empty() )
         neighbors.insert( neighborhood_.begin(), neighborhood_.end() );
   }

   const uint_t process = blockforest_.getProcess();
   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      auto & block = it->second;
      WALBERLA_ASSERT_NOT_NULLPTR( block.get() );
      const auto & neighborhood = block->getNeighborhood();
      for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
         if( n->getProcess() != process )
            neighbors.insert( n->getProcess() );
   }

   neighborhood_.clear();
   for( auto n = neighbors.begin(); n != neighbors.end(); ++n )
      neighborhood_.push_back( *n );
}



void PhantomBlockForest::prepareMigration()
{
   // update process IDs (ranks) of all phantom block neighbors (= every block tells its neighbors where it will migrate to)
   
   WALBERLA_LOG_PROGRESS( "PhantomBlockForest migration: every block tells its neighbors where it will migrate to "
                          "so that after phantom migration the block neighborhood data is in a valid state" );

   std::map< uint_t, std::map< BlockID, uint_t > > blockProcessMap;
   std::set< mpi::MPIRank > ranksToRecvFrom;

   const uint_t process = blockforest_.getProcess();

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      blockProcessMap[ process ][ it->first ] = it->second->getTargetProcess();
      auto & neighbors = it->second->getNeighborhood();
      for( auto n = neighbors.begin(); n != neighbors.end(); ++n )
         if( n->getProcess() != process )
            ranksToRecvFrom.insert( numeric_cast< mpi::MPIRank >( n->getProcess() ) );
   }
   WALBERLA_ASSERT( ( blocks_.empty() && blockProcessMap.empty() ) || ( blockProcessMap.size() == uint_t(1) ) );

   // TODO: use OpenMP buffer system?
   mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 1943 ); // phantomblockforest = 112 104 97 110 116 111 109 98 108 111 99 107 102 111 114 101 115 116 + 2
   bufferSystem.setReceiverInfo( ranksToRecvFrom, true ); // true = the size of a message from A to B varies

   for( auto rank = ranksToRecvFrom.begin(); rank != ranksToRecvFrom.end(); ++rank )
   {
      WALBERLA_ASSERT_UNEQUAL( uint_c(*rank), process );
      bufferSystem.sendBuffer( *rank ) << blockProcessMap[ process ];
   }

   bufferSystem.sendAll();

   for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
   {
      WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process );
      WALBERLA_ASSERT( blockProcessMap.find( uint_c( recvIt.rank() ) ) == blockProcessMap.end() );
      recvIt.buffer() >> blockProcessMap[ uint_c( recvIt.rank() ) ];
   }

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      auto & block = it->second;
      for( uint_t n = uint_t(0); n != block->getNeighborhoodSize(); ++n )
      {
         const uint_t p = block->getNeighborProcess(n);
         WALBERLA_ASSERT( blockProcessMap.find(p) != blockProcessMap.end() );
         WALBERLA_ASSERT( blockProcessMap[p].find( block->getNeighborId(n) ) != blockProcessMap[p].end() );
         block->setNeighborProcess( n, blockProcessMap[p][ block->getNeighborId(n) ] );
      }
   }

   // update "link" between real blocks and phantoms

   WALBERLA_LOG_PROGRESS( "PhantomBlockForest migration: 'links' between real blocks and phantom blocks are updated "
                          "so that after phantom migration these links are in a valid state" );

   ranksToRecvFrom.clear();

   const auto & blocks = blockforest_.getBlockMap();
   for( auto it = blocks.begin(); it != blocks.end(); ++it )
   {
      auto & targets = it->second->getTargetProcess();
      for( auto tp = targets.begin(); tp != targets.end(); ++tp )
         if( *tp != process )
            ranksToRecvFrom.insert( numeric_cast< mpi::MPIRank >( *tp ) );
   }
   
   mpi::BufferSystem linkBufferSystem( MPIManager::instance()->comm(), 1945 ); // phantomblockforest = 112 104 97 110 116 111 109 98 108 111 99 107 102 111 114 101 115 116 + 4
   linkBufferSystem.setReceiverInfo( ranksToRecvFrom, true ); // true = the size of a message from A to B varies
   
   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      auto & block = it->second;
      auto & id = block->getId();
      auto & sources = block->getSourceProcess();
      
      if( block->sourceBlockHasTheSameSize() )
      {
         WALBERLA_ASSERT_EQUAL( sources.size(), uint_t(1) );
         if( sources[0] != process )
         {
            linkBufferSystem.sendBuffer( sources[0] ) << id << uint8_t(0) << block->getTargetProcess();
         }
         else
         {
            WALBERLA_ASSERT( blocks.find(id) != blocks.end() );
            WALBERLA_ASSERT_EQUAL( blocks.find(id)->second->getTargetProcess().size(), uint_t(1) );
            blocks.find(id)->second->setTargetProcess( uint_t(0), block->getTargetProcess() );
         }
      }
      else if( block->sourceBlockIsLarger() )
      {
         WALBERLA_ASSERT_EQUAL( sources.size(), uint_t(1) );
         if( sources[0] != process )
         {
            linkBufferSystem.sendBuffer( sources[0] ) << id.getFatherId() << uint8_c( id.getBranchId() ) << block->getTargetProcess();
         }
         else
         {
            WALBERLA_ASSERT( blocks.find( id.getFatherId() ) != blocks.end() );
            WALBERLA_ASSERT_EQUAL( blocks.find( id.getFatherId() )->second->getTargetProcess().size(), uint_t(8) );
            blocks.find( id.getFatherId() )->second->setTargetProcess( id.getBranchId(), block->getTargetProcess() );
         }
      }
      else
      {
         WALBERLA_ASSERT( block->sourceBlockIsSmaller() );
         WALBERLA_ASSERT_EQUAL( sources.size(), uint_t(8) );
         for( uint_t i = uint_t(0); i != uint_t(8); ++i )
         {
            if( sources[i] != process )
            {
               linkBufferSystem.sendBuffer( sources[i] ) << BlockID( id, i ) << uint8_t(0) << block->getTargetProcess();
            }
            else
            {
               WALBERLA_ASSERT( blocks.find( BlockID( id, i ) ) != blocks.end() );
               WALBERLA_ASSERT_EQUAL( blocks.find( BlockID( id, i ) )->second->getTargetProcess().size(), uint_t(1) );
               blocks.find( BlockID( id, i ) )->second->setTargetProcess( uint_t(0), block->getTargetProcess() );
            }
         }
      }
   }
   
   linkBufferSystem.sendAll();

   for( auto recvIt = linkBufferSystem.begin(); recvIt != linkBufferSystem.end(); ++recvIt )
   {
      WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process );
      
      while( !(recvIt.buffer().isEmpty()) )
      {
         BlockID id;
         uint8_t index;
         uint_t tp;
      
         recvIt.buffer() >> id >> index >> tp;
         
         WALBERLA_ASSERT( blocks.find(id) != blocks.end() );
         WALBERLA_ASSERT_GREATER( blocks.find(id)->second->getTargetProcess().size(), index );
         
         blocks.find(id)->second->setTargetProcess( uint_c( index ), tp );
      }
   }
}



} // namespace blockforest
} // namespace walberla
