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
//! \file BlockForest.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockForest.h"
#include "BlockForestFile.h"
#include "BlockNeighborhoodSection.h"
#include "SetupBlockForest.h"
#include "core/Abort.h"
#include "core/EndianIndependentSerialization.h"
#include "core/debug/CheckFunctions.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/SetReduction.h"
#include "core/mpi/MPIManager.h"

#include <fstream>
#include <memory>
#include <set>
#include <stack>
#include <utility>


namespace walberla {
namespace blockforest {



bool BlockForest::BlockInformation::operator==( const BlockInformation & rhs ) const
{
   if( nodes_.size() != rhs.nodes_.size() )
      return false;

   for( uint_t i = 0; i != nodes_.size(); ++i ) {
      if( !(nodes_[i]) ) {
         if( rhs.nodes_[i] )
            return false;
      }
      else if( !(rhs.nodes_[i]) || *(nodes_[i]) != *(rhs.nodes_[i]) )
         return false;
   }

   return true;
}



void BlockForest::BlockInformation::getAllBlocks( std::vector< shared_ptr< IBlockID > > & blocks ) const
{
   const uint_t treeIdMarker = uintPow2( forest_.getTreeIdDigits() - 1 );

   for( uint_t i = 0; i != nodes_.size(); ++i )
   {
      if( !(nodes_[i]) )
         continue;

      std::stack< std::pair< const Node *, BlockID > > stack;

      stack.push( std::make_pair( nodes_[i].get(), BlockID( i, treeIdMarker ) ) );

      while( !stack.empty() ) {

         std::pair< const Node *, BlockID > node = stack.top();
         stack.pop();

         WALBERLA_ASSERT_NOT_NULLPTR( node.first );

         if( !node.first->children_.empty() ) {
            for( uint_t c = 8; c-- != 0; )
               stack.push( std::make_pair( node.first->children_[c].get(), BlockID( node.second, c ) ) );
         }
         else {
            blocks.push_back( make_shared< BlockID >( node.second ) );
         }
      }
   }
}



bool BlockForest::BlockInformation::getId( BlockID & id, const real_t x, const real_t y, const real_t z ) const
{
   const AABB & domain = forest_.getDomain();

   if( nodes_.empty() || !domain.contains(x,y,z) )
      return false;

   const real_t rootBlockXSize = forest_.getRootBlockXSize();
   const real_t rootBlockYSize = forest_.getRootBlockYSize();
   const real_t rootBlockZSize = forest_.getRootBlockZSize();

   uint_t xi = static_cast< uint_t >( ( x - domain.xMin() ) / rootBlockXSize );
   uint_t yi = static_cast< uint_t >( ( y - domain.yMin() ) / rootBlockYSize );
   uint_t zi = static_cast< uint_t >( ( z - domain.zMin() ) / rootBlockZSize );

   if( xi >= forest_.getXSize() ) xi = forest_.getXSize() - 1; // shouldn't happen, ...
   if( yi >= forest_.getYSize() ) yi = forest_.getYSize() - 1; // ... but might happen due to ...
   if( zi >= forest_.getZSize() ) zi = forest_.getZSize() - 1; // ... floating point inaccuracy?

   const uint_t index = zi * forest_.getYSize() * forest_.getXSize() + yi * forest_.getXSize() + xi;

   id.clear();
   id.init( index, uint_c(1) << ( forest_.getTreeIdDigits() - 1 ) );

   auto node = nodes_[index];

   if( !node )
      return false;

   AABB aabb = AABB::createFromMinMaxCorner( domain.xMin() + static_cast< real_t >( xi ) * rootBlockXSize,
                                             domain.yMin() + static_cast< real_t >( yi ) * rootBlockYSize,
                                             domain.zMin() + static_cast< real_t >( zi ) * rootBlockZSize,
                                             ( xi+1 == forest_.getXSize() ) ? domain.xMax() : domain.xMin() + static_cast< real_t >( xi+1 ) * rootBlockXSize,
                                             ( yi+1 == forest_.getYSize() ) ? domain.yMax() : domain.yMin() + static_cast< real_t >( yi+1 ) * rootBlockYSize,
                                             ( zi+1 == forest_.getZSize() ) ? domain.zMax() : domain.zMin() + static_cast< real_t >( zi+1 ) * rootBlockZSize );
   WALBERLA_ASSERT( aabb.contains(x,y,z) );

   while( !node->children_.empty() ) {

      uint_t branchId = 0;

      const real_t xMid = ( aabb.xMin() + aabb.xMax() ) / real_c(2);
      const real_t yMid = ( aabb.yMin() + aabb.yMax() ) / real_c(2);
      const real_t zMid = ( aabb.zMin() + aabb.zMax() ) / real_c(2);

      if( x >= xMid ) ++branchId;
      if( y >= yMid ) branchId += 2;
      if( z >= zMid ) branchId += 4;

      aabb.initMinMaxCorner( ( branchId & 1 ) ? xMid : aabb.xMin(), ( branchId & 2 ) ? yMid : aabb.yMin(), ( branchId & 4 ) ? zMid : aabb.zMin(),
                             ( branchId & 1 ) ? aabb.xMax() : xMid, ( branchId & 2 ) ? aabb.yMax() : yMid, ( branchId & 4 ) ? aabb.zMax() : zMid );

      WALBERLA_ASSERT( aabb.contains(x,y,z) );
      WALBERLA_ASSERT_LESS( branchId, node->children_.size() );
      WALBERLA_ASSERT_NOT_NULLPTR( node->children_[ branchId ] );

      id.appendBranchId( branchId );

      node = node->children_[ branchId ];
   }

   return true;
}



const BlockForest::BlockInformation::Node * BlockForest::BlockInformation::getNode( const BlockID & id ) const
{
   if( nodes_.empty() )
      return nullptr;

   const uint_t treeIdDigits = forest_.getTreeIdDigits();

   WALBERLA_ASSERT_GREATER_EQUAL( id.getUsedBits(), treeIdDigits );
   WALBERLA_ASSERT_EQUAL( ( id.getUsedBits() - treeIdDigits ) % 3, 0 );

   BlockID blockId( id );

   const uint_t levels = ( blockId.getUsedBits() - treeIdDigits ) / 3;

   std::vector< uint_t > branchId( levels );

   for( uint_t i = levels; i-- != 0; ) {
      branchId[i] = blockId.getBranchId();
      blockId.removeBranchId();
   }

   const uint_t index = blockId.getTreeIndex();

   if( index >= nodes_.size() || !(nodes_[index]) )
      return nullptr;

   auto node = nodes_[index];

   for( uint_t i = 0; i != levels; ++i ) {
      if( node->children_.empty() )
         return nullptr;
      WALBERLA_ASSERT_NOT_NULLPTR( node->children_[ branchId[i] ] );
      node = node->children_[ branchId[i] ];
   }

   return node.get();
}



const BlockForest::BlockInformation::Node * BlockForest::BlockInformation::getNode( const real_t x, const real_t y, const real_t z ) const
{
   const AABB & domain = forest_.getDomain();

   if( nodes_.empty() || !domain.contains(x,y,z) )
      return nullptr;

   const real_t rootBlockXSize =  forest_.getRootBlockXSize();
   const real_t rootBlockYSize =  forest_.getRootBlockYSize();
   const real_t rootBlockZSize =  forest_.getRootBlockZSize();

   uint_t xi = static_cast< uint_t >( ( x - domain.xMin() ) / rootBlockXSize );
   uint_t yi = static_cast< uint_t >( ( y - domain.yMin() ) / rootBlockYSize );
   uint_t zi = static_cast< uint_t >( ( z - domain.zMin() ) / rootBlockZSize );

   if( xi >= forest_.getXSize() ) xi = forest_.getXSize() - 1; // shouldn't happen, ...
   if( yi >= forest_.getYSize() ) yi = forest_.getYSize() - 1; // ... but might happen due to ...
   if( zi >= forest_.getZSize() ) zi = forest_.getZSize() - 1; // ... floating point inaccuracy?

   const uint_t index = zi * forest_.getYSize() * forest_.getXSize() + yi * forest_.getXSize() + xi;

   auto node = nodes_[index];

   if( !node )
      return nullptr;

   AABB aabb = AABB::createFromMinMaxCorner( domain.xMin() + static_cast< real_t >( xi ) * rootBlockXSize,
                                             domain.yMin() + static_cast< real_t >( yi ) * rootBlockYSize,
                                             domain.zMin() + static_cast< real_t >( zi ) * rootBlockZSize,
                                             ( xi+1 == forest_.getXSize() ) ? domain.xMax() : domain.xMin() + static_cast< real_t >( xi+1 ) * rootBlockXSize,
                                             ( yi+1 == forest_.getYSize() ) ? domain.yMax() : domain.yMin() + static_cast< real_t >( yi+1 ) * rootBlockYSize,
                                             ( zi+1 == forest_.getZSize() ) ? domain.zMax() : domain.zMin() + static_cast< real_t >( zi+1 ) * rootBlockZSize );
   WALBERLA_ASSERT( aabb.contains(x,y,z) );

   while( !node->children_.empty() ) {

      uint_t branchId = 0;

      const real_t xMid = ( aabb.xMin() + aabb.xMax() ) / real_c(2);
      const real_t yMid = ( aabb.yMin() + aabb.yMax() ) / real_c(2);
      const real_t zMid = ( aabb.zMin() + aabb.zMax() ) / real_c(2);

      if( x >= xMid ) ++branchId;
      if( y >= yMid ) branchId += 2;
      if( z >= zMid ) branchId += 4;

      aabb.initMinMaxCorner( ( branchId & 1 ) ? xMid : aabb.xMin(), ( branchId & 2 ) ? yMid : aabb.yMin(), ( branchId & 4 ) ? zMid : aabb.zMin(),
                             ( branchId & 1 ) ? aabb.xMax() : xMid, ( branchId & 2 ) ? aabb.yMax() : yMid, ( branchId & 4 ) ? aabb.zMax() : zMid );

      WALBERLA_ASSERT( aabb.contains(x,y,z) );
      WALBERLA_ASSERT_LESS( branchId, node->children_.size() );
      WALBERLA_ASSERT_NOT_NULLPTR( node->children_[ branchId ] );

      node = node->children_[ branchId ];
   }

   return node.get();
}



BlockForest::BlockForest( const uint_t process, const SetupBlockForest& forest, const bool keepGlobalBlockInformation ) :

   BlockStorage( forest.getDomain(), forest.isXPeriodic(), forest.isYPeriodic(), forest.isZPeriodic() ),

   process_( process ), processIdBytes_( forest.getProcessIdBytes() ), depth_( forest.getDepth() ), treeIdDigits_( forest.getTreeIdDigits() ),
   insertBuffersIntoProcessNetwork_( forest.insertBuffersIntoProcessNetwork() ), modificationStamp_( uint_t(0) ),
   recalculateBlockLevelsInRefresh_( true ), alwaysRebalanceInRefresh_( false ), allowMultipleRefreshCycles_( true ),
   reevaluateMinTargetLevelsAfterForcedRefinement_( false ),
   checkForEarlyOutInRefresh_( true ), checkForLateOutInRefresh_( true ), allowChangingDepth_( true ), checkForEarlyOutAfterLoadBalancing_( false ),
   phantomBlockMigrationIterations_( uint_t(0) ),
   nextCallbackBeforeBlockDataIsPackedHandle_( uint_t(0) ), nextCallbackBeforeBlockDataIsUnpackedHandle_( uint_t(0) ),
   nextCallbackAfterBlockDataIsUnpackedHandle_( uint_t(0) ),
   snapshotExists_( false ), snapshotDepth_( uint_t(0) ), snapshotBlockDataItems_( uint_t(0) )
{
   blockInformation_ = make_shared< BlockInformation >( *this );

   size_[0] = forest.getXSize();
   size_[1] = forest.getYSize();
   size_[2] = forest.getZSize();

   if( forest.isWorkerProcess( process ) ) {

      std::vector< const SetupBlock* > blocks;

      forest.getProcessSpecificBlocks( blocks, process );

      WALBERLA_ASSERT( !blocks.empty() );

      std::set< uint_t > neighbors;

      for( uint_t i = 0; i != blocks.size(); ++i ) {

         WALBERLA_ASSERT( blocks_.find( blocks[i]->getId() ) == blocks_.end() );

         blocks_[ blocks[i]->getId() ] = std::make_shared< Block >( *this, blocks[i] );

         for( uint_t j = 0; j != blocks[i]->getNeighborhoodSize(); ++j )
            if( blocks[i]->getNeighbor(j)->getProcess() != process )
               neighbors.insert( blocks[i]->getNeighbor(j)->getProcess() );
      }

      if( insertBuffersIntoProcessNetwork_ && ( process + 1 ) < forest.getNumberOfProcesses() &&
          forest.isBufferProcess( process + 1 ) )
         neighbors.insert( process + 1 );

      for( std::set< uint_t >::iterator it = neighbors.begin(); it != neighbors.end(); ++it )
         neighborhood_.push_back( *it );
   }
   else if( insertBuffersIntoProcessNetwork_ )
   {
      WALBERLA_ASSERT( forest.isBufferProcess( process ) );
      WALBERLA_ASSERT_GREATER( process, 0 );

      neighborhood_.push_back( process - 1 );

      if( ( process + 1 ) < forest.getNumberOfProcesses() )
      {
         WALBERLA_ASSERT( !( forest.isBufferProcess( process - 1 ) && forest.isWorkerProcess( process + 1 ) ) );

         if( forest.isBufferProcess( process + 1 ) )
            neighborhood_.push_back( process + 1 );
      }
   }

   if( keepGlobalBlockInformation ) {
      constructBlockInformation( forest );
#ifndef NDEBUG
      checkBlockInformationConsistency( forest );
#endif
   }

   registerRefreshTimer();
}



BlockForest::BlockForest( const uint_t process, const char* const filename, const bool broadcastFile, const bool keepGlobalBlockInformation ) :

   BlockStorage( AABB(), false, false, false ),

   process_( process ),  modificationStamp_( uint_t(0) ),
   recalculateBlockLevelsInRefresh_( true ), alwaysRebalanceInRefresh_( false ), allowMultipleRefreshCycles_( true ),
   reevaluateMinTargetLevelsAfterForcedRefinement_( false ),
   checkForEarlyOutInRefresh_( true ), checkForLateOutInRefresh_( true ), allowChangingDepth_( true ), checkForEarlyOutAfterLoadBalancing_( false ),
   phantomBlockMigrationIterations_( uint_t(0) ),
   nextCallbackBeforeBlockDataIsPackedHandle_( uint_t(0) ), nextCallbackBeforeBlockDataIsUnpackedHandle_( uint_t(0) ),
   nextCallbackAfterBlockDataIsUnpackedHandle_( uint_t(0) ),
   snapshotExists_( false ), snapshotDepth_( uint_t(0) ), snapshotBlockDataItems_( uint_t(0) )
{
   blockInformation_ = make_shared< BlockInformation >( *this );

   uint_t offset = 0;
   std::vector< uint8_t > buffer;

   if( broadcastFile && (mpi::MPIManager::instance()->numProcesses() > 1) )
   {
      std::ifstream file;
      uint_t length = 0;

      WALBERLA_ROOT_SECTION() {
         file.open( filename, std::ios::binary );
         if( file.fail() )
            WALBERLA_ABORT( "Loading BlockForest from file \'" << filename << "\' failed:\n"
                            "Opening the file failed. Does the file even exist?" );
         file.seekg( 0, std::ios::end );
         length = uint_c( static_cast< std::streamoff >( file.tellg() ) );
         file.seekg( 0, std::ios::beg );
      }
      MPI_Bcast( reinterpret_cast< void* >( &length ), 1, MPITrait< uint_t >::type(), 0, MPI_COMM_WORLD );

      buffer.resize( length );

      WALBERLA_ROOT_SECTION() {
         file.read( reinterpret_cast< char* >( &(buffer[0]) ), numeric_cast< std::streamsize >( length ) );
         file.close();
      }
      MPI_Bcast( reinterpret_cast< void* >( &(buffer[0]) ), numeric_cast< int >( length ), MPITrait< uint8_t >::type(), 0, MPI_COMM_WORLD );
   }
   else
   {
      std::ifstream file;
      file.open( filename, std::ios::binary );
      if( file.fail() )
         WALBERLA_ABORT( "Loading BlockForest from file \'" << filename << "\' failed:\n"
                         "Opening the file failed. Does the file even exist?" );

      file.seekg( 0, std::ios::end );
      const uint_t length = uint_c( static_cast< std::streamoff >( file.tellg() ) );
      file.seekg( 0, std::ios::beg );

      buffer.resize( length );

      file.read( reinterpret_cast< char* >( &(buffer[0]) ), numeric_cast< std::streamsize >( length ) );
      file.close();
   }

   // HEADER

   // domain AABB

   real_t domain[6];

   for( uint_t i = 0; i != 6; ++i ) {
      domain[i] = byteArrayToReal< real_t >( buffer, offset );
      offset += sizeof( real_t ) + 1 + 2;
   }

   domain_.initMinMaxCorner( domain[0], domain[1], domain[2], domain[3], domain[4], domain[5] );

   // number of coarse/root blocks in each direction

   for( uint_t i = 0; i != 3; ++i ) {
      size_[i] = byteArrayToUint( buffer, offset, 4 );
      offset += 4;
   }

   // domain periodicity

   for( uint_t i = 0; i != 3; ++i ) {
      periodic_[i] = ( byteArrayToUint( buffer, offset, 1 ) == uint_c(1) );
      ++offset;
   }

   // block forest depth (= number of levels - 1)

   depth_ = byteArrayToUint( buffer, offset, 1 );
   ++offset;

   // treeIdDigits (= number of bits used for storing the tree ID [tree ID marker + tree index])

   treeIdDigits_ = byteArrayToUint( buffer, offset, 1 );
   ++offset;

   // processIdBytes (= number of bytes required for storing process IDs)

   processIdBytes_ = byteArrayToUint( buffer, offset, 1 );
   ++offset;

   // insertBuffersIntoProcessNetwork?

   insertBuffersIntoProcessNetwork_ = ( byteArrayToUint( buffer, offset, 1 ) == uint_c(1) );
   ++offset;

   // number of processes

   const uint_t numberOfProcesses = byteArrayToUint( buffer, offset, 4 );
   offset += 4;

   if( numberOfProcesses != uint_c( MPIManager::instance()->numProcesses() ) )
      WALBERLA_ABORT( "Loading BlockForest from file \'" << filename << "\' failed:\n"
                      "The number of MPI processes (" << MPIManager::instance()->numProcesses() << ") does not match the number of "
                      "processes stored in this file (" << numberOfProcesses << ")" );

   if( process >= numberOfProcesses )
      WALBERLA_ABORT( "Loading BlockForest from file \'" << filename << "\' failed:\n"
                      "The process index \'" << process << "\' exceeds the number of available processes (" << numberOfProcesses << ")!" );

   // SUID MAPPING

   // number of SUIDs

   const uint_t numberOfSUIDs = byteArrayToUint( buffer, offset, 1 );
   ++offset;

   const uint_t suidBytes = ( ( numberOfSUIDs % 8 == 0 ) ? ( numberOfSUIDs / 8 ) : ( numberOfSUIDs / 8 + 1 ) );

   std::vector< SUID > suidMap;

   // for every SUID ...

   for( uint_t i = 0; i != numberOfSUIDs; ++i ) {

      // length of its identifier string

      const uint_t identifierLength = byteArrayToUint( buffer, offset, 1 );
      ++offset;

      // the identifier string

      const std::string identifier( reinterpret_cast< const char* >( &(buffer[offset]) ), identifierLength );
      offset += identifierLength;

      SUID suid( identifier, false );

      suidMap.push_back( suid );
   }

   // BLOCK DATA

   const uint_t blockIdBytes = getBlockIdBytes();

   // calculate offsets to block and neighborhood data of each process

   std::vector< uint_t > offsetBlocks( numberOfProcesses );
   std::vector< uint_t > offsetNeighbors( numberOfProcesses );

   for( uint_t i = 0; i != numberOfProcesses; ++i )
   {
      offsetBlocks[i] = offset;

      const uint_t numberOfBlocks = byteArrayToUint( buffer, offset, 2 );
      offset += 2 + numberOfBlocks * ( blockIdBytes + suidBytes );

      offsetNeighbors[i] = offset;

      offset += 2 + byteArrayToUint( buffer, offset, 2 ) * processIdBytes_;
   }

   // process neighborhood (= all neighboring processes)

   const uint_t numberOfNeighbors = byteArrayToUint( buffer, offsetNeighbors[ process_ ], 2 );

   for( uint_t i = 0; i != numberOfNeighbors; ++i )
      neighborhood_.push_back( byteArrayToUint( buffer, offsetNeighbors[ process_ ] + 2 + i * processIdBytes_, processIdBytes_ ) );

   // number of blocks associated with this process

   const uint_t numberOfBlocks = byteArrayToUint( buffer, offsetBlocks[ process_ ], 2 );

   if( numberOfBlocks > 0 )
   {
      // construct neighborhood reconstruction blocks (required for initializing the blocks of this process)

      BlockReconstruction::AABBReconstruction aabbReconstruction( domain_, size_[0], size_[1], size_[2], treeIdDigits_ );
      BlockReconstruction::NeighborhoodReconstruction< Block > neighborhoodReconstruction( domain_, periodic_[0], periodic_[1], periodic_[2] );

      std::vector< BlockReconstruction::NeighborhoodReconstructionBlock > neighbors;

      for( uint_t i = 0; i != numberOfBlocks; ++i ) {

         offset = offsetBlocks[ process_ ] + 2 + i * ( blockIdBytes + suidBytes );

         const BlockID id( buffer, offset, blockIdBytes );

         Set<SUID> state;
         std::vector< bool > suidBoolVec = byteArrayToBoolVector( buffer, offset + blockIdBytes, suidBytes );
         for( uint_t j = 0; j != suidBoolVec.size(); ++j ) {
            WALBERLA_ASSERT( !suidBoolVec[ j ] || j < suidMap.size() );
            if( suidBoolVec[j])
               state += suidMap[j];
         }

         neighbors.emplace_back( id, process_, state, aabbReconstruction );
      }

      for( uint_t i = 0; i != numberOfNeighbors; ++i ) {

         const uint_t neighborProcess = neighborhood_[i];
         const uint_t numberOfNeighborBlocks = byteArrayToUint( buffer, offsetBlocks[ neighborProcess ], 2 );

         for( uint_t j = 0; j != numberOfNeighborBlocks; ++j ) {

            offset = offsetBlocks[ neighborProcess ] + 2 + j * ( blockIdBytes + suidBytes );

            const BlockID id( buffer, offset, blockIdBytes );

            Set<SUID> state;
            std::vector< bool > suidBoolVec = byteArrayToBoolVector( buffer, offset + blockIdBytes, suidBytes );
            for( uint_t k = 0; k != suidBoolVec.size(); ++k ) {
               WALBERLA_ASSERT( !suidBoolVec[k] || k < suidMap.size() );
               if( suidBoolVec[k])
                  state += suidMap[k];
            }

            neighbors.emplace_back( id, neighborProcess, state, aabbReconstruction );
         }
      }

      // for each block ...

      for( uint_t i = 0; i != numberOfBlocks; ++i ) {

         offset = offsetBlocks[ process_ ] + 2 + i * ( blockIdBytes + suidBytes );

         // block ID

         const BlockID id( buffer, offset, blockIdBytes );

         // block state (SUID set)

         Set<SUID> state;
         std::vector< bool > suidBoolVec = byteArrayToBoolVector( buffer, offset + blockIdBytes, suidBytes );
         for( uint_t j = 0; j != suidBoolVec.size(); ++j ) {
            WALBERLA_ASSERT( !suidBoolVec[j] || j < suidMap.size() );
            if( suidBoolVec[j])
               state += suidMap[j];
         }

         // create block using the just constructed reconstruction information

         AABB aabb;
         const uint_t level = aabbReconstruction( aabb, id );

         auto block = std::make_shared< Block >( *this, id, aabb, state, level, neighborhoodReconstruction, neighbors );

         blocks_[ id ] = block;
      }
   }

   if( keepGlobalBlockInformation ) {

      std::vector< BlockID > ids;
      std::vector< shared_ptr< BlockInformation::Node > > nodes;

      for( uint_t i = 0; i != numberOfProcesses; ++i ) {

         const uint_t numBlocks = byteArrayToUint( buffer, offsetBlocks[ i ], 2 );

         for( uint_t j = 0; j != numBlocks; ++j ) {

            offset = offsetBlocks[ i ] + 2 + j * ( blockIdBytes + suidBytes );

            ids.emplace_back( buffer, offset, blockIdBytes );

            Set<SUID> state;
            std::vector< bool > suidBoolVec = byteArrayToBoolVector( buffer, offset + blockIdBytes, suidBytes );
            for( uint_t k = 0; k != suidBoolVec.size(); ++k ) {
               WALBERLA_ASSERT( !suidBoolVec[k] || k < suidMap.size() );
               if( suidBoolVec[k])
                  state += suidMap[k];
            }

            nodes.push_back( make_shared< BlockInformation::Node >( i, state ) );
         }
      }

      constructBlockInformation( ids, nodes );
   }

   registerRefreshTimer();
}



void BlockForest::getBlockID( IBlockID& id, const real_t x, const real_t y, const real_t z ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< BlockID* >( &id ), &id );

   if( blockInformation_->active() ) {
      if( !blockInformation_->getId( *static_cast< BlockID* >( &id ), x, y, z ) ) {
         WALBERLA_ABORT( "Getting block ID failed: There exists no block at global location (" << x << "," << y << "," << z <<")!" );
      }
   }
   else {
      const Block* const block = getBlock(x,y,z);
      if( block == nullptr ) {
         WALBERLA_ABORT( "Getting block ID failed: Locally, there exists no block at global location (" << x << "," << y << "," << z <<")!\n"
                         "                         (for simulation global information you have to explicitly construct the block forest to "
                         "contain global knowledge)");
      }
      id = block->getId();
   }
}



void BlockForest::getAABB( AABB& aabb, const IBlockID& id ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &id ), &id );

   if( blockInformation_->active() ) {
      const BlockID& bid = *static_cast< const BlockID* >( &id );
      if( !blockInformation_->exists( bid ) ) {
         WALBERLA_ABORT( "Getting block AABB failed: There exists no block with block ID \'" << id << "\'" );
      }
      getAABBFromBlockId( aabb, bid );
   }
   else {
      const Block* const block = getBlock( id );
      if( block == nullptr )
      {
         const BlockID& bid = *static_cast< const BlockID* >( &id );
         aabb = getAABBFromBlockId( bid );
      }
      else
      {
         aabb = block->getAABB();
      }
   }
}



void BlockForest::getState( Set<SUID>& state, const IBlockID& id ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &id ), &id );

   if( blockInformation_->active() ) {
      if( !blockInformation_->getState( state, *static_cast< const BlockID* >( &id ) ) ) {
         WALBERLA_ABORT( "Getting block state failed: There exists no block with block ID \'" << id << "\'" );
      }
   }
   else {
      const Block* const block = getBlock( id );
      if( block == nullptr ) {
         WALBERLA_ABORT( "Getting block state failed: Locally, there exists no block with block ID \'" << id << "\'\n"
                         "                            (for simulation global information you have to explicitly construct "
                         "the block forest to contain global knowledge)" );
      }
      state = block->getState();
   }
}



void BlockForest::getProcessRank( uint_t& rank, const IBlockID& id ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &id ), &id );

   if( blockInformation_->active() ) {
      if( !blockInformation_->getProcess( rank, *static_cast< const BlockID* >( &id ) ) ) {
         WALBERLA_ABORT( "Getting block process rank failed: There exists no block with block ID \'" << id << "\'" );
      }
   }
   else {
      const Block* const block = getBlock( id );
      if( block == nullptr ) {
         WALBERLA_ABORT( "Getting block process rank failed: Locally, there exists no block with block ID \'" << id << "\'\n"
                         "                                   (for simulation global information you have to explicitly construct "
                         "the block forest to contain global knowledge)" );
      }
      rank = block->getProcess();
   }
}



void BlockForest::getRootBlockAABB( AABB& aabb, const uint_t x, const uint_t y, const uint_t z ) const {

   if( blockInformation_->active() ) {

      if( !blockInformation_->rootBlockExists(x,y,z) ) {
         WALBERLA_ABORT( "Getting root block AABB failed: There exists no root block [" << x << "," << y << "," << z <<"]!" );
      }

      const real_t xw = getRootBlockXSize();
      const real_t yw = getRootBlockYSize();
      const real_t zw = getRootBlockZSize();

      aabb.initMinMaxCorner( domain_.xMin() + static_cast< real_t >(  x  ) * xw,
                             domain_.yMin() + static_cast< real_t >(  y  ) * yw,
                             domain_.zMin() + static_cast< real_t >(  z  ) * zw,
                             ( x+1 == size_[0] ) ? domain_.xMax() : domain_.xMin() + static_cast< real_t >( x+1 ) * xw,
                             ( y+1 == size_[1] ) ? domain_.yMax() : domain_.yMin() + static_cast< real_t >( y+1 ) * yw,
                             ( z+1 == size_[2] ) ? domain_.zMax() : domain_.zMin() + static_cast< real_t >( z+1 ) * zw );
   }
   else {
      const Block* const block = getRootBlock(x,y,z);
      if( block == nullptr ) {
         WALBERLA_ABORT( "Getting root block AABB failed: Locally, there exists no root block [" << x << "," << y << "," << z <<"]!\n"
                         "                                (for simulation global information you have to explicitly construct "
                         "the block forest to contain global knowledge)" );
      }
      aabb = block->getAABB();
   }
}



void BlockForest::getRootBlockState( Set<SUID>& state, const uint_t x, const uint_t y, const uint_t z ) const {

   if( blockInformation_->active() ) {
      if( !blockInformation_->getRootBlockState( state, x, y, z ) ) {
         WALBERLA_ABORT( "Getting root block state failed: There exists no root block [" << x << "," << y << "," << z <<"]!" );
      }
   }
   else {
      const Block* const block = getRootBlock(x,y,z);
      if( block == nullptr ) {
         WALBERLA_ABORT( "Getting root block state failed: Locally, there exists no root block [" << x << "," << y << "," << z <<"]!\n"
                         "                                 (for simulation global information you have to explicitly construct "
                         "the block forest to contain global knowledge)" );
      }
      state = block->getState();
   }
}



void BlockForest::getRootBlockProcessRank( uint_t& rank, const uint_t x, const uint_t y, const uint_t z ) const {

   if( blockInformation_->active() ) {
      if( !blockInformation_->getRootBlockProcess( rank, x, y, z ) ) {
         WALBERLA_ABORT( "Getting root block process rank failed: There exists no root block [" << x << "," << y << "," << z <<"]!" );
      }
   }
   else {
      const Block* const block = getRootBlock(x,y,z);
      if( block == nullptr ) {
         WALBERLA_ABORT( "Getting root block process rank failed: Locally, there exists no root block [" << x << "," << y << "," << z <<"]!\n"
                         "                                        (for simulation global information you have to explicitly construct "
                         "the block forest to contain global knowledge)" );
      }
      rank = block->getProcess();
   }
}



void BlockForest::getForestCoordinates( uint_t& x, uint_t& y, uint_t& z, const BlockID& id ) const {

   WALBERLA_ASSERT_GREATER_EQUAL( id.getUsedBits(), treeIdDigits_ );
   WALBERLA_ASSERT_EQUAL( ( id.getUsedBits() - treeIdDigits_ ) % 3, 0 );

   BlockID blockId( id );

   const uint_t levels = ( blockId.getUsedBits() - treeIdDigits_ ) / 3;
   for( uint_t i = levels; i-- != 0; )
      blockId.removeBranchId();

   uint_t index = blockId.getTreeIndex();

       z = index / ( size_[0] * size_[1] );
   index = index % ( size_[0] * size_[1] );
       y = index / size_[0];
       x = index % size_[0];
}



std::map< uint_t, std::vector< Vector3<real_t> > > BlockForest::getNeighboringProcessOffsets() const
{
   std::map< uint_t, std::vector< Vector3<real_t> > > offsets;

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      const auto & block = it->second;
      const AABB & blockAABB = block->getAABB();

      const real_t eps[] = { real_c( 1.0E-6 ) * ( blockAABB.xMax() - blockAABB.xMin() ),
                             real_c( 1.0E-6 ) * ( blockAABB.yMax() - blockAABB.yMin() ),
                             real_c( 1.0E-6 ) * ( blockAABB.zMax() - blockAABB.zMin() ) };

      uint_t i[] = { uint_c(0), uint_c(0), uint_c(0) };

      for( i[2] = 0; i[2] != 3; ++i[2] ) {
         for( i[1] = 0; i[1] != 3; ++i[1] ) {
            for( i[0] = 0; i[0] != 3; ++i[0] )
            {
               if( i[0] == 1 && i[1] == 1 && i[2] == 1 )
                  continue;

               const auto sectionIndex = getBlockNeighborhoodSectionIndex( i[0], i[1], i[2] );
               const auto neighbors = block->getNeighborhoodSection( sectionIndex );

               for( auto neighbor = neighbors.begin(); neighbor != neighbors.end(); ++neighbor )
               {
                  if( (*neighbor)->getProcess() == getProcess() )
                     continue;

                  AABB aabb;
                  getAABBFromBlockId( aabb, (*neighbor)->getId() );

                  Vector3< real_t > offset;

                  for( uint_t j = 0; j != 3; ++j )
                  {
                     if( i[j] == 0 ) offset[j] = aabb.max(j) - blockAABB.min(j);
                     else if( i[j] == 2 ) offset[j] = aabb.min(j) - blockAABB.max(j);
                     else offset[j] = real_c(0);

                     if( realIsEqual( offset[j], real_c(0), eps[j] ) )
                        offset[j] = real_c(0);

                     WALBERLA_ASSERT( realIsIdentical( offset[j], real_c(0) ) ||
                             ( isPeriodic(j) && realIsEqual( std::fabs( offset[j] ), getDomain().max(j) - getDomain().min(j), eps[j] ) ) );
                  }

                  const auto process = (*neighbor)->getProcess();
                  auto & processOffsets = offsets[ process ];

                  bool newOffset = true;
                  for( auto vec = processOffsets.begin(); vec != processOffsets.end() && newOffset; ++vec )
                  {
                     if( realIsEqual( (*vec)[0], offset[0], eps[0] ) && realIsEqual( (*vec)[1], offset[1], eps[1] ) &&
                         realIsEqual( (*vec)[2], offset[2], eps[2] ) ) newOffset = false;
                  }
                  if( newOffset )
                     processOffsets.push_back( offset );
               }
            }
         }
      }
   }
   return offsets;
}



void BlockForest::refresh()
{
   WALBERLA_LOG_PROGRESS( "BlockForest refresh: starting distributed refresh of the block structure" );

   bool rebalanceAndRedistribute( true );
   bool additionalRefreshCycleRequired( false );

   refreshTiming_[ "block level determination" ].start();

   for( auto block = blocks_.begin(); block != blocks_.end(); ++block )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block->second );
      block->second->setTargetLevel( block->second->getLevel() );
   }

   if( recalculateBlockLevelsInRefresh_ )
   {
      WALBERLA_LOG_PROGRESS( "BlockForest refresh: determining new block levels" );

      bool rerun( true );
      while( rerun )
      {
         rebalanceAndRedistribute = determineBlockTargetLevels( additionalRefreshCycleRequired, rerun );
      }

      if( !rebalanceAndRedistribute )
         WALBERLA_LOG_PROGRESS( "BlockForest refresh: block levels do not change" );
   }
   else
   {
      WALBERLA_LOG_PROGRESS( "BlockForest refresh: block levels do not change" );
   }

   refreshTiming_[ "block level determination" ].end();

   if( alwaysRebalanceInRefresh_ && refreshPhantomBlockMigrationPreparationFunction_ )
      rebalanceAndRedistribute = true;

   while( rebalanceAndRedistribute )
   {
      // create phantom block forest

      WALBERLA_LOG_PROGRESS( "BlockForest refresh: creating phantom forest/blocks" );

      refreshTiming_[ "phantom forest creation" ].start();

      PhantomBlockForest phantomForest( *this );
      phantomForest.initialize( refreshBlockStateDeterminationFunction_, allowChangingDepth_ );

      refreshTiming_[ "phantom forest creation" ].end();

      // move phantom blocks between processes (= dynamic load balancing)

      WALBERLA_MPI_SECTION()
      {
         refreshTiming_[ "phantom block redistribution (= load balancing)" ].start();

         if( refreshPhantomBlockMigrationPreparationFunction_ )
         {
            phantomForest.assignBlockData( refreshPhantomBlockDataAssignmentFunction_ );

            WALBERLA_LOG_PROGRESS( "BlockForest refresh: performing phantom block redistribution/load balancing" );

            uint_t iteration = uint_t(0);
            bool runAgain( true );
            while( runAgain )
            {
               WALBERLA_LOG_PROGRESS( "BlockForest refresh: decide about which phantom blocks need to migrate" );
               runAgain = phantomForest.calculateMigrationInformation( refreshPhantomBlockMigrationPreparationFunction_, iteration );
               WALBERLA_LOG_PROGRESS( "BlockForest refresh: migrate phantom blocks" );
               phantomForest.migrate( refreshPhantomBlockDataPackFunction_, refreshPhantomBlockDataUnpackFunction_ );
               ++iteration;
            }
            phantomBlockMigrationIterations_ = iteration;

            WALBERLA_LOG_PROGRESS( "BlockForest refresh: phantom block redistribution/load balancing finished after " << phantomBlockMigrationIterations_ << " iterations" );
         }

         refreshTiming_[ "phantom block redistribution (= load balancing)" ].end();
      }

      // update block forest: transfer block data (includes migrating, deleting and creating blocks), adapt depth, adapt process neighborhood, etc.

      bool performUpdate( true );
      if( checkForEarlyOutAfterLoadBalancing_ )
      {
         performUpdate = false;
         const auto & phantomBlocks = phantomForest.getBlockMap();
         for( auto phantom = phantomBlocks.begin(); phantom != phantomBlocks.end() && !performUpdate; ++phantom )
         {
            const auto & sourceProcess = phantom->second->getSourceProcess();
            for( auto it = sourceProcess.begin(); it != sourceProcess.end() && !performUpdate; ++it )
               performUpdate = ( *it != process_ );
         }
         mpi::allReduceInplace( performUpdate, mpi::LOGICAL_OR );
      }

      if( performUpdate )
      {
         WALBERLA_LOG_PROGRESS( "BlockForest refresh: update block structure" );

         refreshTiming_[ "block structure update (includes data migration)" ].start();
         update( phantomForest );
         refreshTiming_[ "block structure update (includes data migration)" ].end();

         WALBERLA_LOG_PROGRESS( "BlockForest refresh: updating block structure finished" );
      }
      else
      {
         WALBERLA_LOG_PROGRESS( "BlockForest refresh: block structure update is not required, skipping update" );
      }

      // ----- block forest is again in a valid state -----

      if( additionalRefreshCycleRequired )
      {
         WALBERLA_ASSERT( recalculateBlockLevelsInRefresh_ );
         WALBERLA_ASSERT( allowMultipleRefreshCycles_ );

         refreshTiming_[ "block level determination" ].start();

         for( auto block = blocks_.begin(); block != blocks_.end(); ++block )
         {
            WALBERLA_ASSERT_NOT_NULLPTR( block->second );
            block->second->setTargetLevel( block->second->getLevel() );
         }

         WALBERLA_LOG_PROGRESS( "BlockForest refresh: determining new block levels (again -> more refresh cycles required!)" );

         bool rerun( true );
         while( rerun )
         {
            rebalanceAndRedistribute = determineBlockTargetLevels( additionalRefreshCycleRequired, rerun );
         }

         if( !rebalanceAndRedistribute )
            WALBERLA_LOG_PROGRESS( "BlockForest refresh: block levels do not change" );

         refreshTiming_[ "block level determination" ].end();
      }
      else
         rebalanceAndRedistribute = false;
   }

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: finished!" );
}



void BlockForest::createSnapshot( const std::vector<uint_t> & sendTo, const std::vector<uint_t> & recvFrom )
{
   WALBERLA_LOG_PROGRESS( "BlockForest create snapshot: schedule MPI receives (1)" );

   std::vector< MPI_Request > request( recvFrom.size() + sendTo.size() );

   // schedule recvs for snapshot sizes

   std::vector< uint_t > msgSize( recvFrom.size() );

   for( uint_t i = uint_t(0); i != recvFrom.size(); ++i )
   {
      MPI_Irecv( static_cast< void * >(&(msgSize[i])), 1, MPITrait< uint_t >::type(), int_c(recvFrom[i]), 880, // snapshot = 115 110 097 112 115 104 111 116
                 mpi::MPIManager::instance()->comm(), &(request[i]) );
   }

   // local snapshot

   WALBERLA_LOG_PROGRESS( "BlockForest create snapshot: creating local snapshot" );

   mpi::SendBuffer buffer;

   buffer << neighborhood_ << uint_c( blocks_.size() );
   for( auto block = blocks_.begin(); block != blocks_.end(); ++block )
   {
      WALBERLA_ASSERT_EQUAL( block->first, block->second->getId() );

      buffer << block->first;
      block->second->toBuffer( buffer );

      for( auto dataItem = blockDataItem_.begin(); dataItem != blockDataItem_.end(); ++dataItem )
      {
         auto blockDataHandlingWrapper = dataItem->getDataHandling( block->second.get() );
         if( blockDataHandlingWrapper )
            blockDataHandlingWrapper->serialize( block->second.get(), dataItem->getId(), buffer );
      }
   }

   // send snapshot size

   WALBERLA_LOG_PROGRESS( "BlockForest create snapshot: sending snapshot size" );

   for( uint_t i = uint_t(0); i != sendTo.size(); ++i )
   {
      uint_t bufferSize = buffer.size();
      MPI_Isend( static_cast< void * >( &bufferSize ), 1, MPITrait< uint_t >::type(), int_c(sendTo[i]), 880,
                 mpi::MPIManager::instance()->comm(), &(request[ recvFrom.size() + i ]) );
   }

   // wait for sizes to be exchanged

   if( ! request.empty() )
      MPI_Waitall( int_c( request.size() ), &(request[0]), MPI_STATUSES_IGNORE );

   // schedule recvs

   WALBERLA_LOG_PROGRESS( "BlockForest create snapshot: schedule MPI receives (2)" );

   std::vector< mpi::RecvBuffer > recvBuffer( recvFrom.size() );

   for( uint_t i = uint_t(0); i != recvFrom.size(); ++i )
   {
      recvBuffer[i].resize( msgSize[i] );
      MPI_Irecv( static_cast< void * >( recvBuffer[i].ptr() ), int_c( msgSize[i] ), MPITrait< mpi::RecvBuffer::ElementType >::type(),
                 int_c(recvFrom[i]), 881, mpi::MPIManager::instance()->comm(), &(request[i]) );
   }

   // send data

   WALBERLA_LOG_PROGRESS( "BlockForest create snapshot: sending snapshot" );

   for( uint_t i = uint_t(0); i != sendTo.size(); ++i )
   {
      MPI_Isend( static_cast< void * >( buffer.ptr() ), int_c( buffer.size() ), MPITrait< mpi::SendBuffer::ElementType >::type(),
                 int_c(sendTo[i]), 881, mpi::MPIManager::instance()->comm(), &(request[ recvFrom.size() + i ]) );
   }

   // wait for everything to finish

   if( ! request.empty() )
      MPI_Waitall( int_c( request.size() ), &(request[0]), MPI_STATUSES_IGNORE );

   // make sure everybody is ready to store the new snapshot
   // if a MPI error occurs, we cannot store this snapshot, but must rely on an older snapshot

   WALBERLA_LOG_PROGRESS( "BlockForest create snapshot: waiting for all processes to finish MPI calls" );

   WALBERLA_MPI_BARRIER()

   // store snapshot

   WALBERLA_LOG_PROGRESS( "BlockForest create snapshot: storing snapshot" );

   snapshotDepth_ = depth_;

   WALBERLA_ASSERT_EQUAL( process_, uint_c( mpi::MPIManager::instance()->rank() ) );
   mpi::RecvBuffer processBuffer( buffer );
   snapshot_[ process_ ] = processBuffer;

   for( uint_t i = uint_t(0); i != recvFrom.size(); ++i )
   {
      snapshot_[ recvFrom[i] ] = recvBuffer[i];
      recvBuffer[i].reset();
   }

   snapshotBlockDataItems_ = blockDataItem_.size();
   snapshotExists_ = true;
}



void BlockForest::restoreSnapshot( const SnapshotRestoreFunction & processMapping, const bool rebalance )
{
   WALBERLA_CHECK( snapshotExists_ );
   WALBERLA_CHECK_EQUAL( snapshotBlockDataItems_, blockDataItem_.size() );

   WALBERLA_LOG_PROGRESS( "BlockForest restore snapshot: clearing current data" );

   blocks_.clear();
   neighborhood_.clear();

   WALBERLA_LOG_PROGRESS( "BlockForest restore snapshot: restoring block data from last snapshot" );

   process_ = uint_c( mpi::MPIManager::instance()->rank() );
   depth_ = snapshotDepth_;

   std::set< uint_t > neighborhood;

   for( auto it = snapshot_.begin(); it != snapshot_.end(); ++it )
   {
      if( processMapping( it->first ) == process_ )
      {
         mpi::RecvBuffer & buffer = it->second;

         std::vector< uint_t > pNeighborhood;
         buffer >> pNeighborhood;

         for( auto n = pNeighborhood.begin(); n != pNeighborhood.end(); ++n )
         {
            WALBERLA_ASSERT_LESS( processMapping(*n), uint_c( mpi::MPIManager::instance()->numProcesses() ) );
            neighborhood.insert( processMapping(*n) );
         }

         uint_t blocks( uint_t(0) );
         buffer >> blocks;

         for( uint_t i = uint_t(0); i != blocks; ++i )
         {
            BlockID id;
            buffer >> id;

            AABB aabb;
            const uint_t level = getAABBFromBlockId( aabb, id );

            WALBERLA_ASSERT( blocks_.find( id ) == blocks_.end() );
            blocks_[ id ] = std::make_shared< Block >( *this, id, aabb, level, buffer, processMapping );

            Block * block = blocks_[ id ].get();
            for( auto dataItem = blockDataItem_.begin(); dataItem != blockDataItem_.end(); ++dataItem )
            {
               auto blockDataHandlingWrapper = dataItem->getDataHandling( block );
               if( blockDataHandlingWrapper )
               {
                  addBlockData( block, dataItem->getId(), blockDataHandlingWrapper->deserialize( block ) );
                  blockDataHandlingWrapper->deserialize( block, dataItem->getId(), buffer );
               }
               else
               {
                  addBlockData( block, dataItem->getId(), nullptr );
               }
            }
         }
      }
   }

   WALBERLA_LOG_PROGRESS( "BlockForest restore snapshot: restoring process neighborhood" );

   for( auto n = neighborhood.begin(); n != neighborhood.end(); ++n )
      neighborhood_.push_back( *n );

   rebuildProcessesWithBlocksCommunicator();

   if( containsGlobalBlockInformation() )
   {
      WALBERLA_LOG_PROGRESS( "BlockForest restore snapshot: constructing global block information" );
      constructBlockInformation();
   }

   ++modificationStamp_;

   if( ! callbackAfterBlockDataIsUnpacked_.empty() )
      WALBERLA_LOG_PROGRESS( "BlockForest restore snapshot: executing call back functions after data was restored" );

   for( auto f = callbackAfterBlockDataIsRestored_.begin(); f != callbackAfterBlockDataIsRestored_.end(); ++f )
      f->second();

   WALBERLA_MPI_SECTION()
   {
      if( rebalance && refreshPhantomBlockMigrationPreparationFunction_ )
      {
         // create phantom block forest

         WALBERLA_LOG_PROGRESS( "BlockForest restore snapshot: starting data structure refresh" );
         WALBERLA_LOG_PROGRESS( "BlockForest refresh: creating phantom forest/blocks" );
         PhantomBlockForest phantomForest( *this );
         phantomForest.initialize( PhantomBlockForest::BlockStateDeterminationFunction(), false );

         // move phantom blocks between processes (= dynamic load balancing)

         if( refreshPhantomBlockMigrationPreparationFunction_ )
         {
            phantomForest.assignBlockData( refreshPhantomBlockDataAssignmentFunction_ );

            WALBERLA_LOG_PROGRESS( "BlockForest refresh: performing phantom block redistribution/load balancing" );

            uint_t iteration = uint_t(0);
            bool runAgain( true );
            while( runAgain )
            {
               WALBERLA_LOG_PROGRESS( "BlockForest refresh: decide about which phantom blocks need to migrate" );
               runAgain = phantomForest.calculateMigrationInformation( refreshPhantomBlockMigrationPreparationFunction_, iteration );
               WALBERLA_LOG_PROGRESS( "BlockForest refresh: migrate phantom blocks" );
               phantomForest.migrate( refreshPhantomBlockDataPackFunction_, refreshPhantomBlockDataUnpackFunction_ );
               ++iteration;
            }
            phantomBlockMigrationIterations_ = iteration;

            WALBERLA_LOG_PROGRESS( "BlockForest refresh: phantom block redistribution/load balancing finished after " << phantomBlockMigrationIterations_ << " iterations" );
         }

         // update block forest: transfer block data (includes migrating, deleting and creating blocks), adapt depth, adapt process neighborhood, etc.

         bool performUpdate( true );
         if( checkForEarlyOutAfterLoadBalancing_ )
         {
            performUpdate = false;
            const auto & phantomBlocks = phantomForest.getBlockMap();
            for( auto phantom = phantomBlocks.begin(); phantom != phantomBlocks.end() && !performUpdate; ++phantom )
            {
               const auto & sourceProcess = phantom->second->getSourceProcess();
               for( auto it = sourceProcess.begin(); it != sourceProcess.end(); ++it )
                  performUpdate = ( *it != process_ );
            }
            mpi::allReduceInplace( performUpdate, mpi::LOGICAL_OR );
         }

         if( performUpdate )
         {
            WALBERLA_LOG_PROGRESS( "BlockForest refresh: update block structure" );
            update( phantomForest );
            WALBERLA_LOG_PROGRESS( "BlockForest refresh: updating block structure finished" );
         }
         else
         {
            WALBERLA_LOG_PROGRESS( "BlockForest refresh: block structure update is not required, skipping update" );
         }

         WALBERLA_LOG_PROGRESS( "BlockForest restore snapshot: data structure refresh finished" );
      }
   }
}



void BlockForest::saveToFile( const std::string & filename, FileIOMode fileIOMode ) const
{
   Set<SUID> suids;
   for( auto block = blocks_.begin(); block != blocks_.end(); ++block )
      suids += block->second->getState();

   std::vector<SUID> localSuids( suids.begin(), suids.end() );
   std::vector<SUID> allSuids = mpi::allReduceSet( localSuids, mpi::UNION );

   suids.clear();
   suids.insert( allSuids.begin(), allSuids.end() );

   saveToFile( filename, suids, fileIOMode );
}



/// ATTENTION: 'blockStates' must be identical for every process!
void BlockForest::saveToFile( const std::string & filename, const Set<SUID> & blockStates, FileIOMode fileIOMode ) const
{
   WALBERLA_CHECK_LESS( blockStates.size(), uint_c(256), "When saving the block structure to file, only 255 different SUIDs are allowed!" );

   const uint_t suidBytes = ( ( blockStates.size() % 8 == 0 ) ? ( blockStates.size() / 8 ) : ( blockStates.size() / 8 + 1 ) );

   std::map< SUID, std::vector< bool > > suidMap;

   uint_t i = 0;
   for( Set<SUID>::const_iterator it = blockStates.begin(); it != blockStates.end(); ++it ) {

      std::vector< bool > suidBoolVec( 8 * suidBytes );
      suidBoolVec[i] = true;
      suidMap[ *it ] = suidBoolVec;
      ++i;
   }

   saveToFile( filename, fileIOMode, suidMap, suidBytes );
}



bool BlockForest::equal( const BlockStorage* rhs ) const {

   const BlockForest* forest = dynamic_cast< const BlockForest* >( rhs );

   if( forest != rhs )
      return false;

   if( process_ != forest->process_ || processIdBytes_ != forest->processIdBytes_ ||
       size_[0] != forest->size_[0] || size_[1] != forest->size_[1] || size_[2] != forest->size_[2] ||
       depth_ != forest->depth_ || treeIdDigits_ != forest->treeIdDigits_ || blocks_.size() != forest->blocks_.size() )
      return false;

   auto    it = blocks_.begin();
   auto rhsit = forest->blocks_.begin();

   while( it != blocks_.end() ) {
      if( it->first != rhsit->first || *(it->second) != *(rhsit->second) )
         return false;
      ++it;
      ++rhsit;
   }

   if( insertBuffersIntoProcessNetwork_ != forest->insertBuffersIntoProcessNetwork_ || neighborhood_.size() != forest->neighborhood_.size() )
      return false;

   for( uint_t i = 0; i != neighborhood_.size(); ++i )
      if( neighborhood_[i] != forest->neighborhood_[i] ) return false;

   return *blockInformation_ == *(forest->blockInformation_);
}



void BlockForest::constructBlockInformation( const SetupBlockForest & forest )
{
   std::vector< const SetupBlock * > blocks;
   forest.getBlocks( blocks );

   blockInformation_->clear();

   blockInformation_->nodes_.resize( forest.getNumberOfTrees() );
   auto & nodeForest = blockInformation_->nodes_;

   for( uint_t i = 0; i != blocks.size(); ++i )
   {
      const SetupBlock * const block = blocks[i];
      BlockID id( block->getId() );

      std::stack< uint_t > index;

      for( uint_t k = 0; k != block->getLevel(); ++k )
      {
         index.push( id.getBranchId() );
         id.removeBranchId();
      }
      index.push( id.getTreeIndex() );

      if( !(nodeForest[ index.top() ]) )
      {
         if( index.size() == 1 )
            nodeForest[ index.top() ] = make_shared< BlockInformation::Node >( block->getProcess(), block->getState() );
         else
            nodeForest[ index.top() ] = make_shared< BlockInformation::Node >();
      }

      auto node = nodeForest[ index.top() ];
      index.pop();

      while( !index.empty() )
      {
         //id.appendBranchId( index.top() );

         if( index.size() == 1 ) {
            node->setChild( index.top(), make_shared< BlockInformation::Node >( block->getProcess(), block->getState() ) );
         }
         else if( index.top() >= node->children_.size() || node->children_[ index.top() ] == nullptr ) {
            node->setChild( index.top(), make_shared< BlockInformation::Node >() );
         }

         WALBERLA_ASSERT_LESS( index.top(), node->children_.size() );

         node = node->children_[ index.top() ];
         index.pop();
      }
   }
}



void BlockForest::constructBlockInformation( const std::vector< BlockID > & ids,
                                             const std::vector< shared_ptr< BlockInformation::Node > > & nodes )
{
   WALBERLA_ASSERT_EQUAL( ids.size(), nodes.size() );

   blockInformation_->clear();

   blockInformation_->nodes_.resize( size_[0] * size_[1] * size_[2] );
   std::vector< shared_ptr< BlockInformation::Node > > & nodeForest = blockInformation_->nodes_;

   for( uint_t i = 0; i != nodes.size(); ++i )
   {
      BlockID id( ids[i] );

      std::stack< uint_t > index;

      const uint_t level = getLevelFromBlockId(id);
      for( uint_t k = 0; k != level ; ++k )
      {
         index.push( id.getBranchId() );
         id.removeBranchId();
      }
      index.push( id.getTreeIndex() );

      if( !(nodeForest[ index.top() ]) )
      {
         if( index.size() == 1 )
            nodeForest[ index.top() ] = nodes[i];
         else
            nodeForest[ index.top() ] = make_shared< BlockInformation::Node >();
      }

      auto node = nodeForest[ index.top() ];
      index.pop();

      while( !index.empty() )
      {
         //id.appendBranchId( index.top() );

         if( index.size() == 1 ) {
            node->setChild( index.top(), nodes[i] );
         }
         else if( index.top() >= node->children_.size() || node->children_[ index.top() ] == nullptr ) {
            node->setChild( index.top(), make_shared< BlockInformation::Node >() );
         }

         WALBERLA_ASSERT_LESS( index.top(), node->children_.size() );

         node = node->children_[ index.top() ];
         index.pop();
      }
   }
}



void BlockForest::constructBlockInformation()
{
   std::vector< std::pair< BlockID, std::pair< uint_t, Set<SUID> > > > data;
   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      data.emplace_back( it->first, std::make_pair( process_, it->second->getState() ) );
   }

   mpi::SendBuffer sBuffer;
   mpi::RecvBuffer rBuffer;

   sBuffer << data;
   data.clear();

   mpi::allGathervBuffer( sBuffer, rBuffer );
   sBuffer.reset();

   std::vector< BlockID > ids;
   std::vector< shared_ptr< BlockInformation::Node > > nodes;

   for( int r = 0; r != MPIManager::instance()->numProcesses(); ++r )
   {
      rBuffer >> data;
      for( auto it = data.begin(); it != data.end(); ++it )
      {
         ids.push_back( it->first );
         nodes.push_back( make_shared< BlockInformation::Node >( it->second.first, it->second.second ) );
      }
      data.clear();
   }
   rBuffer.reset();

   constructBlockInformation( ids, nodes );
}



void BlockForest::registerRefreshTimer()
{
   if( ! refreshTiming_.timerExists( "block level determination" ) )
      refreshTiming_.registerTimer( "block level determination" );
   if( ! refreshTiming_.timerExists( "block level determination (callback function)" ) )
      refreshTiming_.registerTimer( "block level determination (callback function)" );
   if( ! refreshTiming_.timerExists( "phantom forest creation" ) )
      refreshTiming_.registerTimer( "phantom forest creation" );
   if( ! refreshTiming_.timerExists( "phantom block redistribution (= load balancing)" ) )
      refreshTiming_.registerTimer( "phantom block redistribution (= load balancing)" );
   if( ! refreshTiming_.timerExists( "block structure update (includes data migration)" ) )
      refreshTiming_.registerTimer( "block structure update (includes data migration)" );
}



bool BlockForest::determineBlockTargetLevels( bool & additionalRefreshCycleRequired, bool & rerun )
{
   // TODO: use OpenMP buffer systems?

   rerun = false;

   if( !refreshMinTargetLevelDeterminationFunction_ )
      return false;

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - calling callback function to determine new block levels" );

   std::vector< std::pair< const Block *, uint_t > > minTargetLevelsCallback;
   std::vector< const Block * > blocksAlreadyMarkedForRefinement;
   std::vector< int > mapping;

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( it->second );
      if( it->second->getTargetLevel() > it->second->getLevel() )
      {
         WALBERLA_ASSERT_EQUAL( it->second->getTargetLevel(), it->second->getLevel() + uint_t(1) );
         blocksAlreadyMarkedForRefinement.push_back( it->second.get() );
         mapping.push_back(1);
      }
      else
      {
         WALBERLA_ASSERT( it->second->getTargetLevel() == it->second->getLevel() ||
                          ( it->second->getTargetLevel() + uint_t(1) ) == it->second->getLevel() );
         minTargetLevelsCallback.emplace_back( it->second.get(), it->second->getTargetLevel() );
         mapping.push_back(0);
      }
   }

   refreshTiming_[ "block level determination" ].end();
   refreshTiming_[ "block level determination (callback function)" ].start();
   refreshMinTargetLevelDeterminationFunction_( minTargetLevelsCallback, blocksAlreadyMarkedForRefinement, *this );
   refreshTiming_[ "block level determination (callback function)" ].end();
   refreshTiming_[ "block level determination" ].start();

   std::vector< std::pair< const Block *, uint_t > > minTargetLevelsAllBlocks;

   auto it0 = minTargetLevelsCallback.begin();
   auto it1 = blocksAlreadyMarkedForRefinement.begin();
   for( auto it = mapping.begin(); it != mapping.end(); ++it )
   {
      if( *it == 0 )
      {
         WALBERLA_CHECK( it0 != minTargetLevelsCallback.end() );
         minTargetLevelsAllBlocks.push_back( *it0 );
         it0++;
      }
      else
      {
         WALBERLA_ASSERT_EQUAL( *it, 1 );
         WALBERLA_CHECK( it1 != blocksAlreadyMarkedForRefinement.end() );
         WALBERLA_CHECK_NOT_NULLPTR( *it1 );
         WALBERLA_CHECK_EQUAL( (*it1)->getTargetLevel(), (*it1)->getLevel() + uint_t(1) );
         minTargetLevelsAllBlocks.emplace_back( *it1, (*it1)->getTargetLevel() );
         it1++;
      }
   }
   WALBERLA_CHECK( it0 == minTargetLevelsCallback.end() );
   WALBERLA_CHECK( it1 == blocksAlreadyMarkedForRefinement.end() );
   WALBERLA_CHECK_EQUAL( minTargetLevelsAllBlocks.size(), blocks_.size() );

   bool levelChangesPossible( false );

   std::map< BlockID, uint_t > minTargetLevels;
   auto bit = blocks_.begin();
   for( auto it = minTargetLevelsAllBlocks.begin(); it != minTargetLevelsAllBlocks.end(); ++it )
   {
      WALBERLA_CHECK_NOT_NULLPTR( it->first );
      WALBERLA_CHECK( bit != blocks_.end() );
      WALBERLA_CHECK_EQUAL( bit->second.get(), it->first );
      if( !allowMultipleRefreshCycles_ )
         WALBERLA_CHECK_LESS_EQUAL( it->second, it->first->getLevel() + uint_t(1) );
      if( !allowChangingDepth_ )
         WALBERLA_CHECK_LESS_EQUAL( it->second, depth_ );

      minTargetLevels[ it->first->getId() ] = it->second;

      if( it->first->getLevel() != it->second )
         levelChangesPossible = true;

      ++bit;
   }

   if( checkForEarlyOutInRefresh_ )
   {
      mpi::allReduceInplace( levelChangesPossible, mpi::LOGICAL_OR );
      if( !levelChangesPossible )
         return false;
   }

   std::map< uint_t, std::map< BlockID, uint_t > > localTargetLevelsForNeighbors;
   std::set< mpi::MPIRank > ranksToRecvFrom;

   additionalRefreshCycleRequired = false;

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( it->second );

      auto & block = it->second;
      const auto & id = block->getId();

      const uint_t level = block->getLevel();
      uint_t minTargetLevel = minTargetLevels[id];

      if( allowMultipleRefreshCycles_ )
      {
         if( minTargetLevel > ( level + uint_t(1) ) )
         {
            additionalRefreshCycleRequired = true;
            minTargetLevel = level + uint_t(1);
            minTargetLevels[id] = minTargetLevel;
         }
      }
#ifndef NDEBUG
      else
      {
         WALBERLA_ASSERT_LESS_EQUAL( minTargetLevel, level + uint_t(1) );
      }
#endif

      const uint_t targetLevel = std::max( level, minTargetLevel );
      block->setTargetLevel( targetLevel );

      WALBERLA_ASSERT_EQUAL( block->getProcess(), process_ );

      for( uint_t n = 0; n != block->getNeighborhoodSize(); ++n )
      {
         const uint_t np = block->getNeighborProcess(n);
         if( np != process_ )
         {
            localTargetLevelsForNeighbors[ np ][ id ] = targetLevel;
            ranksToRecvFrom.insert( numeric_cast< mpi::MPIRank >( np ) );
         }
      }
   }

   if( allowMultipleRefreshCycles_ )
      mpi::allReduceInplace( additionalRefreshCycleRequired, mpi::LOGICAL_OR );

   std::map< uint_t, std::map< BlockID, uint_t > > targetLevelsFromNeighbors;

   mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 1182 ); // blockforest = 98 108 111 99 107 102 111 114 101 115 116
   bufferSystem.setReceiverInfo( ranksToRecvFrom, false ); // ATTENTION: false = the size of a message from A to B is always the same!

   // SPLIT

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - determining blocks that must be split (= block level will increase)" );

   bool forcedRefinement( false );
   for( uint_t i = 0; i != depth_; ++i )
   {
      WALBERLA_LOG_PROGRESS( "BlockForest refresh: - starting split iteration " << (i+uint_t(1)) );
      WALBERLA_LOG_PROGRESS( "BlockForest refresh:   + sending local block target levels to neighbors" );

      for( auto neighbor = localTargetLevelsForNeighbors.begin(); neighbor != localTargetLevelsForNeighbors.end(); ++neighbor )
      {
         WALBERLA_ASSERT_UNEQUAL( neighbor->first, process_ );
         bufferSystem.sendBuffer( neighbor->first ) << neighbor->second;
      }

      bufferSystem.sendAll();

      WALBERLA_LOG_PROGRESS( "BlockForest refresh:   + receiving block target levels from neighbors" );

      for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
      {
         WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process_ );
         recvIt.buffer() >> targetLevelsFromNeighbors[ uint_c( recvIt.rank() ) ];
      }

      WALBERLA_LOG_PROGRESS( "BlockForest refresh:   + determining local blocks that must be split in order to guarantee 2:1 balance" );

      for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      {
         WALBERLA_ASSERT_NOT_NULLPTR( it->second );
         auto & block = it->second;

         // "block->getTargetLevel() == block->getLevel()": we only need to check blocks that are not yet marked for refinement/splitting
         for( uint_t n = 0; block->getTargetLevel() == block->getLevel() && n != block->getNeighborhoodSize(); ++n )
         {
            const uint_t p = block->getNeighborProcess(n);
            const auto & id = block->getNeighborId(n);

            uint_t neighborTargetLevel( uint_t(0) );
            if( p == process_ )
            {
               WALBERLA_ASSERT( blocks_.find(id) != blocks_.end() );
               WALBERLA_ASSERT_NOT_NULLPTR( blocks_[id] );
               neighborTargetLevel = blocks_[id]->getTargetLevel();
            }
            else
            {
               WALBERLA_ASSERT( targetLevelsFromNeighbors.find(p) != targetLevelsFromNeighbors.end() );
               WALBERLA_ASSERT( targetLevelsFromNeighbors[p].find(id) != targetLevelsFromNeighbors[p].end() );
               neighborTargetLevel = targetLevelsFromNeighbors[p][id];
            }

            if( neighborTargetLevel > (block->getTargetLevel() + uint_t(1)) )
            {
               WALBERLA_ASSERT_EQUAL( neighborTargetLevel, block->getTargetLevel() + uint_t(2) );

               const uint_t targetLevel = block->getTargetLevel() + uint_t(1);
               block->setTargetLevel( targetLevel );
               forcedRefinement = true;

               const auto & bid = block->getId();
               for( uint_t nn = 0; nn != block->getNeighborhoodSize(); ++nn )
               {
                  const uint_t np = block->getNeighborProcess(nn);

                  if( np != process_ )
                  {
                     WALBERLA_ASSERT( localTargetLevelsForNeighbors.find(np) != localTargetLevelsForNeighbors.end() );
                     WALBERLA_ASSERT( localTargetLevelsForNeighbors[np].find(bid) != localTargetLevelsForNeighbors[np].end() );

                     localTargetLevelsForNeighbors[ np ][ bid ] = targetLevel;
                  }
               }
            }
         }
      }
   }

   if( reevaluateMinTargetLevelsAfterForcedRefinement_ )
   {
      mpi::allReduceInplace( forcedRefinement, mpi::LOGICAL_OR );
      rerun = forcedRefinement;
   }

   // prepare data structures required for merge

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - determining blocks that must be checked for potential merge" );

   std::map< BlockID, std::map< BlockID, uint_t > > octetMember; // local octet members that are not yet marked for splitting BUT marked for a potential merge

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( it->second );

      auto & block = it->second;
      const auto & id = block->getId();

      WALBERLA_ASSERT( minTargetLevels.find(id) != minTargetLevels.end() );

      if( block->getLevel() > uint_t(0) && block->getTargetLevel() == block->getLevel() && minTargetLevels[id] < block->getLevel() )
      {
         const auto fatherId = id.getFatherId();
         std::map< BlockID, uint_t > octetMembers;
         for( uint_t n = 0; n != block->getNeighborhoodSize(); ++n )
         {
            const auto & nId = block->getNeighborId(n);
            if( getLevelFromBlockId( nId ) > uint_t(0) && nId.getFatherId() == fatherId )
               octetMembers[ nId ] = block->getNeighborProcess(n);
         }
         if( octetMembers.size() == uint_c(7) )
            octetMember[ id ] = octetMembers;
      }
   }

   mpi::BufferSystem mergeBufferSystem( MPIManager::instance()->comm(), 1183 ); // blockforest = 98 108 111 99 107 102 111 114 101 115 116 + 1
   mergeBufferSystem.setReceiverInfo( ranksToRecvFrom, true ); // true = message sizes between A and B vary

   // MERGE

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - determining blocks that must be merged (= block level will decrease)" );

   for( uint_t i = 0; i != depth_; ++i )
   {
      WALBERLA_LOG_PROGRESS( "BlockForest refresh: - starting merge iteration " << (i+uint_t(1)) );
      WALBERLA_LOG_PROGRESS( "BlockForest refresh:   + sending local block target levels to neighbors" );

      for( auto neighbor = localTargetLevelsForNeighbors.begin(); neighbor != localTargetLevelsForNeighbors.end(); ++neighbor )
      {
         WALBERLA_ASSERT_UNEQUAL( neighbor->first, process_ );
         bufferSystem.sendBuffer( neighbor->first ) << neighbor->second;
      }

      bufferSystem.sendAll();

      WALBERLA_LOG_PROGRESS( "BlockForest refresh:   + receiving block target levels from neighbors" );

      for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
      {
         WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process_ );
         recvIt.buffer() >> targetLevelsFromNeighbors[ uint_c( recvIt.rank() ) ];
      }

      // tell neighbors about my intent to merge

      WALBERLA_LOG_PROGRESS( "BlockForest refresh:   + decide which blocks might get merged" );

      std::map< uint_t, std::set< BlockID > > intentToMerge;

      for( auto it = octetMember.begin(); it != octetMember.end(); it++ )
      {
         WALBERLA_ASSERT( blocks_.find(it->first) != blocks_.end() );
         WALBERLA_ASSERT_NOT_NULLPTR( blocks_[it->first] );
         auto & block = blocks_[it->first];
         bool mergeCapable( true );
         for( uint_t n = 0; mergeCapable && n != block->getNeighborhoodSize(); ++n )
         {
            const uint_t np = block->getNeighborProcess(n);
            const auto & nId = block->getNeighborId(n);

            uint_t neighborTargetLevel( uint_t(0) );
            if( np == process_ )
            {
               WALBERLA_ASSERT( blocks_.find(nId) != blocks_.end() );
               WALBERLA_ASSERT_NOT_NULLPTR( blocks_[nId] );
               neighborTargetLevel = blocks_[nId]->getTargetLevel();
            }
            else
            {
               WALBERLA_ASSERT( targetLevelsFromNeighbors.find(np) != targetLevelsFromNeighbors.end() );
               WALBERLA_ASSERT( targetLevelsFromNeighbors[np].find(nId) != targetLevelsFromNeighbors[np].end() );
               neighborTargetLevel = targetLevelsFromNeighbors[np][nId];
            }

            mergeCapable = ( neighborTargetLevel <= block->getTargetLevel() );
         }
         if( mergeCapable )
            intentToMerge[ process_ ].insert( block->getId() );
      }

      WALBERLA_ASSERT_LESS_EQUAL( intentToMerge.size(), uint_t(1) );
      WALBERLA_ASSERT_EQUAL( ranksToRecvFrom.size(), localTargetLevelsForNeighbors.size() );

      WALBERLA_LOG_PROGRESS( "BlockForest refresh:   + tell neighbors about blocks that can get merged" );

      for( auto rank = ranksToRecvFrom.begin(); rank != ranksToRecvFrom.end(); ++rank )
      {
         WALBERLA_ASSERT_UNEQUAL( *rank, process_ );
         mergeBufferSystem.sendBuffer( *rank ) << intentToMerge[ process_ ];
      }

      mergeBufferSystem.sendAll();

      for( auto recvIt = mergeBufferSystem.begin(); recvIt != mergeBufferSystem.end(); ++recvIt )
      {
         WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process_ );
         recvIt.buffer() >> intentToMerge[ uint_c( recvIt.rank() ) ];
      }

      // perform merge if possible

      WALBERLA_LOG_PROGRESS( "BlockForest refresh:   + decide which blocks must get merged" );

      if( !intentToMerge.empty() )
      {
         WALBERLA_ASSERT( intentToMerge.find( process_ ) != intentToMerge.end() );

         const std::set< BlockID > & blocksToCheck = intentToMerge[ process_ ];
         for( auto it = blocksToCheck.begin(); it != blocksToCheck.end(); ++it )
         {
            const BlockID & id = *it;
            WALBERLA_ASSERT( octetMember.find(id) != octetMember.end() );
            const std::map< BlockID, uint_t > & octetMembers = octetMember[id];
            WALBERLA_ASSERT( octetMembers.size() == uint_t(7) );
            bool mergePossible( true );
            for( auto oit = octetMembers.begin(); mergePossible && oit != octetMembers.end(); ++oit )
            {
               WALBERLA_ASSERT( intentToMerge.find(oit->second) != intentToMerge.end() );
               const std::set< BlockID > & blocksWillingToMerge = intentToMerge[ oit->second ];
               mergePossible = ( blocksWillingToMerge.find(oit->first) != blocksWillingToMerge.end() );
            }

            if( mergePossible )
            {
               WALBERLA_ASSERT( blocks_.find(id) != blocks_.end() );
               WALBERLA_ASSERT_NOT_NULLPTR( blocks_[id] );
               auto & block = blocks_[id];

               WALBERLA_ASSERT( block->getTargetLevel() > uint_t(0) );
               WALBERLA_ASSERT( block->getTargetLevel() == block->getLevel() );
               WALBERLA_ASSERT( minTargetLevels.find(id) != minTargetLevels.end() );
               WALBERLA_ASSERT( minTargetLevels[id] < block->getTargetLevel() );

               const uint_t targetLevel = block->getTargetLevel() - uint_t(1);
               block->setTargetLevel( targetLevel );

#ifdef NDEBUG
               if( i < ( depth_ - uint_t(1) ) )
               {
#endif
                  for( uint_t n = 0; n != block->getNeighborhoodSize(); ++n )
                  {
                     const uint_t np = block->getNeighborProcess(n);

                     if( np != process_ )
                     {
                        WALBERLA_ASSERT( localTargetLevelsForNeighbors.find(np) != localTargetLevelsForNeighbors.end() );
                        WALBERLA_ASSERT( localTargetLevelsForNeighbors[np].find(id) != localTargetLevelsForNeighbors[np].end() );

                        localTargetLevelsForNeighbors[ np ][ id ] = targetLevel;
                     }
                  }
#ifdef NDEBUG
               }
#endif
               // remove blocks selected for merge from 'octetMember'
               octetMember.erase(id);
            }
         }
      }
   }

#ifndef NDEBUG
   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - we are in debug mode, so let's perform some consistency checks ..." );
   for( auto neighbor = localTargetLevelsForNeighbors.begin(); neighbor != localTargetLevelsForNeighbors.end(); ++neighbor )
   {
      WALBERLA_ASSERT_UNEQUAL( neighbor->first, process_ );
      bufferSystem.sendBuffer( neighbor->first ) << neighbor->second;
   }
   bufferSystem.sendAll();
   for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
   {
      WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process_ );
      recvIt.buffer() >> targetLevelsFromNeighbors[ uint_c( recvIt.rank() ) ];
   }
   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( it->second );
      auto & block = it->second;
      WALBERLA_ASSERT( block->getTargetLevel() == block->getLevel() ||
                       block->getTargetLevel() == ( block->getLevel() - uint_t(1) ) ||
                       block->getTargetLevel() == ( block->getLevel() + uint_t(1) ) );
      WALBERLA_ASSERT( minTargetLevels.find( block->getId() ) != minTargetLevels.end() );
      WALBERLA_ASSERT( block->getTargetLevel() >= minTargetLevels[ block->getId() ] );
      for( uint_t n = 0; n != block->getNeighborhoodSize(); ++n )
      {
         const uint_t np = block->getNeighborProcess(n);
         const auto & nId = block->getNeighborId(n);
         const uint_t tl = block->getTargetLevel();
         uint_t ntl( tl );
         if( np == process_ )
         {
            WALBERLA_ASSERT( blocks_.find(nId) != blocks_.end() );
            WALBERLA_ASSERT_NOT_NULLPTR( blocks_[nId] );
            ntl = blocks_[nId]->getTargetLevel();
         }
         else
         {
            WALBERLA_ASSERT( targetLevelsFromNeighbors.find(np) != targetLevelsFromNeighbors.end() );
            WALBERLA_ASSERT( targetLevelsFromNeighbors[np].find(nId) != targetLevelsFromNeighbors[np].end() );
            ntl = targetLevelsFromNeighbors[np][nId];
         }
         WALBERLA_ASSERT( tl == ntl || tl == (ntl + uint_t(1)) || tl == (ntl - uint_t(1)) );
      }
   }
#endif

   if( checkForLateOutInRefresh_ && !alwaysRebalanceInRefresh_ ) // if alwaysRebalanceInRefresh_ == true, we do not need this check which is supposed to determine whether or not rebalancing is necessary
   {
      bool levelChanges( false );
      for( auto it = blocks_.begin(); !levelChanges && it != blocks_.end(); ++it )
      {
         WALBERLA_ASSERT_NOT_NULLPTR( it->second );
         auto & block = it->second;
         if( block->getLevel() != block->getTargetLevel() )
            levelChanges = true;
      }

      mpi::allReduceInplace( levelChanges, mpi::LOGICAL_OR );

      return levelChanges;
   }

   return true; // there was no check for level changes -> we have to assume that something changed
}



void BlockForest::update( PhantomBlockForest & phantomForest )
{
   //////////////
   // CALLBACK //
   //////////////

   if( ! callbackBeforeBlockDataIsPacked_.empty() )
      WALBERLA_LOG_PROGRESS( "BlockForest refresh: - executing call back functions before block data is packed" );

   for( auto f = callbackBeforeBlockDataIsPacked_.begin(); f != callbackBeforeBlockDataIsPacked_.end(); ++f )
      f->second( *this, phantomForest );

   //////////////////////////////////////
   // SCHEDULE RECV's FOR BUFFER SIZES //
   //////////////////////////////////////

   std::map< uint_t, uint_t > numberOfBlocksToRecv; // does not include local transfers

   const auto & phantomBlocks = phantomForest.getBlockMap();
   for( auto phantom = phantomBlocks.begin(); phantom != phantomBlocks.end(); ++phantom )
   {
      auto & pBlock = phantom->second;
      auto & sourceProcesses = pBlock->getSourceProcess();
      if( pBlock->getSourceLevel() != pBlock->getLevel() || sourceProcesses[0] != process_ )
      {
         WALBERLA_ASSERT( blocks_.find( pBlock->getId() ) == blocks_.end() );
         for( auto p = sourceProcesses.begin(); p != sourceProcesses.end(); ++p )
         {
            if( *p != process_ )
            {
               if( numberOfBlocksToRecv.find( *p ) == numberOfBlocksToRecv.end() )
                  numberOfBlocksToRecv[*p] = uint_t(1);
               else
                  ++(numberOfBlocksToRecv[*p]);
            }
         }
      }
#ifndef NDEBUG
      else { WALBERLA_ASSERT( blocks_.find( pBlock->getId() ) != blocks_.end() ); }
#endif
   }

   std::map< uint_t, std::vector< uint_t > > recvBufferSizes; // does not include local transfers

   for( auto it = numberOfBlocksToRecv.begin(); it != numberOfBlocksToRecv.end(); ++it )
   {
      WALBERLA_ASSERT_GREATER( it->second, uint_t(0) );
      recvBufferSizes[ it->first ].resize( it->second );
   }

   std::vector< MPI_Request > recvBufferSizesRequests( numberOfBlocksToRecv.size() ); // do not resize this vector!

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - schedule receives for buffer sizes" );

   uint_t i( uint_t(0) );
   for( auto it = recvBufferSizes.begin(); it != recvBufferSizes.end(); ++it )
   {
      MPI_Irecv( static_cast< void * >( &(it->second[0]) ), int_c( it->second.size() ), MPITrait< uint_t >::type(),
                 int_c( it->first ), 0, MPIManager::instance()->comm(), &(recvBufferSizesRequests[i]) );
      ++i;
   }

   ///////////////////////////////
   // FETCH TARGET BLOCK STATES //
   ///////////////////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - fetch target block states" );

   std::set< mpi::MPIRank > ranksToRecvFrom;

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      auto & block = it->second;
      auto & targetProcesses = block->getTargetProcess();
      if( block->getTargetLevel() != block->getLevel() || targetProcesses[0] != process_ )
      {
         WALBERLA_ASSERT( targetProcesses.size() == uint_t(1) || targetProcesses.size() == uint_t(8) );
         for( auto p = targetProcesses.begin(); p != targetProcesses.end(); ++p )
            if( *p != process_ )
               ranksToRecvFrom.insert( numeric_cast< mpi::MPIRank >(*p) );
      }
   }

   mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 1184 ); // blockforest = 98 108 111 99 107 102 111 114 101 115 116 + 2
   bufferSystem.setReceiverInfo( ranksToRecvFrom, true ); // ATTENTION: true = the size of a message from A to B varies

   std::map< uint_t, std::map< BlockID, Set<SUID> > > blockStates;

   for( auto phantom = phantomBlocks.begin(); phantom != phantomBlocks.end(); ++phantom )
   {
      auto & pBlock = phantom->second;
      auto & sourceProcesses = pBlock->getSourceProcess();
      if( pBlock->getSourceLevel() != pBlock->getLevel() || sourceProcesses[0] != process_ )
      {
         for( auto p = sourceProcesses.begin(); p != sourceProcesses.end(); ++p )
            blockStates[ *p ][ pBlock->getId() ] = pBlock->getState();
      }
   }

   for( auto it = blockStates.begin(); it != blockStates.end(); ++it )
   {
      if( it->first != process_ )
         bufferSystem.sendBuffer( numeric_cast< mpi::MPIRank >( it->first ) ) << it->second;
   }

   bufferSystem.sendAll();

   std::map< BlockID, Set<SUID> > targetBlockStates;

   for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
   {
      WALBERLA_ASSERT_UNEQUAL( uint_c( recvIt.rank() ), process_ );
      std::map< BlockID, Set<SUID> > states;
      recvIt.buffer() >> states;
      for( auto state = states.begin(); state != states.end(); ++state )
      {
         WALBERLA_ASSERT( targetBlockStates.find( state->first ) == targetBlockStates.end() );
         targetBlockStates[ state->first ] = state->second;
      }
   }
   if( blockStates.find( process_ ) != blockStates.end() )
   {
      std::map< BlockID, Set<SUID> > & states = blockStates[ process_ ];
      for( auto state = states.begin(); state != states.end(); ++state )
      {
         WALBERLA_ASSERT( targetBlockStates.find( state->first ) == targetBlockStates.end() );
         targetBlockStates[ state->first ] = state->second;
      }
   }

   ///////////////
   // PACK DATA //
   ///////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - determine blocks whose data will be packed" );

   std::vector< std::pair< Block *, std::vector< mpi::SendBuffer > > > blocksToPack; // includes data that is NOT transfered via MPI but copied locally

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
   {
      auto & block = it->second;
      auto & targetProcesses = block->getTargetProcess();
      if( block->getTargetLevel() != block->getLevel() || targetProcesses[0] != process_ )
      {
         WALBERLA_ASSERT( targetProcesses.size() == uint_t(1) || targetProcesses.size() == uint_t(8) );
         blocksToPack.emplace_back( block.get(), std::vector< mpi::SendBuffer >( targetProcesses.size() ) );
      }
   }

   // do not resize 'blocksToPack' or all mpi::SendBuffer vectors within 'blocksToPack' after this point!

   std::map< uint_t, std::vector< mpi::SendBuffer * > > processesToSendTo; // does not include data that is copied locally
   std::vector< mpi::SendBuffer * > localBlocks; // data that is copied locally

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - packing block data into buffers" );

   //#ifdef _OPENMP
   //#pragma omp parallel for schedule(dynamic)
   //#endif
   for( int j = 0; j < int_c( blocksToPack.size() ); ++j )
   {
      Block * block = blocksToPack[uint_c(j)].first;
      auto & buffers = blocksToPack[uint_c(j)].second;

      auto & id = block->getId();

      auto & targetProcesses = block->getTargetProcess();

      std::vector< Set<SUID> > targetState;

      // header = sender block ID + receiver block ID

      if( block->targetBlockHasTheSameSize() )
      {
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(1) );
         WALBERLA_ASSERT( buffers[0].isEmpty() );
         WALBERLA_ASSERT_EQUAL( targetProcesses.size(), uint_t(1) );

         buffers[0] << id << id;

         if( targetProcesses[0] != process_ )
            processesToSendTo[ targetProcesses[0] ].push_back( &(buffers[0]) );
         else
            localBlocks.push_back( &(buffers[0]) );

         WALBERLA_ASSERT( targetBlockStates.find(id) != targetBlockStates.end() );
         targetState.push_back( targetBlockStates[id] );
      }
      else if( block->targetBlockIsLarger() )
      {
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(1) );
         WALBERLA_ASSERT( buffers[0].isEmpty() );
         WALBERLA_ASSERT_EQUAL( targetProcesses.size(), uint_t(1) );

         auto fid = id.getFatherId();
         buffers[0] << id << fid;

         if( targetProcesses[0] != process_ )
            processesToSendTo[ targetProcesses[0] ].push_back( &(buffers[0]) );
         else
            localBlocks.push_back( &(buffers[0]) );

         WALBERLA_ASSERT( targetBlockStates.find(fid) != targetBlockStates.end() );
         targetState.push_back( targetBlockStates[fid] );
      }
      else
      {
         WALBERLA_ASSERT( block->targetBlockIsSmaller() );
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(8) );
         WALBERLA_ASSERT_EQUAL( targetProcesses.size(), uint_t(8) );
         for( uint_t c = 0; c != 8; ++c )
         {
            WALBERLA_ASSERT( buffers[c].isEmpty() );

            BlockID cid( id, c );
            buffers[c] << id << cid;

            if( targetProcesses[c] != process_ )
               processesToSendTo[ targetProcesses[c] ].push_back( &(buffers[c]) );
            else
               localBlocks.push_back( &(buffers[c]) );

            WALBERLA_ASSERT( targetBlockStates.find(cid) != targetBlockStates.end() );
            targetState.push_back( targetBlockStates[cid] );
         }
      }

      // sender block "state"

      for( auto buffer = buffers.begin(); buffer != buffers.end(); ++buffer )
         *buffer << block->getState();

      // block data

      if( block->targetBlockHasTheSameSize() )
      {
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(1) );
         WALBERLA_ASSERT_EQUAL( targetProcesses.size(), uint_t(1) );
         WALBERLA_ASSERT_EQUAL( targetState.size(), uint_t(1) );

         for( auto dataItem = blockDataItem_.begin(); dataItem != blockDataItem_.end(); ++dataItem )
         {
            auto blockDataHandlingWrapper = dataItem->getDataHandling( block, targetState[0] );
            if( blockDataHandlingWrapper )
               blockDataHandlingWrapper->serialize( block, dataItem->getId(), buffers[0] );
         }
      }
      else if( block->targetBlockIsLarger() )
      {
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(1) );
         WALBERLA_ASSERT_EQUAL( targetProcesses.size(), uint_t(1) );
         WALBERLA_ASSERT_EQUAL( targetState.size(), uint_t(1) );

         for( auto dataItem = blockDataItem_.begin(); dataItem != blockDataItem_.end(); ++dataItem )
         {
            auto blockDataHandlingWrapper = dataItem->getDataHandling( block, targetState[0] );
            if( blockDataHandlingWrapper )
            {
               auto downcastBlockDataHandlingWrapper = dynamic_pointer_cast< blockforest::internal::BlockDataHandlingWrapper >( blockDataHandlingWrapper );
               if( !downcastBlockDataHandlingWrapper )
               {
                  WALBERLA_ABORT( "Dynamic data structure refresh failed!\n"
                                  "For the dynamic refresh to work, all registered block data items must implement the blockforest::BlockDataHandling interface\n"
                                  "_AND_ this block data handling object must be registered at an instance of class BlockForest or StructuredBlockForest!\n"
                                  "Registering the block data handling object through a base class [BlockStorage/StructuredBlockStorage] (shared) pointer will not work!\n"
                                  "For block data item '" << dataItem->getIdentifier() << "' a fitting block data handling object is missing or was registered\n"
                                  "through a base class pointer and not directly at an instance of class BlockForest or StructuredBlockForest.\n" );
               }
               downcastBlockDataHandlingWrapper->serializeFineToCoarse( block, dataItem->getId(), buffers[0] );
            }
         }
      }
      else
      {
         WALBERLA_ASSERT( block->targetBlockIsSmaller() );
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(8) );
         WALBERLA_ASSERT_EQUAL( targetProcesses.size(), uint_t(8) );
         WALBERLA_ASSERT_EQUAL( targetState.size(), uint_t(8) );

         for( auto dataItem = blockDataItem_.begin(); dataItem != blockDataItem_.end(); ++dataItem )
         {
            for( uint_t c = uint_t(0); c != uint_t(8); ++c )
            {
               auto blockDataHandlingWrapper = dataItem->getDataHandling( block, targetState[c] );
               if( blockDataHandlingWrapper )
               {
                  auto downcastBlockDataHandlingWrapper = dynamic_pointer_cast< blockforest::internal::BlockDataHandlingWrapper >( blockDataHandlingWrapper );
                  if( !downcastBlockDataHandlingWrapper )
                  {
                     WALBERLA_ABORT( "Dynamic data structure refresh failed!\n"
                                     "For the dynamic refresh to work, all registered block data items must implement the blockforest::BlockDataHandling interface\n"
                                     "_AND_ this block data handling object must be registered at an instance of class BlockForest or StructuredBlockForest!\n"
                                     "Registering the block data handling object through a base class [BlockStorage/StructuredBlockStorage] (shared) pointer will not work!\n"
                                     "For block data item '" << dataItem->getIdentifier() << "' a fitting block data handling object is missing or was registered\n"
                                     "through a base class pointer and not directly at an instance of class BlockForest or StructuredBlockForest.\n" );
                  }
                  downcastBlockDataHandlingWrapper->serializeCoarseToFine( block, dataItem->getId(), buffers[c], c );
               }
            }
         }
      }
   }

   ///////////////////
   // DELETE BLOCKS //
   ///////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - deleting blocks that are split, merged, and/or transfered" );

   for( auto block = blocksToPack.begin(); block != blocksToPack.end(); ++block )
      blocks_.erase( block->first->getId() );

   ///////////////
   // SEND DATA //
   ///////////////

   std::map< uint_t, std::vector< uint_t > > sendBufferSizes; // does not include local transfers

   std::vector< MPI_Request > sendBufferSizesRequests( processesToSendTo.size() ); // do not resize this vector!
   std::map< uint_t, std::vector< MPI_Request > > blockDataSendRequests;

   for( auto it = processesToSendTo.begin(); it != processesToSendTo.end(); ++it )
   {
      WALBERLA_ASSERT( sendBufferSizes.find( it->first ) == sendBufferSizes.end() );

      auto & sizes = sendBufferSizes[ it->first ];
      WALBERLA_ASSERT( sizes.empty() );
      for( auto buffer = it->second.begin(); buffer != it->second.end(); ++buffer )
         sizes.push_back( (*buffer)->size() );

      blockDataSendRequests[ it->first ].resize( it->second.size() ); // do not resize this vector after this point!
   }

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - sending block data" );

   i = uint_t(0);
   for( auto it = processesToSendTo.begin(); it != processesToSendTo.end(); ++it )
   {
      WALBERLA_ASSERT( sendBufferSizes.find( it->first ) != sendBufferSizes.end() );

      auto & sizes = sendBufferSizes[ it->first ];
      WALBERLA_ASSERT_EQUAL( sizes.size(), it->second.size() );
      MPI_Isend( static_cast< void * >( &(sizes[0]) ), int_c( sizes.size() ), MPITrait< uint_t >::type(), int_c( it->first ), 0,
                 MPIManager::instance()->comm(), &(sendBufferSizesRequests[i]) );

      int j(0);
      for( auto buffer = it->second.begin(); buffer != it->second.end(); ++buffer )
      {
         WALBERLA_ASSERT( blockDataSendRequests.find( it->first ) != blockDataSendRequests.end() );
         WALBERLA_ASSERT_GREATER( blockDataSendRequests[ it->first ].size(), j );

         MPI_Isend( static_cast< void * >( (*buffer)->ptr() ), int_c( (*buffer)->size() ), MPITrait< mpi::SendBuffer::ElementType >::type(),
                    int_c( it->first ), j+1, MPIManager::instance()->comm(), &(blockDataSendRequests[it->first][uint_c(j)]) );
         ++j;
      }

      ++i;
   }

   ///////////////////
   // CREATE BLOCKS //
   ///////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - allocating new blocks" );

   for( auto phantom = phantomBlocks.begin(); phantom != phantomBlocks.end(); ++phantom )
   {
      auto & pBlock = phantom->second;
      if( pBlock->getSourceLevel() != pBlock->getLevel() || pBlock->getSourceProcess()[0] != process_ )
      {
         WALBERLA_ASSERT( blocks_.find( pBlock->getId() ) == blocks_.end() );
         blocks_[ pBlock->getId() ] = std::make_shared< Block >( *this, *pBlock );
      }
      else // update neighborhood of existing blocks
      {
         WALBERLA_ASSERT( blocks_.find( pBlock->getId() ) != blocks_.end() );
         blocks_[ pBlock->getId() ]->resetNeighborhood( *pBlock );
      }
   }

   // adapt depth

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - adapting block structure (use current state of phantom forest as reference)" );

   if( allowChangingDepth_ )
      depth_ = phantomForest.getDepth();
#ifndef NDEBUG
   else { WALBERLA_ASSERT_EQUAL( depth_, phantomForest.getDepth() ); }
#endif

   // copy process neighborhood information from phantom forest

   neighborhood_ = phantomForest.getNeighboringProcesses();

   // Some processes might now be empty (= without blocks) _or_ are now not empty anymore -> rebuild communicator that only contains processes with blocks

   rebuildProcessesWithBlocksCommunicator();

   //////////////////////////////////////
   // WAIT FOR RECV's FOR BUFFER SIZES //
   //////////////////////////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - wait for buffer sizes" );

   if( ! recvBufferSizesRequests.empty() )
      MPI_Waitall( int_c( recvBufferSizesRequests.size() ), &(recvBufferSizesRequests[0]), MPI_STATUSES_IGNORE );

   ////////////////////////////////////
   // SCHEDULE RECV's FOR BLOCK DATA //
   ////////////////////////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - schedule block data receive operations" );

   std::map< uint_t, std::vector< mpi::RecvBuffer > > recvBlockData;
   std::map< uint_t, std::vector< MPI_Request > > blockDataRecvRequests;

   for( auto it = recvBufferSizes.begin(); it != recvBufferSizes.end(); ++it )
   {
      WALBERLA_ASSERT_UNEQUAL( it->first, process_ );
      auto & buffers = recvBlockData[ it->first ];
      auto & requests = blockDataRecvRequests[ it->first ];
      auto & sizes = it->second;
      buffers.resize( sizes.size() );
      requests.resize( sizes.size() );
   }
   auto & recvLocalBlocks = recvBlockData[ process_ ];

   for( auto it = recvBufferSizes.begin(); it != recvBufferSizes.end(); ++it )
   {
      WALBERLA_ASSERT_UNEQUAL( it->first, process_ );
      auto & buffers = recvBlockData[ it->first ];
      auto & requests = blockDataRecvRequests[ it->first ];
      auto & sizes = it->second;
      WALBERLA_ASSERT_EQUAL( buffers.size(), sizes.size() );
      WALBERLA_ASSERT_EQUAL( requests.size(), sizes.size() );
      for( i = uint_t(0); i != buffers.size(); ++i )
      {
         buffers[i].resize( sizes[i] );
         MPI_Irecv( static_cast< void * >( buffers[i].ptr() ), int_c( sizes[i] ), MPITrait< mpi::RecvBuffer::ElementType >::type(),
                    int_c( it->first ), int_c(i)+1, MPIManager::instance()->comm(), &(requests[i]) );
      }
   }

   ///////////////////////////
   // COPY LOCAL BLOCK DATA //
   ///////////////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - perform local data transfer" );

   for( auto buffer = localBlocks.begin(); buffer != localBlocks.end(); ++buffer )
      recvLocalBlocks.emplace_back( **buffer );

   ////////////////////////////////////
   // WAIT FOR RECV's FOR BLOCK DATA //
   ////////////////////////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - wait for block data to be received" );

   for( auto it = blockDataRecvRequests.begin(); it != blockDataRecvRequests.end(); ++it )
   {
      auto & requests = it->second;
      MPI_Waitall( int_c( requests.size() ), &(requests[0]), MPI_STATUSES_IGNORE );
   }

   //////////////////////////////////
   // WAIT FOR ALL SENDS TO FINISH //
   //////////////////////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - wait for block data sends to complete" );

   if( ! sendBufferSizesRequests.empty() )
      MPI_Waitall( int_c( sendBufferSizesRequests.size() ), &(sendBufferSizesRequests[0]), MPI_STATUSES_IGNORE );

   for( auto it = blockDataSendRequests.begin(); it != blockDataSendRequests.end(); ++it )
   {
      auto & requests = it->second;
      MPI_Waitall( int_c( requests.size() ), &(requests[0]), MPI_STATUSES_IGNORE );
   }

   ////////////////////////////////////////
   // CLEAR SEND BUFFERS (= FREE MEMORY) //
   ////////////////////////////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - clear send buffers (= free memory)" );

   for( auto it = blocksToPack.begin(); it != blocksToPack.end(); ++it )
   {
      auto & buffers = it->second;
      for( auto buffer = buffers.begin(); buffer != buffers.end(); ++buffer )
         buffer->reset();
   }

   ////////////////////////////////
   // PREPARE DATA FOR UNPACKING //
   ////////////////////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - prepare received data for unpacking" );

   std::map< Block *, std::vector< std::pair< Set<SUID>, mpi::RecvBuffer * > > > blocksToUnpack; // includes data that is NOT transfered via MPI but copied locally

   for( auto it = recvBlockData.begin(); it != recvBlockData.end(); ++it )
   {
      auto & buffers = it->second;
      for( auto buffer = buffers.begin(); buffer != buffers.end(); ++buffer )
      {
         // header = sender block ID + receiver block ID & sender block "state"

         BlockID sId;
         BlockID rId;
         Set<SUID> state;
         (*buffer) >> sId >> rId >> state;

         WALBERLA_ASSERT( blocks_.find( rId ) != blocks_.end() );
         WALBERLA_ASSERT( phantomBlocks.find( rId ) != phantomBlocks.end() );
         Block * block = blocks_[ rId ].get();
         const auto & phantom = phantomBlocks.find(rId)->second;

         if( phantom->sourceBlockHasTheSameSize() || phantom->sourceBlockIsLarger() )
         {
            WALBERLA_ASSERT( blocksToUnpack.find( block ) == blocksToUnpack.end() );
            blocksToUnpack[ block ].push_back( std::make_pair( state, &(*buffer) ) );
         }
         else
         {
            auto & bufferPtrs = blocksToUnpack[ block ];
            if( bufferPtrs.empty() )
               bufferPtrs.resize( uint_t(8), std::make_pair( Set<SUID>::emptySet(), static_cast< mpi::RecvBuffer * >(nullptr) ) );
            WALBERLA_ASSERT_EQUAL( sId.getUsedBits(), rId.getUsedBits() + uint_t(3) );
            bufferPtrs[ sId.getBranchId() ] = std::make_pair( state, &(*buffer) );
         }
      }
   }

#ifndef NDEBUG
   for( auto it = phantomBlocks.begin(); it != phantomBlocks.end(); ++it )
   {
      WALBERLA_ASSERT( blocks_.find( it->first ) != blocks_.end() );
      auto & block = blocks_[ it->first ];
      if( !(it->second->sourceBlockHasTheSameSize()) || it->second->getSourceProcess()[0] != process_ )
      {
         WALBERLA_ASSERT( blocksToUnpack.find( block.get() ) != blocksToUnpack.end() );
         auto & buffers = blocksToUnpack[ block.get() ];
         for( auto buffer = buffers.begin(); buffer != buffers.end(); ++buffer )
            WALBERLA_ASSERT_NOT_NULLPTR( buffer->second );
      }
   }
#endif

   //////////////////////////////
   // GLOBAL BLOCK INFORMATION //
   //////////////////////////////

   // adapt global block information stored on every process (this is only done if the block forest is set-up to store global information!)

   if( containsGlobalBlockInformation() )
   {
      WALBERLA_LOG_PROGRESS( "BlockForest refresh: - update global block information that is stored locally on every process" );
      constructBlockInformation();
   }

   ++modificationStamp_;

   //////////////
   // CALLBACK //
   //////////////

   if( ! callbackBeforeBlockDataIsUnpacked_.empty() )
      WALBERLA_LOG_PROGRESS( "BlockForest refresh: - executing call back functions before block data is unpacked" );

   for( auto f = callbackBeforeBlockDataIsUnpacked_.begin(); f != callbackBeforeBlockDataIsUnpacked_.end(); ++f )
      f->second( *this, phantomForest );

   /////////////////
   // UNPACK DATA //
   /////////////////

   WALBERLA_LOG_PROGRESS( "BlockForest refresh: - unpacking block data from buffers" );

   std::vector< std::pair< Block *, std::vector< std::pair< Set<SUID>, mpi::RecvBuffer * > > > > dataToUnpack;

   for( auto it = blocksToUnpack.begin(); it != blocksToUnpack.end(); ++it )
      dataToUnpack.emplace_back( it->first, it->second );

   //#ifdef _OPENMP
   //#pragma omp parallel for schedule(dynamic)
   //#endif
   for( int j = 0; j < int_c( dataToUnpack.size() ); ++j )
   {
      Block * block = dataToUnpack[uint_c(j)].first;
      auto & buffers = dataToUnpack[uint_c(j)].second;

      WALBERLA_ASSERT( phantomBlocks.find( block->getId() ) != phantomBlocks.end() );
      const auto & phantom = phantomBlocks.find( block->getId() )->second;

      // block data

      if( phantom->sourceBlockHasTheSameSize() )
      {
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(1) );

         for( auto dataItem = blockDataItem_.begin(); dataItem != blockDataItem_.end(); ++dataItem )
         {
            // allocate
            {
               auto blockDataHandlingWrapper = dataItem->getDataHandling( block );
               if( blockDataHandlingWrapper )
                  addBlockData( block, dataItem->getId(), blockDataHandlingWrapper->deserialize( block ) );
               else
                  addBlockData( block, dataItem->getId(), nullptr );
            }
            // fill with sent data
            {
               auto blockDataHandlingWrapper = dataItem->getDataHandling( block, buffers[0].first );
               if( blockDataHandlingWrapper )
               {
                  WALBERLA_ASSERT_NOT_NULLPTR( buffers[0].second );
                  blockDataHandlingWrapper->deserialize( block, dataItem->getId(), *(buffers[0].second) );
               }
            }
         }
      }
      else if( phantom->sourceBlockIsLarger() )
      {
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(1) );

         for( auto dataItem = blockDataItem_.begin(); dataItem != blockDataItem_.end(); ++dataItem )
         {
            // allocate
            {
               auto blockDataHandlingWrapper = dataItem->getDataHandling( block );
               if( blockDataHandlingWrapper )
               {
                  auto downcastBlockDataHandlingWrapper = dynamic_pointer_cast< blockforest::internal::BlockDataHandlingWrapper >( blockDataHandlingWrapper );
                  if( !downcastBlockDataHandlingWrapper )
                  {
                     WALBERLA_ABORT( "Dynamic data structure refresh failed!\n"
                                     "For the dynamic refresh to work, all registered block data items must implement the blockforest::BlockDataHandling interface\n"
                                     "_AND_ this block data handling object must be registered at an instance of class BlockForest or StructuredBlockForest!\n"
                                     "Registering the block data handling object through a base class [BlockStorage/StructuredBlockStorage] (shared) pointer will not work!\n"
                                     "For block data item '" << dataItem->getIdentifier() << "' a fitting block data handling object is missing or was registered\n"
                                     "through a base class pointer and not directly at an instance of class BlockForest or StructuredBlockForest.\n" );
                  }
                  addBlockData( block, dataItem->getId(), downcastBlockDataHandlingWrapper->deserializeCoarseToFine( block ) );
               }
               else
                  addBlockData( block, dataItem->getId(), nullptr );
            }
            // fill with sent data
            {
               auto blockDataHandlingWrapper = dataItem->getDataHandling( block, buffers[0].first );
               if( blockDataHandlingWrapper )
               {
                  auto downcastBlockDataHandlingWrapper = dynamic_pointer_cast< blockforest::internal::BlockDataHandlingWrapper >( blockDataHandlingWrapper );
                  if( !downcastBlockDataHandlingWrapper )
                  {
                     WALBERLA_ABORT( "Dynamic data structure refresh failed!\n"
                                     "For the dynamic refresh to work, all registered block data items must implement the blockforest::BlockDataHandling interface\n"
                                     "_AND_ this block data handling object must be registered at an instance of class BlockForest or StructuredBlockForest!\n"
                                     "Registering the block data handling object through a base class [BlockStorage/StructuredBlockStorage] (shared) pointer will not work!\n"
                                     "For block data item '" << dataItem->getIdentifier() << "' a fitting block data handling object is missing or was registered\n"
                                     "through a base class pointer and not directly at an instance of class BlockForest or StructuredBlockForest.\n" );
                  }
                  WALBERLA_ASSERT_NOT_NULLPTR( buffers[0].second );
                  downcastBlockDataHandlingWrapper->deserializeCoarseToFine( block, dataItem->getId(), *(buffers[0].second) );
               }
            }
         }
      }
      else
      {
         WALBERLA_ASSERT( phantom->sourceBlockIsSmaller() );
         WALBERLA_ASSERT_EQUAL( buffers.size(), uint_t(8) );

         for( auto dataItem = blockDataItem_.begin(); dataItem != blockDataItem_.end(); ++dataItem )
         {
            // allocate
            {
               auto blockDataHandlingWrapper = dataItem->getDataHandling( block );
               if( blockDataHandlingWrapper )
               {
                  auto downcastBlockDataHandlingWrapper = dynamic_pointer_cast< blockforest::internal::BlockDataHandlingWrapper >( blockDataHandlingWrapper );
                  if( !downcastBlockDataHandlingWrapper )
                  {
                     WALBERLA_ABORT( "Dynamic data structure refresh failed!\n"
                                     "For the dynamic refresh to work, all registered block data items must implement the blockforest::BlockDataHandling interface\n"
                                     "_AND_ this block data handling object must be registered at an instance of class BlockForest or StructuredBlockForest!\n"
                                     "Registering the block data handling object through a base class [BlockStorage/StructuredBlockStorage] (shared) pointer will not work!\n"
                                     "For block data item '" << dataItem->getIdentifier() << "' a fitting block data handling object is missing or was registered\n"
                                     "through a base class pointer and not directly at an instance of class BlockForest or StructuredBlockForest.\n" );
                  }
                  addBlockData( block, dataItem->getId(), downcastBlockDataHandlingWrapper->deserializeFineToCoarse( block ) );
               }
               else
                  addBlockData( block, dataItem->getId(), nullptr );
            }
            // fill with sent data
            for( uint_t c = uint_t(0); c != uint_t(8); ++c )
            {
               auto blockDataHandlingWrapper = dataItem->getDataHandling( block, buffers[c].first );
               if( blockDataHandlingWrapper )
               {
                  auto downcastBlockDataHandlingWrapper = dynamic_pointer_cast< blockforest::internal::BlockDataHandlingWrapper >( blockDataHandlingWrapper );
                  if( !downcastBlockDataHandlingWrapper )
                  {
                     WALBERLA_ABORT( "Dynamic data structure refresh failed!\n"
                                     "For the dynamic refresh to work, all registered block data items must implement the blockforest::BlockDataHandling interface\n"
                                     "_AND_ this block data handling object must be registered at an instance of class BlockForest or StructuredBlockForest!\n"
                                     "Registering the block data handling object through a base class [BlockStorage/StructuredBlockStorage] (shared) pointer will not work!\n"
                                     "For block data item '" << dataItem->getIdentifier() << "' a fitting block data handling object is missing or was registered\n"
                                     "through a base class pointer and not directly at an instance of class BlockForest or StructuredBlockForest.\n" );
                  }
                  WALBERLA_ASSERT_NOT_NULLPTR( buffers[c].second );
                  downcastBlockDataHandlingWrapper->deserializeFineToCoarse( block, dataItem->getId(), *(buffers[c].second), c );
               }
            }
         }
      }

      for( auto it = buffers.begin(); it != buffers.end(); ++it )
         it->second->reset();
   }

   //////////////
   // CALLBACK //
   //////////////

   if( ! callbackAfterBlockDataIsUnpacked_.empty() )
      WALBERLA_LOG_PROGRESS( "BlockForest refresh: - executing call back functions after block data was unpacked" );

   for( auto f = callbackAfterBlockDataIsUnpacked_.begin(); f != callbackAfterBlockDataIsUnpacked_.end(); ++f )
      f->second( *this, phantomForest );
}


/// For a description of the file format see BlockForestFile.h
/// \attention 'suidMap' and 'suidBytes' must be identical for every process!
/// \see BlockForestFile.h
void BlockForest::saveToFile( const std::string & filename, FileIOMode fileIOMode,
                              const std::map< SUID, std::vector< bool > > & suidMap, const uint_t suidBytes ) const
{
   // process data

   const uint_t blockIdBytes = getBlockIdBytes();

   uint_t dataSize = uint_t(2) + blocks_.size() * ( blockIdBytes + suidBytes ) + uint_t(2) + neighborhood_.size() * processIdBytes_;
   if( MPIManager::instance()->rank() == 0 )
   {
      dataSize += internal::FILE_HEADER_SIZE; // header
      ++dataSize; // number of SUIDs
      for( auto suid = suidMap.begin(); suid != suidMap.end(); ++suid )
         dataSize += uint_t(1) + uint_c( suid->first.getIdentifier().length() );
   }

   std::vector< uint8_t > processDataBuffer( dataSize );
   uint_t offset = 0;

   if( MPIManager::instance()->rank() == 0 )
   {
      // header
      storeFileHeader( processDataBuffer, offset );

      // SUIDs

      std::vector< SUID > suids( suidMap.size() );

      for( auto suid = suidMap.begin(); suid != suidMap.end(); ++suid )
      {
         for( uint_t i = uint_t(0); i != suids.size(); ++i )
         {
            if( suid->second[i] )
            {
               suids[i] = suid->first;
               break;
            }
         }
      }

      uintToByteArray( suids.size(), processDataBuffer, offset, 1 );
      ++offset;

      // for every SUID ...

      for( auto it = suids.begin(); it != suids.end(); ++it )
      {
         // length of its identifier string

         const uint_t length = it->getIdentifier().length();
         WALBERLA_CHECK_LESS( length, 256, "SUID identifiers are allowed to consist of 255 characters at most when saving the block structure to file!" );

         uintToByteArray( length, processDataBuffer, offset, 1 );
         ++offset;

         // the identifier string

         const char * str = it->getIdentifier().c_str();
         for( uint_t j = 0; j != length; ++j )
            processDataBuffer[offset+j] = *reinterpret_cast< const uint8_t* >( str + j );
         offset += length;
      }
   }

   // number of blocks (can be '0' -> buffer process!)
   uintToByteArray( uint_c( blocks_.size() ), processDataBuffer, offset, uint_t(2) );
   offset += uint_t(2);

   for( auto block = blocks_.begin(); block != blocks_.end(); ++block )
   {
      // block ID
      block->second->getId().toByteArray( processDataBuffer, offset, blockIdBytes );
      offset += blockIdBytes;

      // block state (SUID set)
      if( suidBytes > 0 )
      {
         std::vector< bool > suidBoolVec( 8 * suidBytes );

         const Set<SUID> & state = block->second->getState();
         for( auto suid = state.begin(); suid != state.end(); ++suid )
         {
            WALBERLA_CHECK( suidMap.find( *suid ) != suidMap.end(), "Block state SUID missing from SUID list saved to file."
                                                                    "\n- SUID = " << *suid << "\n- block ID = " << block->first <<
                                                                    "\n- block AABB = " << block->second->getAABB() );
            //Elementwise OR of all elements
            for (uint_t i = 0; i < suidBoolVec.size(); ++i) {
               suidBoolVec[i] = suidBoolVec[i] || suidMap.find( *suid )->second[i];
            }
         }

         boolVectorToByteArray( suidBoolVec, processDataBuffer, offset );
         offset += suidBytes;
      }
   }

   // process neighborhood
   uintToByteArray( uint_c( neighborhood_.size() ), processDataBuffer, offset, uint_t(2) );
   offset += uint_t(2);

   for( auto neighbor = neighborhood_.begin(); neighbor != neighborhood_.end(); ++neighbor )
   {
      uintToByteArray( *neighbor, processDataBuffer, offset, processIdBytes_ );
      offset += processIdBytes_;
   }

   // store data to file

   WALBERLA_NON_MPI_SECTION()
   {
      fileIOMode = MASTER_SLAVE;
   }

   // use serial I/O for versions of OpenMPI that produce segmentation faults when using MPI-IO with a 3D Cartesian MPI
   // communicator (see waLBerla issue #73)
   if (!MPIManager::instance()->isCommMPIIOValid())
   {
      fileIOMode = MASTER_SLAVE;
   }

   if( fileIOMode == MPI_PARALLEL )
   {
      MPI_File mpiFile = MPI_FILE_NULL;
      int result = MPI_SUCCESS;
      result = MPI_File_open( MPIManager::instance()->comm(), const_cast<char*>( filename.c_str() ), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpiFile );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while opening file \"" << filename << "\" for writing. MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      const MPI_Offset filesize = numeric_cast<MPI_Offset>( mpi::allReduce( dataSize, mpi::SUM, MPIManager::instance()->comm() ) );
      MPI_File_set_size( mpiFile, filesize );

      uint_t exscanResult;
      MPI_Exscan( &dataSize, &exscanResult, 1, MPITrait<uint_t>::type(), MPI_SUM, MPIManager::instance()->comm() );
      if( MPIManager::instance()->rank() == 0 )
         exscanResult = uint_t( 0 );

      MPI_Datatype arraytype;
      MPI_Type_contiguous( int_c( dataSize ), MPITrait< uint8_t >::type(), &arraytype );
      MPI_Type_commit( &arraytype );

      result = MPI_File_set_view( mpiFile, numeric_cast<MPI_Offset>( exscanResult ), MPITrait< uint8_t >::type(), arraytype, const_cast<char*>( "native" ), MPI_INFO_NULL );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      result = MPI_File_write_all( mpiFile, reinterpret_cast<char*>( &processDataBuffer[0] ), int_c( processDataBuffer.size() ), MPITrait< uint8_t >::type(), MPI_STATUS_IGNORE );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while writing to file \"" << filename << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      result = MPI_File_close( &mpiFile );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while closing file \"" << filename << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      MPI_Type_free( &arraytype );
   }
   else if( fileIOMode == MASTER_SLAVE )
   {
      std::set< mpi::MPIRank > ranksToRecvFrom;

      if( mpi::MPIManager::instance()->rank() == 0 ) // do _NOT_ use WALBERLA_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
      {
         WALBERLA_ASSERT_EQUAL( mpi::MPIManager::instance()->rank(), 0 );
         for( mpi::MPIRank rank = 1; rank < mpi::MPIManager::instance()->numProcesses(); ++rank )
            ranksToRecvFrom.insert( rank );
      }

      mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 1185 ); // blockforest = 98 108 111 99 107 102 111 114 101 115 116 + 3
      bufferSystem.setReceiverInfo( ranksToRecvFrom, true ); // ATTENTION: true = the size of a message from A to B varies

      if( mpi::MPIManager::instance()->rank() != 0 ) // do _NOT_ use WALBERLA_NON_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
      {
         WALBERLA_ASSERT_UNEQUAL( process_, uint_c(0) );
         bufferSystem.sendBuffer( 0 ) << processDataBuffer;
      }

      bufferSystem.sendAll();

      if( mpi::MPIManager::instance()->rank() == 0 ) // do _NOT_ use WALBERLA_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
      {
         std::ofstream file( filename.c_str(), std::ofstream::binary );

         // data from root process
         file.write( reinterpret_cast< const char* >( &(processDataBuffer[0]) ), numeric_cast< std::streamsize >( processDataBuffer.size() ) );

         std::vector< std::vector< uint8_t > > buffers( uint_c( mpi::MPIManager::instance()->numProcesses() ) );

         // data from all other processes
         for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
         {
            recvIt.buffer() >> buffers[ uint_c( recvIt.rank() ) ];
         }

         auto it = buffers.begin();
         ++it;
         while( it != buffers.end() )
         {
            file.write( reinterpret_cast< const char* >( &((*it)[0]) ), numeric_cast< std::streamsize >( it->size() ) );
            it->clear();
            ++it;
         }

         file.close();
      }
      else
      {
         // begin()/end() must also be called on each slave process in order to
         // properly finalize the communication
         WALBERLA_CHECK( bufferSystem.begin() == bufferSystem.end() );
      }
   }
   else
   {
      for( int r = 0; r < MPIManager::instance()->numProcesses(); ++r )
      {
         if( r == MPIManager::instance()->rank() )
         {
            std::ofstream file;

            if( r == 0 )
               file.open( filename.c_str(), std::ofstream::out | std::ofstream::binary );
            else
               file.open( filename.c_str(), std::ofstream::out | std::ofstream::app | std::ofstream::binary );

            file.write( reinterpret_cast< const char* >( &(processDataBuffer[0]) ), numeric_cast< std::streamsize >( processDataBuffer.size() ) );
            file.close();
         }
         MPI_Barrier( MPIManager::instance()->comm() );
      }
   }
}



void BlockForest::storeFileHeader( std::vector< uint8_t > & data, uint_t & offset ) const
{
   // domain AABB

   for( uint_t i = 0; i != 3; ++i )
      offset += realToByteArray( domain_.min(i), data, offset );

   for( uint_t i = 0; i != 3; ++i )
      offset += realToByteArray( domain_.max(i), data, offset );

   // number of coarse/root blocks in each direction

   for( uint_t i = 0; i != 3; ++i ) {
      uintToByteArray( size_[i], data, offset, 4 );
      offset += 4;
   }

   // domain periodicity

   for( uint_t i = 0; i != 3; ++i ) {
      uintToByteArray( periodic_[i] ? uint_c(1) : uint_c(0), data, offset, 1 );
      ++offset;
   }

   // block forest depth (= number of levels - 1)

   uintToByteArray( depth_, data, offset, 1 );
   ++offset;

   // treeIdDigits (= number of bits used for storing the tree ID [tree ID marker + tree index])

   uintToByteArray( treeIdDigits_, data, offset, 1 );
   ++offset;

   // processIdBytes (= number of bytes required for storing process IDs)

   uintToByteArray( getProcessIdBytes(), data, offset, 1 );
   ++offset;

   // insertBuffersIntoProcessNetwork?

   uintToByteArray( insertBuffersIntoProcessNetwork_ ? uint_c(1) : uint_c(0), data, offset, 1 );
   ++offset;

   // number of processes

   uintToByteArray( uint_c( mpi::MPIManager::instance()->numProcesses() ), data, offset, 4 );
   offset += 4;
}



#ifndef NDEBUG

void BlockForest::checkBlockInformationConsistency( const SetupBlockForest& forest ) const {

   WALBERLA_ASSERT( containsGlobalBlockInformation() );

   std::vector< const SetupBlock* > blocks;
   forest.getBlocks( blocks );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      const SetupBlock* const block = blocks[i];

      uint_t process = 0;
      WALBERLA_ASSERT( getBlockInformation().getProcess( process, block->getId() ) );
      WALBERLA_ASSERT_EQUAL( process, block->getProcess() );
      Set<SUID> state;
      WALBERLA_ASSERT( getBlockInformation().getState( state, block->getId() ) );
      WALBERLA_ASSERT_EQUAL( state, block->getState() );
      WALBERLA_ASSERT_EQUAL( getBlockInformation().existsRemotely( block->getId() ), ( block->getProcess() != process_ ) );
      WALBERLA_ASSERT_EQUAL( blockExistsLocally( block->getId() ), ( block->getProcess() == process_ ) );

      const real_t x = ( block->getAABB().xMax() + block->getAABB().xMin() ) / real_c(2);
      const real_t y = ( block->getAABB().yMax() + block->getAABB().yMin() ) / real_c(2);
      const real_t z = ( block->getAABB().zMax() + block->getAABB().zMin() ) / real_c(2);

      WALBERLA_ASSERT( getBlockInformation().getProcess(process,x,y,z) );
      WALBERLA_ASSERT_EQUAL( process, block->getProcess() );
      WALBERLA_ASSERT( getBlockInformation().getState(state, x,y,z) );
      WALBERLA_ASSERT_EQUAL( state, block->getState() );
      WALBERLA_ASSERT_EQUAL( getBlockInformation().existsRemotely(x,y,z), ( block->getProcess() != process_ ) );
      WALBERLA_ASSERT_EQUAL( blockExistsLocally(x,y,z), ( block->getProcess() == process_ ) );
   }
}

#endif



} // namespace blockforest
} // namespace walberla
