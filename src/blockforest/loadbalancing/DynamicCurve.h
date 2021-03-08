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
//! \file DynamicCurve.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "NoPhantomData.h"
#include "blockforest/BlockForest.h"
#include "blockforest/HilbertCurveConstruction.h"
#include "blockforest/PhantomBlockForest.h"

#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIManager.h"

#include <cmath>
#include <map>
#include <set>
#include <stack>
#include <type_traits>
#include <vector>



namespace walberla {
namespace blockforest {



namespace internal {

template< typename PhantomData_T, typename Value_T, bool weighted = true >
class BlockIDSorter
{
public:
   using weight_t = typename PhantomData_T::weight_t;
   using Blocks_T = std::vector<std::vector<std::pair<BlockID, weight_t>>>;
   BlockIDSorter( const Blocks_T & blocks ) : blocks_( blocks ) {}
   bool operator()( const Value_T & lhs, const Value_T & rhs ) const
   {
      return blocks_[ uint_c(lhs.first) ][ lhs.second ].first < blocks_[ uint_c(rhs.first) ][ rhs.second ].first;
   }
private:
   const Blocks_T & blocks_;
};

template< typename PhantomData_T, typename Value_T >
class BlockIDSorter< PhantomData_T, Value_T, false >
{
public:
   using Blocks_T = std::vector<std::vector<BlockID>>;
   BlockIDSorter( const Blocks_T & blocks ) : blocks_( blocks ) {}
   bool operator()( const Value_T & lhs, const Value_T & rhs ) const
   {
      return blocks_[ uint_c(lhs.first) ][ lhs.second ] < blocks_[ uint_c(rhs.first) ][ rhs.second ];
   }
private:
   const Blocks_T & blocks_;
};

template< typename pid_t, typename idx_t >
struct Node
{
   Node() : index_( std::make_pair( pid_t(0), std::numeric_limits< idx_t >::max() ) ) {}
   Node( const std::pair< pid_t, idx_t > & index ) : index_( index ) {}
   std::pair< pid_t, idx_t > index_;
   std::vector< Node > children_;
};

} // namespace internal


/**
 *  This class implements Hilber and Morton space filling curves for load balancing.
 *
 *  All algorithms are implemented to work levelwise. Load balancing with levels ignored is possible
 *  by specifying levelwise = false in the constructor.
**/
template< typename PhantomData_T >
class DynamicCurveBalance
{
public:

   using weight_t = typename PhantomData_T::weight_t;
   using pid_t = mpi::MPIRank;
   using idx_t = uint16_t; // limits the maximum number of blocks per process to 65536
   using Node = internal::Node<pid_t, idx_t>;

   DynamicCurveBalance( const bool hilbert = true, const bool allGather = true, const bool levelwise = true ) :
      hilbert_( hilbert ), allGather_( allGather ), levelwise_(levelwise)
   {}

   bool operator()( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                    std::set< uint_t > & processesToRecvFrom,
                    const PhantomBlockForest & phantomForest,
                    const uint_t iteration ) const;

   void setMaxBlocksPerProcess(const int maxBlocks) {maxBlocksPerProcess_ = maxBlocks;}
   int  getMaxBlocksPerProcess() const {return maxBlocksPerProcess_;}

private:

   void allGatherWeighted( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                           std::set< uint_t > & processesToRecvFrom,
                           const PhantomBlockForest & phantomForest ) const;
                                 
   void allGatherNoWeight( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                           std::set< uint_t > & processesToRecvFrom,
                           const PhantomBlockForest & phantomForest ) const;

   void masterWeighted( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                        std::set< uint_t > & processesToRecvFrom,
                        const PhantomBlockForest & phantomForest ) const;
                              
   void masterNoWeight( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                        std::set< uint_t > & processesToRecvFrom,
                        const PhantomBlockForest & phantomForest ) const;                              

   void hilbertOrderWeighted( const std::vector< std::vector< std::pair< BlockID, weight_t > > > & allBlocks,
                              std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                              const PhantomBlockForest & phantomForest ) const;
                               
   void hilbertOrderNoWeight( const std::vector< std::vector< BlockID > > & allBlocks,
                              std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                              const PhantomBlockForest & phantomForest ) const;

   void addBlockToForest( std::vector< shared_ptr< Node > > & forest,
                          const std::pair< pid_t, idx_t > & index, BlockID & id, const uint_t level ) const;

   void mortonOrderWeighted( const std::vector< std::vector< std::pair< BlockID, weight_t > > > & allBlocks,
                             std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                             const PhantomBlockForest & phantomForest ) const;
                               
   void mortonOrderNoWeight( const std::vector< std::vector< BlockID > > & allBlocks,
                             std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                             const PhantomBlockForest & phantomForest ) const;
                              
   void balanceWeighted( const std::vector< std::vector< std::pair< BlockID, weight_t > > > & allBlocks,
                         const std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                         std::vector< std::vector<pid_t> > & targets,
                         std::vector< std::set<pid_t> > & sender ) const;
                               
   void balanceNoWeight( const std::vector< std::vector< BlockID > > & allBlocks,
                         const std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                         std::vector< std::vector<pid_t> > & targets,
                         std::vector< std::set<pid_t> > & sender ) const;
                               
   void masterEnd( std::vector< std::vector<pid_t> > & targets,
                   std::vector< std::set<pid_t> > & sender,
                   std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                   std::set< uint_t > & processesToRecvFrom ) const;
                               
   void finalAssignment( const uint_t index, const std::vector< std::vector<pid_t> > & targets,
                         const std::vector< std::set<pid_t> > & sender,
                         std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                         std::set< uint_t > & processesToRecvFrom ) const;

   bool weightedBlocks() const
   {
      return ! std::is_same< PhantomData_T, NoPhantomData >::value;
   }
   
   template< typename T >
   pid_t pid_c( const T & value ) const { return numeric_cast< pid_t >( value ); }
   
   template< typename T >
   idx_t idx_c( const T & value ) const { return numeric_cast< idx_t >( value ); }
   
   

   bool hilbert_;
   bool allGather_;
   /// All gets for levels are wrapped like
   /// \code
   /// levelwise_ ? getCorrectLevel() : 0
   /// \endcode
   ///
   /// This allows to use the same algorithm for levelwise balancing as well as for balancing without levels.
   bool levelwise_;

   int  maxBlocksPerProcess_ = std::numeric_limits<int>::max(); //!< limits the maximum number of blocks per process
};



template< typename PhantomData_T >
bool DynamicCurveBalance< PhantomData_T >::operator()( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                                                                std::set< uint_t > & processesToRecvFrom,
                                                                const PhantomBlockForest & phantomForest, const uint_t ) const
{
   // Do not change or modifiy this check.
   // The Hilbert curve construction relies on "std::numeric_limits< idx_t >::max()" being an invalid number of blocks for a process
   WALBERLA_CHECK_LESS( targetProcess.size(), std::numeric_limits< idx_t >::max() );
   
   if( allGather_ )
   {
      if( weightedBlocks() )
         allGatherWeighted( targetProcess, processesToRecvFrom, phantomForest );
      else
         allGatherNoWeight( targetProcess, processesToRecvFrom, phantomForest );
   }
   else
   {
      if( weightedBlocks() )
         masterWeighted( targetProcess, processesToRecvFrom, phantomForest );
      else
         masterNoWeight( targetProcess, processesToRecvFrom, phantomForest );
   }

   return false;
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::allGatherWeighted( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                                                                       std::set< uint_t > & processesToRecvFrom,
                                                                       const PhantomBlockForest & phantomForest ) const
{  
   std::vector< std::pair< BlockID, weight_t > > localBlocks;

   for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it )
   {
      weight_t weight = it->first->template getData< PhantomData_T >().weight();
      localBlocks.push_back( std::make_pair( it->first->getId(), weight ) );
   }
   
   mpi::SendBuffer sendBuffer;
   mpi::RecvBuffer recvBuffer;
   
   sendBuffer << localBlocks;
   mpi::allGathervBuffer( sendBuffer, recvBuffer );
   sendBuffer.reset();
   
   const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );
      
   std::vector< std::vector< std::pair< BlockID, weight_t > > > allBlocks( processes ); // one vector for every process
   
   for( uint_t p = uint_t(0); p != processes; ++p )
      recvBuffer >> allBlocks[p];
   recvBuffer.reset();
   
   const uint_t numLevels = levelwise_ ? phantomForest.getNumberOfLevels() : uint_t(1);
   std::vector< std::vector< std::pair< pid_t, idx_t > > > blocksPerLevel( numLevels ); // for every level one vector of pair(source process ID, index in 'allBlocks')

   if( hilbert_ )
      hilbertOrderWeighted( allBlocks, blocksPerLevel, phantomForest );
   else
      mortonOrderWeighted( allBlocks, blocksPerLevel, phantomForest );

   std::vector< std::vector<pid_t> > targets( processes ); // for every process targets for all phantoms
   std::vector< std::set<pid_t> > sender( processes ); // for every process 'processesToRecvFrom'
   
   balanceWeighted( allBlocks, blocksPerLevel, targets, sender );

   finalAssignment( phantomForest.getBlockForest().getProcess(), targets, sender, targetProcess, processesToRecvFrom );
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::allGatherNoWeight( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                                                                       std::set< uint_t > & processesToRecvFrom,
                                                                       const PhantomBlockForest & phantomForest ) const
{
   std::vector< BlockID > localBlocks;

   for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it )
      localBlocks.push_back( it->first->getId() );
   
   mpi::SendBuffer sendBuffer;
   mpi::RecvBuffer recvBuffer;
   
   sendBuffer << localBlocks;
   mpi::allGathervBuffer( sendBuffer, recvBuffer );
   sendBuffer.reset();
   
   const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );
   
   std::vector< std::vector< BlockID > > allBlocks( processes ); // one vector for every process
   
   for( uint_t p = uint_t(0); p != processes; ++p )
      recvBuffer >> allBlocks[p];
   recvBuffer.reset();

   const uint_t numLevels = levelwise_ ? phantomForest.getNumberOfLevels() : uint_t(1);
   std::vector< std::vector< std::pair< pid_t, idx_t > > > blocksPerLevel( numLevels ); // for every level one vector of pair(source process ID, index in 'allBlocks')

   if( hilbert_ )
      hilbertOrderNoWeight( allBlocks, blocksPerLevel, phantomForest );
   else
      mortonOrderNoWeight( allBlocks, blocksPerLevel, phantomForest );
   
   std::vector< std::vector<pid_t> > targets( processes ); // for every process targets for all phantoms
   std::vector< std::set<pid_t> > sender( processes ); // for every process 'processesToRecvFrom'
   
   balanceNoWeight( allBlocks, blocksPerLevel, targets, sender );

   finalAssignment( phantomForest.getBlockForest().getProcess(), targets, sender, targetProcess, processesToRecvFrom );
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::masterWeighted( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                                                                    std::set< uint_t > & processesToRecvFrom,
                                                                    const PhantomBlockForest & phantomForest ) const
{
   std::vector< std::pair< BlockID, weight_t > > localBlocks;

   for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it )
   {
      weight_t weight = it->first->template getData< PhantomData_T >().weight();
      localBlocks.push_back( std::make_pair( it->first->getId(), weight ) );
   } 
   
   std::set< mpi::MPIRank > ranksToRecvFrom;

   if( mpi::MPIManager::instance()->rank() == 0 ) // do _NOT_ use WALBERLA_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      for( int rank = 1; rank < mpi::MPIManager::instance()->numProcesses(); ++rank )
         ranksToRecvFrom.insert( ranksToRecvFrom.end(), rank );
   }

   mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 2669 ); // blockforestglobalsortedid = 098 108 111 099 107 102 111 114 101 115 116 103 108 111 098 097 108 115 111 114 116 101 100 105 100
   bufferSystem.setReceiverInfo( ranksToRecvFrom, true );

   if( mpi::MPIManager::instance()->rank() != 0 ) // do _NOT_ use WALBERLA_NON_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      bufferSystem.sendBuffer( 0 ) << localBlocks;
   }

   bufferSystem.sendAll();
   
   std::vector< std::vector< std::pair< BlockID, weight_t > > > allBlocks; // one vector for every process (including root)
   std::vector< std::vector< std::pair< pid_t, idx_t > > > blocksPerLevel; // for every level one vector of pair(source process ID, index in 'allBlocks')

   std::vector< std::vector<pid_t> > targets; // for every process targets for all phantoms
   std::vector< std::set<pid_t> > sender; // for every process 'processesToRecvFrom'
   
   if( mpi::MPIManager::instance()->rank() == 0 ) // do _NOT_ use WALBERLA_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );
      
      allBlocks.resize( processes );
      for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
      {
         const uint_t source = uint_c( recvIt.rank() );
         WALBERLA_ASSERT( allBlocks[ source ].empty() );
         recvIt.buffer() >> allBlocks[ source ];
      }
      WALBERLA_ASSERT_EQUAL( mpi::MPIManager::instance()->rank(), 0 );
      WALBERLA_ASSERT( allBlocks[0].empty() );
      allBlocks[0] = localBlocks;

      const uint_t numLevels = levelwise_ ? phantomForest.getNumberOfLevels() : uint_t(1);
      blocksPerLevel.resize( numLevels );

      if( hilbert_ )
         hilbertOrderWeighted( allBlocks, blocksPerLevel, phantomForest );
      else
         mortonOrderWeighted( allBlocks, blocksPerLevel, phantomForest );

      targets.resize( processes );
      sender.resize( processes );
      
      balanceWeighted( allBlocks, blocksPerLevel, targets, sender );
   }
   else
   {
      // begin()/end() must also be called on each slave process in order to
      // properly finalize the communication
      WALBERLA_CHECK( bufferSystem.begin() == bufferSystem.end() );
   }
   
   masterEnd( targets, sender, targetProcess, processesToRecvFrom );   
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::masterNoWeight( std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                                                                    std::set< uint_t > & processesToRecvFrom,
                                                                    const PhantomBlockForest & phantomForest ) const
{
   std::vector< BlockID > localBlocks;

   for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it )
      localBlocks.push_back( it->first->getId() );
   
   std::set< mpi::MPIRank > ranksToRecvFrom;

   if( mpi::MPIManager::instance()->rank() == 0 ) // do _NOT_ use WALBERLA_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      for( int rank = 1; rank < mpi::MPIManager::instance()->numProcesses(); ++rank )
         ranksToRecvFrom.insert( ranksToRecvFrom.end(), rank );
   }

   mpi::BufferSystem bufferSystem( MPIManager::instance()->comm(), 2669 ); // blockforestglobalsortedid = 098 108 111 099 107 102 111 114 101 115 116 103 108 111 098 097 108 115 111 114 116 101 100 105 100
   bufferSystem.setReceiverInfo( ranksToRecvFrom, true );

   if( mpi::MPIManager::instance()->rank() != 0 ) // do _NOT_ use WALBERLA_NON_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      bufferSystem.sendBuffer( 0 ) << localBlocks;
   }

   bufferSystem.sendAll();
   
   std::vector< std::vector< BlockID > > allBlocks; // one vector for every process (including root)
   std::vector< std::vector< std::pair< pid_t, idx_t > > > blocksPerLevel; // for every level one vector of pair(source process ID, index in 'allBlocks')

   std::vector< std::vector<pid_t> > targets; // for every process targets for all phantoms
   std::vector< std::set<pid_t> > sender; // for every process 'processesToRecvFrom'
   
   if( mpi::MPIManager::instance()->rank() == 0 ) // do _NOT_ use WALBERLA_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );

      allBlocks.resize( processes );
      for( auto recvIt = bufferSystem.begin(); recvIt != bufferSystem.end(); ++recvIt )
      {
         const uint_t source = uint_c( recvIt.rank() );
         WALBERLA_ASSERT( allBlocks[ source ].empty() );
         recvIt.buffer() >> allBlocks[ source ];
      }
      WALBERLA_ASSERT_EQUAL( mpi::MPIManager::instance()->rank(), 0 );
      WALBERLA_ASSERT( allBlocks[0].empty() );
      allBlocks[0] = localBlocks;

      const uint_t numLevels = levelwise_ ? phantomForest.getNumberOfLevels() : uint_t(1);
      blocksPerLevel.resize( numLevels );

      if( hilbert_ )
         hilbertOrderNoWeight( allBlocks, blocksPerLevel, phantomForest );
      else
         mortonOrderNoWeight( allBlocks, blocksPerLevel, phantomForest );

      targets.resize( processes );
      sender.resize( processes );
      
      balanceNoWeight( allBlocks, blocksPerLevel, targets, sender );
   }
   else
   {
      // begin()/end() must also be called on each slave process in order to
      // properly finalize the communication
      WALBERLA_CHECK( bufferSystem.begin() == bufferSystem.end() );
   }
   
   masterEnd( targets, sender, targetProcess, processesToRecvFrom );   
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::hilbertOrderWeighted( const std::vector< std::vector< std::pair< BlockID, typename PhantomData_T::weight_t > > > & allBlocks,
                                                                          std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                                                                          const PhantomBlockForest & phantomForest) const
{
   // construct forest of octrees

   const auto & blockforest = phantomForest.getBlockForest();

   std::vector< shared_ptr< Node > > forest( blockforest.getXSize() * blockforest.getYSize() * blockforest.getZSize() );

   for( pid_t p = pid_t(0); p != pid_c( allBlocks.size() ); ++p )
   {
      for( idx_t i = idx_t(0); i != idx_c( allBlocks[uint_c(p)].size() ); ++i )
      {
         BlockID id( allBlocks[uint_c(p)][i].first );
         addBlockToForest( forest, std::make_pair(p,i), id, blockforest.getLevelFromBlockId(id) );
      }
   }

   // traverse forest in hilbert order

   uint_t y( uint_t(0) );
   uint_t x( uint_t(0) );

   uint_t yLoopEnd( blockforest.getYSize() );
   uint_t xLoopEnd( blockforest.getXSize() );

   for( uint_t z = 0; z != blockforest.getZSize(); ++z )
   {
      while( y != yLoopEnd )
      {
         if( yLoopEnd == 0 ) --y;
         while( x != xLoopEnd )
         {
            if( xLoopEnd == 0 ) --x;

            WALBERLA_ASSERT_LESS( z*blockforest.getXSize()*blockforest.getYSize() + y*blockforest.getXSize() + x, forest.size() );

            auto & root = forest[ z*blockforest.getXSize()*blockforest.getYSize() + y*blockforest.getXSize() + x ];

            if( root )
            {
               std::stack< Node * > stack;
               std::stack< uint_t > orientation;

               stack.push( root.get() );
               orientation.push( uint_t(0) );

               while( !stack.empty() )
               {
                  Node * const node = stack.top();
                  uint_t hilbertIndex = orientation.top();

                  stack.pop();
                  orientation.pop();

                  WALBERLA_ASSERT_NOT_NULLPTR( node );

                  if( ! node->children_.empty() )
                  {
                     for( uint_t c = 8; c-- != 0; )
                     {
                        stack.push( &(node->children_[ hilbertOrder[hilbertIndex][c] ]) );
                        orientation.push( hilbertOrientation[hilbertIndex][c] );
                     }
                  }
                  else
                  {
                     auto & index = node->index_;
                     const uint_t level = levelwise_ ? blockforest.getLevelFromBlockId( allBlocks[ uint_c(index.first) ][ index.second ].first ) : uint_t(0);
                     WALBERLA_ASSERT_LESS( level, levelwise_ ? phantomForest.getNumberOfLevels() : uint_t(1) );
                     blocksPerLevel[ level ].push_back( index );
                  }
               }
            }

            if( xLoopEnd != 0 ) ++x;
         }
         WALBERLA_ASSERT_EQUAL( x, xLoopEnd );
         xLoopEnd = ( xLoopEnd == 0 ) ? blockforest.getXSize() : 0;
         if( yLoopEnd != 0 ) ++y;
      }
      WALBERLA_ASSERT_EQUAL( y, yLoopEnd );
      yLoopEnd = ( yLoopEnd == 0 ) ? blockforest.getYSize() : 0;
   }
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::hilbertOrderNoWeight( const std::vector< std::vector< BlockID > > & allBlocks,
                                                                          std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                                                                          const PhantomBlockForest & phantomForest) const
{
   // construct forest of octrees

   const auto & blockforest = phantomForest.getBlockForest();

   std::vector< shared_ptr< Node > > forest( blockforest.getXSize() * blockforest.getYSize() * blockforest.getZSize() );

   for( pid_t p = pid_t(0); p != pid_c( allBlocks.size() ); ++p )
   {
      for( idx_t i = idx_t(0); i != idx_c( allBlocks[uint_c(p)].size() ); ++i )
      {
         BlockID id( allBlocks[uint_c(p)][i] );
         addBlockToForest( forest, std::make_pair(p,i), id, blockforest.getLevelFromBlockId(id) );
      }
   }

   // traverse forest in hilbert order

   uint_t y( uint_t(0) );
   uint_t x( uint_t(0) );

   uint_t yLoopEnd( blockforest.getYSize() );
   uint_t xLoopEnd( blockforest.getXSize() );

   for( uint_t z = 0; z != blockforest.getZSize(); ++z )
   {
      while( y != yLoopEnd )
      {
         if( yLoopEnd == 0 ) --y;
         while( x != xLoopEnd )
         {
            if( xLoopEnd == 0 ) --x;

            WALBERLA_ASSERT_LESS( z*blockforest.getXSize()*blockforest.getYSize() + y*blockforest.getXSize() + x, forest.size() );

            auto & root = forest[ z*blockforest.getXSize()*blockforest.getYSize() + y*blockforest.getXSize() + x ];

            if( root )
            {
               std::stack< Node * > stack;
               std::stack< uint_t > orientation;

               stack.push( root.get() );
               orientation.push( uint_t(0) );

               while( !stack.empty() )
               {
                  Node * const node = stack.top();
                  uint_t hilbertIndex = orientation.top();

                  stack.pop();
                  orientation.pop();

                  WALBERLA_ASSERT_NOT_NULLPTR( node );

                  if( ! node->children_.empty() )
                  {
                     for( uint_t c = 8; c-- != 0; )
                     {
                        stack.push( &(node->children_[ hilbertOrder[hilbertIndex][c] ]) );
                        orientation.push( hilbertOrientation[hilbertIndex][c] );
                     }
                  }
                  else
                  {
                     auto & index = node->index_;
                     const uint_t level = levelwise_ ? blockforest.getLevelFromBlockId( allBlocks[ uint_c(index.first) ][ index.second ] ) : uint_t(0);
                     WALBERLA_ASSERT_LESS( level , levelwise_ ? phantomForest.getNumberOfLevels() : uint_t(1) );
                     blocksPerLevel[ level ].push_back( index );
                  }
               }
            }

            if( xLoopEnd != 0 ) ++x;
         }
         WALBERLA_ASSERT_EQUAL( x, xLoopEnd );
         xLoopEnd = ( xLoopEnd == 0 ) ? blockforest.getXSize() : 0;
         if( yLoopEnd != 0 ) ++y;
      }
      WALBERLA_ASSERT_EQUAL( y, yLoopEnd );
      yLoopEnd = ( yLoopEnd == 0 ) ? blockforest.getYSize() : 0;
   }
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::addBlockToForest( std::vector< shared_ptr< Node > > & forest,
                                                                      const std::pair< pid_t, idx_t > & index, BlockID & id, const uint_t level ) const
{
   std::stack< uint_t > path;

   for( uint_t l = uint_t(0); l != level; ++l )
   {
      path.push( id.getBranchId() );
      id.removeBranchId();
   }
   path.push( id.getTreeIndex() );

   if( ! (forest[ path.top() ]) )
   {
      if( path.size() == uint_t(1) )
         forest[ path.top() ] = walberla::make_shared< Node >( index );
      else
      {
         forest[ path.top() ] = walberla::make_shared< Node >();
         forest[ path.top() ]->children_.resize( uint_t(8) );
      }
   }

   Node * node = forest[ path.top() ].get();
   path.pop();

   while( !path.empty() )
   {
      WALBERLA_ASSERT_EQUAL( node->children_.size(), uint_t(8) );
      WALBERLA_ASSERT_LESS( path.top(), uint_t(8) );

      if( path.size() == uint_t(1) )
      {
         WALBERLA_ASSERT_EQUAL( node->children_[ path.top() ].index_.first, pid_t(0) );
         WALBERLA_ASSERT_EQUAL( node->children_[ path.top() ].index_.second, std::numeric_limits< idx_t >::max() );
         WALBERLA_ASSERT( node->children_[ path.top() ].children_.empty() );
         node->children_[ path.top() ].index_ = index;
      }
      else if( node->children_[ path.top() ].children_.empty() )
      {
         WALBERLA_ASSERT_EQUAL( node->index_.first, pid_t(0) );
         WALBERLA_ASSERT_EQUAL( node->index_.second, std::numeric_limits< idx_t >::max() );
         node->children_[ path.top() ].children_.resize( uint_t(8) );
      }

      node = &(node->children_[ path.top() ]);
      path.pop();
   }
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::mortonOrderWeighted( const std::vector< std::vector< std::pair< BlockID, typename PhantomData_T::weight_t > > > & allBlocks,
                                                                         std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                                                                         const PhantomBlockForest & phantomForest) const
{
   const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );
   
   for( uint_t p = uint_t(0); p != processes; ++p )
   {
      for( uint_t i = uint_t(0); i != allBlocks[p].size(); ++i )
      {
         uint_t level = levelwise_ ? phantomForest.getBlockForest().getLevelFromBlockId( allBlocks[p][i].first ) : uint_t(0);
         WALBERLA_ASSERT_LESS( level, blocksPerLevel.size() );
         blocksPerLevel[ level ].push_back( std::make_pair( pid_c(p), idx_c(i) ) );
      }
   }
      
#if defined(_OPENMP) && ((__INTEL_COMPILER < 1700) || (__INTEL_COMPILER > 1900)) // Disable OpenMP for Intel 2018/2019 due to a bug
   #pragma omp parallel for schedule(static)
#endif
   for( int i = 0; i < int_c( blocksPerLevel.size() ); ++i )
   {
      std::vector< std::pair< pid_t, idx_t > > & blocks = blocksPerLevel[uint_c(i)];
      internal::BlockIDSorter< PhantomData_T, std::pair< pid_t, idx_t >, true > sorter( allBlocks );
      std::sort( blocks.begin(), blocks.end(), sorter );
   }
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::mortonOrderNoWeight( const std::vector< std::vector< BlockID > > & allBlocks,
                                                                         std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                                                                         const PhantomBlockForest & phantomForest) const
{
   const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );
   
   for( uint_t p = uint_t(0); p != processes; ++p )
   {
      for( uint_t i = uint_t(0); i != allBlocks[p].size(); ++i )
      {
         uint_t level = levelwise_ ? phantomForest.getBlockForest().getLevelFromBlockId( allBlocks[p][i] ) : uint_t(0);
         WALBERLA_ASSERT_LESS( level, blocksPerLevel.size() );
         blocksPerLevel[ level ].push_back( std::make_pair( pid_c(p), idx_c(i) ) );
      }
   }
   
#ifdef _OPENMP
   #pragma omp parallel for schedule(static)
#endif
   for( int i = 0; i < int_c( blocksPerLevel.size() ); ++i )
   {
      std::vector< std::pair< pid_t, idx_t > > & blocks = blocksPerLevel[uint_c(i)];
      internal::BlockIDSorter< PhantomData_T, std::pair< pid_t, idx_t >, false > sorter( allBlocks );
      std::sort( blocks.begin(), blocks.end(), sorter );
   }
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::balanceWeighted( const std::vector< std::vector< std::pair< BlockID, typename PhantomData_T::weight_t > > > & allBlocks,
                                                                     const std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                                                                     std::vector< std::vector<pid_t> > & targets,
                                                                     std::vector< std::set<pid_t> > & sender ) const
{
   const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );
   
   for( uint_t p = uint_t(0); p != processes; ++p )
      targets[p].resize( allBlocks[p].size() );
   
   for( uint_t i = 0; i < blocksPerLevel.size(); ++i )
   {
      const std::vector< std::pair< pid_t, idx_t > > & blocks = blocksPerLevel[i];
      
      long double totalWeight( 0 );
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
      {
         WALBERLA_ASSERT_LESS( block->first, allBlocks.size() );
         WALBERLA_ASSERT_LESS( block->second, allBlocks[ uint_c( block->first ) ].size() );
         totalWeight += numeric_cast< long double >( allBlocks[ uint_c( block->first ) ][ block->second ].second );
      }
      
      uint_t c( uint_t(0) );
      for( uint_t p = uint_t(0); p != processes; ++p )
      {
         const long double pWeight = totalWeight / numeric_cast< long double >( processes - p );
         long double weight( 0 );
         int numBlocks( 0 );
         while( c < blocks.size() &&
                ( isIdentical(weight, 0.0l) || 
                  std::abs( pWeight - weight - numeric_cast< long double >( allBlocks[ uint_c( blocks[c].first ) ][ blocks[c].second ].second ) ) <=
                  std::abs( pWeight - weight ) ) &&
                numBlocks < maxBlocksPerProcess_ )
         {
            targets[ uint_c( blocks[c].first ) ][ blocks[c].second ] = pid_c(p);
            sender[p].insert( blocks[c].first );
            const long double addedWeight = numeric_cast< long double >( allBlocks[ uint_c( blocks[c].first ) ][ blocks[c].second ].second );
            weight += addedWeight;
            totalWeight -= addedWeight;
            ++c;
            ++numBlocks;
         }
      }
      while( c < blocks.size() )
      {
         targets[ uint_c( blocks[c].first ) ][ blocks[c].second ] = pid_c( processes - uint_t(1) );
         sender[ processes - uint_t(1) ].insert( blocks[c].first );
         ++c;
      }
   }
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::balanceNoWeight( const std::vector< std::vector< BlockID > > & allBlocks,
                                                                     const std::vector< std::vector< std::pair< pid_t, idx_t > > > & blocksPerLevel,
                                                                     std::vector< std::vector<pid_t> > & targets,
                                                                     std::vector< std::set<pid_t> > & sender ) const
{
   const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );   
   
   for( uint_t p = uint_t(0); p != processes; ++p )
      targets[p].resize( allBlocks[p].size() );
   
   for( uint_t i = 0; i < blocksPerLevel.size(); ++i )
   {
      const std::vector< std::pair< pid_t, idx_t > > & blocks = blocksPerLevel[i];

      const uint_t base = uint_c( blocks.size() ) / processes;
      const uint_t rest = uint_c( blocks.size() ) % processes;

      uint_t c( uint_t(0) );
      for( uint_t p = uint_t(0); p != processes; ++p )
      {
         const uint_t numBlocks = ( p < rest ) ? ( base + uint_t(1) ) : base;
         for( uint_t j = uint_t(0); j < numBlocks; ++j )
         {
            WALBERLA_ASSERT_LESS( c, blocks.size() );
            targets[ uint_c( blocks[c].first ) ][ blocks[c].second ] = pid_c(p);
            sender[p].insert( blocks[c].first );
            ++c;
         }
      }
   }
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::masterEnd( std::vector< std::vector<pid_t> > & targets,
                                                               std::vector< std::set<pid_t> > & sender,
                                                               std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                                                               std::set< uint_t > & processesToRecvFrom ) const
{
#ifndef NDEBUG
   for( uint_t p = uint_t(0); p != targets.size(); ++p )
   {
      for( auto t = targets[p].begin(); t != targets[p].end(); ++t )
      {
         WALBERLA_ASSERT_LESS( *t, sender.size() );
         WALBERLA_ASSERT( sender[uint_c(*t)].find( int_c(p) ) != sender[uint_c(*t)].end() );
      }
   }
   for( uint_t p = uint_t(0); p != sender.size(); ++p )
   {
      for( auto s = sender[p].begin(); s != sender[p].end(); ++s )
      {
         WALBERLA_ASSERT_LESS( *s, targets.size() );
         WALBERLA_ASSERT( std::find( targets[uint_c(*s)].begin(), targets[uint_c(*s)].end(), int_c(p) ) != targets[uint_c(*s)].end() );
      }
   }
#endif

   std::set< mpi::MPIRank > ranksToRecvFrom;

   if( mpi::MPIManager::instance()->rank() != 0 ) // do _NOT_ use WALBERLA_NON_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      ranksToRecvFrom.insert(0);
   }

   mpi::BufferSystem resultsBufferSystem( MPIManager::instance()->comm(), 2670 ); // blockforestglobalsortedid = 098 108 111 099 107 102 111 114 101 115 116 103 108 111 098 097 108 115 111 114 116 101 100 105 100 + 1
   resultsBufferSystem.setReceiverInfo( ranksToRecvFrom, true );

   if( mpi::MPIManager::instance()->rank() == 0 ) // do _NOT_ use WALBERLA_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      for( int rank = 1; rank < mpi::MPIManager::instance()->numProcesses(); ++rank )
         resultsBufferSystem.sendBuffer( rank ) << targets[uint_c(rank)] << sender[uint_c(rank)];
   }

   resultsBufferSystem.sendAll();

   if( mpi::MPIManager::instance()->rank() != 0 ) // do _NOT_ use WALBERLA_NON_ROOT_SECTION ! (-> buffer system must use 'comm' which corresponds to 'rank' / block structure communicator + ranks)
   {
      for( auto recvIt = resultsBufferSystem.begin(); recvIt != resultsBufferSystem.end(); ++recvIt )
      {
         WALBERLA_ASSERT_EQUAL( recvIt.rank(), 0 );
         WALBERLA_ASSERT_EQUAL( targets.size(), uint_t(0) );
         WALBERLA_ASSERT_EQUAL( sender.size(), uint_t(0) );

         targets.resize(1);
         sender.resize(1);
         recvIt.buffer() >> targets[0] >> sender[0];
      }
   }
   else
   {
      // begin()/end() must also be called on each slave process in order to
      // properly finalize the communication
      WALBERLA_CHECK( resultsBufferSystem.begin() == resultsBufferSystem.end() );
   }
   
   finalAssignment( uint_t(0), targets, sender, targetProcess, processesToRecvFrom );
}



template< typename PhantomData_T >
void DynamicCurveBalance< PhantomData_T >::finalAssignment( const uint_t index, const std::vector< std::vector<pid_t> > & targets,
                                                                     const std::vector< std::set<pid_t> > & sender,
                                                                     std::vector< std::pair< const PhantomBlock *, uint_t > > & targetProcess,
                                                                     std::set< uint_t > & processesToRecvFrom ) const
{
   WALBERLA_ASSERT_GREATER( targets.size(), index );
   WALBERLA_ASSERT_EQUAL( targetProcess.size(), targets[index].size() );
   for( uint_t i = uint_t(0); i != targetProcess.size(); ++i )
      targetProcess[i].second = uint_c( targets[index][i] );

   for( auto s = sender[index].begin(); s != sender[index].end(); ++s )
      processesToRecvFrom.insert( uint_c(*s) ) ;
}

template< typename PhantomData_T >
using DynamicLevelwiseCurveBalance [[deprecated("Use DynamicCurveBalance with the corresponding argument!")]] = DynamicCurveBalance<PhantomData_T> ;

} // namespace blockforest
} // namespace walberla
