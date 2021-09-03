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
//! \file UniformDirectScheme.impl.h
//! \ingroup blockforest
//! \author MartinBauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "UniformDirectScheme.h"
#include "blockforest/BlockNeighborhoodSection.h"

namespace walberla {
namespace blockforest {
namespace communication {

//===================================================================================================================
//
//  BufferlessScheme::CommInfo
//
//===================================================================================================================

template< typename Stencil >
bool UniformDirectScheme<Stencil>::CommInfo::sortByLocal( const CommInfo & lhs, const CommInfo & rhs )
{
   // lexicographical compare in order:  dataIdx, localBlockId, remoteBlockId
   if ( lhs.dataIdx < rhs.dataIdx ) return true;
   if ( lhs.dataIdx > rhs.dataIdx ) return false;
   WALBERLA_ASSERT_EQUAL( lhs.dataIdx, rhs.dataIdx );

   if( lhs.localBlockId < rhs.localBlockId ) return true;
   if( lhs.localBlockId > rhs.localBlockId ) return false;
   WALBERLA_ASSERT_EQUAL( lhs.localBlockId, rhs.localBlockId );

   if( lhs.remoteBlockId < rhs.remoteBlockId ) return true;
   if( lhs.remoteBlockId > rhs.remoteBlockId ) return false;
   WALBERLA_ASSERT_EQUAL( lhs.remoteBlockId, rhs.remoteBlockId );

   return lhs.dir < rhs.dir;
}


template< typename Stencil >
bool UniformDirectScheme<Stencil>::CommInfo::sortByRemote( const CommInfo & lhs, const CommInfo & rhs )
{
   // lexicographical compare in order:  dataIdx, remoteBlockId, localBlockId
   if ( lhs.dataIdx < rhs.dataIdx ) return true;
   if ( lhs.dataIdx > rhs.dataIdx ) return false;
   WALBERLA_ASSERT_EQUAL( lhs.dataIdx, rhs.dataIdx );

   if( lhs.remoteBlockId < rhs.remoteBlockId ) return true;
   if( lhs.remoteBlockId > rhs.remoteBlockId ) return false;
   WALBERLA_ASSERT_EQUAL( lhs.remoteBlockId, rhs.remoteBlockId );

   if( lhs.localBlockId < rhs.localBlockId ) return true;
   if( lhs.localBlockId > rhs.localBlockId ) return false;
   WALBERLA_ASSERT_EQUAL( lhs.localBlockId, rhs.localBlockId );

   return stencil::inverseDir[lhs.dir] < stencil::inverseDir[rhs.dir];
}



//===================================================================================================================
//
//  BufferlessScheme
//
//===================================================================================================================

template< typename Stencil >
inline void UniformDirectScheme<Stencil>::communicate()
{
   startCommunication();
   wait();
}


template< typename Stencil >
void UniformDirectScheme<Stencil>::setup()
{
   if ( ! setupRequired_ )
      return;

   WALBERLA_ASSERT( !communicationRunning_ );

   WALBERLA_DEBUG_SECTION() {
      for( auto it = mpiRequests_.begin(); it != mpiRequests_.end(); ++it )
      {
         WALBERLA_ASSERT_EQUAL( *it, MPI_REQUEST_NULL );
      }
   }

   mpiDatatypes_.clear();
   mpiRequests_.clear();
   sendInfos_.clear();
   recvInfos_.clear();

   auto forest = blockForest_.lock();
   for( auto it = forest->begin(); it != forest->end(); ++it )
   {
      Block * block = dynamic_cast< Block * >( it.get() );

      if( !selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
         continue;

      for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
      {
         const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );

         if( block->getNeighborhoodSectionSize(neighborIdx) == uint_t(0) )
            continue;

         WALBERLA_ASSERT( block->neighborhoodSectionHasEquallySizedBlock(neighborIdx) );
         WALBERLA_ASSERT_EQUAL( block->getNeighborhoodSectionSize(neighborIdx), uint_t(1) );

         const BlockID & nBlockId = block->getNeighborId( neighborIdx, uint_t(0) );

         if( !selectable::isSetSelected( block->getNeighborState( neighborIdx, uint_t(0) ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
            continue;

         auto nProcess = block->getNeighborProcess( neighborIdx, uint_t( 0 ) );

         for( uint_t dataIdx = 0; dataIdx < dataInfos_.size(); ++dataIdx )
         {
            // These two infos do not belong to the same data exchange
            // here we just say for example: "recv from west", "send to west"
            // the sorting according to communication partners is done in a second step
            CommInfo info = { dataIdx, block->getId(), nBlockId, nProcess, *dir };
            sendInfos_.push_back( info );
            recvInfos_.push_back( info );
         }
      }
   }

   // Sort sends and receives so they match correctly with remote sends and receives
   // Otherwise the following problem could occur:
   // Assume two processes P1 and P2 with two blocks each P1_A, P1_B and P2_A and P2_B
   // to make sure that on both processes the Send/Recvs for block A are scheduled first
   // this sorting has to be done.
   std::sort( sendInfos_.begin(), sendInfos_.end(), CommInfo::sortByLocal );
   std::sort( recvInfos_.begin(), recvInfos_.end(), CommInfo::sortByRemote );

   for( auto it = sendInfos_.begin(); it != sendInfos_.end(); ++it )
   {
      auto block = forest->getBlock( it->localBlockId );
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      mpiDatatypes_.push_back( dataInfos_[ it->dataIdx ]->getSendDatatype( block, it->dir )  );
   }

   for( auto it = recvInfos_.begin(); it != recvInfos_.end(); ++it )
   {
      auto block = forest->getBlock( it->localBlockId );
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      mpiDatatypes_.push_back( dataInfos_[ it->dataIdx ]->getRecvDatatype( block, it->dir )  );
   }

   mpiRequests_.resize( sendInfos_.size() + recvInfos_.size(), MPI_REQUEST_NULL );

   setupRequired_ = false;
}


template< typename Stencil >
void UniformDirectScheme<Stencil>::startCommunication()
{
   WALBERLA_ASSERT( !communicationRunning_ );

   setup();

   if( mpiRequests_.empty() )
      return;

   communicationRunning_ = true;

   const MPI_Comm comm = MPIManager::instance()->comm();

   auto forest = blockForest_.lock();

   auto requestIt  = mpiRequests_.begin();
   auto datatypeIt = mpiDatatypes_.begin();
   for( auto it = sendInfos_.begin(); it != sendInfos_.end(); ++it, ++requestIt, ++datatypeIt )
   {
      WALBERLA_ASSERT_UNEQUAL( requestIt,  mpiRequests_.end()  );
      WALBERLA_ASSERT_UNEQUAL( datatypeIt, mpiDatatypes_.end() );

      domain_decomposition::IBlock * const block = forest->getBlock( it->localBlockId );
      WALBERLA_ASSERT_NOT_NULLPTR( block );

      WALBERLA_ASSERT_EQUAL( *requestIt, MPI_REQUEST_NULL );
      WALBERLA_ASSERT_UNEQUAL( **datatypeIt, MPI_DATATYPE_NULL );

      const shared_ptr<UniformMPIDatatypeInfo> & dataInfo = dataInfos_ [ it->dataIdx ];
      MPI_Isend( dataInfo->getSendPointer( block, it->dir ),
                 dataInfo->getNumberOfItemsToCommunicate( block, it->dir ),
                 **datatypeIt,
                 int_c( it->remoteProcess ),
                 tag_,
                 comm,
                 &( *requestIt ) );

      WALBERLA_ASSERT_UNEQUAL( *requestIt, MPI_REQUEST_NULL );
   }

   for( auto it = recvInfos_.begin(); it != recvInfos_.end(); ++it, ++requestIt, ++datatypeIt )
   {
      WALBERLA_ASSERT_UNEQUAL( requestIt,  mpiRequests_.end()  );
      WALBERLA_ASSERT_UNEQUAL( datatypeIt, mpiDatatypes_.end() );

      domain_decomposition::IBlock * const block = forest->getBlock( it->localBlockId );
      WALBERLA_ASSERT_NOT_NULLPTR( block );

      WALBERLA_ASSERT_EQUAL(   *requestIt,  MPI_REQUEST_NULL );
      WALBERLA_ASSERT_UNEQUAL( **datatypeIt, MPI_DATATYPE_NULL );

      const shared_ptr<UniformMPIDatatypeInfo> & dataInfo = dataInfos_ [ it->dataIdx ];

      MPI_Irecv( dataInfo->getRecvPointer( block, it->dir ),
                 dataInfo->getNumberOfItemsToCommunicate( block, it->dir ),
                 **datatypeIt,
                 int_c( it->remoteProcess ),
                 tag_,
                 comm,
                 &( *requestIt ) );

      WALBERLA_ASSERT_UNEQUAL( *requestIt, MPI_REQUEST_NULL );
   }

   WALBERLA_ASSERT_EQUAL( requestIt, mpiRequests_.end() );
   WALBERLA_ASSERT_EQUAL( datatypeIt, mpiDatatypes_.end() );
}



template< typename Stencil >
void UniformDirectScheme<Stencil>::wait()
{
   if( mpiRequests_.empty() )
      return;

   WALBERLA_ASSERT( communicationRunning_ );

   MPI_Waitall( int_c( mpiRequests_.size() ), &(mpiRequests_.front()), MPI_STATUSES_IGNORE );

   communicationRunning_ = false;
}


template< typename Stencil >
void UniformDirectScheme<Stencil>::addDataToCommunicate(  const shared_ptr<UniformMPIDatatypeInfo> & dataInfo )
{
   WALBERLA_ASSERT( !communicationRunning_ );

   setupRequired_ = true;
   dataInfos_.push_back( dataInfo );
}



} // namespace communication
} // namespace blockforest
} // namespace walberla


