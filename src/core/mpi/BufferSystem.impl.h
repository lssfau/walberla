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
//! \file GenericBufferSystem.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"
#include "core/debug/CheckFunctions.h"


namespace walberla {
namespace mpi {


template< typename Rb, typename Sb>
std::set<int> GenericBufferSystem<Rb, Sb>::activeTags_;


//======================================================================================================================
//
//  Iterator
//
//======================================================================================================================


template< typename Rb, typename Sb>
GenericBufferSystem<Rb, Sb>::iterator::iterator( GenericBufferSystem<Rb, Sb> & bufferSystem, bool begin )
    : bufferSystem_( bufferSystem), currentRecvBuffer_( nullptr ), currentSenderRank_( -1 )
{
   if ( begin ) // init iterator
      ++(*this);
}

template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::iterator::operator++()
{
   currentRecvBuffer_ = bufferSystem_.waitForNext( currentSenderRank_ );
   if ( ! currentRecvBuffer_ ) {
      WALBERLA_ASSERT_EQUAL( currentSenderRank_, -1 );
   } else
   {
      bufferSystem_.bytesReceived_ += currentRecvBuffer_->size() * sizeof(RecvBuffer::ElementType);
      bufferSystem_.numberOfReceives_ += 1;
   }
}

template< typename Rb, typename Sb>
bool GenericBufferSystem<Rb, Sb>::iterator::operator==( const typename GenericBufferSystem<Rb, Sb>::iterator & other )
{
   // only equality checks with end iterators are allowed
   WALBERLA_ASSERT( other.currentSenderRank_ == -1 || currentSenderRank_ == -1 );

   return ( currentSenderRank_ == other.currentSenderRank_ );
}

template< typename Rb, typename Sb>
bool GenericBufferSystem<Rb, Sb>::iterator::operator!=( const typename GenericBufferSystem<Rb, Sb>::iterator & other )
{
   // only equality checks with end iterators are allowed
   WALBERLA_ASSERT( other.currentSenderRank_ == -1 || currentSenderRank_ == -1 );

   return ( currentSenderRank_ != other.currentSenderRank_ );
}


template< typename Rb, typename Sb>
template<typename Range>
void GenericBufferSystem<Rb, Sb>::setReceiverInfo( const Range & range, bool changingSize )
{
   setReceiverInfo( range.begin(), range.end(), changingSize );
}

template< typename Rb, typename Sb>
template<typename RankIter>
void GenericBufferSystem<Rb, Sb>::setReceiverInfo( RankIter rankBegin, RankIter rankEnd, bool changingSize )
{
   WALBERLA_ASSERT( ! communicationRunning_ );

   recvInfos_.clear();
   for ( auto it = rankBegin; it != rankEnd; ++it )
   {
      const MPIRank sender = *it;
      recvInfos_[ sender ].size = INVALID_SIZE;
   }

   sizeChangesEverytime_ = changingSize;
   setCommunicationType( false );
}

template< typename Rb, typename Sb>
inline size_t GenericBufferSystem<Rb, Sb>::size() const
{
   size_t sum = 0;
   for( auto iter = sendInfos_.begin(); iter != sendInfos_.end(); ++iter )
   {
      sum += iter->second.buffer.size();
   }
   return sum;
}




//======================================================================================================================
//
//  Constructors
//
//======================================================================================================================

template< typename Rb, typename Sb>
GenericBufferSystem<Rb, Sb>::GenericBufferSystem( const MPI_Comm & communicator, int tag )
   : knownSizeComm_  ( communicator, tag ),
     unknownSizeComm_( communicator, tag ),
     unknownSizeCommIProbe_( communicator, tag ),
     noMPIComm_( communicator, tag ),
     currentComm_    ( nullptr ),
     sizeChangesEverytime_( true ),
     communicationRunning_( false )
{
}

template< typename Rb, typename Sb>
GenericBufferSystem<Rb, Sb>::GenericBufferSystem( const GenericBufferSystem &other )
   : knownSizeComm_  ( other.knownSizeComm_.getCommunicator(), other.knownSizeComm_.getTag() ),
     unknownSizeComm_( other.knownSizeComm_.getCommunicator(), other.knownSizeComm_.getTag() ),
     unknownSizeCommIProbe_( other.knownSizeComm_.getCommunicator(), other.knownSizeComm_.getTag() ),
     noMPIComm_      ( other.knownSizeComm_.getCommunicator(), other.knownSizeComm_.getTag() ),
     currentComm_ ( nullptr ),
     sizeChangesEverytime_( other.sizeChangesEverytime_ ),
     communicationRunning_( other.communicationRunning_ ),
     recvInfos_( other.recvInfos_ ),
     sendInfos_( other.sendInfos_ )
{
   WALBERLA_ASSERT( !communicationRunning_, "Can't copy GenericBufferSystem while communication is running" );
   if( other.currentComm_ == &other.knownSizeComm_ )
      currentComm_ = &knownSizeComm_;
   else if ( other.currentComm_ == &other.unknownSizeComm_ )
      currentComm_ = &unknownSizeComm_;
   else if ( other.currentComm_ == &other.unknownSizeCommIProbe_ )
      currentComm_ = &unknownSizeCommIProbe_;
   else if ( other.currentComm_ == &other.noMPIComm_ )
      currentComm_ = &noMPIComm_;
   else
      currentComm_ = nullptr; // receiver information not yet set
}

template< typename Rb, typename Sb>
GenericBufferSystem<Rb, Sb> & GenericBufferSystem<Rb, Sb>::operator=( const GenericBufferSystem<Rb, Sb> & other )
{
   WALBERLA_ASSERT( !communicationRunning_, "Can't copy GenericBufferSystem while communication is running" );

   sizeChangesEverytime_ = other.sizeChangesEverytime_;
   communicationRunning_ = other.communicationRunning_;
   recvInfos_ = other.recvInfos_;
   sendInfos_ = other.sendInfos_;

   if( other.currentComm_ == &other.knownSizeComm_ )
      currentComm_ = &knownSizeComm_;
   else if ( other.currentComm_ == &other.unknownSizeComm_ )
      currentComm_ = &unknownSizeComm_;
   else if ( other.currentComm_ == &other.unknownSizeCommIProbe_ )
      currentComm_ = &unknownSizeCommIProbe_;
   else if ( other.currentComm_ == &other.noMPIComm_ )
      currentComm_ = &noMPIComm_;
   else
      currentComm_ = nullptr; // receiver information not yet set

   return *this;
}

//======================================================================================================================
//
//  Receive Information Setup
//
//======================================================================================================================


//**********************************************************************************************************************
/*! Sets receiver information, when message sizes are unknown
*
* \param ranksToRecvFrom  set of all ranks where messages are received
* \param changingSize     true if the message size is different in each communication step.
*                         If false the message size is exchanged once and is expected to be constant
*                         If true the message size is exchanged before each communication step.
*                         The behavior can be changed later one using setReceiverInfo() or sizeHasChanged().
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::setReceiverInfo( const std::set<MPIRank> & ranksToRecvFrom, bool changingSize )
{
   WALBERLA_ASSERT( ! communicationRunning_ );

   recvInfos_.clear();
   for ( auto it = ranksToRecvFrom.begin(); it != ranksToRecvFrom.end(); ++it )
   {
      const MPIRank sender = *it;
      recvInfos_[ sender ].size = INVALID_SIZE;
   }

   sizeChangesEverytime_ = changingSize;
   setCommunicationType( false ); // no size information on first run -> UnknownSizeCommunication
}



//**********************************************************************************************************************
/*! Sets receiver information, when message sizes are known
*
* \param ranksToRecvFrom  Map containing all ranks, where messages are received from, as keys
*                         and the message sizes as values.
*                         The message sizes are expected to be constant for all communication step until
*                         behavior is changed with setReceiverInfo*() or sizeHasChanged()
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::setReceiverInfo( const std::map<MPIRank,MPISize> & ranksToRecvFrom )
{
   WALBERLA_ASSERT( ! communicationRunning_ );

   recvInfos_.clear();
   for ( auto it = ranksToRecvFrom.begin(); it != ranksToRecvFrom.end(); ++it )
   {
      const MPIRank sender       = it->first;
      const MPISize senderSize   = it->second;
      WALBERLA_ASSERT_GREATER( senderSize, 0 );
      recvInfos_[ sender ].size   = senderSize;
   }

   sizeChangesEverytime_ = false;
   setCommunicationType( true );
}



//**********************************************************************************************************************
/*! Sets receiver information, using SendBuffers (symmetric communication)
*
* Gives the GenericBufferSystem the information that messages are received from the same processes that we
* send to (i.e. from all ranks where SendBuffers were already filled )
* sendBuffer() has to be called before, and corresponding SendBuffers have to be filled.
*
*
* \param useSizeFromSendBuffers  If true, the sizes are expected to be known and equal to the size
*                                of the SendBuffers. SendBuffers with zero size are ignored.
*                                If false, all SendBuffers (also with zero size) are registered as ranks where
*                                messages are received from. The size is unknown, and communicated before.
*
* \param changingSize            if true the size is communicated before every communication step.
*                                if useSizeFromSendBuffer==true and changingSize==true, the size is not
*                                communicated in the first step but in all following steps.
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::setReceiverInfoFromSendBufferState( bool useSizeFromSendBuffers, bool changingSize )
{
   WALBERLA_ASSERT( ! communicationRunning_ );

   recvInfos_.clear();
   for ( auto it = sendInfos_.begin(); it != sendInfos_.end(); ++it )
   {
      const MPIRank sender = it->first;
      const Sb & buffer = it->second.buffer;

      if ( buffer.size() == 0 && useSizeFromSendBuffers )
         continue;

      recvInfos_[ sender ].size  = useSizeFromSendBuffers ? int_c( buffer.size() )  : INVALID_SIZE;
   }

   sizeChangesEverytime_ = changingSize;

   setCommunicationType( useSizeFromSendBuffers );
}



//**********************************************************************************************************************
/*! Notifies that GenericBufferSystem that message sizes have changed ( and optionally are changing in all following steps)
*
* Useful when setReceiverInfo was set such that GenericBufferSystem assumes constant message sizes for all steps.
* Can only be called if no communication is currently running.
*
* \param alwaysChangingSize  if true the message sizes is communicated in all following steps, if false
*                            only in the next step.
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::sizeHasChanged( bool alwaysChangingSize )
{
   WALBERLA_ASSERT( ! communicationRunning_ );

   sizeChangesEverytime_ = alwaysChangingSize;
   setCommunicationType( false );
}



//======================================================================================================================
//
//  Step 1: Schedule Receives and ISends
//
//======================================================================================================================



//**********************************************************************************************************************
/*! Returns an existing SendBuffer, or creates a new one (only if !isCommunicationRunning() )
*
* \param rank  the rank where the buffer should be sent to
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
Sb & GenericBufferSystem<Rb, Sb>::sendBuffer( MPIRank rank )
{
   return sendInfos_[rank].buffer;
}





//**********************************************************************************************************************
/*! Sends filled ( nonzero length) SendBuffers to their destination ranks
*
* If some of the SendBuffers have been sent manually before using send(int rank) they are skipped,
* only the remaining buffers are sent.
*
* If communication was not started before, it is started with this function.
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::sendAll()
{
   WALBERLA_ASSERT_NOT_NULLPTR( currentComm_ ); // call setReceiverInfo first!

   if ( !communicationRunning_ )
      startCommunication();

   for( auto iter = sendInfos_.begin(); iter != sendInfos_.end(); ++iter )
   {
      if ( ! iter->second.alreadySent )
      {
         if ( iter->second.buffer.size() > 0 )
         {
            bytesSent_     += iter->second.buffer.size() * sizeof(SendBuffer::ElementType);
            numberOfSends_ += 1;
            currentComm_->send( iter->first, iter->second.buffer );
         }

         iter->second.alreadySent = true;
      }
   }
}


//**********************************************************************************************************************
/*! Sends a single SendBuffer to its destination rank.
*
* If SendBuffer is empty no message is sent.
*
* If communication was not started before, it is started with this function.
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::send( MPIRank rank )
{
   WALBERLA_ASSERT_NOT_NULLPTR( currentComm_ ); // call setReceiverInfo first!

   if ( !communicationRunning_ )
      startCommunication();

   auto iter = sendInfos_.find(rank);
   WALBERLA_ASSERT( iter != sendInfos_.end() );   // no send buffer was created for this rank
   WALBERLA_ASSERT( ! iter->second.alreadySent ); // this buffer has already been sent

   if ( iter->second.buffer.size() > 0 )
   {
      bytesSent_     += iter->second.buffer.size() * sizeof(SendBuffer::ElementType);
      numberOfSends_ += 1;
      currentComm_->send( rank, iter->second.buffer );
   }

   iter->second.alreadySent = true;
}



//======================================================================================================================
//
//  Private Helper Functions
//
//======================================================================================================================


//**********************************************************************************************************************
/*! Starts communication
*
* - schedules receives and reserves space for MPI_Request vectors in the currentComm_ member
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::startCommunication()
{
   // Clear receive buffers
   for( auto iter = recvInfos_.begin(); iter != recvInfos_.end(); ++iter )  {
      iter->second.buffer.clear();
   }

   const auto tag = currentComm_->getTag();
   WALBERLA_CHECK_EQUAL(activeTags_.find(tag), activeTags_.end(),
                        "Another communication with the same MPI tag is currently in progress.");
   activeTags_.insert(tag);

   WALBERLA_CHECK( ! communicationRunning_ );

   currentComm_->scheduleReceives( recvInfos_ );
   communicationRunning_ = true;

   bytesSent_        = 0;
   bytesReceived_    = 0;

   numberOfSends_    = 0;
   numberOfReceives_ = 0;
}



//**********************************************************************************************************************
/*! Cleanup after communication
*
* - wait for sends to complete
* - clear buffers
* - manage sizeChangesEverytime
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::endCommunication()
{
   WALBERLA_CHECK( communicationRunning_ );
   currentComm_->waitForSends();

   // Clear send buffers
   for( auto iter = sendInfos_.begin(); iter != sendInfos_.end(); ++iter )
   {
      iter->second.alreadySent = false;
      iter->second.buffer.clear();
   }


   if( !sizeChangesEverytime_ )
      setCommunicationType( true );

   communicationRunning_ = false;

   activeTags_.erase( activeTags_.find( currentComm_->getTag() ) );
}



//**********************************************************************************************************************
/*! Helper function for iterator
*
*  See documentation of AbstractCommunication::waitForNext()
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
Rb * GenericBufferSystem<Rb, Sb>::waitForNext( MPIRank & fromRank )
{
   WALBERLA_ASSERT( communicationRunning_ );

   fromRank = currentComm_->waitForNextReceive( recvInfos_ );

   if( fromRank >= 0 )
      return & ( recvInfos_[fromRank].buffer );
   else
   {
      endCommunication();
      return nullptr;
   }

}



//**********************************************************************************************************************
/*! Sets the communication type to known size, unknown size or NoMPI comm
*/
//**********************************************************************************************************************
template< typename Rb, typename Sb>
void GenericBufferSystem<Rb, Sb>::setCommunicationType( const bool knownSize )
{
   WALBERLA_NON_MPI_SECTION()
   {
      currentComm_ = &noMPIComm_;
   }

   WALBERLA_MPI_SECTION()
   {
      if( knownSize )
         currentComm_ = &knownSizeComm_;
      else if ( useIProbe_ )
         currentComm_ = &unknownSizeCommIProbe_;
      else
         currentComm_ = &unknownSizeComm_;
   }
}


//======================================================================================================================
//
//  Rank Ranges
//
//======================================================================================================================

template< typename Rb, typename Sb>
typename GenericBufferSystem<Rb, Sb>::RankRange GenericBufferSystem<Rb,Sb>::noRanks()
{
   return RankRange();
}
template< typename Rb, typename Sb>
typename GenericBufferSystem<Rb, Sb>::RankRange GenericBufferSystem<Rb,Sb>::allRanks()
{
   int numProcesses = MPIManager::instance()->numProcesses();
   RankRange r;
   std::generate_n( std::inserter(r, r.end()), numProcesses, [&]{ return MPIRank(r.size()); } );
   return r;
}
template< typename Rb, typename Sb>
typename GenericBufferSystem<Rb, Sb>::RankRange GenericBufferSystem<Rb,Sb>::allRanksButRoot()
{
   int numProcesses = MPIManager::instance()->numProcesses();
   RankRange r;
   std::generate_n( std::inserter(r, r.end()), numProcesses-1, [&]{ return MPIRank(r.size()+size_t(1)); } );
   return r;
}
template< typename Rb, typename Sb>
typename GenericBufferSystem<Rb, Sb>::RankRange GenericBufferSystem<Rb,Sb>::onlyRank( int includedRank )
{
   WALBERLA_ASSERT_LESS( includedRank, MPIManager::instance()->numProcesses() );
   return RankRange { includedRank };
}

template< typename Rb, typename Sb>
typename GenericBufferSystem<Rb, Sb>::RankRange GenericBufferSystem<Rb,Sb>::onlyRoot()
{
   return RankRange { 0 };
}



} // namespace mpi
} // namespace walberla







