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
//! \file BufferSystem.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "BufferSystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"
#include "core/debug/CheckFunctions.h"


namespace walberla {
namespace mpi {



std::set<int> BufferSystem::activeTags_;

//======================================================================================================================
//
//  Iterator
//
//======================================================================================================================



BufferSystem::iterator::iterator( BufferSystem & bufferSystem, bool begin )
    : bufferSystem_( bufferSystem), currentRecvBuffer_( nullptr ), currentSenderRank_( -1 )
{
   if ( begin ) // init iterator
      ++(*this);
}


void BufferSystem::iterator::operator++()
{
   currentRecvBuffer_ = bufferSystem_.waitForNext( currentSenderRank_ );
   if ( ! currentRecvBuffer_ ) {
      WALBERLA_ASSERT_EQUAL( currentSenderRank_, -1 );
   }
}

bool BufferSystem::iterator::operator==( const BufferSystem::iterator & other )
{
   // only equality checks with end iterators are allowed
   WALBERLA_ASSERT( other.currentSenderRank_ == -1 || currentSenderRank_ == -1 );

   return ( currentSenderRank_ == other.currentSenderRank_ );
}

bool BufferSystem::iterator::operator!=( const BufferSystem::iterator & other )
{
   // only equality checks with end iterators are allowed
   WALBERLA_ASSERT( other.currentSenderRank_ == -1 || currentSenderRank_ == -1 );

   return ( currentSenderRank_ != other.currentSenderRank_ );
}




//======================================================================================================================
//
//  Constructors
//
//======================================================================================================================


BufferSystem::BufferSystem( const MPI_Comm & communicator, int tag )
   : knownSizeComm_  ( communicator, tag ),
     unknownSizeComm_( communicator, tag ),
     noMPIComm_( communicator, tag ),
     currentComm_    ( nullptr ),
     sizeChangesEverytime_( true ),
     communicationRunning_( false )
{
}


BufferSystem::BufferSystem( const BufferSystem &other )
   : knownSizeComm_  ( other.knownSizeComm_.getCommunicator(), other.knownSizeComm_.getTag() ),
     unknownSizeComm_( other.knownSizeComm_.getCommunicator(), other.knownSizeComm_.getTag() ),
     noMPIComm_      ( other.knownSizeComm_.getCommunicator(), other.knownSizeComm_.getTag() ),
     currentComm_ ( nullptr ),
     sizeChangesEverytime_( other.sizeChangesEverytime_ ),
     communicationRunning_( other.communicationRunning_ ),
     recvInfos_( other.recvInfos_ ),
     sendInfos_( other.sendInfos_ )
{
   WALBERLA_ASSERT( !communicationRunning_, "Can't copy BufferSystem while communication is running" );
   if( other.currentComm_ == &other.knownSizeComm_ )
      currentComm_ = &knownSizeComm_;
   else if ( other.currentComm_ == &other.unknownSizeComm_ )
      currentComm_ = &unknownSizeComm_;
   else if ( other.currentComm_ == &other.noMPIComm_ )
      currentComm_ = &noMPIComm_;
   else
      currentComm_ = nullptr; // receiver information not yet set
}


BufferSystem & BufferSystem::operator=( const BufferSystem & other )
{
   WALBERLA_ASSERT( !communicationRunning_, "Can't copy BufferSystem while communication is running" );

   sizeChangesEverytime_ = other.sizeChangesEverytime_;
   communicationRunning_ = other.communicationRunning_;
   recvInfos_ = other.recvInfos_;
   sendInfos_ = other.sendInfos_;

   if( other.currentComm_ == &other.knownSizeComm_ )
      currentComm_ = &knownSizeComm_;
   else if ( other.currentComm_ == &other.unknownSizeComm_ )
      currentComm_ = &unknownSizeComm_;
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
void BufferSystem::setReceiverInfo( const std::set<MPIRank> & ranksToRecvFrom, bool changingSize )
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
void BufferSystem::setReceiverInfo( const std::map<MPIRank,MPISize> & ranksToRecvFrom )
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
* Gives the BufferSystem the information that messages are received from the same processes that we
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
void BufferSystem::setReceiverInfoFromSendBufferState( bool useSizeFromSendBuffers, bool changingSize )
{
   WALBERLA_ASSERT( ! communicationRunning_ );

   recvInfos_.clear();
   for ( auto it = sendInfos_.begin(); it != sendInfos_.end(); ++it )
   {
      const MPIRank sender = it->first;
      const SendBuffer & buffer = it->second.buffer;

      if ( buffer.size() == 0 && useSizeFromSendBuffers )
         continue;

      recvInfos_[ sender ].size  = useSizeFromSendBuffers ? int_c( buffer.size() )  : INVALID_SIZE;
   }

   sizeChangesEverytime_ = changingSize;

   setCommunicationType( useSizeFromSendBuffers );
}



//**********************************************************************************************************************
/*! Notifies that BufferSystem that message sizes have changed ( and optionally are changing in all following steps)
*
* Useful when setReceiverInfo was set such that BufferSystem assumes constant message sizes for all steps.
* Can only be called if no communication is currently running.
*
* \param alwaysChangingSize  if true the message sizes is communicated in all following steps, if false
*                            only in the next step.
*/
//**********************************************************************************************************************
void BufferSystem::sizeHasChanged( bool alwaysChangingSize )
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
SendBuffer & BufferSystem::sendBuffer( MPIRank rank )
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
void BufferSystem::sendAll()
{
   WALBERLA_ASSERT_NOT_NULLPTR( currentComm_ ); // call setReceiverInfo first!

   if ( !communicationRunning_ )
      startCommunication();

   for( auto iter = sendInfos_.begin(); iter != sendInfos_.end(); ++iter )
   {
      if ( ! iter->second.alreadySent )
      {
         if ( iter->second.buffer.size() > 0 )
            currentComm_->send( iter->first, iter->second.buffer );

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
void BufferSystem::send( MPIRank rank )
{
   WALBERLA_ASSERT_NOT_NULLPTR( currentComm_ ); // call setReceiverInfo first!

   if ( !communicationRunning_ )
      startCommunication();

   auto iter = sendInfos_.find(rank);
   WALBERLA_ASSERT( iter != sendInfos_.end() );   // no send buffer was created for this rank
   WALBERLA_ASSERT( ! iter->second.alreadySent ); // this buffer has already been sent

   if ( iter->second.buffer.size() > 0 )
      currentComm_->send( rank, iter->second.buffer );

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
void BufferSystem::startCommunication()
{
   const auto tag = currentComm_->getTag();
   WALBERLA_CHECK_EQUAL(activeTags_.find(tag), activeTags_.end(),
                        "Another communication with the same MPI tag is currently in progress.");
   activeTags_.insert(tag);

   WALBERLA_CHECK( ! communicationRunning_ );

   currentComm_->scheduleReceives( recvInfos_ );
   communicationRunning_ = true;
}



//**********************************************************************************************************************
/*! Cleanup after communication
*
* - wait for sends to complete
* - clear buffers
* - manage sizeChangesEverytime
*/
//**********************************************************************************************************************
void BufferSystem::endCommunication()
{
   WALBERLA_CHECK( communicationRunning_ );
   currentComm_->waitForSends();

   // Clear send buffers
   for( auto iter = sendInfos_.begin(); iter != sendInfos_.end(); ++iter )
   {
      iter->second.alreadySent = false;
      iter->second.buffer.clear();
   }

   // Clear receive buffers
   for( auto iter = recvInfos_.begin(); iter != recvInfos_.end(); ++iter )  {
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
RecvBuffer * BufferSystem::waitForNext( MPIRank & fromRank )
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
void BufferSystem::setCommunicationType( const bool knownSize )
{
   WALBERLA_NON_MPI_SECTION()
   {
      currentComm_ = &noMPIComm_;
   }

   WALBERLA_MPI_SECTION()
   {
      if( knownSize )
         currentComm_ = &knownSizeComm_;
      else
         currentComm_ = &unknownSizeComm_;
   }
}


//======================================================================================================================
//
//  Rank Ranges
//
//======================================================================================================================

// using boost::counting_range didn't work on all supported compilers
// so the range is created explicitly

BufferSystem::RankRange BufferSystem::noRanks()
{
   return RankRange ( RankCountIter( 0 ),
                      RankCountIter( 0 ) );
}

BufferSystem::RankRange BufferSystem::allRanks()
{
   int numProcesses = MPIManager::instance()->numProcesses();
   return RankRange ( RankCountIter( 0            ),
                      RankCountIter( numProcesses ) );
}

BufferSystem::RankRange BufferSystem::allRanksButRoot()
{
   int numProcesses = MPIManager::instance()->numProcesses();
   return RankRange ( RankCountIter( 1            ),
                      RankCountIter( numProcesses ) );
}

BufferSystem::RankRange BufferSystem::onlyRank( int includedRank )
{
   WALBERLA_ASSERT_LESS( includedRank, MPIManager::instance()->numProcesses() );
   return RankRange ( RankCountIter( includedRank   ),
                      RankCountIter( includedRank+1 ) );
}


BufferSystem::RankRange BufferSystem::onlyRoot()
{
   return RankRange ( RankCountIter( 0 ),
                      RankCountIter( 1 ) );
}



} // namespace mpi
} // namespace walberla







