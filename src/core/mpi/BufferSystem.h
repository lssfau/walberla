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
//! \file BufferSystem.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Header file for BufferSystem
//
//======================================================================================================================

#pragma once

#include "BufferSystemHelper.h"
#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include <map>
#include <set>
#include <vector>


namespace walberla {
namespace mpi  {

template< typename RecvBuffer_T, typename SendBuffer_T>
class GenericOpenMPBufferSystem;

//**********************************************************************************************************************
/*! Manages MPI Communication with a set of known communication partners.
*
*   Features / Restrictions:
*    - usable in every case where normal MPI_Send, MPI_Recv is needed
*    - communication partners have to be known, message size not necessarily
*    - unknown message sizes possible ( -> automatic extra message to exchange sizes )
*    - Implemented with non-blocking MPI calls, and MPI_Waitany to process a message as soon
*      as it was received while still waiting for other messages
*
*
* \ingroup mpi
*
*
*  Example:
*  --------
*
   \code
      BufferSystem bs ( MPI_COMM_WORLD );
      // we assume every process has exactly two neighbors ( periodic )
      // and the neighbor ranks are stored in 'neighborRank1' and 'neighborRank2'
      bs.sendBuffer( neighborRank1 ) << 42;
      bs.sendBuffer( neighborRank2 ) << 4242;

      // We expect to receive the same amount as we send:
      bs.setReceiverInfoFromSendBufferState( true, false );

      bs.sendAll();

      for ( auto i = bs.begin(); i != bs.end(); ++i ) {
         cout << "Message received from " << i.rank() << endl;
         int messageContent;
         i.buffer() >> messageContent;
         cout << "Received " << messageContent << endl;
      }
   \endcode
*
*
* Usage:
* ------
*
*    1. Setup:
*        - define receive information using setReceiverInfo() or setReceiverInfoFromSendBufferState()
*          With these functions one defines the processes that we receive from and optionally the size
*          of the received messages. If the message sizes are unknown, they have to be communicated first.
*          One also defines if the sizes stay constant or if they change in each communication step
*           (size message is then sent before every content message)
*        - The receiver and send information can be changed, if no communication is currently running.
*
*    2. Communication Step:
*        - Optionally call scheduleReceives() -> starts communication step and causes MPI_IRecv's to be called.
*          This is also automatically called on first send operation.
*        - fill send buffers
*        - call send() for each buffer after filling the buffer, or call sendAll() after filling all buffers.
*          The send*() functions return immediately, internally MPI_ISend's are called.
*          The send*() functions start the communication step if it was not already started by scheduleReceives()
*        - Receiving: iterate over incoming messages using begin() and end(). Internally a MPI_Waitany is executed that
*          returns as soon as a single message was received. Then this message can be processed while waiting
*          for the other messages.
*          \attention Even if the current process does not receive anything, call begin() and end()
*                     to finish the communication step
*        - When iteration has reached the end the communication step is finished.
*        - when communication has finished all Send- and RecvBuffers are automatically cleared
*
* When running multiple BufferSystems concurrently different MPI tags have to be used
* for the systems: the tag can be passed in the constructor.
*
*/
//**********************************************************************************************************************
template< typename RecvBuffer_T = RecvBuffer, typename SendBuffer_T = SendBuffer>
class GenericBufferSystem
{
public:
   class iterator;

   //**Construction and Destruction*************************************************************************************
   /*!\name Constructors */
   //@{
   explicit GenericBufferSystem( const MPI_Comm & communicator, int tag = 0 );
   GenericBufferSystem( const GenericBufferSystem & other );
   GenericBufferSystem & operator=( const GenericBufferSystem & other );
   ~GenericBufferSystem() = default;
   //@}
   //*******************************************************************************************************************


   //** Receiver Registration        ***********************************************************************************
   /*! \name Receiver Registration  */
   //@{
   template<typename RankIter> void setReceiverInfo( RankIter begin, RankIter end,              bool changingSize );
   template<typename Range>    void setReceiverInfo( const Range & range,                       bool changingSize );
              void setReceiverInfo( const std::set<MPIRank> & ranksToRecvFrom, bool changingSize );

   void setReceiverInfo( const std::map<MPIRank,MPISize> & ranksToRecvFrom );
   void setReceiverInfoFromSendBufferState( bool useSizeFromSendBuffers, bool changingSize );

   void sizeHasChanged( bool alwaysChangingSize = false );
   //@}
   //*******************************************************************************************************************


   //** Executing Communication Step   *********************************************************************************
   /*! \name Executing Communication Step */
   //@{
   void scheduleReceives() { startCommunication(); }

   SendBuffer_T & sendBuffer ( MPIRank rank );
   SendBuffer_T & sendBuffer ( uint_t  rank ) { return sendBuffer( int_c( rank ) ); }
   inline size_t size() const;


   void sendAll();
   void send( MPIRank rank );

   iterator begin() { WALBERLA_ASSERT( communicationRunning_); return iterator( *this, true ); }
   iterator end()   {                                          return iterator( *this, false); }
   //@}
   //*******************************************************************************************************************


   //** Iterator        ************************************************************************************************
   /*! \name Iterator  */
   //@{
   class iterator
   {
   public:
      MPIRank        rank()   { return currentSenderRank_;  }
      RecvBuffer_T & buffer() { return *currentRecvBuffer_; }

      void operator++();
      bool operator==( const iterator & other );
      bool operator!=( const iterator & other );

   private:
      iterator( GenericBufferSystem & bufferSystem, bool begin );

      GenericBufferSystem & bufferSystem_;

      RecvBuffer_T * currentRecvBuffer_;
      MPIRank        currentSenderRank_;

      friend class GenericBufferSystem;
   };
   friend class iterator;
   //@}
   //*******************************************************************************************************************


   //** Status Queries        ******************************************************************************************
   /*! \name Status Queries  */
   //@{
   bool isSizeCommunicatedInNextStep() const { return (currentComm_ == &unknownSizeComm_); }
   bool isCommunicationRunning() const       { return communicationRunning_;               }
   bool isReceiverInformationSet() const     { return currentComm_ != NULL;                }
   //@}
   //*******************************************************************************************************************

   void useIProbe(const bool use) { useIProbe_ = use; }
   bool isIProbedUsed() const { return useIProbe_; }

   ///Bytes sent during the current or last communication
   int64_t getBytesSent() const { return bytesSent_; }
   ///Bytes received during the current or last communication
   int64_t getBytesReceived() const { return bytesReceived_; }

   ///Communication partners during current or last send operation
   int64_t getNumberOfSends() const { return numberOfSends_; }
   ///Communication partners during current or last receive operation
   int64_t getNumberOfReceives() const { return numberOfReceives_; }


   //* Rank Ranges     *************************************************************************************************
   /*! \name Rank Ranges  */
   //@{
   using RankRange = std::set<MPIRank>;
   static RankRange noRanks();
   static RankRange allRanks();
   static RankRange allRanksButRoot();
   static RankRange onlyRoot();
   static RankRange onlyRank( MPIRank includedRank );
   //@}
   //*******************************************************************************************************************

protected:

   friend class GenericOpenMPBufferSystem<RecvBuffer_T, SendBuffer_T>;

   void startCommunication();
   void endCommunication();
   RecvBuffer_T * waitForNext( MPIRank & fromRank );
   void setCommunicationType( const bool knownSize );

   internal::KnownSizeCommunication<RecvBuffer_T, SendBuffer_T>         knownSizeComm_;
   internal::UnknownSizeCommunication<RecvBuffer_T, SendBuffer_T>       unknownSizeComm_;
   internal::UnknownSizeCommunicationIProbe<RecvBuffer_T, SendBuffer_T> unknownSizeCommIProbe_;
   internal::NoMPICommunication<RecvBuffer_T, SendBuffer_T>             noMPIComm_;
   internal::AbstractCommunication<RecvBuffer_T, SendBuffer_T> *        currentComm_;  //< after receiver setup, this points to unknown- or knownSizeComm_

   bool sizeChangesEverytime_; //< if set to true, the receiveSizeUnknown_ is set to true before communicating
   bool communicationRunning_; //< indicates if a communication step is currently running


   /// Info about the message to be received from a certain rank:
   /// information holds the buffer and, if known, the message size
   std::map<MPIRank, typename internal::AbstractCommunication<RecvBuffer_T, SendBuffer_T>::ReceiveInfo> recvInfos_;


   struct SendInfo {
      SendInfo() : alreadySent(false) {}
      SendBuffer_T buffer;
      bool alreadySent;
   };
   std::map<MPIRank, SendInfo> sendInfos_;

   //stores tags of running communications in debug mode to ensure that
   //each concurrently running communication uses different tags
   static std::set<int> activeTags_;

   bool useIProbe_ = false; ///< switch between IProbe and two message communication for varying size communication

   int64_t bytesSent_        = 0; ///< number of bytes sent during last communication
   int64_t bytesReceived_    = 0; ///< number of bytes received during last communication

   int64_t numberOfSends_    = 0; ///< number of communication partners during last send
   int64_t numberOfReceives_ = 0; ///< number of communication partners during last receive
};

using BufferSystem = GenericBufferSystem<RecvBuffer, SendBuffer>;

} // namespace mpi
} // namespace walberla

#include "BufferSystem.impl.h"

