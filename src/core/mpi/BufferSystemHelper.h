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
//! \file BufferSystemHelper.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "MPIWrapper.h"
#include "RecvBuffer.h"
#include "SendBuffer.h"
#include "core/DataTypes.h"

#include <list>
#include <map>
#include <vector>


namespace walberla {
namespace mpi {
namespace internal {

   template< typename RecvBuffer_T, typename SendBuffer_T>
   class AbstractCommunication
   {
   public:
      AbstractCommunication( const MPI_Comm & communicator, int tag = 0 )
         : communicator_( communicator), tag_(tag)
      {}

      virtual ~AbstractCommunication() {}


      struct ReceiveInfo {
         RecvBuffer_T buffer;
         MPISize size;
      };


      /*************************************************************************************************************//**
      * Send buffer content to receiver using MPI_ISend, request is stored internally -> see waitForSends()
      *****************************************************************************************************************/
      virtual void send( MPIRank rank, const SendBuffer_T & sendBuffer ) = 0;


      /*************************************************************************************************************//**
      * Wait for all sends to complete.
      *****************************************************************************************************************/
      virtual void waitForSends() = 0;




      /*************************************************************************************************************//**
      * Executes MPI_Irecv calls, with the recvInfos as target memory.
      *
      * \param recvInfos  Map that has entries for all ranks where messages are expected.
      *                   The KnownSizeCommunication also expects the size field to be valid ( defining the expected
      *                   message size)
      *                   Do not change/read the recvInfos after scheduleReceive() call
      *                   and before waitForNextReceive() is completed
      *
      *****************************************************************************************************************/
      virtual void    scheduleReceives  ( std::map<MPIRank, ReceiveInfo> & recvInfos ) = 0;



      /*************************************************************************************************************//**
      * Waits for the next receive to complete and returns.
      *
      * \param recvInfo  the same receive info that was passed to scheduleReceives()
      *
      * \return The rank where the data was received -> recvInfos[rank] is now valid
      *         INVALID_RANK if all messages were received.
      *****************************************************************************************************************/
      virtual MPIRank waitForNextReceive( std::map<MPIRank, ReceiveInfo> & recvInfos ) = 0;


      virtual int getTag() const { return tag_; }

      virtual MPI_Comm getCommunicator() const { return communicator_; }

   protected:
      MPI_Comm communicator_;
      int tag_;
   };


   template< typename RecvBuffer_T, typename SendBuffer_T>
   class KnownSizeCommunication : public AbstractCommunication<RecvBuffer_T, SendBuffer_T>
   {
   public:
      using typename AbstractCommunication<RecvBuffer_T, SendBuffer_T>::ReceiveInfo;

      KnownSizeCommunication( const MPI_Comm & communicator, int tag = 0 )
           : AbstractCommunication<RecvBuffer_T, SendBuffer_T>( communicator, tag ), sending_(false), receiving_(false) {}

      virtual ~KnownSizeCommunication() {}

      virtual void send( MPIRank receiver, const SendBuffer_T & sendBuffer );
      virtual void waitForSends();

      virtual void    scheduleReceives  ( std::map<MPIRank, ReceiveInfo> & recvInfos );

      /// size field of recvInfos is expected to be valid
      virtual MPIRank waitForNextReceive( std::map<MPIRank, ReceiveInfo> & recvInfos );

   private:
      bool sending_;
      bool receiving_;

      std::vector<MPI_Request> sendRequests_;
      std::vector<MPI_Request> recvRequests_;
   };


   template< typename RecvBuffer_T, typename SendBuffer_T>
   class UnknownSizeCommunication : public AbstractCommunication<RecvBuffer_T, SendBuffer_T>
   {
   public:
      using typename AbstractCommunication<RecvBuffer_T, SendBuffer_T>::ReceiveInfo;

      UnknownSizeCommunication( const MPI_Comm & communicator, int tag = 0 )
           :  AbstractCommunication<RecvBuffer_T, SendBuffer_T>( communicator, tag ), sending_(false), receiving_(false) {}

      virtual ~UnknownSizeCommunication() {}

      virtual void send( MPIRank receiver, const SendBuffer_T & sendBuffer );
      virtual void waitForSends();

      virtual void scheduleReceives( std::map<MPIRank, ReceiveInfo> & recvInfos );

      /// size field of recvInfos can be invalid, is filled in with the actual message size
      virtual MPIRank waitForNextReceive( std::map<MPIRank, ReceiveInfo> & recvInfos );

   private:
      bool sending_;
      bool receiving_;

      std::vector<MPI_Request> sendRequests_;
      std::list<MPISize>       outgoingBufferForSizes_;


      std::vector<MPI_Request> recvRequests_;
      std::vector<bool>        sizeAlreadyReceived_;
   };


   template< typename RecvBuffer_T, typename SendBuffer_T>
   class UnknownSizeCommunicationIProbe : public AbstractCommunication<RecvBuffer_T, SendBuffer_T>
   {
   public:
      using typename AbstractCommunication<RecvBuffer_T, SendBuffer_T>::ReceiveInfo;

      UnknownSizeCommunicationIProbe( const MPI_Comm & communicator, int tag = 0 )
           :  AbstractCommunication<RecvBuffer_T, SendBuffer_T>( communicator, tag ), sending_(false), receiving_(false) {}

      virtual ~UnknownSizeCommunicationIProbe() {}

      virtual void send( MPIRank receiver, const SendBuffer_T & sendBuffer );
      virtual void waitForSends();

      virtual void scheduleReceives( std::map<MPIRank, ReceiveInfo> & recvInfos );

      /// size field of recvInfos can be invalid, is filled in with the actual message size
      virtual MPIRank waitForNextReceive( std::map<MPIRank, ReceiveInfo> & recvInfos );

   private:
      bool sending_;
      bool receiving_;
      int  pendingReceives_;

      std::vector<MPI_Request> sendRequests_;
   };


   template< typename RecvBuffer_T, typename SendBuffer_T>
   class NoMPICommunication : public AbstractCommunication<RecvBuffer_T, SendBuffer_T>
   {
   public:
      using typename AbstractCommunication<RecvBuffer_T, SendBuffer_T>::ReceiveInfo;

      NoMPICommunication( const MPI_Comm & communicator, int tag = 0 )
         : AbstractCommunication<RecvBuffer_T, SendBuffer_T>( communicator, tag ), received_( false ) {}

      virtual ~NoMPICommunication() {}

      virtual void send( MPIRank receiver, const SendBuffer_T & sendBuffer );
      virtual void waitForSends();

      virtual void scheduleReceives( std::map<MPIRank, ReceiveInfo> & recvInfos );

      /// size field of recvInfos can be invalid, is filled in with the actual message size
      virtual MPIRank waitForNextReceive( std::map<MPIRank, ReceiveInfo> & recvInfos );

   private:
      bool         received_;
      RecvBuffer_T tmpBuffer_;
   };



} // namespace internal
} // namespace mpi
} // namespace walberla

#include "BufferSystemHelper.impl.h"
