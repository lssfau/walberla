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
//! \file Broadcast.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "RecvBuffer.h"
#include "SendBuffer.h"

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPIWrapper.h"

#include <vector>


namespace walberla {
namespace mpi {


//======================================================================================================================
/*!
 *  \brief Broadcasts an arbitrary sized object from one process to all other processes
 *
 *  T has to be able to be packed into SendBuffer and unpacked from RecvBuffer.
 *
 *  \param object     The object to be broadcasted
 *  \param senderRank The rank of the process sending the object
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
template< typename T >
void broadcastObject( T & object, int senderRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( senderRank, 0 );
      return;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );
   int numProcesses;
   MPI_Comm_size( comm, &numProcesses );

   if( numProcesses == 1 )
      return;

   if( myRank == senderRank )
   {
      // Send object to all other processes

      SendBuffer sendBuffer;
      sendBuffer << object;

      size_t size =  sendBuffer.size();

      MPI_Bcast( &size, 1, MPITrait<size_t>::type(), senderRank, comm );
      MPI_Bcast( sendBuffer.ptr(), int_c(size), MPITrait<SendBuffer::ElementType>::type(), senderRank, comm );
   }
   else // => myRank != senderRank
   {
      // Receive object from process with rank 'senderRank'

      size_t size;
      MPI_Bcast( &size, 1, MPITrait<size_t>::type(), senderRank, comm );

      RecvBuffer recvBuffer;
      recvBuffer.resize( size );
      MPI_Bcast( recvBuffer.ptr(), int_c(size), MPITrait<SendBuffer::ElementType>::type(), senderRank, comm );

      recvBuffer >> object;

      WALBERLA_ASSERT( recvBuffer.isEmpty() );
   }
}

} // namespace mpi
} // namespace walberla
