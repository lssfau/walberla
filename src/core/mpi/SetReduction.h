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
//! \file SetReduction.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include <vector>
#include <algorithm>

#include "core/debug/Debug.h"

#include "core/math/Uint.h"

#include "core/mpi/MPIWrapper.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/BufferDataTypeExtensions.h"

namespace walberla {
namespace mpi {

enum SetOperation { INTERSECTION, UNION };

//======================================================================================================================
/*!
*  \brief Reduces a set of values on all processes without using global mpi communication
*
*  The algorithm performs log(n) communication steps, where n is the number of processes in mpiCommunicator. The returned
*  vector is a sorted set of unique values.
*
*  \param values           The local input values. Duplicates will internally be removed. Values have to be buffer packable
*                          and sortable.
*  \param op               The operation to be performed on the set: intersection or union.
*  \param mpiCommunicator  MPI communicator used for the reduction.
*  \param mpiTag           MPI tag used for the reduction.
*/
//======================================================================================================================
template< typename T >
std::vector<T> allReduceSet( std::vector<T> values, SetOperation op, MPI_Comm mpiCommunicator = MPI_COMM_WORLD, int mpiTag = 0 )
{
   // create a sorted set from the vector
   std::sort( values.begin(), values.end() );
   values.erase( std::unique( values.begin(), values.end() ), values.end() );

   WALBERLA_NON_MPI_SECTION()
   {
      return values;
   }

   int rank = 0;
   MPI_Comm_rank( mpiCommunicator, &rank );
   int numberOfRanks = 0;
   MPI_Comm_size( mpiCommunicator, &numberOfRanks );

   int virtualNumberOfRanks = math::uintIsPowerOfTwo( uint_c( numberOfRanks ) ) ? numberOfRanks : int_c( math::uintPow2( math::uintMSBPosition( uint_c( numberOfRanks ) ) ) );
   int virtualRank     = rank + numberOfRanks;
   bool hasVirtualRank = virtualRank < virtualNumberOfRanks;

   WALBERLA_ASSERT_GREATER_EQUAL( virtualNumberOfRanks, numberOfRanks );
   WALBERLA_ASSERT_LESS( virtualNumberOfRanks / 2, numberOfRanks );
   
   std::vector<T> tmp0, tmp1;
   int lvl = 0;
   int lvlRank = rank >> lvl;
   int virtLvlRank = virtualRank >> lvl;
   while( ( virtualNumberOfRanks - 1 ) >> lvl > 0 )
   {
      int commPartner        = ( ( lvlRank % 2 )     == 0 ) ? ( rank        + ( 1 << lvl ) ) : (        rank - ( 1 << lvl ) );
      int virtualCommPartner = ( ( virtLvlRank % 2 ) == 0 ) ? ( virtualRank + ( 1 << lvl ) ) : ( virtualRank - ( 1 << lvl ) );

      SendBuffer sendBuffer;
      sendBuffer << values;

      int sendBufferSize = int_c( sendBuffer.size() );
      MPI_Request sendRequests[] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL };
      MPI_Isend( &sendBufferSize, 1, MPI_INT, commPartner % numberOfRanks, mpiTag, mpiCommunicator, &sendRequests[0] );
      MPI_Isend( sendBuffer.ptr(), sendBufferSize, MPI_CHAR, commPartner % numberOfRanks, mpiTag, mpiCommunicator, &sendRequests[1] );
      if( hasVirtualRank )
      {
         MPI_Isend( &sendBufferSize, 1, MPI_INT, virtualCommPartner % numberOfRanks, mpiTag, mpiCommunicator, &sendRequests[2] );
         MPI_Isend( sendBuffer.ptr(), sendBufferSize, MPI_CHAR, virtualCommPartner % numberOfRanks, mpiTag, mpiCommunicator, &sendRequests[3] );
      }

      int recvBufferSize = 0;
      MPI_Recv( &recvBufferSize, 1, MPI_INT, commPartner % numberOfRanks, mpiTag, mpiCommunicator, MPI_STATUS_IGNORE );
      RecvBuffer recvBuffer, virtualRecvBuffer;
      recvBuffer.resize( numeric_cast<size_t>( recvBufferSize ) );
      MPI_Recv( recvBuffer.ptr(), recvBufferSize, MPI_CHAR, commPartner % numberOfRanks, mpiTag, mpiCommunicator, MPI_STATUS_IGNORE );
      recvBuffer >> tmp0;
      tmp1.clear();
      switch( op )
      {
      case INTERSECTION:
         std::set_intersection( values.begin(), values.end(), tmp0.begin(), tmp0.end(), std::back_inserter( tmp1 ) );
         break;
      case UNION:
         std::set_union( values.begin(), values.end(), tmp0.begin(), tmp0.end(), std::back_inserter( tmp1 ) );
         break;
      }
      swap( values, tmp1 );

      if( hasVirtualRank )
      {
         int virtualRecvBufferSize = 0;
         MPI_Recv( &virtualRecvBufferSize, 1, MPI_INT, virtualCommPartner % numberOfRanks, mpiTag, mpiCommunicator, MPI_STATUS_IGNORE );
         virtualRecvBuffer.resize( numeric_cast<size_t>( virtualRecvBufferSize ) );
         MPI_Recv( virtualRecvBuffer.ptr(), virtualRecvBufferSize, MPI_CHAR, virtualCommPartner % numberOfRanks, mpiTag, mpiCommunicator, MPI_STATUS_IGNORE );
         virtualRecvBuffer >> tmp0;
         tmp1.clear();
         switch( op )
         {
         case INTERSECTION:
            std::set_intersection( values.begin(), values.end(), tmp0.begin(), tmp0.end(), std::back_inserter( tmp1 ) );
            break;
         case UNION:
            std::set_union( values.begin(), values.end(), tmp0.begin(), tmp0.end(), std::back_inserter( tmp1 ) );
            break;
         }
         swap( values, tmp1 );
      }

      MPI_Waitall( hasVirtualRank ? 4 : 2, sendRequests, MPI_STATUSES_IGNORE );

      ++lvl;
      lvlRank = rank >> lvl;
      virtLvlRank = virtualRank >> lvl;
   }

   return values;
}

} // namespace mpi
} // namespace walberla
