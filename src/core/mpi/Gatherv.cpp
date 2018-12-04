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
//! \file Gatherv.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Gatherv.h"

#include "Broadcast.h"
#include "BufferDataTypeExtensions.h"



namespace walberla {
namespace mpi {



/////////////////////////////////////////////
// gatherv specializations for std::string //
// TODO: improve! Especially allGatherv .. //
/////////////////////////////////////////////

template<>
std::vector< std::string > gatherv( const std::vector< std::string > & values, int recvRank, MPI_Comm comm )
{
   WALBERLA_NON_MPI_SECTION()
   {
      return std::vector< std::string >( values );
   }

   mpi::SendBuffer sb;
   mpi::RecvBuffer rb;
   sb << values;

   gathervBuffer( sb, rb, recvRank, comm );

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   if( myRank == recvRank )
   {
      int numProcesses;
      MPI_Comm_size( comm, &numProcesses );

      std::vector< std::string > result;

      for( int i = 0; i < numProcesses; ++i )
      {
         std::vector< std::string > tmp;
         rb >> tmp;
         for( auto it = tmp.begin(); it != tmp.end(); ++it )
            result.push_back( *it );
      }

      return result;
   }

   return std::vector< std::string >();
}

template<>
std::vector< std::string > allGatherv( const std::vector< std::string > & values, MPI_Comm comm )
{
   auto result = gatherv( values, 0, comm );

   broadcastObject( result, 0, comm );

   return result;
}



void gathervBuffer( const mpi::SendBuffer & sendBuffer, mpi::RecvBuffer & recvBuffer, int targetRank, MPI_Comm comm )
{
   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( targetRank, 0 );
      recvBuffer = sendBuffer;
      return;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   const bool isGatherProcess = ( myRank == targetRank );

   int numProcesses;
   MPI_Comm_size( comm, &numProcesses );
   WALBERLA_ASSERT_GREATER( numProcesses, 0 );

   std::vector<int> sizes;
   if ( isGatherProcess )
      sizes.resize( uint_c(numProcesses) );

   int sendBufferSize = int_c( sendBuffer.size() );

   // Gather the message sizes on root process
   MPI_Gather( &sendBufferSize,                   1, MPITrait<int>::type(),
               isGatherProcess? &sizes[0] : nullptr, 1, MPITrait<int>::type(),
               targetRank, comm );

   int totalSize = 0;
   std::vector<int> displacements ( sizes.size(), 0 );
   for( uint_t i=0; i < sizes.size(); ++i  )
   {
      displacements[i] = totalSize;
      totalSize += sizes[i];
   }

   if ( isGatherProcess )
      recvBuffer.resize( numeric_cast< size_t >( totalSize ) );

   MPI_Gatherv( sendBuffer.ptr(), int_c( sendBuffer.size() ), MPITrait< mpi::SendBuffer::ElementType >::type(),
                recvBuffer.ptr(),
                isGatherProcess? &sizes[0] : nullptr,
                isGatherProcess? &displacements[0] : nullptr,
                MPITrait< mpi::RecvBuffer::ElementType >::type(), targetRank, comm );
}



void allGathervBuffer( const mpi::SendBuffer & sendBuffer, mpi::RecvBuffer & recvBuffer, MPI_Comm comm )
{
   WALBERLA_NON_MPI_SECTION()
   {
      recvBuffer = sendBuffer;
      return;
   }

   int numProcesses;
   MPI_Comm_size( comm, &numProcesses );
   WALBERLA_ASSERT_GREATER( numProcesses, 0 );

   std::vector<int> sizes;
   int sendBufferSize = int_c( sendBuffer.size() );
   sizes.resize( uint_c(numProcesses) );

   MPI_Allgather( &sendBufferSize, 1, MPITrait<int>::type(),
                  &sizes[0],       1, MPITrait<int>::type(), comm );

   int totalSize = 0;
   std::vector<int> displacements( sizes.size(), 0 );
   for( uint_t i = 0; i < sizes.size(); ++i )
   {
      displacements[i] = totalSize;
      totalSize += sizes[i];
   }
   recvBuffer.resize( numeric_cast< size_t >( totalSize ) );

   MPI_Allgatherv( sendBuffer.ptr(), int_c( sendBuffer.size() ), MPITrait< mpi::SendBuffer::ElementType >::type(),
                   recvBuffer.ptr(), &sizes[0], &displacements[0],
                   MPITrait< mpi::RecvBuffer::ElementType >::type(), comm );
}



} // namespace mpi
} // namespace walberla
