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
//! \file GathervTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Tests for core/mpi/Gatherv.h
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIManager.h"

#include <limits>
#include <vector>


void runTestGatherv( int recvRank )
{
   using namespace walberla;

   const int rank         = walberla::MPIManager::instance()->rank();
   const int numProcesses = walberla::MPIManager::instance()->numProcesses();

   // Static size of values

   std::vector<int> values( 1, rank );

   std::vector<int> result = walberla::mpi::gatherv( values, recvRank );

   if( rank == recvRank )
   {
      WALBERLA_CHECK_EQUAL( result.size(), numeric_cast<size_t>( numProcesses ) );
      for( int i = 0; i < numProcesses; ++i )
        WALBERLA_CHECK_EQUAL( result[ uint_c(i) ], i );
   }
   else
   {
      WALBERLA_CHECK( result.empty() );
   }

   // Dynamic size of values

   values.assign( numeric_cast< std::vector<int>::size_type >(rank), rank );

   result = walberla::mpi::gatherv( values, recvRank );

   if( rank == recvRank )
   {
      const auto expected_size = numProcesses * ( numProcesses + 1 ) / 2 - numProcesses; // Gaussian sum 0 + 1 + ... + n - 1

      WALBERLA_CHECK_EQUAL( result.size(), numeric_cast<size_t>( expected_size ) );
      uint_t pos = 0;
      for( int i = 0; i < numProcesses; ++i )
         for( int j = 0; j < i; ++j )
            WALBERLA_CHECK_EQUAL( result[pos++], i );
   }
   else
   {
      WALBERLA_CHECK( result.empty() );
   }

   // gathervBuffer

   std::vector<int> vSend;
   vSend.push_back( 2 * rank );
   vSend.push_back( 2 * rank + 1 );
   if( rank % 2 == 0 )
      vSend.push_back( 23 );

   walberla::mpi::SendBuffer sBuffer;
   walberla::mpi::RecvBuffer rBuffer;
   sBuffer << vSend;

   walberla::mpi::gathervBuffer( sBuffer, rBuffer, recvRank );

   if( rank == recvRank )
   {
      std::vector<int> vRecv;
      for( int i = 0; i != numProcesses; ++i )
      {
         rBuffer >> vRecv;
         WALBERLA_CHECK_EQUAL( vRecv.size(), ( i % 2 == 0 ) ? uint_t(3) : uint_t(2) );
         WALBERLA_CHECK_EQUAL( vRecv[0], 2 * i );
         WALBERLA_CHECK_EQUAL( vRecv[1], 2 * i + 1 );
         if( i % 2 == 0 )
            WALBERLA_CHECK_EQUAL( vRecv[2], 23 );
      }
   }
   else
   {
      WALBERLA_CHECK( rBuffer.isEmpty() );
   }
}

void runTestAllGatherv()
{
   using namespace walberla;

   const int rank         = walberla::MPIManager::instance()->rank();
   const int numProcesses = walberla::MPIManager::instance()->numProcesses();

   // Static size of values

   std::vector<int> values( 1, rank );

   std::vector<int> result = walberla::mpi::allGatherv( values );

   WALBERLA_CHECK_EQUAL( result.size(), numeric_cast<size_t>( numProcesses ) );
   for( int i = 0; i < numProcesses; ++i )
      WALBERLA_CHECK_EQUAL( result[ uint_c(i) ], i );

   // Dynamic size of values

   values.assign( numeric_cast< std::vector<int>::size_type >(rank), rank );

   result = walberla::mpi::allGatherv( values );

   const auto expected_size = numProcesses * ( numProcesses + 1 ) / 2 - numProcesses; // Gaussian sum 0 + 1 + ... + n - 1

   WALBERLA_CHECK_EQUAL( result.size(), numeric_cast<size_t>( expected_size ) );
   uint_t pos = 0;
   for( int i = 0; i < numProcesses; ++i )
      for( int j = 0; j < i; ++j )
         WALBERLA_CHECK_EQUAL( result[pos++], i );

   // allGathervBuffer

   std::vector<int> vSend;
   vSend.push_back( 2 * rank );
   vSend.push_back( 2 * rank + 1 );
   if( rank % 2 == 0 )
      vSend.push_back( 23 );

   walberla::mpi::SendBuffer sBuffer;
   walberla::mpi::RecvBuffer rBuffer;
   sBuffer << vSend;

   walberla::mpi::allGathervBuffer( sBuffer, rBuffer );

   std::vector<int> vRecv;
   for( int i = 0; i != numProcesses; ++i )
   {
      rBuffer >> vRecv;
      WALBERLA_CHECK_EQUAL( vRecv.size(), ( i % 2 == 0 ) ? uint_t(3) : uint_t(2) );
      WALBERLA_CHECK_EQUAL( vRecv[0], 2 * i );
      WALBERLA_CHECK_EQUAL( vRecv[1], 2 * i + 1 );
      if( i % 2 == 0 )
         WALBERLA_CHECK_EQUAL( vRecv[2], 23 );
   }
}

int main( int argc, char * argv[] )
{
   using namespace walberla;

   debug::enterTestMode();

   WALBERLA_MPI_SECTION()
   {
      MPIManager::instance()->initializeMPI( &argc, &argv );
      MPIManager::instance()->useWorldComm();
   }

   runTestAllGatherv();

   const int numProcesses = walberla::MPIManager::instance()->numProcesses();
   for( int rank = 0; rank <numProcesses; ++rank )
   {
      runTestGatherv(rank);
   }

   return EXIT_SUCCESS;
}
