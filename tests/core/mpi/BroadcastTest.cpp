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
//! \file BroadcastTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Tests for core/mpi/Broadcast.h
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Broadcast.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/MPIManager.h"

#include <limits>
#include <cstdlib>
#include <vector>


template< typename T >
std::vector<T> fibonacci()
{
   std::vector<T> fib;
   fib.push_back( T(0) );
   fib.push_back( T(1) );

   while( std::numeric_limits<T>::max() - *(fib.end() - 1) >= *(fib.end() - 2) )
      fib.push_back( *(fib.end() - 1) + *(fib.end() - 2) );

   return fib;
}

void runTest( int senderRank )
{
   using namespace walberla;

   const std::vector<unsigned long long> fib( fibonacci<unsigned long long>() );
   const int  testInt = 42;

   std::vector< unsigned long long > allFib;
   int allTestInt = 0;

   WALBERLA_EXCLUSIVE_WORLD_SECTION( senderRank )
   {
      allFib = fib;
      allTestInt = testInt;
   }

   walberla::mpi::broadcastObject(allFib, senderRank);
   walberla::mpi::broadcastObject(allTestInt, senderRank);

   WALBERLA_CHECK_EQUAL( allFib, fib );
   WALBERLA_CHECK_EQUAL( allTestInt, testInt );
}

int main( int argc, char * argv[] )
{
   using namespace walberla;

   debug::enterTestMode();

   MPIManager::instance()->initializeMPI( &argc, &argv );

   for( int rank = 0; rank < MPIManager::instance()->numProcesses(); ++rank )
      runTest(rank);

   return EXIT_SUCCESS;
}
