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
//! \file GatherTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Tests for core/mpi/Gather.h
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Gather.h"
#include "core/mpi/MPIManager.h"

#include <limits>
#include <vector>


void runTestGather( int recvRank )
{
   using namespace walberla;

   const int rank         = walberla::MPIManager::instance()->rank();
   const int numProcesses = walberla::MPIManager::instance()->numProcesses();

   int value = rank;

   std::vector<int> result = walberla::mpi::gather( value, recvRank );

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
}

void runTestAllGather()
{
   using namespace walberla;

   const int rank         = walberla::MPIManager::instance()->rank();
   const int numProcesses = walberla::MPIManager::instance()->numProcesses();

   int value = rank;

   std::vector<int> result = walberla::mpi::allGather( value );

   WALBERLA_CHECK_EQUAL( result.size(), numeric_cast<size_t>( numProcesses ) );
   for( int i = 0; i < numProcesses; ++i )
      WALBERLA_CHECK_EQUAL( result[ uint_c(i) ], i );

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

   runTestAllGather();

   const int numProcesses = walberla::MPIManager::instance()->numProcesses();
   for( int rank = 0; rank <numProcesses; ++rank )
   {
      runTestGather(rank);
   }

   return EXIT_SUCCESS;
}
