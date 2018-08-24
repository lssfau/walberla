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
//! \file SetReductionTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"

#include "core/debug/TestSubsystem.h"

#include "core/logging/Logging.h"

#include "core/mpi/MPIManager.h"
#include "core/mpi/SetReduction.h"
#include "core/mpi/BufferDataTypeExtensions.h"

#include <algorithm>
#include <string>
#include <vector>

namespace set_reduction_test {

using namespace walberla;

void testRankSet()
{
   int rank = mpi::MPIManager::instance()->rank();
   int numberOfRanks = mpi::MPIManager::instance()->numProcesses();

   std::vector<int> value;
   value.push_back( rank );

   std::vector<int> reducedSetUnion = mpi::allReduceSet( value, mpi::UNION );
   std::vector<int> reducedSetIntersection = mpi::allReduceSet( value, mpi::INTERSECTION );

   WALBERLA_CHECK( reducedSetIntersection.empty() || ( numberOfRanks == 1 && reducedSetIntersection.size() == 1 ) );

   //std::ostringstream oss;
   //for( auto it = reducedSetUnion.begin(); it != reducedSetUnion.end(); ++it )
   //{
   //   oss << *it << ' ';
   //}
   //WALBERLA_LOG_DEVEL( oss.str() );

   WALBERLA_CHECK_EQUAL( reducedSetUnion.size(), numberOfRanks );

   int ctr = 0;
   for( auto it = reducedSetUnion.begin(); it != reducedSetUnion.end(); ++it )
   {
      WALBERLA_CHECK_EQUAL( *it, ctr );
      ++ctr;
   }
   

   value.clear();
   for( int i = 0; i < numberOfRanks; ++i )
   {
      value.push_back( i );
      value.push_back( i );
      value.push_back( i );
   }

   reducedSetUnion = mpi::allReduceSet( value, mpi::UNION );
   reducedSetIntersection = mpi::allReduceSet( value, mpi::INTERSECTION );

   WALBERLA_CHECK_EQUAL( reducedSetUnion.size(), numberOfRanks );

   ctr = 0;
   auto it2 = reducedSetIntersection.begin();
   for( auto it = reducedSetUnion.begin(); it != reducedSetUnion.end(); ++it, ++it2 )
   {
      WALBERLA_CHECK_EQUAL( *it, ctr );
      WALBERLA_CHECK_EQUAL( *it2, ctr );
      ++ctr;
   }   
}

static const int         NUM_FRUITS = 5;
static const std::string FRUITS[] = { "apple", "banana", "pear", "melon", "grapefruit" };

void testStrings()
{
   int rank = mpi::MPIManager::instance()->rank();
   int numProcesses = mpi::MPIManager::instance()->numProcesses();
   
   std::vector< std::string > values;
   values.push_back( FRUITS[rank % NUM_FRUITS] );

   std::vector< std::string > reducedValuesUnion = mpi::allReduceSet( values, mpi::UNION );

   values.emplace_back("GRAPES" );
   std::vector< std::string > reducedValuesIntersection = mpi::allReduceSet( values, mpi::INTERSECTION );

   if( numProcesses == 1 )
   {
      WALBERLA_CHECK_EQUAL( reducedValuesIntersection.size(), values.size() );
   }
   else
   {
      WALBERLA_CHECK_EQUAL( reducedValuesIntersection.size(), 1 );
      WALBERLA_CHECK_EQUAL( reducedValuesIntersection.front(), "GRAPES" );
   }

   //std::ostringstream oss;
   //for( auto it = reducedValuesUnion.begin(); it != reducedValuesUnion.end(); ++it )
   //{
   //   oss << *it << ' ';
   //}
   //WALBERLA_LOG_DEVEL( oss.str() );

   int numberOfProcesses = mpi::MPIManager::instance()->numProcesses();
   std::vector< std::string > expectedSet( FRUITS, FRUITS + std::min( numberOfProcesses, NUM_FRUITS ) );
   std::sort( expectedSet.begin(), expectedSet.end() );

   WALBERLA_CHECK_EQUAL( reducedValuesUnion.size(), expectedSet.size() );
   for( size_t i = 0; i < expectedSet.size(); ++i )
   {
      WALBERLA_CHECK_EQUAL( reducedValuesUnion[i], expectedSet[i] );
   }
}


} // namespace set_reduction_test

int main( int argc, char * argv[] )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   WALBERLA_MPI_SECTION() { walberla::MPIManager::instance()->useWorldComm(); }

   set_reduction_test::testRankSet();
   set_reduction_test::testStrings();
}
