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
//! \file ReduceTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Tests for core/mpi/Reduce.h
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"

#include <limits>
#include <vector>


void runTestReduce( int recvRank )
{
   using namespace walberla;

   const int rank = walberla::MPIManager::instance()->rank();
   const int numProcesses = walberla::MPIManager::instance()->numProcesses();


   int value = rank;

   int result = walberla::mpi::reduce( value, walberla::mpi::SUM, recvRank );
   walberla::mpi::reduceInplace( value, walberla::mpi::SUM, recvRank );

   if( rank == recvRank )
   {
      int sum = 0;
      for( int i = 0; i < numProcesses; ++i )
         sum += i;
      WALBERLA_CHECK_EQUAL( value, sum );
      WALBERLA_CHECK_EQUAL( result, sum );
   }
   else
   {
      WALBERLA_CHECK_EQUAL( value, rank );
      WALBERLA_CHECK_EQUAL( result, 0 );
   }

   value = rank;

   walberla::mpi::reduceInplace( value, walberla::mpi::PRODUCT, recvRank );

   if( rank == recvRank )
   {
      WALBERLA_CHECK_EQUAL( value, 0 );
   }
   else
   {
      WALBERLA_CHECK_EQUAL( value, rank );
   }


   Vector3<int> intVec3{rank-1, rank, rank+1};
   auto resultIntVec3 = walberla::mpi::reduce( intVec3, walberla::mpi::SUM, recvRank );
   walberla::mpi::reduceInplace( intVec3, walberla::mpi::SUM, recvRank );
   if( rank == recvRank )
   {
      int sum = 0;
      for( int i = 0; i < numProcesses; ++i )
         sum += i;

      math::Vector3<int> refResultIntVec3 {sum - numProcesses, sum, sum + numProcesses};

      for(uint_t i = 0; i < 3; ++i)
      {
         WALBERLA_CHECK_EQUAL( intVec3[i], refResultIntVec3[i], "Reduce Vector3, component " << i);
         WALBERLA_CHECK_EQUAL( resultIntVec3[i], refResultIntVec3[i], "ReduceInPlace Vector3, component " << i);
      }

   }
   else
   {
      Vector3<int> refResultIntVec3{rank-1, rank, rank+1};
      for(uint_t i = 0; i < 3; ++i)
      {
         WALBERLA_CHECK_EQUAL( intVec3[i], refResultIntVec3[i], "Reduce Vector3, component " << i);
         WALBERLA_CHECK_EQUAL( resultIntVec3[i], 0, "ReduceInPlace Vector3, component " << i);
      }
   }


   std::vector<int> some_ints( 100 );
   for( int i = 0; i < 100; ++i )
      some_ints[ uint_c(i) ] = rank + i;

   walberla::mpi::reduceInplace( some_ints, walberla::mpi::SUM, recvRank );

   if( rank == recvRank )
   {
      int sum = 0;
      for( int i = 0; i < numProcesses; ++i )
         sum += i;

      for( int i = 0; i < 100; ++i )
         WALBERLA_CHECK_EQUAL( some_ints[ uint_c(i) ], sum + numProcesses * i );
   }

   {
      std::vector<int> empty;
      mpi::reduceInplace( empty,  walberla::mpi::SUM, recvRank );
      WALBERLA_CHECK_EQUAL( empty.size(), 0 );
   }
}

void runTestAllReduce()
{
   using namespace walberla;

   const int rank = walberla::MPIManager::instance()->rank();
   const int numProcesses = walberla::MPIManager::instance()->numProcesses();


   int value = rank;

   walberla::mpi::allReduceInplace( value, walberla::mpi::SUM );

   int sum = 0;
   for( int i = 0; i < numProcesses; ++i )
      sum += i;
   WALBERLA_CHECK_EQUAL( value, sum );

   value = rank;

   walberla::mpi::allReduceInplace( value, walberla::mpi::PRODUCT );

   WALBERLA_CHECK_EQUAL( value, 0 );


   Vector3<int> intVec3{rank-1, rank, rank+1};
   walberla::mpi::allReduceInplace( intVec3, walberla::mpi::SUM );
   math::Vector3<int> refResultIntVec3 {sum - numProcesses, sum, sum + numProcesses};

   for(uint_t i = 0; i < 3; ++i)
   {
      WALBERLA_CHECK_EQUAL( intVec3[i], refResultIntVec3[i], "AllReduce Vector3, component " << i);
   }


   std::vector<int> some_ints( 100 );
   for( int i = 0; i < 100; ++i )
      some_ints[ uint_c(i) ] = rank + i;

   walberla::mpi::allReduceInplace( some_ints, walberla::mpi::SUM );
   for( int i = 0; i < 100; ++i )
      WALBERLA_CHECK_EQUAL( some_ints[ uint_c(i) ], sum + numProcesses * i );

   {
      std::vector<int> empty;
      mpi::allReduceInplace( empty,  walberla::mpi::SUM );
      WALBERLA_CHECK_EQUAL( empty.size(), 0 );
   }
}

void runTestAllReduceBool()
{
   using namespace walberla;

   const int rank = walberla::MPIManager::instance()->rank();
   const int numProcesses = walberla::MPIManager::instance()->numProcesses();

   const bool allTrue  = true;
   const bool allFalse = false;
   const bool oneTrue  = rank == 0;
   const bool mixed    = ( rank % 2 ) != 0;

   WALBERLA_CHECK_EQUAL( mpi::allReduce( allTrue,  mpi::LOGICAL_AND ), true  );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( allFalse, mpi::LOGICAL_AND ), false );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( oneTrue,  mpi::LOGICAL_AND ), numProcesses == 1 );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( mixed,    mpi::LOGICAL_AND ), false );

   WALBERLA_CHECK_EQUAL( mpi::allReduce( allTrue,  mpi::LOGICAL_OR ), true  );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( allFalse, mpi::LOGICAL_OR ), false );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( oneTrue,  mpi::LOGICAL_OR ), true  );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( mixed,    mpi::LOGICAL_OR ), numProcesses >= 2  );

   WALBERLA_CHECK_EQUAL( mpi::allReduce( allTrue,  mpi::LOGICAL_XOR ), numProcesses == 1 );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( allFalse, mpi::LOGICAL_XOR ), false );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( oneTrue,  mpi::LOGICAL_XOR ), true  );
   WALBERLA_CHECK_EQUAL( mpi::allReduce( mixed,    mpi::LOGICAL_XOR ), numProcesses == 2 );

   {
      bool copyAllTrue  = allTrue;
      bool copyAllFalse = allFalse;
      bool copyOneTrue  = oneTrue;
      bool copyMixed    = mixed;

      mpi::allReduceInplace( copyAllTrue,  mpi::LOGICAL_AND );
      mpi::allReduceInplace( copyAllFalse, mpi::LOGICAL_AND );
      mpi::allReduceInplace( copyOneTrue,  mpi::LOGICAL_AND );
      mpi::allReduceInplace( copyMixed,    mpi::LOGICAL_AND );

      WALBERLA_CHECK_EQUAL( copyAllTrue,  true  );
      WALBERLA_CHECK_EQUAL( copyAllFalse, false );
      WALBERLA_CHECK_EQUAL( copyOneTrue,  numProcesses == 1 );
      WALBERLA_CHECK_EQUAL( copyMixed,    false );
   }

   {
      bool copyAllTrue  = allTrue;
      bool copyAllFalse = allFalse;
      bool copyOneTrue  = oneTrue;
      bool copyMixed    = mixed;

      mpi::allReduceInplace( copyAllTrue,  mpi::LOGICAL_OR );
      mpi::allReduceInplace( copyAllFalse, mpi::LOGICAL_OR );
      mpi::allReduceInplace( copyOneTrue,  mpi::LOGICAL_OR );
      mpi::allReduceInplace( copyMixed,    mpi::LOGICAL_OR );

      WALBERLA_CHECK_EQUAL( copyAllTrue,  true  );
      WALBERLA_CHECK_EQUAL( copyAllFalse, false );
      WALBERLA_CHECK_EQUAL( copyOneTrue,  true  );
      WALBERLA_CHECK_EQUAL( copyMixed,    numProcesses >= 2  );
   }

   {
      bool copyAllTrue  = allTrue;
      bool copyAllFalse = allFalse;
      bool copyOneTrue  = oneTrue;
      bool copyMixed    = mixed;

      mpi::allReduceInplace( copyAllTrue,  mpi::LOGICAL_XOR );
      mpi::allReduceInplace( copyAllFalse, mpi::LOGICAL_XOR );
      mpi::allReduceInplace( copyOneTrue,  mpi::LOGICAL_XOR );
      mpi::allReduceInplace( copyMixed,    mpi::LOGICAL_XOR );

      WALBERLA_CHECK_EQUAL( copyAllTrue,  numProcesses == 1 );
      WALBERLA_CHECK_EQUAL( copyAllFalse, false );
      WALBERLA_CHECK_EQUAL( copyOneTrue,  true  );
      WALBERLA_CHECK_EQUAL( copyMixed,    numProcesses == 2 );
   }


   math::Vector3<bool> boolVec;
   boolVec[0] = true;
   boolVec[1] = false;
   boolVec[2] = rank == 0;

   {
      math::Vector3<bool> result( boolVec );
      mpi::allReduceInplace( result,  mpi::LOGICAL_AND );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), numProcesses == 1 );
   }
   {
      math::Vector3<bool> result( boolVec );
      mpi::allReduceInplace( result,  mpi::LOGICAL_OR );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
   }
   {
      math::Vector3<bool> result( boolVec );
      mpi::allReduceInplace( result,  mpi::LOGICAL_XOR );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), numProcesses == 1 );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
   }


   std::vector<bool> bools( 4u );
   bools[0] = true;
   bools[1] = false;
   bools[2] = rank == 0;
   bools[3] = ( rank % 2 ) != 0;

   {
      std::vector<bool> result( bools );
      mpi::allReduceInplace( result,  mpi::BITWISE_AND );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), numProcesses == 1 );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[3] ), false );
   }
   {
      std::vector<bool> result( bools );
      mpi::allReduceInplace( result,  mpi::BITWISE_OR );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[3] ), numProcesses >= 2  );
   }
   {
      std::vector<bool> result( bools );
      mpi::allReduceInplace( result,  mpi::BITWISE_XOR );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), numProcesses == 1 );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
      WALBERLA_CHECK_EQUAL( static_cast<bool>( result[3] ), numProcesses == 2 );
   }
   {
      std::vector<bool> empty;
      mpi::allReduceInplace( empty,  mpi::BITWISE_XOR );
      WALBERLA_CHECK_EQUAL( empty.size(), 0 );
   }
}

void runTestReduceBool( int recvRank )
{
   using namespace walberla;

   const int rank = walberla::MPIManager::instance()->rank();
   const int numProcesses = walberla::MPIManager::instance()->numProcesses();

   const bool allTrue  = true;
   const bool allFalse = false;
   const bool oneTrue  = rank == 0;
   const bool mixed    = ( rank % 2 ) != 0;

   {
      bool result0 = mpi::reduce( allTrue,  mpi::LOGICAL_AND, recvRank );
      bool result1 = mpi::reduce( allFalse, mpi::LOGICAL_AND, recvRank );
      bool result2 = mpi::reduce( oneTrue,  mpi::LOGICAL_AND, recvRank );
      bool result3 = mpi::reduce( mixed,    mpi::LOGICAL_AND, recvRank );

      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( result0, true  );
         WALBERLA_CHECK_EQUAL( result1, false );
         WALBERLA_CHECK_EQUAL( result2, numProcesses == 1 );
         WALBERLA_CHECK_EQUAL( result3, false );
      }
      else
      {
         WALBERLA_CHECK_EQUAL( result0, false );
         WALBERLA_CHECK_EQUAL( result1, false );
         WALBERLA_CHECK_EQUAL( result2, false );
         WALBERLA_CHECK_EQUAL( result3, false );
      }
   }

   {
      bool result0 = mpi::reduce( allTrue,  mpi::LOGICAL_OR, recvRank);
      bool result1 = mpi::reduce( allFalse, mpi::LOGICAL_OR, recvRank);
      bool result2 = mpi::reduce( oneTrue,  mpi::LOGICAL_OR, recvRank);
      bool result3 = mpi::reduce( mixed,    mpi::LOGICAL_OR, recvRank);

      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( result0, true  );
         WALBERLA_CHECK_EQUAL( result1, false );
         WALBERLA_CHECK_EQUAL( result2, true  );
         WALBERLA_CHECK_EQUAL( result3, numProcesses >= 2 );
      }
      else
      {
         WALBERLA_CHECK_EQUAL( result0, false );
         WALBERLA_CHECK_EQUAL( result1, false );
         WALBERLA_CHECK_EQUAL( result2, false );
         WALBERLA_CHECK_EQUAL( result3, false );
      }
   }

   {
      bool result0 = mpi::reduce( allTrue,  mpi::LOGICAL_XOR, recvRank);
      bool result1 = mpi::reduce( allFalse, mpi::LOGICAL_XOR, recvRank);
      bool result2 = mpi::reduce( oneTrue,  mpi::LOGICAL_XOR, recvRank);
      bool result3 = mpi::reduce( mixed,    mpi::LOGICAL_XOR, recvRank);

      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( result0, numProcesses == 1 );
         WALBERLA_CHECK_EQUAL( result1, false );
         WALBERLA_CHECK_EQUAL( result2, true  );
         WALBERLA_CHECK_EQUAL( result3, numProcesses == 2 );
      }
      else
      {
         WALBERLA_CHECK_EQUAL( result0, false );
         WALBERLA_CHECK_EQUAL( result1, false );
         WALBERLA_CHECK_EQUAL( result2, false );
         WALBERLA_CHECK_EQUAL( result3, false );
      }
   }

   {
      bool copyAllTrue  = allTrue;
      bool copyAllFalse = allFalse;
      bool copyOneTrue  = oneTrue;
      bool copyMixed    = mixed;

      mpi::reduceInplace( copyAllTrue,  mpi::LOGICAL_AND, recvRank );
      mpi::reduceInplace( copyAllFalse, mpi::LOGICAL_AND, recvRank );
      mpi::reduceInplace( copyOneTrue,  mpi::LOGICAL_AND, recvRank );
      mpi::reduceInplace( copyMixed,    mpi::LOGICAL_AND, recvRank );

      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( copyAllTrue,  true  );
         WALBERLA_CHECK_EQUAL( copyAllFalse, false );
         WALBERLA_CHECK_EQUAL( copyOneTrue,  numProcesses == 1 );
         WALBERLA_CHECK_EQUAL( copyMixed,    false );
      }
   }

   {
      bool copyAllTrue  = allTrue;
      bool copyAllFalse = allFalse;
      bool copyOneTrue  = oneTrue;
      bool copyMixed    = mixed;

      mpi::reduceInplace( copyAllTrue,  mpi::LOGICAL_OR, recvRank );
      mpi::reduceInplace( copyAllFalse, mpi::LOGICAL_OR, recvRank );
      mpi::reduceInplace( copyOneTrue,  mpi::LOGICAL_OR, recvRank );
      mpi::reduceInplace( copyMixed,    mpi::LOGICAL_OR, recvRank );

      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( copyAllTrue,  true  );
         WALBERLA_CHECK_EQUAL( copyAllFalse, false );
         WALBERLA_CHECK_EQUAL( copyOneTrue,  true  );
         WALBERLA_CHECK_EQUAL( copyMixed,    numProcesses >= 2  );
      }
   }

   {
      bool copyAllTrue  = allTrue;
      bool copyAllFalse = allFalse;
      bool copyOneTrue  = oneTrue;
      bool copyMixed    = mixed;

      mpi::reduceInplace( copyAllTrue,  mpi::LOGICAL_XOR, recvRank );
      mpi::reduceInplace( copyAllFalse, mpi::LOGICAL_XOR, recvRank );
      mpi::reduceInplace( copyOneTrue,  mpi::LOGICAL_XOR, recvRank );
      mpi::reduceInplace( copyMixed,    mpi::LOGICAL_XOR, recvRank );

      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( copyAllTrue,  numProcesses == 1 );
         WALBERLA_CHECK_EQUAL( copyAllFalse, false );
         WALBERLA_CHECK_EQUAL( copyOneTrue,  true  );
         WALBERLA_CHECK_EQUAL( copyMixed,    numProcesses == 2 );
      }
   }

   math::Vector3<bool> boolVec;
   boolVec[0] = true;
   boolVec[1] = false;
   boolVec[2] = rank == 0;

   {
      math::Vector3<bool> result = mpi::reduce( boolVec,  mpi::LOGICAL_AND, recvRank );
      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), numProcesses == 1 );
      } else
      {
         for(uint_t i = 0; i < 3; ++i) WALBERLA_CHECK_EQUAL( static_cast<bool>( result[i] ), false  );
      }
   }
   {
      math::Vector3<bool> result = mpi::reduce( boolVec,  mpi::LOGICAL_OR, recvRank );
      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
      } else
      {
         for(uint_t i = 0; i < 3; ++i) WALBERLA_CHECK_EQUAL( static_cast<bool>( result[i] ), false  );
      }
   }
   {
      math::Vector3< bool > result = mpi::reduce(boolVec, mpi::LOGICAL_XOR, recvRank);
      if (rank == recvRank)
      {
         WALBERLA_CHECK_EQUAL(static_cast< bool >(result[0]), numProcesses == 1);
         WALBERLA_CHECK_EQUAL(static_cast< bool >(result[1]), false);
         WALBERLA_CHECK_EQUAL(static_cast< bool >(result[2]), true);
      } else
      {
         for(uint_t i = 0; i < 3; ++i) WALBERLA_CHECK_EQUAL( static_cast<bool>( result[i] ), false  );
      }
   }

   {
      math::Vector3<bool> result( boolVec );
      mpi::reduceInplace( result,  mpi::LOGICAL_AND, recvRank );
      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), numProcesses == 1 );
      }
   }
   {
      math::Vector3<bool> result( boolVec );
      mpi::reduceInplace( result,  mpi::LOGICAL_OR, recvRank );

      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
      }
   }
   {
      math::Vector3<bool> result( boolVec );
      mpi::reduceInplace( result,  mpi::LOGICAL_XOR, recvRank );
      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), numProcesses == 1 );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
      }
   }


   std::vector<bool> bools( 4u );
   bools[0] = true;
   bools[1] = false;
   bools[2] = rank == 0;
   bools[3] = ( rank % 2 ) != 0;

   {
      std::vector<bool> result( bools );
      mpi::reduceInplace( result,  mpi::BITWISE_AND, recvRank );
      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), numProcesses == 1 );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[3] ), false );
      }
   }
   {
      std::vector<bool> result( bools );
      mpi::reduceInplace( result,  mpi::BITWISE_OR, recvRank );

      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), true  );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[3] ), numProcesses >= 2  );
      }
   }
   {
      std::vector<bool> result( bools );
      mpi::reduceInplace( result,  mpi::BITWISE_XOR, recvRank );
      if( rank == recvRank )
      {
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[0] ), numProcesses == 1 );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[1] ), false );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[2] ), true  );
         WALBERLA_CHECK_EQUAL( static_cast<bool>( result[3] ), numProcesses == 2 );
      }
   }

   {
      std::vector<bool> empty;
      mpi::reduceInplace( empty,  mpi::BITWISE_XOR, recvRank );
      WALBERLA_CHECK_EQUAL( empty.size(), 0 );
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

   runTestAllReduce();
   runTestAllReduceBool();

   for( int rank = 0; rank < MPIManager::instance()->numProcesses(); ++rank )
   {
      runTestReduce( rank );
      runTestReduceBool( rank );
   }

   return EXIT_SUCCESS;
}
