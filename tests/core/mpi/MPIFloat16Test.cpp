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
//! \file MPIFloat16.cpp
//! \ingroup core
//! \author Michael Zikeli <michael.zikeli@fau.de>
//! \brief This test is to check whether the self defined MPI_Datatype and the self defined Operators are working.
//!    To verify the type, some parts of the BufferSystemTest are just copied.
//!    To verify the operations, a simple AllReduce is executed for each operation.
//!        For now only { SUM, MIN, MAX } are implemented, thus only those are tested.
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/Environment.h"
#include "core/mpi/Reduce.h"


namespace walberla::mpifloat16test
{

using mpi::BufferSystem;
using namespace std::literals::chrono_literals;

#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
void symmetricCommunication()
{
   const int MSG_SIZE = 10;

   auto mpiManager = MPIManager::instance();

   const int numProcesses  = mpiManager->numProcesses();
   const int rank          = mpiManager->worldRank();
   const int leftNeighbor  = (rank-1+numProcesses)  % numProcesses;
   const int rightNeighbor = (rank+1) % numProcesses;

   WALBERLA_CHECK_GREATER_EQUAL( numProcesses, 3 );


   BufferSystem bs ( MPI_COMM_WORLD );

   // Pack Message to left neighbor containing own rank
   for( int i=0; i< MSG_SIZE; ++i )
      bs.sendBuffer( leftNeighbor )  << numeric_cast<float16>(rank) + float16(0.3);

   // Pack Message to right neighbor containing own rank
   for( int i=0; i< MSG_SIZE; ++i )
      bs.sendBuffer( rightNeighbor ) << numeric_cast<float16>(rank) - float16(0.3);

   bs.setReceiverInfoFromSendBufferState( true, false );
   bs.sendAll();

   for( auto it = bs.begin(); it != bs.end(); ++it )
   {
      WALBERLA_CHECK ( it.rank() == leftNeighbor || it.rank() == rightNeighbor );
      WALBERLA_CHECK_EQUAL( it.buffer().size(), MSG_SIZE * sizeof(float16) + MSG_SIZE * mpi::BUFFER_DEBUG_OVERHEAD );

      auto receivedVal = float16(-1);
      it.buffer() >> receivedVal;

      WALBERLA_CHECK_EQUAL( typeid(receivedVal), typeid(float16) );

      if ( it.rank() == leftNeighbor )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( receivedVal, numeric_cast<float16>( it.rank() ) - float16(0.3) );
         WALBERLA_CHECK_FLOAT_UNEQUAL( receivedVal, numeric_cast<float16>( it.rank() ) + float16(0.3), 0.5);
      } else {
         WALBERLA_CHECK_FLOAT_EQUAL( receivedVal, numeric_cast<float16>( it.rank() ) + float16(0.3) );
         WALBERLA_CHECK_FLOAT_UNEQUAL( receivedVal, numeric_cast<float16>( it.rank() ) - float16(0.3), 0.5);
      }
   }

   WALBERLA_CHECK_EQUAL( bs.getBytesSent(), (MSG_SIZE * sizeof(float16) + MSG_SIZE * mpi::BUFFER_DEBUG_OVERHEAD) * 2 );
   WALBERLA_CHECK_EQUAL( bs.getBytesReceived(), (MSG_SIZE * sizeof(float16) + MSG_SIZE * mpi::BUFFER_DEBUG_OVERHEAD) * 2 );
}

void reduce( )
{

   auto mpiManager = MPIManager::instance();

   const int numProcesses  = mpiManager->numProcesses();
   const int rank          = mpiManager->worldRank();

   const auto startValue = numeric_cast<float16>(rank) + float16(0.3);

   // SUM
   auto value = startValue;

   walberla::mpi::allReduceInplace( value, walberla::mpi::SUM );

   auto sum = float16( 0.0 );
   for( int i = 0; i < numProcesses; ++i )
      sum += numeric_cast<float16>(i) + float16(0.3);
   WALBERLA_CHECK_FLOAT_EQUAL( value, sum );
   WALBERLA_CHECK_FLOAT_UNEQUAL( value, ((numProcesses*(numProcesses-1)/2.)+0.3*numProcesses), 1e-10 );

   // MIN
   value = startValue;

   walberla::mpi::allReduceInplace( value, walberla::mpi::MIN );
   WALBERLA_CHECK_FLOAT_EQUAL( value, numeric_cast<float16>(0.3) );

   // MAX
   value = startValue;

   walberla::mpi::allReduceInplace( value, walberla::mpi::MAX );
   WALBERLA_CHECK_FLOAT_EQUAL( value, numeric_cast<float16>(numProcesses-1)+numeric_cast<float16>(0.3) );

}
#endif // WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

int main( int argc, char** argv )
{
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
   mpi::Environment mpiEnv( argc, argv );
   debug::enterTestMode();
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );

   auto mpiManager   = MPIManager::instance();
   auto numProcesses = mpiManager->numProcesses();

   if(numProcesses <= 2)
   {
      WALBERLA_ABORT("This test has to be executed on at least 3 processes. Executed on " <<  numProcesses);
      return EXIT_FAILURE;
   }

   WALBERLA_LOG_INFO_ON_ROOT("Testing Symmetric Communication...");
   symmetricCommunication();

   WALBERLA_LOG_INFO_ON_ROOT("Testing Reduce operations...");
   reduce( );
#else
   WALBERLA_LOG_WARNING_ON_ROOT( "\nWALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT is not enabled. So this test cannot be run!\n" )
#endif // WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

   return EXIT_SUCCESS;

} // mpifloat16test::main()

} // namespace walberla::mpifloat16test

int main( int argc, char** argv )
{
   return walberla::mpifloat16test::main( argc, argv );
} // main()
