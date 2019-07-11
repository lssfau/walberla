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
//! \file BufferSystemTest.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Tests for BufferSystem: symmetric and asymmetric MPI communication tests
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/Environment.h"

#include <random>

#include <cmath>
#include <iostream>
#include <set>
#include <thread>
#include <chrono>


using namespace walberla;
using mpi::BufferSystem;
using namespace std::literals::chrono_literals;



using base_generator_type = std::mt19937;

/**
 * Utility function for sleeping a random time
 * used to simulate a variable process load
 */
void randomSleep( int maxTimeInMs = 20 )
{
   static base_generator_type generator(42u);
   static unsigned  int counter =0;

   counter += 100;

   int rank          = MPIManager::instance()->worldRank();
   unsigned int seed = static_cast<unsigned int>(std::time(nullptr)) + static_cast<unsigned int>(rank*1000) + counter;
   generator.seed(seed);

   std::uniform_int_distribution<> uni_dist(0,maxTimeInMs);

   int sleepTime = uni_dist(generator);
   std::this_thread::sleep_for( sleepTime * 1ms );
}


/**
 * Every process sends a message containing its own rank
 * to the neighboring processes (1D , periodic boundary)
 */
void symmetricCommunication()
{
   const int MSG_SIZE = 10;

   auto mpiManager = MPIManager::instance();

   int numProcesses  = mpiManager->numProcesses();
   int rank          = mpiManager->worldRank();
   int leftNeighbor  = (rank-1+numProcesses)  % numProcesses;
   int rightNeighbor = (rank+1) % numProcesses;

   WALBERLA_CHECK_GREATER_EQUAL( numProcesses, 3 );


   BufferSystem bs ( MPI_COMM_WORLD );

   // Pack Message to left neighbor containing own rank
   for( int i=0; i< MSG_SIZE; ++i )
      bs.sendBuffer( leftNeighbor )  <<  rank;

   // Pack Message to right neighbor containing own rank
   for( int i=0; i< MSG_SIZE; ++i )
      bs.sendBuffer( rightNeighbor ) << rank;

   bs.setReceiverInfoFromSendBufferState( true, false );
   randomSleep();
   bs.sendAll();

   // In between we could do some computation
   randomSleep();

   for( auto it = bs.begin(); it != bs.end(); ++it )
   {
      WALBERLA_CHECK ( it.rank() == leftNeighbor || it.rank() == rightNeighbor );
      WALBERLA_CHECK_EQUAL( it.buffer().size(), MSG_SIZE * sizeof(int) + MSG_SIZE * mpi::BUFFER_DEBUG_OVERHEAD );

      int receivedVal = -1;
      it.buffer() >> receivedVal;

      WALBERLA_CHECK_EQUAL( receivedVal, it.rank() );
   }

   WALBERLA_CHECK_EQUAL( bs.getBytesSent(), (MSG_SIZE * sizeof(int) + MSG_SIZE * mpi::BUFFER_DEBUG_OVERHEAD) * 2 );
   WALBERLA_CHECK_EQUAL( bs.getBytesReceived(), (MSG_SIZE * sizeof(int) + MSG_SIZE * mpi::BUFFER_DEBUG_OVERHEAD) * 2 );
}

/**
 * Every process sends a message as big as his rank number
 * to the neighboring processes (1D , periodic boundary)
 */
void asymmetricCommunication(const bool useIProbe)
{
   auto mpiManager = MPIManager::instance();

   int numProcesses  = mpiManager->numProcesses();
   int rank          = mpiManager->worldRank();
   int leftNeighbor  = (rank-1+numProcesses)  % numProcesses;
   int rightNeighbor = (rank+1) % numProcesses;

   WALBERLA_CHECK_GREATER_EQUAL( numProcesses, 3 );


   BufferSystem bs ( MPI_COMM_WORLD );
   bs.useIProbe(useIProbe);

   // Set receiver information
   std::set<int> receiveFrom;
   if ( leftNeighbor  > 0 ) receiveFrom.insert( leftNeighbor );
   if ( rightNeighbor > 0 ) receiveFrom.insert( rightNeighbor );
   bs.setReceiverInfo( receiveFrom, false );


   const uint_t NUM_STEPS = 3;

   for ( uint_t step = 0; step < NUM_STEPS; ++step )
   {
      // Pack Messages to neighbors containing rank times rank value
      for( int i=0; i< rank; ++i )  bs.sendBuffer( leftNeighbor )  <<  rank;
      for( int i=0; i< rank; ++i )  bs.sendBuffer( rightNeighbor ) << rank;

      randomSleep();
      bs.sendAll();

      // In between we could do some computation
      randomSleep();

      for( auto it = bs.begin(); it != bs.end(); ++it )
      {
         if ( it.rank() == leftNeighbor )
         {
            for( int i=0; i < leftNeighbor; ++i ) {
               int value = -1;
               it.buffer() >> value;
               WALBERLA_CHECK_EQUAL( value, leftNeighbor );
            }
         }
         else if ( it.rank() == rightNeighbor )
         {
            for( int i=0; i < rightNeighbor; ++i ) {
               int value = -1;
               it.buffer() >> value;
               WALBERLA_CHECK_EQUAL( value, rightNeighbor );
            }
         }
         else
            WALBERLA_CHECK( false ); // unexpected sender

         WALBERLA_CHECK( it.buffer().isEmpty() );
      }
   }

   WALBERLA_CHECK_EQUAL( bs.getBytesSent(), int64_c(sizeof(int) + mpi::BUFFER_DEBUG_OVERHEAD) * int64_c(rank + rank) );
   WALBERLA_CHECK_EQUAL( bs.getBytesReceived(), int64_c(sizeof(int) + mpi::BUFFER_DEBUG_OVERHEAD) * int64_c(leftNeighbor + rightNeighbor) );
}


// like asymmetricCommunication, but the message size is a random value
// that changes every communication step
void timeVaryingCommunication(const bool useIProbe)
{
   auto mpiManager = MPIManager::instance();

   int numProcesses  = mpiManager->numProcesses();
   int rank          = mpiManager->worldRank();
   int leftNeighbor  = (rank-1+numProcesses)  % numProcesses;
   int rightNeighbor = (rank+1) % numProcesses;

   WALBERLA_CHECK_GREATER_EQUAL( numProcesses, 3 );

   BufferSystem bs ( MPI_COMM_WORLD );
   bs.useIProbe(useIProbe);

   // artificial special case: no message from root
   bs.sendBuffer( rightNeighbor );
   bs.sendBuffer( leftNeighbor );
   bs.setReceiverInfoFromSendBufferState( false, true );


   const uint_t NUM_STEPS = 5;
   for ( uint_t step = 1; step <= NUM_STEPS; ++step )
   {
      for( uint_t i=0; i < std::max<uint_t>( uint_c(rank * leftNeighbor) * step % 17, 1ul); ++i )
         bs.sendBuffer( leftNeighbor ) << i;
      bs.send( leftNeighbor );

      for( uint_t i=0; i < std::max<uint_t>( uint_c(rank * rightNeighbor) * step % 17, 1ul); ++i )
         bs.sendBuffer( rightNeighbor ) << i;
      bs.send( rightNeighbor );


      WALBERLA_CHECK( bs.isCommunicationRunning()  );

      for( auto it = bs.begin(); it != bs.end(); ++it )
      {
         if ( it.rank() == leftNeighbor )
         {
            for( uint_t i=0; i < std::max<uint_t>( uint_c(rank * leftNeighbor) * step % 17, 1ul); ++i ) {
               uint_t value = 0;
               it.buffer() >> value;
               WALBERLA_CHECK_EQUAL( value, i );
            }
         }
         else if ( it.rank() == rightNeighbor )
         {
            for( uint_t i=0; i < std::max<uint_t>( uint_c(rank * rightNeighbor) * step % 17,1ul); ++i ) {
               uint_t value = 0;
               it.buffer() >> value;
               WALBERLA_CHECK_EQUAL( value, i );
            }
         }
         else
            WALBERLA_CHECK( false ); // unexpected sender

         WALBERLA_CHECK( it.buffer().isEmpty() );

      }
      WALBERLA_CHECK( ! bs.isCommunicationRunning()  );
   }
}





/**
 * Gathering using asymmetric communication
 *    every process sends a message of size rank*sizeof(int) containing only its own rank to root process
 *    i.e. rank 1 sends a "1" once, rank 2 sends a message containing two "2"'s ...
 */

void gatherUsingAsymmetricCommunication(const bool useIProbe)
{
   int rank          = MPIManager::instance()->worldRank();
   int numProcesses  = MPIManager::instance()->numProcesses();

   WALBERLA_CHECK_GREATER_EQUAL( numProcesses, 3 );

   const int TAG=42;

   BufferSystem bs (MPI_COMM_WORLD, TAG );
   bs.useIProbe(useIProbe);


   if ( rank ==0 )
      bs.setReceiverInfo( BufferSystem::allRanksButRoot(), true );
   else
      bs.setReceiverInfo( std::set<mpi::MPIRank>(), true );

   if(rank > 0)
   {
      for( int i=0; i < rank; ++i )
         bs.sendBuffer(0) << rank;
   }

   bs.sendAll();
   randomSleep();

   for( auto it = bs.begin(); it != bs.end(); ++it )
   {
      WALBERLA_CHECK( rank == 0); // only root should receive something

      for( int i=0; i < it.rank(); ++i )
      {
         int received = -1;
         it.buffer() >> received;
         WALBERLA_CHECK_EQUAL( received, it.rank() );
      }
   }

}


void selfSend()
{
   int rank          = MPIManager::instance()->worldRank();
   int numProcesses  = MPIManager::instance()->numProcesses();

   WALBERLA_CHECK_GREATER_EQUAL( numProcesses, 3 );

   const int TAG=42;

   BufferSystem bs (MPI_COMM_WORLD, TAG );


   if ( rank ==0 )
      bs.setReceiverInfo( BufferSystem::allRanks(), true );
   else
      bs.setReceiverInfo( std::set<mpi::MPIRank>(), true );

   bs.sendBuffer(0) << rank;

   bs.sendAll();
   randomSleep();

   for( auto it = bs.begin(); it != bs.end(); ++it )
   {
      WALBERLA_CHECK( rank == 0); // only root should receive something

      int received = -1;
      it.buffer() >> received;
      WALBERLA_CHECK_EQUAL( received, it.rank() );
   }
}

void copyTest()
{
   int rank = MPIManager::instance()->worldRank();

   BufferSystem bs1( MPI_COMM_WORLD, 3 );
   {
      BufferSystem bs2( MPI_COMM_WORLD, 7 );
      bs2.sendBuffer(rank) << int(42);
      bs2.setReceiverInfoFromSendBufferState( true, false );
      bs2.sendAll();

      for ( auto i = bs2.begin(); i != bs2.end(); ++i )
      {
         int messageContent;
         i.buffer() >> messageContent;
         WALBERLA_CHECK_EQUAL(messageContent, 42);
      }

      bs1 = bs2;

   }

   bs1.sendBuffer(rank) << int(42);
   bs1.sendAll();
   for ( auto i = bs1.begin(); i != bs1.end(); ++i )
   {
      int messageContent;
      i.buffer() >> messageContent;
      WALBERLA_CHECK_EQUAL(messageContent, 42);
   }

}

int main(int argc, char**argv)
{
   mpi::Environment mpiEnv( argc, argv );
   debug::enterTestMode();

   auto mpiManager = MPIManager::instance();
   int numProcesses  = mpiManager->numProcesses();

   if(numProcesses <= 2)
   {
      WALBERLA_ABORT("This test has to be executed on at least 3 processes. Executed on " <<  numProcesses);
      return 1;
   }

   WALBERLA_LOG_INFO_ON_ROOT("Testing Symmetric Communication...");
   symmetricCommunication();

   WALBERLA_LOG_INFO_ON_ROOT("Testing Asymmetric Communication...");
   asymmetricCommunication(false);
   asymmetricCommunication(true);

   WALBERLA_LOG_INFO_ON_ROOT("Testing time-varying Communication...");
   timeVaryingCommunication(false);
   timeVaryingCommunication(true);

   WALBERLA_LOG_INFO_ON_ROOT("Testing Gather Operation...");
   gatherUsingAsymmetricCommunication(false);
   gatherUsingAsymmetricCommunication(true);

   WALBERLA_LOG_INFO_ON_ROOT("Testing self-send...");
   selfSend();

   WALBERLA_LOG_INFO_ON_ROOT("Testing Buffer System copy...");
   copyTest();

   return EXIT_SUCCESS;
}
