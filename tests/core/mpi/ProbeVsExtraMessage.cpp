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
//! \file ProbeVsExtraMessage.h
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Micro Benchmark, measuring time for different variable sized communications
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/math/Random.h"
#include "core/timing/TimingPool.h"

#include <iostream>
#include <sstream>

using namespace walberla;

int getProcToReceiveFrom()
{
   auto mpi = MPIManager::instance();
   int res = mpi->worldRank() -1;
   if ( res < 0 )
      res = mpi->numProcesses() -1;
   return res;
}

int getProcToSendTo()
{
   auto mpi = MPIManager::instance();
   int res = mpi->worldRank() +1;
   if ( res == mpi->numProcesses() )
      res = 0;
   return res;
}

int getRandomMessageSize( uint_t maxMessageSize )
{
   return int_c( math::intRandom( maxMessageSize / 3, maxMessageSize-1 ) );
}

void iprobeVersion( uint_t iterations, uint_t maxMessageSize, WcTimingPool & timingPool )
{
   char * sendBuf = new char[ maxMessageSize ];
   char * recvBuf = new char[ maxMessageSize ];

   auto & timer = timingPool["iprobe"];

   for( uint_t i =0; i < iterations; ++i )
   {
      int messageSize = getRandomMessageSize( maxMessageSize );

      timer.start();
      MPI_Request sendRequest;
      MPI_Isend( sendBuf, messageSize, MPI_BYTE, getProcToSendTo(), 0, MPI_COMM_WORLD, &sendRequest );

      MPI_Status probeStatus;
      MPI_Probe( getProcToReceiveFrom(), 0, MPI_COMM_WORLD, &probeStatus );

      int count = 0;
      MPI_Get_count( &probeStatus, MPI_BYTE, &count );
      MPI_Recv( recvBuf, count, MPI_BYTE, getProcToReceiveFrom(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      MPI_Waitall( 1, &sendRequest, MPI_STATUSES_IGNORE );
      timer.end();
   }

   delete [] sendBuf;
   delete [] recvBuf;
}


void twoMessageVersion( uint_t iterations, uint_t maxMessageSize, WcTimingPool & timingPool )
{
   char * sendBuf = new char[ maxMessageSize ];
   char * recvBuf = new char[ maxMessageSize ];

   auto & timer = timingPool["twoMessages"];

   int recvSize;
   MPI_Request sendRequests[2];

   const int TAG_SIZE_MSG    = 1;
   const int TAG_CONTENT_MSG = 0;

   for( uint_t i =0; i < iterations; ++i )
   {
      int sendSize = getRandomMessageSize( maxMessageSize );

      timer.start();

      MPI_Request recvSizeMsgRequest;
      MPI_Irecv( &recvSize, 1      , MPI_INT,  getProcToReceiveFrom(), TAG_SIZE_MSG,    MPI_COMM_WORLD, &recvSizeMsgRequest );
      MPI_Isend( &sendSize, 1      , MPI_INT,  getProcToSendTo(),      TAG_SIZE_MSG,    MPI_COMM_WORLD, &sendRequests[0]    );
      MPI_Isend( sendBuf , sendSize, MPI_BYTE, getProcToSendTo(),      TAG_CONTENT_MSG, MPI_COMM_WORLD, &sendRequests[1]    );

      // wait for size message to arrive
      MPI_Waitall( 1, &recvSizeMsgRequest, MPI_STATUSES_IGNORE );
      // receive content
      MPI_Recv( recvBuf, recvSize, MPI_BYTE, getProcToReceiveFrom(), TAG_CONTENT_MSG, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

      // wait for sends to finish
      MPI_Waitall( 2, sendRequests, MPI_STATUSES_IGNORE );

      timer.end();
   }

   delete [] sendBuf;
   delete [] recvBuf;

}

void maxMessageSizeVersion( uint_t iterations, uint_t maxMessageSize, WcTimingPool & timingPool )
{
   char * sendBuf = new char[ maxMessageSize ];
   char * recvBuf = new char[ maxMessageSize ];

   auto & timer = timingPool["maxMessageSizeKnown"];

   for( uint_t i =0; i < iterations; ++i )
   {
      int messageSize = getRandomMessageSize( maxMessageSize );
      timer.start();
      MPI_Request sendRequest;
      MPI_Request recvRequest;
      MPI_Irecv( recvBuf, int_c(maxMessageSize), MPI_BYTE,  getProcToReceiveFrom(), 0, MPI_COMM_WORLD, &recvRequest );
      MPI_Isend( sendBuf, messageSize,           MPI_BYTE,  getProcToSendTo(),      0, MPI_COMM_WORLD, &sendRequest );

      MPI_Status status;
      MPI_Waitall( 1, &recvRequest, &status );

      int count = 0;
      MPI_Get_count( &status, MPI_BYTE, &count );
      MPI_Waitall( 1, &sendRequest, MPI_STATUSES_IGNORE );
      timer.end();
   }

   delete [] sendBuf;
   delete [] recvBuf;
}





int main( int argc, char ** argv )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   MPIManager::instance()->useWorldComm();

   using namespace std;

   WALBERLA_ROOT_SECTION() {
      if ( argc != 3 ) {
         cerr << "Wrong number of arguments " << endl;
         cerr << "Usage ./probeVsExtraMessage <iterations> <messageSize> " << endl;
      }
   }
   uint_t iterations;
   uint_t messageSize;
   std::string arg1( argv[1] );
   std::string arg2( argv[2] );
   std::stringstream streamIterations  ( arg1 );
   std::stringstream streamMessageSize ( arg2 );
   streamIterations >> iterations;
   streamMessageSize >> messageSize;

   WcTimingPool tp;
   iprobeVersion         ( iterations, messageSize, tp );
   twoMessageVersion     ( iterations, messageSize, tp );
   maxMessageSizeVersion ( iterations, messageSize, tp );

   WALBERLA_ROOT_SECTION() {
      cout << tp << endl;
   }

   return 0;
}
