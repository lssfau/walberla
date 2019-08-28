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
//! \file MPIGatherScheme.cpp
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Communicating gathered data using MPI
//
//======================================================================================================================

#include "CommonSchemeFunctions.h"
#include "MPIGatherScheme.h"

#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "domain_decomposition/BlockStorage.h"


namespace walberla {
namespace gather {


MPIGatherScheme::MPIGatherScheme( BlockStorage & blockStorage, int gatherRank, uint_t everyNTimestep )
   : blocks_ ( blockStorage ), gatherRankInGlobalComm_( gatherRank),
     setupPhaseDone_( false ), gatherCommunicator_( MPI_COMM_NULL ), gatherRank_ ( -1 ),
     gatherMsgSize_(0), bytesToSend_(0), everyNTimestep_( everyNTimestep )
{
   WALBERLA_ASSERT_LESS         ( gatherRankInGlobalComm_, MPIManager::instance()->numProcesses() );
   WALBERLA_ASSERT_GREATER_EQUAL( gatherRankInGlobalComm_, 0 );
}

MPIGatherScheme::~MPIGatherScheme()
{
   if ( gatherCommunicator_ != MPI_COMM_NULL )
      MPI_Comm_free( &gatherCommunicator_ );
}


void MPIGatherScheme::addPackInfo( const shared_ptr<GatherPackInfo> & pi )
{
   packInfos_.push_back(pi);

   // Setup phase has to be repeated, if a new GatherPackInfo was added
   if ( setupPhaseDone_ )
   {
      if ( gatherCommunicator_ != MPI_COMM_NULL )
         MPI_Comm_free( &gatherCommunicator_ );

      setupPhaseDone_ = false;
   }
}


void MPIGatherScheme::setupGatherCommunicator( bool thisProcessParticipates, MPI_Comm & commOut, int & newRank )
{
   auto mpiManager = MPIManager::instance();

   WALBERLA_ASSERT( commOut == MPI_COMM_NULL );

   // gather a boolean for each rank, true if this process has something to send
   unsigned char participateInGather = thisProcessParticipates ? 1 : 0;

   std::vector<unsigned char> recvBuffer(uint_c( mpiManager->numProcesses() ), 0);
   MPI_Allgather ( &participateInGather, 1, MPITrait<unsigned char>::type (),
                   &(recvBuffer[0]),     1, MPITrait<unsigned char>::type (),
                    mpiManager->comm () );

   // Create a communicator out of all participating processes
   // ( i.e. processes that send more than 0 bytes )
   std::vector<int> sendingProcesses;
   for ( int process = 0; process < mpiManager->numProcesses (); ++process )
   {
      bool participates = recvBuffer[ uint_c(process) ] == 1;

      if (process == gatherRankInGlobalComm_ ) // gather rank always participates
         participates = true;

      if (participates )
         sendingProcesses.push_back ( process );
   }

   MPI_Group origGroup;
   MPI_Comm_group ( mpiManager->comm (), &origGroup );

   MPI_Group newGroup; // create newGroup consisting of sendingProcesses
   MPI_Group_incl ( origGroup, int_c ( sendingProcesses.size () ), &sendingProcesses[0], &newGroup );

   MPI_Comm_create ( mpiManager->comm (), newGroup, &commOut );

   MPI_Group_translate_ranks ( origGroup, 1, &gatherRankInGlobalComm_, newGroup, &newRank );

   MPI_Group_free ( &origGroup );
   MPI_Group_free ( &newGroup );
}


void MPIGatherScheme::runSetupPhase()
{
   auto mpiManager = MPIManager::instance();

   WALBERLA_ASSERT( ! setupPhaseDone_ );

   mpi::GenericSendBuffer<unsigned char>  sendBuffer;
   internal::packData( blocks_, packInfos_, sendBuffer );

   // Assumption: bytesToSend_ stays the same for all timesteps
   //             for enforcement of this assumption the value is stored as member
   bytesToSend_ = int_c( sendBuffer.size() );

   setupGatherCommunicator( bytesToSend_ > 0, gatherCommunicator_, gatherRank_ );

   // Gather the number of bytes to be sent to gatherRank_
   // and create displacement vector needed for MPI_Gatherv
   if ( bytesToSend_ > 0 ||  mpiManager->rank() == gatherRankInGlobalComm_ )
   {
      std::vector<decltype(bytesToSend_)> recvBuffer;

      int nrOfGatherProcesses;
      MPI_Comm_size( gatherCommunicator_, & nrOfGatherProcesses );

      if ( mpiManager->rank() == gatherRankInGlobalComm_ )
         recvBuffer.resize( uint_c( nrOfGatherProcesses ) );

      MPI_Gather( & bytesToSend_,   1, MPITrait<decltype(bytesToSend_)>::type(),
                  recvBuffer.empty()? nullptr : & recvBuffer[0],  1, MPITrait<decltype(bytesToSend_)>::type(),
                  gatherRank_, gatherCommunicator_  );

      WALBERLA_ASSERT_EQUAL( displacementVector_.size(), 0 );
      sendBytesPerProcess_.resize( uint_c( nrOfGatherProcesses ) );
      if ( mpiManager->rank() == gatherRankInGlobalComm_ )
      {
         int sum = 0;
         for( uint_t i = 0; i < uint_c( nrOfGatherProcesses ); ++i )
         {
            displacementVector_.push_back( sum );
            sendBytesPerProcess_[i] = recvBuffer[i];
            sum += int_c ( recvBuffer[i] );
         }
         gatherMsgSize_ = sum;
      }
   }

   setupPhaseDone_ = true;
}




void MPIGatherScheme::communicate()
{
   WALBERLA_NON_MPI_SECTION() {
      return; // do nothing when build without MPI
   }

   auto mpiManager = MPIManager::instance();

   if ( packInfos_.empty() ) {
      WALBERLA_LOG_WARNING( "Communicating MPIGatherScheme without registered PackInfos");
      return;
   }

   if ( !setupPhaseDone_ )
      runSetupPhase( );

   // Process does not participate in gather communication
   if ( gatherCommunicator_ == MPI_COMM_NULL )
      return;

   mpi::GenericSendBuffer<unsigned char>  sendBuffer;
   internal::packData( blocks_, packInfos_, sendBuffer );

   int myRankInGatherComm;
   MPI_Comm_rank( gatherCommunicator_, &myRankInGatherComm );

   // This assert ensures that all registered PackInfos pack the same amount of data in each timestep
   // If this is not the case this Scheme cannot be used!
   WALBERLA_ASSERT_EQUAL( sendBuffer.size(), uint_c( bytesToSend_ ) );

   mpi::GenericRecvBuffer<unsigned char>  recvBuffer;

   int * displacementVectorPtr = nullptr;
   int * sendBytesPerProcessPtr = nullptr;
   if ( mpiManager->rank() == gatherRankInGlobalComm_  ) {
      recvBuffer.resize( uint_c( gatherMsgSize_ ) );
      displacementVectorPtr  = &displacementVector_[0];
      sendBytesPerProcessPtr = &sendBytesPerProcess_[0];
   }

   MPI_Gatherv( sendBuffer.ptr(), int_c( sendBuffer.size() ), MPITrait<unsigned char>::type(),
                recvBuffer.ptr(), sendBytesPerProcessPtr , displacementVectorPtr, MPITrait<unsigned char>::type() ,
                gatherRank_, gatherCommunicator_ );


   // Unpacking
   if ( mpiManager->rank() == gatherRankInGlobalComm_  )
   {
      internal::unpackData( packInfos_, recvBuffer );

      for( size_t  s =0; s < packInfos_.size() ; ++s )
         packInfos_[s]->gatherFinished( );

   }
}


} // namespace gather
} // namespace walberla

