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
//! \file OpenMPBufferSystem.impl.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


namespace walberla {
namespace mpi {



//======================================================================================================================
//
//  Setup
//
//======================================================================================================================

template<typename Rb, typename Sb>
GenericOpenMPBufferSystem<Rb, Sb>::GenericOpenMPBufferSystem( const MPI_Comm & communicator, int tag,
                                                              bool _serialSends, bool _serialRecvs  )
   : bs_( communicator, tag),
     dirty_( true ),
     serialSends_( _serialSends ),
     serialRecvs_( _serialRecvs ),
     sizeChangesEverytime_( true )
{
}


template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::addReceivingFunction( MPIRank rank, const std::function<void ( Rb & buf ) >& recvFunction )
{
   dirty_ = true;
   recvFunctions_[rank] = recvFunction;
}


template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::addSendingFunction( MPIRank rank, const std::function<void ( Sb & buf ) >& sendFunction )
{
   dirty_ = true;
   sendRanks_.push_back( rank );
   sendFunctions_.push_back( sendFunction );
   WALBERLA_ASSERT_EQUAL( sendRanks_.size(), sendFunctions_.size() )
}


template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::setupBufferSystem()
{
   if ( ! dirty_ )
      return;

   std::set<MPIRank> recvRanks;
   std::transform( recvFunctions_.begin(), recvFunctions_.end(), std::inserter(recvRanks, recvRanks.end()),
                   []( typename decltype(recvFunctions_)::value_type r ){ return r.first; });
   bs_.setReceiverInfo( recvRanks, sizeChangesEverytime_ );

   for(auto sendRank : sendRanks_) // Do NOT delete this for loop! This loop is needed ...
      bs_.sendBuffer( sendRank ); // ... so that the "sendBuffer(rank)" call in startCommunicationOpenMP is thread-safe!

   dirty_ = false;
}


//======================================================================================================================
//
//  Sending
//
//======================================================================================================================


template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::startCommunication()
{
   setupBufferSystem();
   if( serialSends_ )
      startCommunicationSerial();
   else
      startCommunicationOpenMP();
}

template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::startCommunicationSerial()
{
   bs_.scheduleReceives();

   WALBERLA_ASSERT_EQUAL( sendRanks_.size(), sendFunctions_.size() )

   const uint_t nrOfSendFunctions = sendFunctions_.size();

   for( uint_t i = 0; i < nrOfSendFunctions; ++i )
   {
      const MPIRank rank = sendRanks_[i];
      SendBuffer & sendBuffer = bs_.sendBuffer( rank );

      sendFunctions_[i]( sendBuffer );
      bs_.send( rank );
   }

   bs_.sendAll(); // for the case where sendFunctions_ is empty
}

template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::startCommunicationOpenMP()
{
   bs_.scheduleReceives();

   WALBERLA_ASSERT_EQUAL( sendRanks_.size(), sendFunctions_.size() )

   const int nrOfSendFunctions = int_c( sendFunctions_.size() );

   #ifdef _OPENMP
   #pragma omp parallel for schedule(dynamic)
   #endif
   for( int i=0; i < nrOfSendFunctions; ++i )
   {
      const MPIRank rank = sendRanks_[ uint_c(i) ];
      SendBuffer & sendBuffer = bs_.sendBuffer( rank ); // This call is thread-safe since all send buffers
                                                        // are already allocated in setupBufferSystem().
      sendFunctions_[ uint_c(i) ]( sendBuffer );

      #ifdef _OPENMP
      #pragma omp critical
      #endif
      bs_.send( rank );
   }

   bs_.sendAll(); // for the case where sendFunctions_ is empty
}




//======================================================================================================================
//
//  Receiving
//
//======================================================================================================================

template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::wait()
{
   if ( serialRecvs_ )
      waitSerial();
   else
      waitOpenMP();
}


template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::waitSerial()
{
   for( auto recvIt = bs_.begin(); recvIt != bs_.end(); ++recvIt )
   {
      const MPIRank rank = recvIt.rank();
      RecvBuffer & recvBuffer = recvIt.buffer();

      // call unpacking
      WALBERLA_ASSERT( recvFunctions_.find( rank ) != recvFunctions_.end() )
      recvFunctions_[rank] ( recvBuffer );
   }
}


template<typename Rb, typename Sb>
void GenericOpenMPBufferSystem<Rb, Sb>::waitOpenMP()
{
   const int numReceives = int_c( bs_.recvInfos_.size() );

   #ifdef _OPENMP
   #pragma omp parallel for schedule(dynamic)
   #endif
   for( int i = 0; i < numReceives; ++i )
   {
      MPIRank recvRank = INVALID_RANK;
      RecvBuffer * recvBuffer = nullptr;

      #ifdef _OPENMP
      #pragma omp critical
      #endif
      recvBuffer = bs_.waitForNext( recvRank );


      WALBERLA_ASSERT_GREATER_EQUAL( recvRank, 0 )
      WALBERLA_ASSERT_NOT_NULLPTR( recvBuffer )

      WALBERLA_ASSERT( recvFunctions_.find( recvRank ) != recvFunctions_.end() )
      recvFunctions_[recvRank] ( *recvBuffer );
   }

   MPIRank rank;
   RecvBuffer * ret = bs_.waitForNext( rank );
   WALBERLA_ASSERT_NULLPTR( ret ) // call last time to finish communication
   WALBERLA_UNUSED( ret );

   WALBERLA_ASSERT( ! bs_.isCommunicationRunning() )
}



} // namespace mpi
} // namespace walberla
