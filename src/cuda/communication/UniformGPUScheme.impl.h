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
//! \file UniformGPUScheme.impl.h
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


namespace walberla {
namespace cuda {
namespace communication {


template<typename Stencil>
UniformGPUScheme<Stencil>::UniformGPUScheme( weak_ptr_wrapper <StructuredBlockForest> bf,
                                             const shared_ptr <cuda::EventRAII> &startWaitEvent,
                                             bool sendDirectlyFromGPU,
                                             const int tag )
        : blockForest_( bf ),
          startWaitEvent_( startWaitEvent ),
          setupBeforeNextCommunication_( true ),
          communicationInProgress_( false ),
          sendFromGPU_( sendDirectlyFromGPU ),
          bufferSystemCPU_( mpi::MPIManager::instance()->comm(), tag ),
          bufferSystemGPU_( mpi::MPIManager::instance()->comm(), tag ) {}


   template<typename Stencil>
   void UniformGPUScheme<Stencil>::startCommunication()
   {
      WALBERLA_ASSERT( !communicationInProgress_ );
      auto forest = blockForest_.lock();

      if( setupBeforeNextCommunication_ ||
          forest->getBlockForest().getModificationStamp() != forestModificationStamp_ )
         setupCommunication();

      // Schedule Receives
      if( sendFromGPU_ )
         bufferSystemGPU_.scheduleReceives();
      else
         bufferSystemCPU_.scheduleReceives();


      if( !sendFromGPU_ )
         for( auto it : headers_ )
            bufferSystemGPU_.sendBuffer( it.first ).clear();

      // Start filling send buffers
      for( auto &iBlock : *forest )
      {
         auto block = dynamic_cast< Block * >( &iBlock );
         for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
         {
            const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );
            if( block->getNeighborhoodSectionSize( neighborIdx ) == uint_t( 0 ))
               continue;
            auto nProcess = mpi::MPIRank( block->getNeighborProcess( neighborIdx, uint_t( 0 )));

            if( streams_.find( *dir ) == streams_.end() )
            {
               streams_.emplace( *dir, StreamRAII::newPriorityStream( -1 ));
            }

            auto &ci = streams_.at( *dir );

            for( auto &pi : packInfos_ )
            {
               auto size = pi->size( *dir, block );
               auto gpuDataPtr = bufferSystemGPU_.sendBuffer( nProcess ).advanceNoResize( size );
               WALBERLA_ASSERT_NOT_NULLPTR( gpuDataPtr );
               WALBERLA_CUDA_CHECK( cudaStreamWaitEvent( ci, *startWaitEvent_, 0 ));
               pi->pack( *dir, gpuDataPtr, block, ci );

               if( !sendFromGPU_ )
               {
                  auto cpuDataPtr = bufferSystemCPU_.sendBuffer( nProcess ).advanceNoResize( size );
                  WALBERLA_ASSERT_NOT_NULLPTR( cpuDataPtr );
                  WALBERLA_CUDA_CHECK(
                          cudaMemcpyAsync( cpuDataPtr, gpuDataPtr, size, cudaMemcpyDeviceToHost, ci ));
               }
            }
         }
      }

      // Busy waiting for packing to finish - then send
      for( auto &ci : streams_ ) WALBERLA_CUDA_CHECK( cudaStreamSynchronize( ci.second ));

      if( sendFromGPU_ )
         bufferSystemGPU_.sendAll();
      else
         bufferSystemCPU_.sendAll();

      communicationInProgress_ = true;
   }


   template<typename Stencil>
   void UniformGPUScheme<Stencil>::wait()
   {
      WALBERLA_ASSERT( communicationInProgress_ );

      auto forest = blockForest_.lock();

      if( sendFromGPU_ )
      {
         for( auto recvInfo = bufferSystemGPU_.begin(); recvInfo != bufferSystemGPU_.end(); ++recvInfo )
         {
            for( auto &header : headers_[recvInfo.rank()] )
            {
               auto &ci = streams_.at( header.dir );
               auto block = dynamic_cast< Block * >( forest->getBlock( header.blockId ));

               for( auto &pi : packInfos_ )
               {
                  auto size = pi->size( header.dir, block );
                  auto gpuDataPtr = recvInfo.buffer().advanceNoResize( size );
                  WALBERLA_ASSERT_NOT_NULLPTR( gpuDataPtr );
                  pi->unpack( stencil::inverseDir[header.dir], gpuDataPtr, block, ci );
               }
            }
         }
      }
      else
      {
         for( auto recvInfo = bufferSystemCPU_.begin(); recvInfo != bufferSystemCPU_.end(); ++recvInfo )
         {
            using namespace std::chrono_literals;
            std::this_thread::sleep_for( 1ms );
            auto &gpuBuffer = bufferSystemGPU_.sendBuffer( recvInfo.rank());

            gpuBuffer.clear();
            for( auto &header : headers_[recvInfo.rank()] ) {
               auto &ci = streams_.at( header.dir );
               auto block = dynamic_cast< Block * >( forest->getBlock( header.blockId ));

               for( auto &pi : packInfos_ )
               {
                  auto size = pi->size( header.dir, block );
                  auto cpuDataPtr = recvInfo.buffer().advanceNoResize( size );
                  auto gpuDataPtr = gpuBuffer.advanceNoResize( size );
                  WALBERLA_ASSERT_NOT_NULLPTR( cpuDataPtr );
                  WALBERLA_ASSERT_NOT_NULLPTR( gpuDataPtr );

                  WALBERLA_CUDA_CHECK( cudaMemcpyAsync( gpuDataPtr, cpuDataPtr, size,
                                                        cudaMemcpyHostToDevice, ci ));
                  pi->unpack( stencil::inverseDir[header.dir], gpuDataPtr, block, ci );
               }
            }
         }
      }

      for( auto &ci : streams_ )
         WALBERLA_CUDA_CHECK( cudaStreamSynchronize( ci.second ));

      communicationInProgress_ = false;
   }


   template<typename Stencil>
   void UniformGPUScheme<Stencil>::setupCommunication()
   {
      auto forest = blockForest_.lock();

      headers_.clear();

      std::map<mpi::MPIRank, mpi::MPISize> receiverInfo; // how many bytes to send to each neighbor

      mpi::BufferSystem headerExchangeBs( mpi::MPIManager::instance()->comm(), 123 );

      for( auto &iBlock : *forest ) {
         auto block = dynamic_cast< Block * >( &iBlock );

         for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir ) {
            // skip if block has no neighbors in this direction
            const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );
            if( block->getNeighborhoodSectionSize( neighborIdx ) == uint_t( 0 ))
               continue;

            WALBERLA_ASSERT( block->neighborhoodSectionHasEquallySizedBlock( neighborIdx ),
                             "Works for uniform setups only" );
            WALBERLA_ASSERT_EQUAL( block->getNeighborhoodSectionSize( neighborIdx ), uint_t( 1 ),
                                   "Works for uniform setups only" );

            const BlockID &nBlockId = block->getNeighborId( neighborIdx, uint_t( 0 ));
            auto nProcess = mpi::MPIRank( block->getNeighborProcess( neighborIdx, uint_t( 0 )));

            for( auto &pi : packInfos_ )
               receiverInfo[nProcess] += mpi::MPISize( pi->size( *dir, block ));

            auto &headerBuffer = headerExchangeBs.sendBuffer( nProcess );
            nBlockId.toBuffer( headerBuffer );
            headerBuffer << *dir;
         }
      }

      headerExchangeBs.setReceiverInfoFromSendBufferState( false, true );
      headerExchangeBs.sendAll();
      for( auto recvIter = headerExchangeBs.begin(); recvIter != headerExchangeBs.end(); ++recvIter ) {
         auto &headerVector = headers_[recvIter.rank()];
         auto &buffer = recvIter.buffer();
         while ( buffer.size()) {
            Header header;
            header.blockId.fromBuffer( buffer );
            buffer >> header.dir;
            headerVector.push_back( header );
         }
      }

      bufferSystemCPU_.setReceiverInfo( receiverInfo );
      bufferSystemGPU_.setReceiverInfo( receiverInfo );

      for( auto it : receiverInfo ) {
         bufferSystemCPU_.sendBuffer( it.first ).resize( it.second );
         bufferSystemGPU_.sendBuffer( it.first ).resize( it.second );
      }

      forestModificationStamp_ = forest->getBlockForest().getModificationStamp();
      setupBeforeNextCommunication_ = false;
   }


   template<typename Stencil>
   void UniformGPUScheme<Stencil>::addPackInfo( const shared_ptr<GeneratedGPUPackInfo> &pi )
   {
      WALBERLA_ASSERT( !communicationInProgress_, "Cannot add pack info while communication is in progress" );
      packInfos_.push_back( pi );
      setupBeforeNextCommunication_ = true;
   }


} // namespace communication
} // namespace cuda
} // namespace walberla
