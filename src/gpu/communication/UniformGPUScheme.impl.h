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
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "gpu/ParallelStreams.h"

namespace walberla {
namespace gpu
{
namespace communication {


   template<typename Stencil>
   UniformGPUScheme<Stencil>::UniformGPUScheme( weak_ptr <StructuredBlockForest> bf,
                                                bool sendDirectlyFromGPU,
                                                bool useLocalCommunication,
                                                const int tag )
        : blockForest_( bf ),
          setupBeforeNextCommunication_( true ),
          communicationInProgress_( false ),
          sendFromGPU_( sendDirectlyFromGPU ),
          useLocalCommunication_(useLocalCommunication),
          bufferSystemCPU_( mpi::MPIManager::instance()->comm(), tag ),
          bufferSystemGPU_( mpi::MPIManager::instance()->comm(), tag ),
          parallelSectionManager_( -1 ),
          requiredBlockSelectors_( Set<SUID>::emptySet() ),
          incompatibleBlockSelectors_( Set<SUID>::emptySet() )
   {
      WALBERLA_MPI_SECTION()
      {
#if !(defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT)
         WALBERLA_CHECK(!sendDirectlyFromGPU)
#endif
      }
   }

   template<typename Stencil>
   UniformGPUScheme<Stencil>::UniformGPUScheme( weak_ptr <StructuredBlockForest> bf,
                                                const Set<SUID> & requiredBlockSelectors,
                                                const Set<SUID> & incompatibleBlockSelectors,
                                                bool sendDirectlyFromGPU,
                                                bool useLocalCommunication,
                                                const int tag )
      : blockForest_( bf ),
        setupBeforeNextCommunication_( true ),
        communicationInProgress_( false ),
        sendFromGPU_( sendDirectlyFromGPU ),
        useLocalCommunication_(useLocalCommunication),
        bufferSystemCPU_( mpi::MPIManager::instance()->comm(), tag ),
        bufferSystemGPU_( mpi::MPIManager::instance()->comm(), tag ),
        parallelSectionManager_( -1 ),
        requiredBlockSelectors_( requiredBlockSelectors ),
        incompatibleBlockSelectors_( incompatibleBlockSelectors )
   {
      WALBERLA_MPI_SECTION()
      {
#if !(defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT)
         WALBERLA_CHECK(!sendDirectlyFromGPU)
#endif
      }
   }


   template<typename Stencil>
   void UniformGPUScheme<Stencil>::startCommunication( gpuStream_t stream )
   {
      WALBERLA_ASSERT( !communicationInProgress_ )
      auto forest = blockForest_.lock();

      auto currentBlockForestStamp = forest->getBlockForest().getModificationStamp();
      if( setupBeforeNextCommunication_ || currentBlockForestStamp != forestModificationStamp_ )
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
      {
         auto parallelSection = parallelSectionManager_.parallelSection( stream );
         for( auto &iBlock : *forest )
         {
            auto senderBlock = dynamic_cast< Block * >( &iBlock );

            if( !selectable::isSetSelected( senderBlock->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
               continue;

            for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
            {
               const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );

               if( senderBlock->getNeighborhoodSectionSize( neighborIdx ) == uint_t( 0 ))
                  continue;

               if( !selectable::isSetSelected( senderBlock->getNeighborState( neighborIdx, uint_t(0) ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
                  continue;

               if( senderBlock->neighborExistsLocally( neighborIdx, uint_t(0) ) && useLocalCommunication_ )
               {
                  auto receiverBlock = dynamic_cast< Block * >( forest->getBlock( senderBlock->getNeighborId( neighborIdx, uint_t(0) )) );
                  for (auto& pi : packInfos_)
                  {
                     pi->communicateLocal(*dir, senderBlock, receiverBlock, stream);
                  }
               }
               else
               {
                  auto nProcess = mpi::MPIRank( senderBlock->getNeighborProcess( neighborIdx, uint_t( 0 )));

                  for( auto &pi : packInfos_ )
                  {
                     parallelSection.run([&](auto s) {
                     auto size = pi->size( *dir, senderBlock );
                     auto gpuDataPtr = bufferSystemGPU_.sendBuffer( nProcess ).advanceNoResize( size );
                     WALBERLA_ASSERT_NOT_NULLPTR( gpuDataPtr )
                     pi->pack( *dir, gpuDataPtr, senderBlock, s );

                     if( !sendFromGPU_ )
                     {
                        auto cpuDataPtr = bufferSystemCPU_.sendBuffer( nProcess ).advanceNoResize( size );
                        WALBERLA_ASSERT_NOT_NULLPTR( cpuDataPtr )
                        WALBERLA_GPU_CHECK( gpuMemcpyAsync( cpuDataPtr, gpuDataPtr, size, gpuMemcpyDeviceToHost, s ))
                     }
                     });
                  }
               }
            }
         }
      }
      // wait for packing to finish
      WALBERLA_GPU_CHECK( gpuStreamSynchronize( stream ) );

      if( sendFromGPU_ )
         bufferSystemGPU_.sendAll();
      else
         bufferSystemCPU_.sendAll();

      communicationInProgress_ = true;
   }


   template<typename Stencil>
   void UniformGPUScheme<Stencil>::wait( gpuStream_t stream )
   {
      WALBERLA_ASSERT( communicationInProgress_ )

      auto forest = blockForest_.lock();

      if( sendFromGPU_ )
      {
         auto parallelSection = parallelSectionManager_.parallelSection( stream );
         for( auto recvInfo = bufferSystemGPU_.begin(); recvInfo != bufferSystemGPU_.end(); ++recvInfo )
         {
            recvInfo.buffer().clear();
            for( auto &header : headers_[recvInfo.rank()] )
            {
               auto block = dynamic_cast< Block * >( forest->getBlock( header.blockId ));

               for( auto &pi : packInfos_ )
               {
                  auto size = pi->size( header.dir, block );
                  auto gpuDataPtr = recvInfo.buffer().advanceNoResize( size );
                  WALBERLA_ASSERT_NOT_NULLPTR( gpuDataPtr )
                  parallelSection.run([&](auto s) {
                     pi->unpack( stencil::inverseDir[header.dir], gpuDataPtr, block, s );
                  });
               }
            }
         }
      }
      else
      {
         auto parallelSection = parallelSectionManager_.parallelSection( stream );
         for( auto recvInfo = bufferSystemCPU_.begin(); recvInfo != bufferSystemCPU_.end(); ++recvInfo )
         {
            auto &gpuBuffer = bufferSystemGPU_.sendBuffer( recvInfo.rank());

            recvInfo.buffer().clear();
            gpuBuffer.clear();
            for( auto &header : headers_[recvInfo.rank()] ) {
               auto block = dynamic_cast< Block * >( forest->getBlock( header.blockId ));

               for( auto &pi : packInfos_ )
               {
                  auto size = pi->size( header.dir, block );
                  auto cpuDataPtr = recvInfo.buffer().advanceNoResize( size );
                  auto gpuDataPtr = gpuBuffer.advanceNoResize( size );
                  WALBERLA_ASSERT_NOT_NULLPTR( cpuDataPtr )
                  WALBERLA_ASSERT_NOT_NULLPTR( gpuDataPtr )
                  parallelSection.run([&](auto s) {
                     WALBERLA_GPU_CHECK( gpuMemcpyAsync( gpuDataPtr, cpuDataPtr, size,
                                                           gpuMemcpyHostToDevice, s ))
                     pi->unpack( stencil::inverseDir[header.dir], gpuDataPtr, block, s );
                  });
               }
            }
         }
      }

      WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
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

         if( !selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
            continue;

         for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir ) {
            // skip if block has no neighbors in this direction
            const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );

            if( block->getNeighborhoodSectionSize( neighborIdx ) == uint_t( 0 ))
               continue;

            WALBERLA_ASSERT( block->neighborhoodSectionHasEquallySizedBlock( neighborIdx ),
                             "Works for uniform setups only" )
            WALBERLA_ASSERT_EQUAL( block->getNeighborhoodSectionSize( neighborIdx ), uint_t( 1 ),
                                   "Works for uniform setups only" )

            const BlockID &nBlockId = block->getNeighborId( neighborIdx, uint_t( 0 ));

            if( !selectable::isSetSelected( block->getNeighborState( neighborIdx, uint_t(0) ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
               continue;

            if( block->neighborExistsLocally( neighborIdx, uint_t(0) ) && useLocalCommunication_ )
               continue;

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
         bufferSystemCPU_.sendBuffer( it.first ).resize( size_t(it.second) );
         bufferSystemGPU_.sendBuffer( it.first ).resize( size_t(it.second) );
      }

      forestModificationStamp_ = forest->getBlockForest().getModificationStamp();
      setupBeforeNextCommunication_ = false;
   }


   template<typename Stencil>
   void UniformGPUScheme<Stencil>::addPackInfo( const shared_ptr<GeneratedGPUPackInfo> &pi )
   {
      WALBERLA_ASSERT( !communicationInProgress_, "Cannot add pack info while communication is in progress" )
      packInfos_.push_back( pi );
      setupBeforeNextCommunication_ = true;
   }

   template< typename Stencil >
   std::function<void()> UniformGPUScheme<Stencil>::getCommunicateFunctor(gpuStream_t stream)
   {
      return [this, stream]() { communicate( stream ); };
   }

   template< typename Stencil >
   std::function<void()> UniformGPUScheme<Stencil>::getStartCommunicateFunctor(gpuStream_t stream)
   {
      return [this, stream]() { startCommunication( stream ); };
   }

   template< typename Stencil >
   std::function<void()> UniformGPUScheme<Stencil>::getWaitFunctor(gpuStream_t stream)
   {
      return [this, stream]() { wait( stream ); };
   }

} // namespace communication
} // namespace gpu
} // namespace walberla
