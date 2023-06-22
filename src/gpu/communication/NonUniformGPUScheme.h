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
//! \file NonUniformGPUScheme.h
//! \ingroup gpu
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIWrapper.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include <thread>

#include "gpu/ErrorChecking.h"
#include "gpu/GPURAII.h"
#include "gpu/GPUWrapper.h"
#include "gpu/ParallelStreams.h"
#include "gpu/communication/CustomMemoryBuffer.h"
#include "gpu/communication/GeneratedNonUniformGPUPackInfo.h"

namespace walberla::gpu::communication
{

template< typename Stencil >
class NonUniformGPUScheme
{
 public:
   enum INDEX { EQUAL_LEVEL = 0, COARSE_TO_FINE = 1, FINE_TO_COARSE = 2 };

   using CpuBuffer_T = walberla::gpu::communication::PinnedMemoryBuffer;
   using GpuBuffer_T = walberla::gpu::communication::GPUMemoryBuffer;

   explicit NonUniformGPUScheme(const weak_ptr< StructuredBlockForest >& bf, bool sendDirectlyFromGPU = false,
                                const int tag = 5432);

   explicit NonUniformGPUScheme(const weak_ptr< StructuredBlockForest >& bf, const Set< SUID >& requiredBlockSelectors,
                                const Set< SUID >& incompatibleBlockSelectors, bool sendDirectlyFromGPU = false,
                                const int tag = 5432);

   ~NonUniformGPUScheme();

   //** Pack Info Registration *****************************************************************************************
   /*! \name Pack Info Registration */
   //@{
   void addPackInfo(const shared_ptr< GeneratedNonUniformGPUPackInfo >& pi);
   //@}
   //*******************************************************************************************************************

   inline void communicateEqualLevel(const uint_t level);
   inline void communicateCoarseToFine(const uint_t fineLevel);
   inline void communicateFineToCoarse(const uint_t fineLevel);

   std::function<void()>  communicateEqualLevelFunctor(const uint_t level) {
      return [level, this](){ NonUniformGPUScheme::communicateEqualLevel(level);};
   }
   std::function<void()>  communicateCoarseToFineFunctor(const uint_t fineLevel) {
      return [fineLevel, this](){ NonUniformGPUScheme::communicateCoarseToFine(fineLevel);};
   }
   std::function<void()>  communicateFineToCoarseFunctor(const uint_t fineLevel) {
      return [fineLevel, this](){ NonUniformGPUScheme::communicateFineToCoarse(fineLevel);};
   }

   inline void startCommunicateEqualLevel(const uint_t level);
   inline void startCommunicateCoarseToFine(const uint_t fineLevel);
   inline void startCommunicateFineToCoarse(const uint_t fineLevel);

   inline void waitCommunicateEqualLevel(const uint_t level);
   inline void waitCommunicateCoarseToFine(const uint_t fineLevel);
   inline void waitCommunicateFineToCoarse(const uint_t fineLevel);

 private:
   void setupCommunication();

   void init();
   void refresh();

   bool isAnyCommunicationInProgress() const;

   void startCommunicationEqualLevel(const uint_t index, std::set< uint_t >& participatingLevels);
   void startCommunicationCoarseToFine(const uint_t index, const uint_t coarsestLevel);
   void startCommunicationFineToCoarse(const uint_t index, const uint_t finestLevel);

   weak_ptr< StructuredBlockForest > blockForest_;
   uint_t forestModificationStamp_;

   std::vector< std::vector< bool > > communicationInProgress_;
   bool sendFromGPU_;
   int baseTag_;

   std::vector< std::vector< mpi::GenericBufferSystem< CpuBuffer_T, CpuBuffer_T > > > bufferSystemCPU_;
   std::vector< std::vector< mpi::GenericBufferSystem< GpuBuffer_T, GpuBuffer_T > > > bufferSystemGPU_;
   std::vector< std::vector< GpuBuffer_T > > localBuffer_;

   std::vector< shared_ptr< GeneratedNonUniformGPUPackInfo > > packInfos_;

   ParallelStreams parallelSectionManager_;

   struct Header
   {
      BlockID receiverId;
      BlockID senderId;
      stencil::Direction dir;
   };
   std::vector< std::vector< std::map< mpi::MPIRank, std::vector< Header > > > > headers_;

   Set< SUID > requiredBlockSelectors_;
   Set< SUID > incompatibleBlockSelectors_;
};

template< typename Stencil >
NonUniformGPUScheme< Stencil >::NonUniformGPUScheme(const weak_ptr< StructuredBlockForest >& bf, bool sendDirectlyFromGPU,
                                                    const int tag)
   : blockForest_(bf), sendFromGPU_(sendDirectlyFromGPU), baseTag_(tag), parallelSectionManager_(-1),
     requiredBlockSelectors_(Set< SUID >::emptySet()), incompatibleBlockSelectors_(Set< SUID >::emptySet())
{
   WALBERLA_MPI_SECTION()
   {
#if !(defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT)
      WALBERLA_CHECK(!sendDirectlyFromGPU)
#endif
   }
   init();
}

template< typename Stencil >
NonUniformGPUScheme< Stencil >::NonUniformGPUScheme(const weak_ptr< StructuredBlockForest >& bf,
                                                    const Set< SUID >& requiredBlockSelectors,
                                                    const Set< SUID >& incompatibleBlockSelectors,
                                                    bool sendDirectlyFromGPU, const int tag)
   : blockForest_(bf), requiredBlockSelectors_(requiredBlockSelectors),
     incompatibleBlockSelectors_(incompatibleBlockSelectors), sendFromGPU_(sendDirectlyFromGPU), baseTag_(tag),
     parallelSectionManager_(-1)
{
   WALBERLA_MPI_SECTION()
   {
#if !(defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT)
      WALBERLA_CHECK(!sendDirectlyFromGPU)
#endif
   }
   init();
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::init()
{
   bufferSystemCPU_.resize(3);
   bufferSystemGPU_.resize(3);
   localBuffer_.resize(3);
   headers_.resize(3);

   communicationInProgress_.resize(3);

   refresh();
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::refresh()
{
   WALBERLA_ASSERT(!isAnyCommunicationInProgress())

   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest,
                              "Trying to access communication for a block storage object that doesn't exist anymore")
   const uint_t levels = forest->getNumberOfLevels();

   for (uint_t i = 0; i != 3; ++i)
   {
      bufferSystemCPU_[i].clear();
      bufferSystemGPU_[i].clear();
      localBuffer_[i].clear();
      headers_[i].clear();
      headers_[i].resize(size_t(levels + uint_t(1)));

      for (uint_t j = 0; j <= levels; ++j)
      {
         headers_[i][j].clear();
         bufferSystemCPU_[i].emplace_back(
            mpi::MPIManager::instance()->comm(), baseTag_ + int_c(i * levels + j));
         bufferSystemGPU_[i].emplace_back(
            mpi::MPIManager::instance()->comm(), baseTag_ + int_c(i * levels + j));
         localBuffer_[i].emplace_back();
      }

      communicationInProgress_[i].resize(size_t(levels + uint_t(1)), false);
   }

#ifndef NDEBUG
   for (auto p = packInfos_.begin(); p != packInfos_.end(); ++p)
      (*p)->clearBufferSizeCheckMap();
#endif

   forestModificationStamp_ = forest->getBlockForest().getModificationStamp();
}

template< typename Stencil >
inline void NonUniformGPUScheme< Stencil >::communicateEqualLevel(const uint_t level)
{
   startCommunicateEqualLevel(level);
   waitCommunicateEqualLevel(level);
}

template< typename Stencil >
inline void NonUniformGPUScheme< Stencil >::communicateCoarseToFine(const uint_t fineLevel)
{
   startCommunicateCoarseToFine(fineLevel);
   waitCommunicateCoarseToFine(fineLevel);
}

template< typename Stencil >
inline void NonUniformGPUScheme< Stencil >::communicateFineToCoarse(const uint_t fineLevel)
{
   startCommunicateFineToCoarse(fineLevel);
   waitCommunicateFineToCoarse(fineLevel);
}

template< typename Stencil >
inline void NonUniformGPUScheme< Stencil >::startCommunicateEqualLevel(const uint_t level)
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest,
                              "Trying to access communication for a block storage object that doesn't exist anymore")
   WALBERLA_ASSERT_LESS(level, forest->getNumberOfLevels())

   if (forestModificationStamp_ != forest->getBlockForest().getModificationStamp()) refresh();

   std::set< uint_t > participatingLevels;
   participatingLevels.insert(level);

   startCommunicationEqualLevel(level, participatingLevels);
}

template< typename Stencil >
inline void NonUniformGPUScheme< Stencil >::startCommunicateCoarseToFine(const uint_t fineLevel)
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest,
                              "Trying to access communication for a block storage object that doesn't exist anymore")
   WALBERLA_ASSERT_GREATER(fineLevel, uint_t(0))
   WALBERLA_ASSERT_LESS(fineLevel, forest->getNumberOfLevels())

   if (forestModificationStamp_ != forest->getBlockForest().getModificationStamp()) refresh();

   const uint_t coarsestLevel = fineLevel - uint_t(1);

   startCommunicationCoarseToFine(fineLevel, coarsestLevel);
}

template< typename Stencil >
inline void NonUniformGPUScheme< Stencil >::startCommunicateFineToCoarse(const uint_t fineLevel)
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest,
                              "Trying to access communication for a block storage object that doesn't exist anymore")
   WALBERLA_ASSERT_GREATER(fineLevel, uint_t(0))
   WALBERLA_ASSERT_LESS(fineLevel, forest->getNumberOfLevels())

   if (forestModificationStamp_ != forest->getBlockForest().getModificationStamp()) refresh();

   const uint_t finestLevel   = fineLevel;

   startCommunicationFineToCoarse(fineLevel, finestLevel);
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::startCommunicationEqualLevel(const uint_t index,
                                                                  std::set< uint_t >& participatingLevels)
{
   if (packInfos_.empty()) return;

   WALBERLA_ASSERT(!communicationInProgress_[EQUAL_LEVEL][index])
   communicationInProgress_[EQUAL_LEVEL][index] = true;

   auto forest = blockForest_.lock();

   // Schedule Receives
   if (sendFromGPU_)
      bufferSystemGPU_[EQUAL_LEVEL][index].scheduleReceives();
   else
      bufferSystemCPU_[EQUAL_LEVEL][index].scheduleReceives();

   if (!sendFromGPU_)
      for (auto it : headers_[EQUAL_LEVEL][index])
         bufferSystemGPU_[EQUAL_LEVEL][index].sendBuffer(it.first).clear();

   // Start filling send buffers
   {
      for (auto& iBlock : *forest)
      {
         auto senderBlock = dynamic_cast< Block* >(&iBlock);

         if (!selectable::isSetSelected(senderBlock->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_))
            continue;

         if (participatingLevels.find(senderBlock->getLevel()) == participatingLevels.end())
            continue;

         for (auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir)
         {
            const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex(*dir);

            if (!(senderBlock->neighborhoodSectionHasEquallySizedBlock(neighborIdx)))
               continue;
            WALBERLA_ASSERT_EQUAL(senderBlock->getNeighborhoodSectionSize(neighborIdx), uint_t(1))
            if (!selectable::isSetSelected(senderBlock->getNeighborState(neighborIdx, uint_t(0)),requiredBlockSelectors_, incompatibleBlockSelectors_))
               continue;

            if( senderBlock->neighborExistsLocally( neighborIdx, uint_t(0) ) )
            {
               auto receiverBlock = dynamic_cast< Block * >( forest->getBlock( senderBlock->getNeighborId( neighborIdx, uint_t(0) )) );
               for (auto& pi : packInfos_)
               {
                  pi->communicateLocalEqualLevel(senderBlock, receiverBlock, *dir, nullptr);
               }
            }
            else
            {
               auto nProcess              = mpi::MPIRank(senderBlock->getNeighborProcess(neighborIdx, uint_t(0)));
               GpuBuffer_T& gpuDataBuffer = bufferSystemGPU_[EQUAL_LEVEL][index].sendBuffer(nProcess);

               for (auto& pi : packInfos_)
               {
                  WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
                  WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.allocSize() - gpuDataBuffer.size(), pi->sizeEqualLevelSend(senderBlock, *dir))

                  pi->packDataEqualLevel(senderBlock, *dir, gpuDataBuffer);

                  if (!sendFromGPU_)
                  {
                     auto gpuDataPtr = gpuDataBuffer.cur();
                     auto size = pi->sizeEqualLevelSend(senderBlock, *dir);
                     auto cpuDataPtr = bufferSystemCPU_[EQUAL_LEVEL][index].sendBuffer(nProcess).advanceNoResize(size);
                     WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
                     WALBERLA_GPU_CHECK(gpuMemcpyAsync(cpuDataPtr, gpuDataPtr, size, gpuMemcpyDeviceToHost))
                  }
               }
            }
         }
      }
   }

   // wait for packing to finish
   WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )


   if (sendFromGPU_)
      bufferSystemGPU_[EQUAL_LEVEL][index].sendAll();
   else
      bufferSystemCPU_[EQUAL_LEVEL][index].sendAll();

   communicationInProgress_[EQUAL_LEVEL][index] = true;
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::startCommunicationCoarseToFine(const uint_t index, const uint_t coarsestLevel)
{
   if (packInfos_.empty()) return;
   WALBERLA_ASSERT(!communicationInProgress_[COARSE_TO_FINE][index])
   communicationInProgress_[COARSE_TO_FINE][index] = true;

   auto forest = blockForest_.lock();

   // Schedule Receives
   if (sendFromGPU_)
      bufferSystemGPU_[COARSE_TO_FINE][index].scheduleReceives();
   else
      bufferSystemCPU_[COARSE_TO_FINE][index].scheduleReceives();

   if (!sendFromGPU_)
      for (auto it : headers_[COARSE_TO_FINE][index])
         bufferSystemGPU_[COARSE_TO_FINE][index].sendBuffer(it.first).clear();

   // Start filling send buffers
   {
      for (auto& iBlock : *forest)
      {
         auto coarseBlock = dynamic_cast< Block* >(&iBlock);
         auto nLevel      = coarseBlock->getLevel();

         if (!selectable::isSetSelected(coarseBlock->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_))
            continue;

         if (nLevel != coarsestLevel) continue;

         for (auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir)
         {
            const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex(*dir);

            if (coarseBlock->getNeighborhoodSectionSize(neighborIdx) == uint_t(0)) continue;
            if (!(coarseBlock->neighborhoodSectionHasSmallerBlocks(neighborIdx))) continue;

            for (uint_t n = 0; n != coarseBlock->getNeighborhoodSectionSize(neighborIdx); ++n)
            {
               const BlockID& fineReceiverId = coarseBlock->getNeighborId(neighborIdx, n);
               if (!selectable::isSetSelected(coarseBlock->getNeighborState(neighborIdx, n), requiredBlockSelectors_,
                                              incompatibleBlockSelectors_))
                  continue;

               if( coarseBlock->neighborExistsLocally( neighborIdx, n ) )
               {
                  auto fineReceiverBlock = dynamic_cast< Block * >( forest->getBlock( fineReceiverId ) );
                  //                  for (auto& pi : packInfos_)
                  //                  {
                  //                     pi->communicateLocalCoarseToFine(coarseBlock, fineReceiverBlock, *dir);
                  //                  }

                  GpuBuffer_T& gpuDataBuffer = localBuffer_[COARSE_TO_FINE][index];
                  WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())

                  for (auto& pi : packInfos_)
                  {
                     WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.allocSize() - gpuDataBuffer.size(), pi->sizeCoarseToFineSend(coarseBlock, fineReceiverId, *dir))
                     pi->communicateLocalCoarseToFine(coarseBlock, fineReceiverBlock, *dir, gpuDataBuffer, nullptr);
                  }
               }
               else
               {
                  auto nProcess              = mpi::MPIRank(coarseBlock->getNeighborProcess(neighborIdx, n));
                  GpuBuffer_T& gpuDataBuffer = bufferSystemGPU_[COARSE_TO_FINE][index].sendBuffer(nProcess);
                  for (auto& pi : packInfos_)
                  {
                     WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
                     WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.allocSize() - gpuDataBuffer.size(), pi->sizeCoarseToFineSend(coarseBlock, fineReceiverId, *dir))

                     pi->packDataCoarseToFine(coarseBlock, fineReceiverId, *dir, gpuDataBuffer);

                     if (!sendFromGPU_)
                     {
                        auto gpuDataPtr = gpuDataBuffer.cur();
                        auto size = pi->sizeCoarseToFineSend(coarseBlock, fineReceiverId, *dir);
                        auto cpuDataPtr =
                           bufferSystemCPU_[COARSE_TO_FINE][index].sendBuffer(nProcess).advanceNoResize(size);
                        WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
                        WALBERLA_GPU_CHECK(gpuMemcpyAsync(cpuDataPtr, gpuDataPtr, size, gpuMemcpyDeviceToHost))
                     }
                  }
               }
            }
         }
         localBuffer_[COARSE_TO_FINE][index].clear();
      }
   }

   // wait for packing to finish
   WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )

   if (sendFromGPU_)
      bufferSystemGPU_[COARSE_TO_FINE][index].sendAll();
   else
      bufferSystemCPU_[COARSE_TO_FINE][index].sendAll();

   communicationInProgress_[COARSE_TO_FINE][index] = true;
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::startCommunicationFineToCoarse(const uint_t index, const uint_t finestLevel)
{
   if (packInfos_.empty()) return;

   WALBERLA_ASSERT(!communicationInProgress_[FINE_TO_COARSE][index])

   communicationInProgress_[FINE_TO_COARSE][index] = true;

   auto forest = blockForest_.lock();

   // Schedule Receives
   if (sendFromGPU_)
      bufferSystemGPU_[FINE_TO_COARSE][index].scheduleReceives();
   else
      bufferSystemCPU_[FINE_TO_COARSE][index].scheduleReceives();

   if (!sendFromGPU_)
      for (auto it : headers_[FINE_TO_COARSE][index])
         bufferSystemGPU_[FINE_TO_COARSE][index].sendBuffer(it.first).clear();

   // Start filling send buffers
   {
      for (auto& iBlock : *forest)
      {
         auto fineBlock = dynamic_cast< Block* >(&iBlock);
         auto nLevel    = fineBlock->getLevel();

         if (!selectable::isSetSelected(fineBlock->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_))
            continue;

         if (nLevel != finestLevel) continue;

         for (auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir)
         {
            const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex(*dir);

            if (fineBlock->getNeighborhoodSectionSize(neighborIdx) == uint_t(0)) continue;
            if (!(fineBlock->neighborhoodSectionHasLargerBlock(neighborIdx))) continue;
            WALBERLA_ASSERT_EQUAL(fineBlock->getNeighborhoodSectionSize(neighborIdx), uint_t(1))

            const BlockID& coarseReceiverId = fineBlock->getNeighborId(neighborIdx, uint_t(0));
            if (!selectable::isSetSelected(fineBlock->getNeighborState(neighborIdx, uint_t(0)), requiredBlockSelectors_,
                                           incompatibleBlockSelectors_))
               continue;
            if( fineBlock->neighborExistsLocally( neighborIdx, uint_t(0) ) )
            {
               auto coarseReceiverBlock = dynamic_cast< Block * >( forest->getBlock( coarseReceiverId ) );
               //               for (auto& pi : packInfos_)
               //               {
               //                  pi->communicateLocalFineToCoarse(fineBlock, coarseReceiverBlock, *dir);
               //               }

               GpuBuffer_T& gpuDataBuffer = localBuffer_[FINE_TO_COARSE][index];
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
               for (auto& pi : packInfos_)
               {
                  WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.allocSize() - gpuDataBuffer.size(), pi->sizeFineToCoarseSend(fineBlock, *dir))
                  pi->communicateLocalFineToCoarse(fineBlock, coarseReceiverBlock, *dir, gpuDataBuffer, nullptr);
               }
            }
            else
            {
               auto nProcess              = mpi::MPIRank(fineBlock->getNeighborProcess(neighborIdx, uint_t(0)));
               GpuBuffer_T& gpuDataBuffer = bufferSystemGPU_[FINE_TO_COARSE][index].sendBuffer(nProcess);

               for (auto& pi : packInfos_)
               {
                  WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
                  WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.allocSize() - gpuDataBuffer.size(), pi->sizeFineToCoarseSend(fineBlock, *dir))

                  pi->packDataFineToCoarse(fineBlock, coarseReceiverId, *dir, gpuDataBuffer);

                  if (!sendFromGPU_)
                  {
                     auto gpuDataPtr = gpuDataBuffer.cur();
                     auto size = pi->sizeFineToCoarseSend(fineBlock, *dir);
                     auto cpuDataPtr = bufferSystemCPU_[FINE_TO_COARSE][index].sendBuffer(nProcess).advanceNoResize(size);
                     WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
                     WALBERLA_GPU_CHECK(gpuMemcpyAsync(cpuDataPtr, gpuDataPtr, size, gpuMemcpyDeviceToHost))
                  }
               }
            }
         }
         localBuffer_[FINE_TO_COARSE][index].clear();
      }
   }

   // wait for packing to finish
   WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )

   if (sendFromGPU_)
      bufferSystemGPU_[FINE_TO_COARSE][index].sendAll();
   else
      bufferSystemCPU_[FINE_TO_COARSE][index].sendAll();

   communicationInProgress_[FINE_TO_COARSE][index] = true;
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::waitCommunicateEqualLevel(const uint_t level)
{
   if (!communicationInProgress_[EQUAL_LEVEL][level] || packInfos_.empty()) return;

   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest,
                              "Trying to access communication for a block storage object that doesn't exist anymore")
   WALBERLA_ASSERT_LESS(level, forest->getNumberOfLevels())

   if (sendFromGPU_)
   {
      // auto parallelSection = parallelSectionManager_.parallelSection( nullptr );
      for (auto recvInfo = bufferSystemGPU_[EQUAL_LEVEL][level].begin();
           recvInfo != bufferSystemGPU_[EQUAL_LEVEL][level].end(); ++recvInfo)
      {
         recvInfo.buffer().clear();
         for (auto& header : headers_[EQUAL_LEVEL][level][recvInfo.rank()])
         {
            auto block = dynamic_cast< Block* >(forest->getBlock(header.receiverId));

            for (auto& pi : packInfos_)
            {
               GpuBuffer_T& gpuDataBuffer = recvInfo.buffer();
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
               // parallelSection.run([&](auto s) {
               pi->unpackDataEqualLevel(block, stencil::inverseDir[header.dir], gpuDataBuffer);
               // });
            }
         }
      }
   }
   else
   {
      for (auto recvInfo = bufferSystemCPU_[EQUAL_LEVEL][level].begin();
           recvInfo != bufferSystemCPU_[EQUAL_LEVEL][level].end(); ++recvInfo)
      {
         auto& gpuBuffer = bufferSystemGPU_[EQUAL_LEVEL][level].sendBuffer(recvInfo.rank());

         recvInfo.buffer().clear();
         gpuBuffer.clear();
         for (auto& header : headers_[EQUAL_LEVEL][level][recvInfo.rank()])
         {
            auto block       = dynamic_cast< Block* >(forest->getBlock(header.receiverId));
            auto senderBlock = dynamic_cast< Block* >(forest->getBlock(header.senderId));

            for (auto& pi : packInfos_)
            {
               auto size       = pi->sizeEqualLevelSend(senderBlock, header.dir);
               auto cpuDataPtr = recvInfo.buffer().advanceNoResize(size);
               auto gpuDataPtr = gpuBuffer.cur(); // advanceNoResize( size );
               WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataPtr)

               WALBERLA_GPU_CHECK(gpuMemcpyAsync(gpuDataPtr, cpuDataPtr, size, gpuMemcpyHostToDevice, nullptr))
               pi->unpackDataEqualLevel(block, stencil::inverseDir[header.dir], gpuBuffer);
            }
         }
      }
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
   }
   communicationInProgress_[EQUAL_LEVEL][level] = false;
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::waitCommunicateCoarseToFine(const uint_t fineLevel)
{
   if (!communicationInProgress_[COARSE_TO_FINE][fineLevel] || packInfos_.empty()) return;

   WALBERLA_ASSERT_GREATER(fineLevel, uint_t(0))

   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest,
                              "Trying to access communication for a block storage object that doesn't exist anymore")
   WALBERLA_ASSERT_LESS(fineLevel, forest->getNumberOfLevels())

   if (sendFromGPU_)
   {
      // auto parallelSection = parallelSectionManager_.parallelSection( nullptr );
      for (auto recvInfo = bufferSystemGPU_[COARSE_TO_FINE][fineLevel].begin();
           recvInfo != bufferSystemGPU_[COARSE_TO_FINE][fineLevel].end(); ++recvInfo)
      {
         recvInfo.buffer().clear();
         for (auto& header : headers_[COARSE_TO_FINE][fineLevel][recvInfo.rank()])
         {
            auto block       = dynamic_cast< Block* >(forest->getBlock(header.receiverId));
            auto senderBlock = dynamic_cast< Block* >(forest->getBlock(header.senderId));

            for (auto& pi : packInfos_)
            {
               // auto size = pi->sizeCoarseToFineSend( senderBlock, block->getId(), header.dir );
               GpuBuffer_T& gpuDataBuffer = recvInfo.buffer();
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
               // parallelSection.run([&](auto s) {
               pi->unpackDataCoarseToFine(block, senderBlock->getId(), stencil::inverseDir[header.dir], gpuDataBuffer);
               // });
            }
         }
      }
   }
   else
   {
      auto parallelSection = parallelSectionManager_.parallelSection(nullptr);
      for (auto recvInfo = bufferSystemCPU_[COARSE_TO_FINE][fineLevel].begin();
           recvInfo != bufferSystemCPU_[COARSE_TO_FINE][fineLevel].end(); ++recvInfo)
      {
         auto& gpuBuffer = bufferSystemGPU_[COARSE_TO_FINE][fineLevel].sendBuffer(recvInfo.rank());

         recvInfo.buffer().clear();
         gpuBuffer.clear();
         for (auto& header : headers_[COARSE_TO_FINE][fineLevel][recvInfo.rank()])
         {
            auto block       = dynamic_cast< Block* >(forest->getBlock(header.receiverId));
            auto senderBlock = dynamic_cast< Block* >(forest->getBlock(header.senderId));

            for (auto& pi : packInfos_)
            {
               auto size       = pi->sizeCoarseToFineSend(senderBlock, block->getId(), header.dir);
               auto cpuDataPtr = recvInfo.buffer().advanceNoResize(size);
               auto gpuDataPtr = gpuBuffer.cur(); // advanceNoResize( size );
               WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataPtr)

               parallelSection.run([&](auto s) {
                  WALBERLA_GPU_CHECK(gpuMemcpyAsync(gpuDataPtr, cpuDataPtr, size, gpuMemcpyHostToDevice, s))
                  pi->unpackDataCoarseToFine(block, senderBlock->getId(), stencil::inverseDir[header.dir], gpuBuffer);
               });
            }
         }
      }
   }
   communicationInProgress_[COARSE_TO_FINE][fineLevel] = false;
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::waitCommunicateFineToCoarse(const uint_t fineLevel)
{
   if (!communicationInProgress_[FINE_TO_COARSE][fineLevel] || packInfos_.empty()) return;

   WALBERLA_ASSERT_GREATER(fineLevel, uint_t(0))

   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest,
                              "Trying to access communication for a block storage object that doesn't exist anymore")
   WALBERLA_ASSERT_LESS(fineLevel, forest->getNumberOfLevels())
   // WALBERLA_ASSERT_EQUAL( forestModificationStamp_, forest->getBlockForest().getModificationStamp() );

   if (sendFromGPU_)
   {
      // auto parallelSection = parallelSectionManager_.parallelSection( nullptr );
      for (auto recvInfo = bufferSystemGPU_[FINE_TO_COARSE][fineLevel].begin();
           recvInfo != bufferSystemGPU_[FINE_TO_COARSE][fineLevel].end(); ++recvInfo)
      {
         recvInfo.buffer().clear();
         for (auto& header : headers_[FINE_TO_COARSE][fineLevel][recvInfo.rank()])
         {
            auto block       = dynamic_cast< Block* >(forest->getBlock(header.receiverId));
            auto senderBlock = dynamic_cast< Block* >(forest->getBlock(header.senderId));

            for (auto& pi : packInfos_)
            {
               GpuBuffer_T& gpuDataBuffer = recvInfo.buffer();
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
               // parallelSection.run([&](auto s) {
               pi->unpackDataFineToCoarse(block, senderBlock->getId(), stencil::inverseDir[header.dir], gpuDataBuffer);
               // });
            }
         }
      }
   }
   else
   {
      auto parallelSection = parallelSectionManager_.parallelSection(nullptr);
      for (auto recvInfo = bufferSystemCPU_[FINE_TO_COARSE][fineLevel].begin();
           recvInfo != bufferSystemCPU_[FINE_TO_COARSE][fineLevel].end(); ++recvInfo)
      {
         auto& gpuBuffer = bufferSystemGPU_[FINE_TO_COARSE][fineLevel].sendBuffer(recvInfo.rank());

         recvInfo.buffer().clear();
         gpuBuffer.clear();
         for (auto& header : headers_[FINE_TO_COARSE][fineLevel][recvInfo.rank()])
         {
            auto block       = dynamic_cast< Block* >(forest->getBlock(header.receiverId));
            auto senderBlock = dynamic_cast< Block* >(forest->getBlock(header.senderId));

            for (auto& pi : packInfos_)
            {
               auto size       = pi->sizeFineToCoarseSend(senderBlock, header.dir);
               auto cpuDataPtr = recvInfo.buffer().advanceNoResize(size);
               auto gpuDataPtr = gpuBuffer.cur(); // advanceNoResize( size );
               WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataPtr)

               parallelSection.run([&](auto s) {
                  WALBERLA_GPU_CHECK(gpuMemcpyAsync(gpuDataPtr, cpuDataPtr, size, gpuMemcpyHostToDevice, s))
                  pi->unpackDataFineToCoarse(block, senderBlock->getId(), stencil::inverseDir[header.dir], gpuBuffer);
               });
            }
         }
      }
   }
   communicationInProgress_[FINE_TO_COARSE][fineLevel] = false;
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::setupCommunication()
{
   WALBERLA_ASSERT_GREATER(packInfos_.size(), uint_c(0),
                           "You have not registered a packInfo yet, thus setupCommunication does not work yet.")
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest,
                              "Trying to access communication for a block storage object that doesn't exist anymore")
   const uint_t levels = forest->getNumberOfLevels();

   std::vector< std::vector< std::map< mpi::MPIRank, mpi::MPISize > > >
      receiverInfo; // how many bytes to send to each neighbor
   std::vector< std::vector< mpi::BufferSystem > > headerExchangeBs;

   receiverInfo.resize(3);
   receiverInfo[EQUAL_LEVEL].resize(levels + uint_c(1));
   receiverInfo[COARSE_TO_FINE].resize(levels + uint_c(1));
   receiverInfo[FINE_TO_COARSE].resize(levels + uint_c(1));

   std::vector< std::vector< mpi::MPISize > > localBufferSize;

   headerExchangeBs.resize(3);
   localBufferSize.resize(3);

   for (uint_t j = 0; j <= levels; ++j)
   {
      headerExchangeBs[EQUAL_LEVEL].push_back(mpi::BufferSystem(mpi::MPIManager::instance()->comm(), 123));
      headerExchangeBs[COARSE_TO_FINE].push_back(mpi::BufferSystem(mpi::MPIManager::instance()->comm(), 123));
      headerExchangeBs[FINE_TO_COARSE].push_back(mpi::BufferSystem(mpi::MPIManager::instance()->comm(), 123));

      localBufferSize[EQUAL_LEVEL].push_back(mpi::MPISize(0));
      localBufferSize[COARSE_TO_FINE].push_back(mpi::MPISize(0));
      localBufferSize[FINE_TO_COARSE].push_back(mpi::MPISize(0));
   }

   for (auto& iBlock : *forest)
   {
      auto block = dynamic_cast< Block* >(&iBlock);
      if (!selectable::isSetSelected(block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_)) continue;

      const BlockID& senderId = block->getId();
      auto nLevel             = block->getLevel();

      for (auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir)
      {
         // skip if block has no neighbors in this direction
         const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex(*dir);
         if (block->getNeighborhoodSectionSize(neighborIdx) == uint_t(0)) continue;

         if (block->neighborhoodSectionHasEquallySizedBlock(neighborIdx))
         {
            WALBERLA_ASSERT_EQUAL(block->getNeighborhoodSectionSize(neighborIdx), uint_t(1))
            if (!selectable::isSetSelected(block->getNeighborState(neighborIdx, uint_t(0)), requiredBlockSelectors_,
                                           incompatibleBlockSelectors_))
               continue;
            if( block->neighborExistsLocally( neighborIdx, uint_t(0) ) )
               continue;

            const BlockID& receiverId = block->getNeighborId(neighborIdx, uint_t(0));
            auto nProcess             = mpi::MPIRank(block->getNeighborProcess(neighborIdx, uint_t(0)));

            for (auto& pi : packInfos_)
            {
               receiverInfo[EQUAL_LEVEL][nLevel][nProcess] += mpi::MPISize(pi->sizeEqualLevelSend(block, *dir));
            }

            auto& headerBuffer = headerExchangeBs[EQUAL_LEVEL][nLevel].sendBuffer(nProcess);
            receiverId.toBuffer(headerBuffer);
            senderId.toBuffer(headerBuffer);
            headerBuffer << *dir;
         }
         else if (block->neighborhoodSectionHasSmallerBlocks(neighborIdx))
         {
            auto fineLevel = nLevel + uint_c(1); // For indexing always the fineLevel is taken to be consistent.
            WALBERLA_ASSERT_LESS(fineLevel, levels)

            for (uint_t n = 0; n != block->getNeighborhoodSectionSize(neighborIdx); ++n)
            {
               const BlockID& receiverId = block->getNeighborId(neighborIdx, n);
               if (!selectable::isSetSelected(block->getNeighborState(neighborIdx, n), requiredBlockSelectors_,
                                              incompatibleBlockSelectors_))
                  continue;
               if( block->neighborExistsLocally( neighborIdx, n ) )
               {
                  for (auto& pi : packInfos_)
                     localBufferSize[COARSE_TO_FINE][fineLevel] += mpi::MPISize(pi->sizeCoarseToFineSend(block, receiverId, *dir));
                  continue;
               }

               auto nProcess = mpi::MPIRank(block->getNeighborProcess(neighborIdx, n));
               for (auto& pi : packInfos_)
                  receiverInfo[COARSE_TO_FINE][fineLevel][nProcess] +=
                     mpi::MPISize(pi->sizeCoarseToFineSend(block, receiverId, *dir));
               auto& headerBuffer = headerExchangeBs[COARSE_TO_FINE][fineLevel].sendBuffer(nProcess);
               receiverId.toBuffer(headerBuffer);
               senderId.toBuffer(headerBuffer);
               headerBuffer << *dir;
            }
         }
         else if (block->neighborhoodSectionHasLargerBlock(neighborIdx))
         {
            WALBERLA_ASSERT_EQUAL(block->getNeighborhoodSectionSize(neighborIdx), uint_t(1))

            const BlockID& receiverId = block->getNeighborId(neighborIdx, uint_t(0));
            if (!selectable::isSetSelected(block->getNeighborState(neighborIdx, uint_t(0)), requiredBlockSelectors_,
                                           incompatibleBlockSelectors_))
               continue;

            if( block->neighborExistsLocally( neighborIdx, uint_t(0) ) )
            {
               for (auto& pi : packInfos_)
                  localBufferSize[FINE_TO_COARSE][nLevel] += mpi::MPISize(pi->sizeFineToCoarseSend(block, *dir));
               continue;
            }

            auto nProcess = mpi::MPIRank(block->getNeighborProcess(neighborIdx, uint_t(0)));
            for (auto& pi : packInfos_)
               receiverInfo[FINE_TO_COARSE][nLevel][nProcess] += mpi::MPISize(pi->sizeFineToCoarseSend(block, *dir));

            auto& headerBuffer = headerExchangeBs[FINE_TO_COARSE][nLevel].sendBuffer(nProcess);
            receiverId.toBuffer(headerBuffer);
            senderId.toBuffer(headerBuffer);
            headerBuffer << *dir;
         }
      }
   }

   for (uint_t i = 0; i != 3; ++i)
   {
      for (uint_t j = 0; j <= levels; ++j)
      {
         headerExchangeBs[i][j].setReceiverInfoFromSendBufferState(false, true);
         headerExchangeBs[i][j].sendAll();
         for (auto recvIter = headerExchangeBs[i][j].begin(); recvIter != headerExchangeBs[i][j].end(); ++recvIter)
         {
            auto& headerVector = headers_[i][j][recvIter.rank()];
            auto& buffer       = recvIter.buffer();
            while (buffer.size())
            {
               Header header;
               header.receiverId.fromBuffer(buffer);
               header.senderId.fromBuffer(buffer);
               buffer >> header.dir;
               headerVector.push_back(header);
            }
         }

         bufferSystemCPU_[i][j].setReceiverInfo(receiverInfo[i][j]);
         bufferSystemGPU_[i][j].setReceiverInfo(receiverInfo[i][j]);

         for (auto it : receiverInfo[i][j])
         {
            bufferSystemCPU_[i][j].sendBuffer(it.first).resize(size_t(it.second));
            bufferSystemGPU_[i][j].sendBuffer(it.first).resize(size_t(it.second));
         }
         if (localBufferSize[i][j] > 0)
            localBuffer_[i][j].resize(size_t(localBufferSize[i][j]));
      }
   }

   forestModificationStamp_      = forest->getBlockForest().getModificationStamp();
}

template< typename Stencil >
bool NonUniformGPUScheme< Stencil >::isAnyCommunicationInProgress() const
{
   for (auto caseIt = communicationInProgress_.begin(); caseIt != communicationInProgress_.end(); ++caseIt)
      for (auto levelIt = caseIt->begin(); levelIt != caseIt->end(); ++levelIt)
         if (*levelIt) return true;

   return false;
}

template< typename Stencil >
NonUniformGPUScheme< Stencil >::~NonUniformGPUScheme()
{
   for (uint_t i = 0; i != bufferSystemGPU_[EQUAL_LEVEL].size(); ++i)
   {
      waitCommunicateEqualLevel(i);
      waitCommunicateCoarseToFine(i);
      waitCommunicateFineToCoarse(i);
   }
}

template< typename Stencil >
void NonUniformGPUScheme< Stencil >::addPackInfo(const shared_ptr< GeneratedNonUniformGPUPackInfo >& pi)
{
   if (isAnyCommunicationInProgress())
   {
      WALBERLA_ABORT("You may not add a PackInfo to a NonUniformBufferedScheme if any communication is in progress!")
   }
   packInfos_.push_back(pi);
   setupCommunication();
}

} // namespace walberla::gpu::communication
