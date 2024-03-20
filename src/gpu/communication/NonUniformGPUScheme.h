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

#include "core/DataTypes.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPIWrapper.h"

#include "domain_decomposition/IBlock.h"

#include "gpu/ErrorChecking.h"
#include "gpu/GPUWrapper.h"
#include "gpu/communication/CustomMemoryBuffer.h"
#include "gpu/communication/GeneratedNonUniformGPUPackInfo.h"

#include "stencil/Directions.h"

#include <memory>
#include <thread>

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
                                int tag = 5432);

   explicit NonUniformGPUScheme(const weak_ptr< StructuredBlockForest >& bf, const Set< SUID >& requiredBlockSelectors,
                                const Set< SUID >& incompatibleBlockSelectors, bool sendDirectlyFromGPU = false,
                                int tag = 5432);

   ~NonUniformGPUScheme();

   //** Pack Info Registration *****************************************************************************************
   /*! \name Pack Info Registration */
   //@{
   void addPackInfo(const shared_ptr< GeneratedNonUniformGPUPackInfo >& pi);
   //@}
   //*******************************************************************************************************************

   inline void communicateEqualLevel(uint_t level);
   inline void communicateCoarseToFine(uint_t fineLevel);
   inline void communicateFineToCoarse(uint_t fineLevel);

   std::function<void()>  communicateEqualLevelFunctor(const uint_t level) {
      return [level, this](){ NonUniformGPUScheme::communicateEqualLevel(level);};
   }
   std::function<void()>  communicateCoarseToFineFunctor(const uint_t fineLevel) {
      return [fineLevel, this](){ NonUniformGPUScheme::communicateCoarseToFine(fineLevel);};
   }
   std::function<void()>  communicateFineToCoarseFunctor(const uint_t fineLevel) {
      return [fineLevel, this](){ NonUniformGPUScheme::communicateFineToCoarse(fineLevel);};
   }

   inline void startCommunicateEqualLevel(uint_t level);
   inline void startCommunicateCoarseToFine(uint_t fineLevel);
   inline void startCommunicateFineToCoarse(uint_t fineLevel);

   inline void waitCommunicateEqualLevel(uint_t level);
   inline void waitCommunicateCoarseToFine(uint_t fineLevel);
   inline void waitCommunicateFineToCoarse(uint_t fineLevel);

 private:
   void setupCommunication();

   void init();
   void refresh();

   [[nodiscard]] bool isAnyCommunicationInProgress() const;

   void startCommunicationEqualLevel(uint_t index, std::set< uint_t >& participatingLevels);
   void startCommunicationCoarseToFine(uint_t index, uint_t coarsestLevel);
   void startCommunicationFineToCoarse(uint_t index, uint_t finestLevel);

   weak_ptr< StructuredBlockForest > blockForest_;
   uint_t forestModificationStamp_{uint_c(0)};

   std::vector< std::vector< bool > > communicationInProgress_;
   bool sendFromGPU_;
   int baseTag_;

   std::vector< std::vector< mpi::GenericBufferSystem< CpuBuffer_T, CpuBuffer_T > > > bufferSystemCPU_;
   std::vector< std::vector< mpi::GenericBufferSystem< GpuBuffer_T, GpuBuffer_T > > > bufferSystemGPU_;
   std::vector< std::vector< GpuBuffer_T > > localBuffer_;
   GpuBuffer_T adaptiveGPUBuffer;

   std::vector< shared_ptr< GeneratedNonUniformGPUPackInfo > > packInfos_;

   struct Header
   {
      BlockID receiverId;
      BlockID senderId;
      stencil::Direction dir;
   };
   std::vector< std::vector< std::map< mpi::MPIRank, std::vector< Header > > > > headers_;

   Set< SUID > requiredBlockSelectors_;
   Set< SUID > incompatibleBlockSelectors_;

   gpuStream_t streams_[Stencil::Q];
};

template< typename Stencil >
NonUniformGPUScheme< Stencil >::NonUniformGPUScheme(const weak_ptr< StructuredBlockForest >& bf, bool sendDirectlyFromGPU,
                                                    const int tag)
   : blockForest_(bf), sendFromGPU_(sendDirectlyFromGPU), baseTag_(tag),
     requiredBlockSelectors_(Set< SUID >::emptySet()), incompatibleBlockSelectors_(Set< SUID >::emptySet())
{
   WALBERLA_MPI_SECTION()
      {
// Open MPI supports compile time CUDA-aware support check
#if (defined(OPEN_MPI) && OPEN_MPI) && !(defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT)
         WALBERLA_CHECK(!sendDirectlyFromGPU)
#endif
      }
   init();

   if(sendFromGPU_){WALBERLA_LOG_DETAIL_ON_ROOT("Using GPU-Direct Communication in NonUniformGPUScheme")}
   else{WALBERLA_LOG_DETAIL_ON_ROOT("Using Communication via CPU Memory")}

}

template< typename Stencil >
NonUniformGPUScheme< Stencil >::NonUniformGPUScheme(const weak_ptr< StructuredBlockForest >& bf,
                                                    const Set< SUID >& requiredBlockSelectors,
                                                    const Set< SUID >& incompatibleBlockSelectors,
                                                    bool sendDirectlyFromGPU, const int tag)
   : blockForest_(bf), requiredBlockSelectors_(requiredBlockSelectors),
     incompatibleBlockSelectors_(incompatibleBlockSelectors), sendFromGPU_(sendDirectlyFromGPU), baseTag_(tag)
{
   WALBERLA_MPI_SECTION()
      {
#if !(defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT)
         WALBERLA_CHECK(!sendDirectlyFromGPU)
#endif
      }
   init();
   if(sendFromGPU_){WALBERLA_LOG_DETAIL_ON_ROOT("Using GPU-Direct Communication in NonUniformGPUScheme")}
   else{WALBERLA_LOG_DETAIL_ON_ROOT("Using Communication via CPU Memory")}
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

   for (uint_t i = 0; i < Stencil::Q; ++i)
   {
      WALBERLA_GPU_CHECK(gpuStreamCreate(&streams_[i]))
   }
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
         bufferSystemCPU_[i].emplace_back(mpi::MPIManager::instance()->comm(), baseTag_ + int_c(i * levels + j));
         bufferSystemGPU_[i].emplace_back(mpi::MPIManager::instance()->comm(), baseTag_ + int_c(i * levels + j));
         localBuffer_[i].emplace_back();
      }

      communicationInProgress_[i].resize(size_t(levels + uint_t(1)), false);
   }

#ifndef NDEBUG
   for (auto & packInfo : packInfos_)
      packInfo->clearBufferSizeCheckMap();
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
               pi->communicateLocalEqualLevel(senderBlock, receiverBlock, *dir, streams_[*dir]);
            }
         }
         else
         {
            auto nProcess              = mpi::MPIRank(senderBlock->getNeighborProcess(neighborIdx, uint_t(0)));
            GpuBuffer_T& gpuDataBuffer = bufferSystemGPU_[EQUAL_LEVEL][index].sendBuffer(nProcess);

            for (auto& pi : packInfos_)
            {
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
               WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.remainingSize(), pi->sizeEqualLevelSend(senderBlock, *dir))
               if(sendFromGPU_)
               {
                  pi->packDataEqualLevel(senderBlock, *dir, gpuDataBuffer, streams_[*dir]);
               }
               else
               {
                  auto gpuDataPtr = gpuDataBuffer.cur();
                  // packDataEqualLevel moves the pointer with advanceNoResize
                  pi->packDataEqualLevel(senderBlock, *dir, gpuDataBuffer, streams_[*dir]);
                  auto size = pi->sizeEqualLevelSend(senderBlock, *dir);
                  auto cpuDataPtr = bufferSystemCPU_[EQUAL_LEVEL][index].sendBuffer(nProcess).advanceNoResize(size);
                  WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
                  WALBERLA_GPU_CHECK(gpuMemcpyAsync(cpuDataPtr, gpuDataPtr, size, gpuMemcpyDeviceToHost, streams_[*dir]))
               }
            }
         }
      }
   }
   // wait for packing to finish
   for (uint_t i = 0; i < Stencil::Q; ++i)
   {
      WALBERLA_GPU_CHECK(gpuStreamSynchronize(streams_[i]))
   }


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

   for (auto it : headers_[COARSE_TO_FINE][index]){
      bufferSystemGPU_[COARSE_TO_FINE][index].sendBuffer(it.first).clear();
   }
   // wait until communication dependent kernels are finished
   WALBERLA_GPU_CHECK(gpuDeviceSynchronize())

   // Start filling send buffers
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
               GpuBuffer_T& gpuDataBuffer = localBuffer_[COARSE_TO_FINE][index];

               for (auto& pi : packInfos_)
               {
                  WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
                  WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.remainingSize(), pi->sizeCoarseToFineSend(coarseBlock, fineReceiverId, *dir))
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
                  WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.remainingSize(), pi->sizeCoarseToFineSend(coarseBlock, fineReceiverId, *dir))
                  if (sendFromGPU_)
                  {
                     pi->packDataCoarseToFine(coarseBlock, fineReceiverId, *dir, gpuDataBuffer, streams_[0]);
                  }
                  else
                  {
                     gpuDataBuffer.clear();
                     auto gpuDataPtr = gpuDataBuffer.cur();
                     // packDataCoarseToFine moves the pointer with advanceNoResize
                     pi->packDataCoarseToFine(coarseBlock, fineReceiverId, *dir, gpuDataBuffer, streams_[0]);
                     auto size = pi->sizeCoarseToFineSend(coarseBlock, fineReceiverId, *dir);
                     auto cpuDataPtr = bufferSystemCPU_[COARSE_TO_FINE][index].sendBuffer(nProcess).advanceNoResize(size);
                     WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
                     WALBERLA_GPU_CHECK(gpuMemcpyAsync(cpuDataPtr, gpuDataPtr, size, gpuMemcpyDeviceToHost, streams_[0]))
                  }
               }
            }
         }
      }
      localBuffer_[COARSE_TO_FINE][index].clear();
   }


   // wait for packing to finish
   for (uint_t i = 0; i < Stencil::Q; ++i)
   {
      WALBERLA_GPU_CHECK(gpuStreamSynchronize(streams_[i]))
   }

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

   for (auto it : headers_[FINE_TO_COARSE][index])
      bufferSystemGPU_[FINE_TO_COARSE][index].sendBuffer(it.first).clear();

   // wait until communication dependent kernels are finished
   WALBERLA_GPU_CHECK(gpuDeviceSynchronize())

   // Start filling send buffers
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
            GpuBuffer_T& gpuDataBuffer = localBuffer_[FINE_TO_COARSE][index];

            for (auto& pi : packInfos_)
            {
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
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
               WALBERLA_ASSERT_GREATER_EQUAL(gpuDataBuffer.remainingSize(), pi->sizeFineToCoarseSend(fineBlock, *dir))
               if (sendFromGPU_)
               {
                  pi->packDataFineToCoarse(fineBlock, coarseReceiverId, *dir, gpuDataBuffer, streams_[0]);
               }
               else
               {
                  gpuDataBuffer.clear();
                  auto gpuDataPtr = gpuDataBuffer.cur();
                  // packDataFineToCoarse moves the pointer with advanceNoResize
                  pi->packDataFineToCoarse(fineBlock, coarseReceiverId, *dir, gpuDataBuffer, streams_[0]);
                  auto size = pi->sizeFineToCoarseSend(fineBlock, *dir);
                  auto cpuDataPtr = bufferSystemCPU_[FINE_TO_COARSE][index].sendBuffer(nProcess).advanceNoResize(size);
                  WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
                  WALBERLA_GPU_CHECK(gpuMemcpyAsync(cpuDataPtr, gpuDataPtr, size, gpuMemcpyDeviceToHost, streams_[0]))
               }
            }
         }
      }
      localBuffer_[FINE_TO_COARSE][index].clear();
   }
   // wait for packing to finish
   for (uint_t i = 0; i < Stencil::Q; ++i)
   {
      WALBERLA_GPU_CHECK(gpuStreamSynchronize(streams_[i]))
   }

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
               pi->unpackDataEqualLevel(block, stencil::inverseDir[header.dir], gpuDataBuffer, streams_[stencil::inverseDir[header.dir]]);
            }
         }
      }
   }
   else
   {
      for (auto recvInfo = bufferSystemCPU_[EQUAL_LEVEL][level].begin();
           recvInfo != bufferSystemCPU_[EQUAL_LEVEL][level].end(); ++recvInfo)
      {
         auto &gpuBuffer = bufferSystemGPU_[EQUAL_LEVEL][level].sendBuffer(recvInfo.rank());

         recvInfo.buffer().clear();
         gpuBuffer.clear();

         for (auto &header : headers_[EQUAL_LEVEL][level][recvInfo.rank()])
         {
            auto block       = dynamic_cast< Block* >(forest->getBlock(header.receiverId));
            for (auto& pi : packInfos_)
            {
               auto size       = pi->sizeEqualLevelSend(block, stencil::inverseDir[header.dir]);
               auto cpuDataPtr = recvInfo.buffer().advanceNoResize(size);
               auto gpuDataPtr = gpuBuffer.cur(); // advanceNoResize( size );
               WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataPtr)

               WALBERLA_GPU_CHECK(gpuMemcpyAsync(gpuDataPtr, cpuDataPtr, size, gpuMemcpyHostToDevice, streams_[stencil::inverseDir[header.dir]]))
               pi->unpackDataEqualLevel(block, stencil::inverseDir[header.dir], gpuBuffer, streams_[stencil::inverseDir[header.dir]]);
            }
         }
      }
   }
   for (uint_t i = 0; i < Stencil::Q; ++i)
   {
      WALBERLA_GPU_CHECK(gpuStreamSynchronize(streams_[i]))
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

   if (sendFromGPU_) {
      for (auto recvInfo = bufferSystemGPU_[COARSE_TO_FINE][fineLevel].begin();
           recvInfo != bufferSystemGPU_[COARSE_TO_FINE][fineLevel].end(); ++recvInfo) {
         recvInfo.buffer().clear();
         for (auto &header: headers_[COARSE_TO_FINE][fineLevel][recvInfo.rank()]) {
            auto fineReceiver = dynamic_cast< Block * >(forest->getBlock(header.receiverId));
            for (auto &pi: packInfos_) {
               GpuBuffer_T &gpuDataBuffer = recvInfo.buffer();
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
               pi->unpackDataCoarseToFine(fineReceiver, header.senderId, stencil::inverseDir[header.dir],
                                          gpuDataBuffer, streams_[0]);
            }
         }
      }
   } else
   {
      for (auto recvInfo = bufferSystemCPU_[COARSE_TO_FINE][fineLevel].begin();
           recvInfo != bufferSystemCPU_[COARSE_TO_FINE][fineLevel].end(); ++recvInfo) {

         adaptiveGPUBuffer.clear();
         adaptiveGPUBuffer.resize(recvInfo.buffer().allocSize());
         recvInfo.buffer().clear();

         for (auto &header: headers_[COARSE_TO_FINE][fineLevel][recvInfo.rank()]) {
            auto fineReceiver = dynamic_cast< Block * >(forest->getBlock(header.receiverId));
            // WALBERLA_ASSERT_NOT_NULLPTR(fineReceiver)
            for (auto &pi: packInfos_) {
               auto size = pi->sizeCoarseToFineReceive(fineReceiver, stencil::inverseDir[header.dir]);

               auto cpuDataPtr = recvInfo.buffer().advanceNoResize(size);
               auto gpuDataPtr = adaptiveGPUBuffer.cur(); // advanceNoResize( size );
               WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataPtr)

               WALBERLA_GPU_CHECK(gpuMemcpyAsync(gpuDataPtr, cpuDataPtr, size, gpuMemcpyHostToDevice, streams_[0]))
               pi->unpackDataCoarseToFine(fineReceiver, header.senderId, stencil::inverseDir[header.dir], adaptiveGPUBuffer, streams_[0]);
            }
         }
      }
   }

   for (uint_t i = 0; i < Stencil::Q; ++i)
   {
      WALBERLA_GPU_CHECK(gpuStreamSynchronize(streams_[i]))
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

   if (sendFromGPU_)
   {
      for (auto recvInfo = bufferSystemGPU_[FINE_TO_COARSE][fineLevel].begin();
           recvInfo != bufferSystemGPU_[FINE_TO_COARSE][fineLevel].end(); ++recvInfo)
      {
         recvInfo.buffer().clear();
         for (auto& header : headers_[FINE_TO_COARSE][fineLevel][recvInfo.rank()])
         {
            auto block       = dynamic_cast< Block* >(forest->getBlock(header.receiverId));
            for (auto& pi : packInfos_)
            {
               GpuBuffer_T& gpuDataBuffer = recvInfo.buffer();
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataBuffer.cur())
               pi->unpackDataFineToCoarse(block, header.senderId, stencil::inverseDir[header.dir], gpuDataBuffer, streams_[0]);
            }
         }
      }
   }
   else
   {
      for (auto recvInfo = bufferSystemCPU_[FINE_TO_COARSE][fineLevel].begin();
           recvInfo != bufferSystemCPU_[FINE_TO_COARSE][fineLevel].end(); ++recvInfo)
      {
         recvInfo.buffer().clear();
         adaptiveGPUBuffer.clear();
         adaptiveGPUBuffer.resize(recvInfo.buffer().allocSize());
         for (auto& header : headers_[FINE_TO_COARSE][fineLevel][recvInfo.rank()])
         {
            auto block       = dynamic_cast< Block* >(forest->getBlock(header.receiverId));
            for (auto& pi : packInfos_)
            {
               auto size       = pi->sizeFineToCoarseSend(block, stencil::inverseDir[header.dir]);
               auto cpuDataPtr = recvInfo.buffer().advanceNoResize(size);
               auto gpuDataPtr = adaptiveGPUBuffer.cur(); // advanceNoResize( size );
               WALBERLA_ASSERT_NOT_NULLPTR(cpuDataPtr)
               WALBERLA_ASSERT_NOT_NULLPTR(gpuDataPtr)

               WALBERLA_GPU_CHECK(gpuMemcpyAsync(gpuDataPtr, cpuDataPtr, size, gpuMemcpyHostToDevice, streams_[0]))
               pi->unpackDataFineToCoarse(block, header.senderId, stencil::inverseDir[header.dir], adaptiveGPUBuffer, streams_[0]);
            }
         }
      }
   }
   for (uint_t i = 0; i < Stencil::Q; ++i)
   {
      WALBERLA_GPU_CHECK(gpuStreamSynchronize(streams_[i]))
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

   std::vector< std::vector< std::map< mpi::MPIRank, mpi::MPISize > > > senderInfo; // how many bytes to send to each neighbor
   std::vector< std::vector< std::set< mpi::MPIRank > > > receiverInfo; // how many bytes to receive from each neighbor

   std::vector< std::vector< shared_ptr< mpi::BufferSystem > > > headerExchangeBs;

   std::vector< std::vector< mpi::MPISize > > localBufferSize;

   localBufferSize.resize(3);
   senderInfo.resize(3);
   receiverInfo.resize(3);

   senderInfo[EQUAL_LEVEL].resize(levels + uint_c(1));
   senderInfo[COARSE_TO_FINE].resize(levels + uint_c(1));
   senderInfo[FINE_TO_COARSE].resize(levels + uint_c(1));

   receiverInfo[EQUAL_LEVEL].resize(levels + uint_c(1));
   receiverInfo[COARSE_TO_FINE].resize(levels + uint_c(1));
   receiverInfo[FINE_TO_COARSE].resize(levels + uint_c(1));

   headerExchangeBs.resize(3);

   for (uint_t j = 0; j <= levels; ++j)
   {
      headerExchangeBs[EQUAL_LEVEL].push_back(make_shared< mpi::BufferSystem >(mpi::MPIManager::instance()->comm(), 123));
      headerExchangeBs[COARSE_TO_FINE].push_back(make_shared< mpi::BufferSystem >(mpi::MPIManager::instance()->comm(), 123));
      headerExchangeBs[FINE_TO_COARSE].push_back(make_shared< mpi::BufferSystem >(mpi::MPIManager::instance()->comm(), 123));

      localBufferSize[EQUAL_LEVEL].push_back(mpi::MPISize(0));
      localBufferSize[COARSE_TO_FINE].push_back(mpi::MPISize(0));
      localBufferSize[FINE_TO_COARSE].push_back(mpi::MPISize(0));
   }

   for (auto& iBlock : *forest)
   {
      auto block = dynamic_cast< Block* >(&iBlock);
      if (!selectable::isSetSelected(block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_)) continue;

      const BlockID& senderId = block->getId();
      auto level       = block->getLevel();

      for (auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir)
      {
         // skip if block has no neighbors in this direction
         const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex(*dir);
         if (block->getNeighborhoodSectionSize(neighborIdx) == uint_t(0)) continue;

         // EQUAL_LEVEL communication
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
               senderInfo[EQUAL_LEVEL][level][nProcess] += mpi::MPISize(pi->sizeEqualLevelSend(block, *dir));
            }

            auto& headerBuffer = headerExchangeBs[EQUAL_LEVEL][level]->sendBuffer(nProcess);
            receiverId.toBuffer(headerBuffer);
            senderId.toBuffer(headerBuffer);
            headerBuffer << *dir;

            receiverInfo[EQUAL_LEVEL][level].insert( nProcess );
         }
         else if (block->neighborhoodSectionHasSmallerBlocks(neighborIdx))
         {
            auto fineLevel = level + uint_c(1); // For indexing always the fineLevel is taken to be consistent.
            WALBERLA_ASSERT_LESS(fineLevel, levels)

            for (uint_t n = 0; n != block->getNeighborhoodSectionSize(neighborIdx); ++n)
            {
               const BlockID& receiverId = block->getNeighborId(neighborIdx, n);

               if (!selectable::isSetSelected(block->getNeighborState(neighborIdx, n), requiredBlockSelectors_, incompatibleBlockSelectors_))
                  continue;

               if( block->neighborExistsLocally( neighborIdx, n ) )
               {
                  for (auto& pi : packInfos_)
                     localBufferSize[COARSE_TO_FINE][fineLevel] += mpi::MPISize(pi->sizeCoarseToFineSend(block, receiverId, *dir));
                  continue;
               }

               auto nProcess = mpi::MPIRank(block->getNeighborProcess(neighborIdx, n));
               for (auto& pi : packInfos_)
                  senderInfo[COARSE_TO_FINE][fineLevel][nProcess] += mpi::MPISize(pi->sizeCoarseToFineSend(block, receiverId, *dir));

               auto& headerBuffer = headerExchangeBs[COARSE_TO_FINE][fineLevel]->sendBuffer(nProcess);
               receiverId.toBuffer(headerBuffer);
               senderId.toBuffer(headerBuffer);
               headerBuffer << *dir;

               receiverInfo[FINE_TO_COARSE][fineLevel].insert( nProcess );
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
                  localBufferSize[FINE_TO_COARSE][level] += mpi::MPISize(pi->sizeFineToCoarseSend(block, *dir));
               continue;
            }

            auto nProcess = mpi::MPIRank(block->getNeighborProcess(neighborIdx, uint_t(0)));
            for (auto& pi : packInfos_)
               senderInfo[FINE_TO_COARSE][level][nProcess] += mpi::MPISize(pi->sizeFineToCoarseSend(block, *dir));

            auto& headerBuffer = headerExchangeBs[FINE_TO_COARSE][level]->sendBuffer(nProcess);
            receiverId.toBuffer(headerBuffer);
            senderId.toBuffer(headerBuffer);
            headerBuffer << *dir;

            receiverInfo[COARSE_TO_FINE][level].insert( nProcess );
         }
      }
   }

   for (uint_t i = 0; i != 3; ++i)
   {
      for (uint_t j = 0; j <= levels; ++j)
      {
         headerExchangeBs[i][j]->setReceiverInfo(receiverInfo[i][j], false);
         headerExchangeBs[i][j]->sendAll();
         for (auto recvIter = headerExchangeBs[i][j]->begin(); recvIter != headerExchangeBs[i][j]->end(); ++recvIter) {
            auto &headerVector = headers_[i][j][recvIter.rank()];
            auto &buffer = recvIter.buffer();
            while (buffer.size()) {
               Header header;
               header.receiverId.fromBuffer(buffer);
               header.senderId.fromBuffer(buffer);
               buffer >> header.dir;
               headerVector.push_back(header);
            }
         }
         bufferSystemCPU_[i][j].setReceiverInfo(receiverInfo[i][j], false);
         bufferSystemGPU_[i][j].setReceiverInfo(receiverInfo[i][j], false);
         for (auto it : senderInfo[i][j])
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
   const uint_t levels = uint_c(communicationInProgress_[0].size());
   for (uint_t i = 0; i != 3; ++i)
      for (uint_t j = 0; j != levels; ++j)
         if (communicationInProgress_[i][j]) return true;

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

   for (uint_t i = 0; i < Stencil::Q; ++i)
   {
      WALBERLA_GPU_CHECK(gpuStreamDestroy(streams_[i]))
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
