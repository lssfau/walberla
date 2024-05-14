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
//! \file NonuniformGeneratedGPUPdfPackInfo.impl.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "NonuniformGeneratedGPUPdfPackInfo.h"

using namespace walberla::lbm_generated::util;

namespace walberla::lbm_generated {

/***********************************************************************************************************************
 *                                                  Factory Functions                                                  *
 **********************************************************************************************************************/


/**
 * Sets up a NonuniformGeneratedPdfPackInfo.
 *
 * @tparam LatticeStorageSpecification_T
 * @tparam PackingKernels_T
 * @param blocks
 * @param pdfFieldID
 * @param dataIdentifier
 * @return
 */
template< typename PdfField_T>
std::shared_ptr< NonuniformGeneratedGPUPdfPackInfo< PdfField_T > >
   setupNonuniformGPUPdfCommunication( const std::weak_ptr< StructuredBlockForest > & blocks,
                                 const BlockDataID pdfFieldID,
                                 const std::string & dataIdentifier)
{
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;

   auto sbf = blocks.lock();
   WALBERLA_CHECK_NOT_NULLPTR(sbf, "Trying to create Nonuniform GPU Packinfo for a block storage object that doesn't exist anymore" );

   auto handling = std::make_shared<NonuniformGPUCommDataHandling< LatticeStorageSpecification_T > >(blocks);
   BlockDataID commDataID = sbf->addBlockData(handling, dataIdentifier);

   return std::make_shared<NonuniformGeneratedGPUPdfPackInfo< PdfField_T > >(pdfFieldID, commDataID);
}


/***********************************************************************************************************************
 *                                          Equal Level Communication                                                  *
 **********************************************************************************************************************/

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::unpackDataEqualLevel(Block* receiver,
                                                                           Direction dir,
                                                                           GpuBuffer_T & buffer, gpuStream_t stream)
{
   auto field = receiver->getData< PdfField_T >(pdfFieldID_);
   CellInterval ci;
   cell_idx_t gls = skipsThroughCoarseBlock(receiver, dir) ? 2 : 1;
   field->getGhostRegion(dir, ci, gls, false);
   uint_t size              = kernels_.size(ci, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.unpackDirection(field, ci, bufferPtr, dir, stream);
}

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::communicateLocalEqualLevel(
   const Block* sender, Block* receiver, stencil::Direction dir, gpuStream_t stream)
{
   auto srcField = const_cast< Block* >(sender)->getData< PdfField_T >(pdfFieldID_);
   auto dstField = receiver->getData< PdfField_T >(pdfFieldID_);

   CellInterval srcRegion;
   CellInterval dstRegion;
   cell_idx_t gls = skipsThroughCoarseBlock(sender, dir) ? 2 : 1;
   srcField->getSliceBeforeGhostLayer(dir, srcRegion, gls, false);
   dstField->getGhostRegion(stencil::inverseDir[dir], dstRegion, gls, false);
   kernels_.localCopyDirection(srcField, srcRegion, dstField, dstRegion, dir, stream);
}


template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::packDataEqualLevelImpl(
   const Block* sender, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream) const
{
   auto field = const_cast< Block* >(sender)->getData< PdfField_T >(pdfFieldID_);
   CellInterval ci;
   cell_idx_t gls = skipsThroughCoarseBlock(sender, dir) ? 2 : 1;
   field->getSliceBeforeGhostLayer(dir, ci, gls, false);
   uint_t size              = kernels_.size(ci, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.packDirection(field, ci, bufferPtr, dir, stream);
}

/***********************************************************************************************************************
 *                                          Coarse to Fine Communication                                               *
 **********************************************************************************************************************/

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::packDataCoarseToFineImpl(
   const Block* coarseSender, const BlockID& fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream) const
{
   auto field = const_cast< Block* >(coarseSender)->getData< PdfField_T >(pdfFieldID_);

   std::vector< std::pair< Direction, CellInterval > > intervals;
   getCoarseBlockCommIntervals(fineReceiver, dir, field, intervals);

   for (auto t : intervals)
   {
      CellInterval ci          = t.second;
      auto bufferPtr = buffer.advanceNoResize(kernels_.size(ci));
      kernels_.packAll(field, ci, bufferPtr, stream);
   }
}

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::unpackDataCoarseToFine(
   Block* fineReceiver, const BlockID& /*coarseSender*/, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream)
{
   auto field = fineReceiver->getData< PdfField_T >(pdfFieldID_);

   std::vector< std::pair< Direction, CellInterval > > intervals;
   getFineBlockCommIntervals(fineReceiver->getId(), dir, field, intervals);

   for (auto t : intervals)
   {
      Direction d              = t.first;
      CellInterval ci          = t.second;
      uint_t size              = kernels_.redistributeSize(ci);
      auto bufferPtr = buffer.advanceNoResize(size);
      kernels_.unpackRedistribute(field, ci, bufferPtr, d, stream);
   }
}

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::communicateLocalCoarseToFine(
   const Block* coarseSender, Block* fineReceiver, stencil::Direction dir, gpuStream_t stream)
{
   auto srcField = const_cast< Block* >(coarseSender)->getData< PdfField_T >(pdfFieldID_);
   auto dstField = fineReceiver->getData< PdfField_T >(pdfFieldID_);

   std::vector< std::pair< Direction, CellInterval > > srcIntervals;
   getCoarseBlockCommIntervals(fineReceiver->getId(), dir, srcField, srcIntervals);

   std::vector< std::pair< Direction, CellInterval > > dstIntervals;
   getFineBlockCommIntervals(fineReceiver->getId(), stencil::inverseDir[dir], dstField, dstIntervals);

   WALBERLA_ASSERT_EQUAL(srcIntervals.size(), dstIntervals.size())

   for(size_t index = 0; index < srcIntervals.size(); index++)
   {
      CellInterval srcInterval = srcIntervals[index].second;

      Direction const unpackDir      = dstIntervals[index].first;
      CellInterval dstInterval = dstIntervals[index].second;

      uint_t packSize      = kernels_.size(srcInterval);

#ifndef NDEBUG
      Direction const packDir        = srcIntervals[index].first;
      WALBERLA_ASSERT_EQUAL(packDir, stencil::inverseDir[unpackDir])
      uint_t unpackSize = kernels_.redistributeSize(dstInterval);
      WALBERLA_ASSERT_EQUAL(packSize, unpackSize)
#endif

      // TODO: This is a dirty workaround. Code-generate direct redistribution!
      unsigned char *buffer;
      WALBERLA_GPU_CHECK( gpuMalloc( &buffer, packSize))
      kernels_.packAll(srcField, srcInterval, buffer, stream);
      kernels_.unpackRedistribute(dstField, dstInterval, buffer, unpackDir, stream);
      WALBERLA_GPU_CHECK(gpuFree(buffer))
   }
}

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::communicateLocalCoarseToFine(
   const Block* coarseSender, Block* fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream)
{
   auto srcField = const_cast< Block* >(coarseSender)->getData< PdfField_T >(pdfFieldID_);
   auto dstField = fineReceiver->getData< PdfField_T >(pdfFieldID_);

   std::vector< std::pair< Direction, CellInterval > > srcIntervals;
   getCoarseBlockCommIntervals(fineReceiver->getId(), dir, srcField, srcIntervals);

   std::vector< std::pair< Direction, CellInterval > > dstIntervals;
   getFineBlockCommIntervals(fineReceiver->getId(), stencil::inverseDir[dir], dstField, dstIntervals);

   WALBERLA_ASSERT_EQUAL(srcIntervals.size(), dstIntervals.size())

   for(size_t index = 0; index < srcIntervals.size(); index++)
   {
      CellInterval srcInterval = srcIntervals[index].second;

      Direction const unpackDir      = dstIntervals[index].first;
      CellInterval dstInterval = dstIntervals[index].second;

      uint_t packSize      = kernels_.size(srcInterval);

#ifndef NDEBUG
      Direction const packDir        = srcIntervals[index].first;
      WALBERLA_ASSERT_EQUAL(packDir, stencil::inverseDir[unpackDir])
      uint_t unpackSize = kernels_.redistributeSize(dstInterval);
      WALBERLA_ASSERT_EQUAL(packSize, unpackSize)
#endif

      auto bufferPtr = buffer.advanceNoResize(packSize);
      kernels_.packAll(srcField, srcInterval, bufferPtr, stream);
      kernels_.unpackRedistribute(dstField, dstInterval, bufferPtr, unpackDir, stream);
   }
}


/***********************************************************************************************************************
 *                                          Fine to Coarse Communication                                               *
 **********************************************************************************************************************/

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::prepareCoalescence(Block* coarseReceiver, gpuStream_t gpuStream)
{
   auto dstField = coarseReceiver->getData<PdfField_T>(pdfFieldID_);

   for(auto it = CommunicationStencil::beginNoCenter(); it != CommunicationStencil::end(); ++it){
      uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(*it);
      if(coarseReceiver->neighborhoodSectionHasSmallerBlocks(nSecIdx)){
         CellInterval ci;
         dstField->getSliceBeforeGhostLayer(*it, ci, 1);
         kernels_.zeroCoalescenceRegion(dstField, ci, *it, gpuStream);
      }
   }
}

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::unpackDataFineToCoarse(
   Block* coarseReceiver, const walberla::BlockID& fineSender, walberla::stencil::Direction dir,
   GpuBuffer_T & buffer, gpuStream_t stream)
{
   auto dstField = coarseReceiver->getData<PdfField_T>(pdfFieldID_);

   CellInterval ci = getCoarseBlockCoalescenceInterval(coarseReceiver, fineSender, dir, dstField);
   uint_t size = kernels_.size(ci, dir);
   unsigned char* bufferPtr = buffer.advanceNoResize(size);
   kernels_.unpackCoalescence(dstField, ci, bufferPtr, dir, stream);
}

template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::communicateLocalFineToCoarse(
   const Block* fineSender, Block* coarseReceiver, walberla::stencil::Direction dir, gpuStream_t stream)
{
   auto varFineSender = const_cast< Block * >(fineSender);
   auto srcField   = varFineSender->getData< PdfField_T >(pdfFieldID_);
   auto srcCommData   = varFineSender->getData< CommData_T >(commDataID_);
   PartialCoalescenceMaskFieldGPU * maskField = &(srcCommData->getMaskFieldGPU());
   auto dstField = coarseReceiver->getData<PdfField_T>(pdfFieldID_);
   Direction invDir = stencil::inverseDir[dir];

   CellInterval srcInterval;
   srcField->getGhostRegion(dir, srcInterval, 2);
   uint_t packSize = kernels_.partialCoalescenceSize(srcInterval, dir);

   CellInterval dstInterval = getCoarseBlockCoalescenceInterval(coarseReceiver, fineSender->getId(),
                                                                invDir, dstField);

#ifndef NDEBUG
   uint_t unpackSize = kernels_.size(dstInterval, invDir);
   WALBERLA_ASSERT_EQUAL(packSize, unpackSize)
#endif

   // TODO: This is a dirty workaround. Code-generate direct redistribution!
   unsigned char *buffer;
   WALBERLA_GPU_CHECK( gpuMalloc( &buffer, packSize))
   kernels_.packPartialCoalescence(srcField, maskField, srcInterval, buffer, dir, stream);
   kernels_.unpackCoalescence(dstField, dstInterval, buffer, invDir, stream);
   WALBERLA_GPU_CHECK(gpuFree(buffer))
}


template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::communicateLocalFineToCoarse(
   const Block* fineSender, Block* coarseReceiver, walberla::stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream)
{
   auto varFineSender = const_cast< Block * >(fineSender);
   auto srcField   = varFineSender->getData< PdfField_T >(pdfFieldID_);
   auto srcCommData   = varFineSender->getData< CommData_T >(commDataID_);
   PartialCoalescenceMaskFieldGPU * maskField = &(srcCommData->getMaskFieldGPU());
   auto dstField = coarseReceiver->getData<PdfField_T>(pdfFieldID_);
   Direction invDir = stencil::inverseDir[dir];

   CellInterval srcInterval;
   srcField->getGhostRegion(dir, srcInterval, 2);
   uint_t packSize = kernels_.partialCoalescenceSize(srcInterval, dir);

   CellInterval dstInterval = getCoarseBlockCoalescenceInterval(coarseReceiver, fineSender->getId(),
                                                                invDir, dstField);

#ifndef NDEBUG
   uint_t unpackSize = kernels_.size(dstInterval, invDir);
   WALBERLA_ASSERT_EQUAL(packSize, unpackSize)
#endif

   auto bufferPtr = buffer.advanceNoResize(packSize);
   kernels_.packPartialCoalescence(srcField, maskField, srcInterval, bufferPtr, dir, stream);
   kernels_.unpackCoalescence(dstField, dstInterval, bufferPtr, invDir, stream);
}


template< typename PdfField_T>
uint_t NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::sizeEqualLevelSend( const Block * sender, stencil::Direction dir)
{
   auto field = const_cast< Block* >(sender)->getData< PdfField_T >(pdfFieldID_);
   CellInterval ci;
   cell_idx_t gls = skipsThroughCoarseBlock(sender, dir) ? 2 : 1;
   field->getSliceBeforeGhostLayer(dir, ci, gls, false);
   return kernels_.size(ci, dir);
}



template< typename PdfField_T>
uint_t NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::sizeCoarseToFineSend ( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir)
{
   auto field = const_cast< Block* >(coarseSender)->getData< PdfField_T >(pdfFieldID_);

   std::vector< std::pair< Direction, CellInterval > > intervals;
   getCoarseBlockCommIntervals(fineReceiver, dir, field, intervals);

   uint_t size = 0;

   for (auto t : intervals)
   {
      CellInterval ci          = t.second;
      size += kernels_.size(ci);
   }
   WALBERLA_ASSERT_GREATER(size, 0)
   return size;
}

template< typename PdfField_T>
uint_t NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::sizeCoarseToFineReceive ( Block* fineReceiver, stencil::Direction dir)
{

   auto field = fineReceiver->getData< PdfField_T >(pdfFieldID_);

   std::vector< std::pair< Direction, CellInterval > > intervals;
   getFineBlockCommIntervals(fineReceiver->getId(), dir, field, intervals);

   uint_t size = 0;
   for (auto t : intervals)
   {
      size += kernels_.redistributeSize(t.second);
   }
   WALBERLA_ASSERT_GREATER(size, 0)
   return size;
}



template< typename PdfField_T>
uint_t NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::sizeFineToCoarseSend ( const Block * sender, stencil::Direction dir)
{
    auto block      = const_cast< Block* >(sender);
    auto srcField   = block->getData< PdfField_T >(pdfFieldID_);

    CellInterval ci;
    srcField->getGhostRegion(dir, ci, 2);
    return kernels_.partialCoalescenceSize(ci, dir);
}



template< typename PdfField_T>
void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::packDataFineToCoarseImpl(
   const Block* fineSender, const walberla::BlockID& /*coarseReceiver*/, walberla::stencil::Direction dir,
   GpuBuffer_T & buffer, gpuStream_t stream) const
{
   auto varBlock = const_cast< Block* >(fineSender);
   auto srcField   = varBlock->getData< PdfField_T >(pdfFieldID_);
   auto commData  = varBlock->getData< CommData_T >(commDataID_);
   PartialCoalescenceMaskFieldGPU * maskField = &(commData->getMaskFieldGPU());

   CellInterval ci;
   srcField->getGhostRegion(dir, ci, 2);
   uint_t size = kernels_.partialCoalescenceSize(ci, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.packPartialCoalescence(srcField, maskField, ci, bufferPtr, dir, stream);
}

/***********************************************************************************************************************
 *                                                  Helper Functions                                                   *
 **********************************************************************************************************************/

template< typename PdfField_T>
inline Vector3< cell_idx_t >
   NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::getNeighborShift(const BlockID& fineBlock,
                                                                     stencil::Direction dir) const
{
   // dir: direction from coarse to fine block, or vice versa
   Vector3< cell_idx_t > shift;

   uint_t const branchId = fineBlock.getBranchId();

   shift[0] = (stencil::cx[dir] == 0) ? (((branchId & uint_t(1)) == uint_t(0)) ? cell_idx_t(-1) : cell_idx_t(1)) :
              cell_idx_t(0);
   shift[1] = (stencil::cy[dir] == 0) ? (((branchId & uint_t(2)) == uint_t(0)) ? cell_idx_t(-1) : cell_idx_t(1)) :
              cell_idx_t(0);
   shift[2] = (Stencil::D == uint_t(3)) ?
              ((stencil::cz[dir] == 0) ? (((branchId & uint_t(4)) == uint_t(0)) ? cell_idx_t(-1) : cell_idx_t(1)) :
               cell_idx_t(0)) :
              cell_idx_t(0);

   return shift;
}

/**
 * Returns the part of a cell interval's hull of given \p width in direction \p dirVec.
 * @param ci        The original cell interval
 * @param dirVec    Direction Vector
 * @param width     Width of the hull
 * @return          Interval forming the part of the hull
 */
template< typename PdfField_T>
inline CellInterval NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::intervalHullInDirection(
   const CellInterval& ci, const Vector3< cell_idx_t > dirVec, cell_idx_t width) const
{
   CellInterval result(ci);
   for (uint_t i = 0; i < Stencil::D; i++)
   {
      if (dirVec[i] == 1)
      {
         result.min()[i] = result.max()[i] + cell_idx_t(1);
         result.max()[i] += width;
      }
      if (dirVec[i] == -1)
      {
         result.max()[i] = result.min()[i] - cell_idx_t(1);
         result.min()[i] -= width;
      }
   }

   return result;
}

/**
 * For edge or corner directions, checks if a coarser block is part of the respective edge or corner intersection.
 * @param block The local block
 * @param dir   The direction to check
 * @return      `true`  if dir is an edge or corner direction skipping through a coarser block.
 */
template< typename PdfField_T>
inline bool NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::skipsThroughCoarseBlock(
   const Block* block, const Direction dir) const
{
   Vector3< cell_idx_t > dirVec(stencil::cx[dir], stencil::cy[dir], stencil::cz[dir]);
   bool coarseBlockFound = false;
   forEachSubdirectionCancel(dirVec, [&](Vector3< cell_idx_t > subdir) {
     coarseBlockFound =
        coarseBlockFound || block->neighborhoodSectionHasLargerBlock(
           blockforest::getBlockNeighborhoodSectionIndex(subdir[0], subdir[1], subdir[2]));
     return !coarseBlockFound;
   });

   return coarseBlockFound;
}

/**
 * For coarse-to-fine and fine-to-coarse communication, returns a list of pairs (Direction, CellInterval)
 * mapping sub-directions of the communication direction to cell intervals on the coarse block interior
 * whose data must be communicated <i>as if</i> communicating in those sub-directions.
 * @param fineBlockID   ID of the fine block
 * @param dir           Direction from the coarse to the fine block
 * @param field         Pointer to the PDF field on the coarse block
 * @param intervals     Vector that will be filled with the computed intervals
 */
template< typename PdfField_T>
inline void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::getCoarseBlockCommIntervals(
   const BlockID& fineBlockID, const Direction dir, const PdfField_T* field,
   std::vector< std::pair< Direction, CellInterval > >& intervals) const
{
   Vector3< cell_idx_t > shift = getNeighborShift(fineBlockID, dir);

   CellInterval mainSlice;
   field->getSliceBeforeGhostLayer(dir, mainSlice, 1, false);

   // In all directions, restrict the slice to the lower or upper half, depending on neighbor shift
   for (uint_t i = 0; i != Stencil::D; ++i)
   {
      if (shift[i] == cell_idx_t(-1))
      {
         WALBERLA_ASSERT_EQUAL(mainSlice.size(i) & 1, 0)
         mainSlice.max()[i] = mainSlice.min()[i] + cell_idx_c(mainSlice.size(i) / uint_t(2)) - cell_idx_t(1);
      }
      if (shift[i] == cell_idx_t(1))
      {
         WALBERLA_ASSERT_EQUAL(mainSlice.size(i) & 1, 0)
         mainSlice.min()[i] = mainSlice.min()[i] + cell_idx_c(mainSlice.size(i) / uint_t(2));
      }
   }

   intervals.emplace_back(dir, mainSlice);

   Vector3< cell_idx_t > const commDirVec{ stencil::cx[dir], stencil::cy[dir], stencil::cz[dir] };

   // Get extended slices in all tangential directions for the diagonal part of communication
   forEachSubdirection(-shift, [&](Vector3< cell_idx_t > t) {
     CellInterval hullInterval = intervalHullInDirection(mainSlice, t, cell_idx_t(1));
     Direction subCommDir      = stencil::vectorToDirection(commDirVec - t);
     if(CommunicationStencil::containsDir(subCommDir)){
        intervals.emplace_back(subCommDir, hullInterval);
     }
   });
}

/**
 * For coarse-to-fine and fine-to-coarse communication, returns a list of pairs (Direction, CellInterval)
 * mapping sub-directions of the communication direction to cell intervals on the fine block whose data must
 * be communicated <i>as if</i> communicating in those sub-directions.
 * @param fineBlockID   ID of the fine block
 * @param dir           Direction from the fine to the coarse block
 * @param field         Pointer to the PDF Field on the fine block
 * @param intervals     Vector that will be filled with the computed intervals
 */
template< typename PdfField_T>
inline void NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::getFineBlockCommIntervals(
   const BlockID& fineBlockID, const Direction dir, const PdfField_T* field,
   std::vector< std::pair< Direction, CellInterval > >& intervals) const
{
   Vector3< cell_idx_t > shift = getNeighborShift(fineBlockID, dir);

   CellInterval mainSlice;
   field->getGhostRegion(dir, mainSlice, 2, false);
   intervals.emplace_back(dir, mainSlice);

   Vector3< cell_idx_t > const commDirVec{ stencil::cx[dir], stencil::cy[dir], stencil::cz[dir] };

   forEachSubdirection(-shift, [&](Vector3< cell_idx_t > t) {
     CellInterval hullInterval = intervalHullInDirection(mainSlice, t, cell_idx_t(2));
     Direction subCommDir      = stencil::vectorToDirection(commDirVec + t);
     if(CommunicationStencil::containsDir(subCommDir)){
        intervals.emplace_back(subCommDir, hullInterval);
     }
   });
}
/**
 * Checks whether or not the block with ID `neighborID` is a neighbor of `block` in direction `dir`.
 */
template< typename PdfField_T>
bool NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::areNeighborsInDirection(
   const Block* block, const BlockID& neighborID, const Vector3< cell_idx_t> dirVec) const
{
   uint_t const nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(dirVec[0], dirVec[1], dirVec[2]);
   uint_t const nSecSize = block->getNeighborhoodSectionSize(nSecIdx);

   for(uint_t i = 0; i < nSecSize; i++){
      if(block->getNeighborId(nSecIdx, i) == neighborID){
         return true;
      }
   }
   return false;
}

template< typename PdfField_T>
CellInterval NonuniformGeneratedGPUPdfPackInfo< PdfField_T >::getCoarseBlockCoalescenceInterval(
   const Block* coarseBlock, const BlockID& fineBlockID, Direction dir, const PdfField_T* field) const
{
   Direction mainDir(dir);
   Vector3< cell_idx_t > commDirVec(stencil::cx[dir], stencil::cy[dir], stencil::cz[dir]);
   Vector3< cell_idx_t > mainDirVec(commDirVec);
   bool isAsymmetric = !areNeighborsInDirection(coarseBlock, fineBlockID, commDirVec);

   // If asymmetric, find the main subdirection
   if(isAsymmetric){
      mainDirVec = Vector3< cell_idx_t >(0);
      forEachSubdirection(commDirVec, [&](Vector3< cell_idx_t > subdirVec){
         if(areNeighborsInDirection(coarseBlock, fineBlockID, subdirVec)){
            // -dir is one main communication direction from F to C, but, due to periodicity,
            // it might not be the only one. Find the main comm direction from the subdirections
            // that is largest in the 1-norm.
            if(subdirVec.sqrLength() > mainDirVec.sqrLength()) mainDirVec = subdirVec;
         }
      });
      mainDir = stencil::vectorToDirection(mainDirVec);
   }

   Vector3< cell_idx_t > shift = getNeighborShift(fineBlockID, mainDir);

   CellInterval mainSlice;
   field->getSliceBeforeGhostLayer(mainDir, mainSlice, 1, false);

   // In all directions, restrict the slice to the lower or upper half, depending on neighbor shift
   for (uint_t i = 0; i != Stencil::D; ++i)
   {
      if (shift[i] == cell_idx_t(-1))
      {
         WALBERLA_ASSERT_EQUAL(mainSlice.size(i) & 1, 0)
         mainSlice.max()[i] = mainSlice.min()[i] + cell_idx_c(mainSlice.size(i) / uint_t(2)) - cell_idx_t(1);
      }
      if (shift[i] == cell_idx_t(1))
      {
         WALBERLA_ASSERT_EQUAL(mainSlice.size(i) & 1, 0)
         mainSlice.min()[i] = mainSlice.min()[i] + cell_idx_c(mainSlice.size(i) / uint_t(2));
      }
   }

   CellInterval commSlice(mainSlice);

   // If asymmetric, find coalescence slice as hull of main slice
   if(isAsymmetric){
      commSlice = intervalHullInDirection(mainSlice, mainDirVec - commDirVec, 1);
   }

   return commSlice;
}

} // walberla::lbm_generated
