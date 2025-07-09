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
//! \file GeneratedNonUniformGPUFieldPackInfo.impl.h
//! \ingroup gpu
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#pragma once

#include "GeneratedNonUniformGPUFieldPackInfo.h"




namespace walberla {

using stencil::Direction;

/***********************************************************************************************************************
 *                                          Equal Level Communication                                                  *
 **********************************************************************************************************************/

template< typename Field_T, typename PackingKernels_T >
void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::packDataEqualLevelImpl(
   const Block* sender, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream) const
{
   auto field = const_cast< Block* >(sender)->getData< Field_T >(fieldID_);
   CellInterval packInterval = field::refinement::equalLevelPackInterval( dir, field->xyzSize(), uint_t(1) );
   uint_t size = sizeEqualLevelSend(sender, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.packEqual(field, packInterval, bufferPtr, stream);
}

template< typename Field_T, typename PackingKernels_T >
void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::unpackDataEqualLevel(Block* receiver, Direction dir, GpuBuffer_T & buffer, gpuStream_t stream)
{
   auto field = receiver->getData< Field_T >(fieldID_);
   CellInterval unpackInterval = field::refinement::equalLevelUnpackInterval(dir, field->xyzSize(), uint_t(1) );
   uint_t size = sizeEqualLevelSend(receiver, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.unpackEqual(field, unpackInterval, bufferPtr, stream);
}

template< typename Field_T, typename PackingKernels_T >
void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::communicateLocalEqualLevel(
   const Block* sender, Block* receiver, stencil::Direction dir, gpuStream_t stream)
{
   auto sendField = const_cast< Block* >(sender)->getData< Field_T >(fieldID_);
   auto receiveField = receiver->getData< Field_T >(fieldID_);
   CellInterval packInterval = field::refinement::equalLevelPackInterval( dir, sendField->xyzSize(), uint_t(1) );
   CellInterval unpackInterval = field::refinement::equalLevelUnpackInterval(stencil::inverseDir[dir], receiveField->xyzSize(), uint_t(1) );
   kernels_.localCopyEqual(sendField, packInterval, receiveField, unpackInterval, stream);
}


/***********************************************************************************************************************
 *                                          Coarse to Fine Communication                                               *
 **********************************************************************************************************************/

template< typename Field_T, typename PackingKernels_T >
void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::packDataCoarseToFineImpl(
   const Block* coarseSender, const BlockID& fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream) const
{
   auto field = const_cast< Block* >(coarseSender)->getData< Field_T >(fieldID_);
   CellInterval packInterval = field::refinement::coarseToFinePackInterval( dir, field->xyzSize(), fineReceiver );
   auto size = sizeCoarseToFineSend(coarseSender, fineReceiver, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.packEqual(field, packInterval, bufferPtr, stream);
}

template< typename Field_T, typename PackingKernels_T >
void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::unpackDataCoarseToFine(
   Block* fineReceiver, const BlockID& /*coarseSender*/, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream)
{
   auto field = fineReceiver->getData< Field_T >( fieldID_ );
   CellInterval unpackInterval = field::refinement::coarseToFineUnpackInterval( dir, field->xyzSize(), fineReceiver->getId() );
   auto size = sizeCoarseToFineReceive(fineReceiver, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.unpackCoarseToFine(field, unpackInterval, bufferPtr, stream);
}

template< typename Field_T, typename PackingKernels_T >
void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::communicateLocalCoarseToFine(
   const Block* coarseSender, Block* fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream)
{
   auto srcField = const_cast< Block* >(coarseSender)->getData< Field_T >(fieldID_);
   CellInterval packInterval = field::refinement::coarseToFinePackInterval( dir, srcField->xyzSize(),  fineReceiver->getId() );
   auto dstField = fineReceiver->getData< Field_T >( fieldID_ );
   CellInterval unpackInterval = field::refinement::coarseToFineUnpackInterval( stencil::inverseDir[dir], dstField->xyzSize(), fineReceiver->getId() );
   auto size = sizeCoarseToFineSend(coarseSender, fineReceiver->getId(), dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.packEqual(srcField, packInterval, bufferPtr, stream);
   kernels_.unpackCoarseToFine(dstField, unpackInterval, bufferPtr, stream);
}


/***********************************************************************************************************************
 *                                          Fine to Coarse Communication                                               *
 **********************************************************************************************************************/

template< typename Field_T, typename PackingKernels_T >
void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::packDataFineToCoarseImpl(
   const Block* fineSender, const walberla::BlockID& coarseReceiver, walberla::stencil::Direction dir,
   GpuBuffer_T & buffer, gpuStream_t stream) const
{
   if( ( ( field::refinement::isEdgeDirection(dir) || field::refinement::isCornerDirection(dir) ) &&
        field::refinement::blocksConnectedByFaces( fineSender, coarseReceiver ) ) ||
       ( field::refinement::isCornerDirection(dir) && field::refinement::blocksConnectedByEdges( fineSender, coarseReceiver ) ) )
      return;

   auto field = const_cast< Block* >(fineSender)->getData< Field_T >(fieldID_);
   CellInterval packInterval = field::refinement::fineToCoarsePackInterval( dir, field->xyzSize() );
   auto size = sizeFineToCoarseSend(fineSender, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.packFineToCoarse(field, packInterval, bufferPtr, stream);
}


template< typename Field_T, typename PackingKernels_T >
   void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::unpackDataFineToCoarse(
      Block* coarseReceiver, const walberla::BlockID& fineSender, walberla::stencil::Direction dir,
      GpuBuffer_T & buffer, gpuStream_t stream)
{
   if( ( ( field::refinement::isEdgeDirection(dir) || field::refinement::isCornerDirection(dir) ) &&
        field::refinement::blocksConnectedByFaces( coarseReceiver, fineSender ) ) ||
       ( field::refinement::isCornerDirection(dir) && field::refinement::blocksConnectedByEdges( coarseReceiver, fineSender ) ) )
      return;

   auto field = coarseReceiver->getData< Field_T >( fieldID_ );
   CellInterval unpackInterval = field::refinement::fineToCoarseUnpackInterval( dir, field->xyzSize(), fineSender );
   auto size = sizeFineToCoarseSend(coarseReceiver, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.unpackEqual(field, unpackInterval, bufferPtr, stream);
}


template< typename Field_T, typename PackingKernels_T >
void GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::communicateLocalFineToCoarse(
   const Block* fineSender, Block* coarseReceiver, walberla::stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream)
{
   if( ( ( field::refinement::isEdgeDirection(dir) || field::refinement::isCornerDirection(dir) ) &&
        field::refinement::blocksConnectedByFaces( fineSender, coarseReceiver->getId() ) ) ||
       ( field::refinement::isCornerDirection(dir) && field::refinement::blocksConnectedByEdges( fineSender, coarseReceiver->getId() ) ) )
      return;

   auto srcField = const_cast< Block* >(fineSender)->getData< Field_T >(fieldID_);
   CellInterval packInterval = field::refinement::fineToCoarsePackInterval( dir, srcField->xyzSize() );
   auto dstField = coarseReceiver->getData< Field_T >( fieldID_ );
   CellInterval unpackInterval = field::refinement::fineToCoarseUnpackInterval( stencil::inverseDir[dir], dstField->xyzSize(), fineSender->getId() );
   auto size = sizeFineToCoarseSend(fineSender, dir);
   auto bufferPtr = buffer.advanceNoResize(size);
   kernels_.packFineToCoarse(srcField, packInterval, bufferPtr, stream);
   kernels_.unpackEqual(dstField, unpackInterval, bufferPtr, stream);
}


template< typename Field_T, typename PackingKernels_T >
uint_t GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::sizeEqualLevelSend( const Block * sender, stencil::Direction dir) const
{
   auto field = sender->getData< Field_T >(fieldID_);
   CellInterval ci = field::refinement::equalLevelPackInterval( dir, field->xyzSize(), uint_t(1) );
   return ci.numCells() * field->fSize() * sizeof(value_type);
}



template< typename Field_T, typename PackingKernels_T >
uint_t GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::sizeCoarseToFineSend ( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir) const
{
   auto field = coarseSender->getData< Field_T >(fieldID_);
   CellInterval ci = field::refinement::coarseToFinePackInterval( dir, field->xyzSize(), fineReceiver );
   return ci.numCells() * field->fSize() * sizeof(value_type);
}

template< typename Field_T, typename PackingKernels_T >
uint_t GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::sizeCoarseToFineReceive ( Block* fineReceiver, stencil::Direction dir) const
{
   auto field = fineReceiver->getData< Field_T >(fieldID_);
   CellInterval ci = field::refinement::coarseToFineUnpackInterval( dir, field->xyzSize(), uint_t(1) );
   return (ci.numCells() >> 2) * field->fSize() * sizeof(value_type);
}


template< typename Field_T, typename PackingKernels_T >
uint_t GeneratedNonUniformGPUFieldPackInfo< Field_T, PackingKernels_T >::sizeFineToCoarseSend ( const Block * sender, stencil::Direction dir) const
{
   auto field = sender->getData< Field_T >(fieldID_);
   CellInterval ci = field::refinement::fineToCoarsePackInterval( dir, field->xyzSize() );
   return (ci.numCells() >> 2) * field->fSize() * sizeof(value_type);
}

} //namespace walberla