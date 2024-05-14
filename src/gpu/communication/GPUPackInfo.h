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
//! \file GPUPackInfo.h
//! \ingroup gpu
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \author João Victor Tozatti Risso <jvtrisso@inf.ufpr.br>
//======================================================================================================================

#pragma once

#include "blockforest/Block.h"

#include "communication/UniformPackInfo.h"

#include "core/debug/Debug.h"
#include "core/math/Vector3.h"
#include "core/mpi/BufferSizeTrait.h"

#include "field/GhostRegions.h"
#include "field/Layout.h"

#include "stencil/Directions.h"

#include <map>
#include <tuple>
#include <vector>

#include "gpu/ErrorChecking.h"
#include "gpu/GPUCopy.h"
#include "gpu/GPUWrapper.h"
#include "gpu/communication/CustomMemoryBuffer.h"

namespace walberla::gpu::communication {


/**
 * \brief Data packing/unpacking for ghost layer based communication of a \ref GPUField.
 *
 * Encapsulate information on how to extract data from blocks that should be
 * communicated to neighboring blocks (see \ref packDataImpl())
 * and how to inject this data in a receiving block (see \ref unpackData()).
 * This involves a host memory buffer and two device-to-host memory copy operations.
 *
 * A special method exists for communication between two blocks which are
 * allocated on the same process (see \ref communicateLocal()).
 * In this case the data does not have be communicated via a host buffer,
 * but can be sent directly. This involves a single device-to-device memory
 * copy operation.
 *
 * Data that is packed in direction "dir" at one block is unpacked in
 * direction "stencil::inverseDir[dir]" at the neighboring block.
 * This behavior must be implemented in \ref communicateLocal()!
 *
 * See \ref MemcpyPackInfo for a more efficient packing/unpacking method
 * where the buffer is stored in device memory rather than in host memory.
 *
 * \ingroup gpu
 * \tparam GPUField_T   A fully qualified \ref GPUField.
 */
template<typename GPUField_T>
class GPUPackInfo : public walberla::communication::UniformPackInfo
{
public:
   typedef typename GPUField_T::value_type FieldType;

   GPUPackInfo( const BlockDataID & bdId )
   : bdId_( bdId ), communicateAllGhostLayers_( true ), numberOfGhostLayers_( 0 ),
     copyAsync_( false ), communicationStream_( 0 )
   {
   }

   GPUPackInfo( const BlockDataID & bdId, const uint_t numberOfGHostLayers )
   : bdId_( bdId ), communicateAllGhostLayers_( false ), numberOfGhostLayers_(  numberOfGHostLayers ),
     copyAsync_( false ), communicationStream_( 0 )
   {
   }

   virtual ~GPUPackInfo() {}

   bool constantDataExchange() const { return mpi::BufferSizeTrait<FieldType>::constantSize; }
   bool threadsafeReceiving()  const { return true; }

   void unpackData(IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer);

   void communicateLocal(const IBlock * sender, IBlock * receiver, stencil::Direction dir);

   void setCommunicationStream( gpuStream_t stream )
   {
      if ( stream != 0 )
      {
         copyAsync_ = true;
         communicationStream_ = stream;
      }
   }

protected:
   void packDataImpl(const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer) const;

   uint_t numberOfGhostLayersToCommunicate( const GPUField_T * const field ) const;

   const BlockDataID bdId_;
   bool   communicateAllGhostLayers_;
   uint_t numberOfGhostLayers_;
   bool copyAsync_;
   gpuStream_t communicationStream_;
   std::map< stencil::Direction, PinnedMemoryBuffer > pinnedRecvBuffers_;
   mutable std::map< stencil::Direction, PinnedMemoryBuffer > pinnedSendBuffers_;
};


template<typename GPUField_T>
void GPUPackInfo<GPUField_T>::unpackData(IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer)
{
   GPUField_T * fieldPtr = receiver->getData< GPUField_T >( bdId_ );
   WALBERLA_ASSERT_NOT_NULLPTR(fieldPtr)

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( fieldPtr ) );

   CellInterval fieldCi = field::getGhostRegion( *fieldPtr, dir, nrOfGhostLayers, false );

   uint_t nrOfBytesToRead = fieldCi.numCells() * fieldPtr->fSize() * sizeof(FieldType);

   unsigned char * bufPtr = buffer.skip(nrOfBytesToRead);

   unsigned char * copyBufferPtr = bufPtr;
   if ( copyAsync_ )
   {
      PinnedMemoryBuffer & pinnedBuffer = pinnedRecvBuffers_[dir];
      pinnedBuffer.clear();
      copyBufferPtr = pinnedBuffer.advance( nrOfBytesToRead );
      // Copy data into pinned memory buffer, in order to transfer it asynchronously to the GPU
      std::copy( bufPtr, static_cast< unsigned char * >( bufPtr + nrOfBytesToRead ), copyBufferPtr );
   }

   gpuStream_t & unpackStream = communicationStream_;

   auto dstOffset = std::make_tuple( uint_c(fieldCi.xMin() + nrOfGhostLayers),
                                     uint_c(fieldCi.yMin() + nrOfGhostLayers),
                                     uint_c(fieldCi.zMin() + nrOfGhostLayers),
                                     uint_c(0) );
   auto srcOffset = std::make_tuple( uint_c(0), uint_c(0), uint_c(0), uint_c(0) );

   auto intervalSize = std::make_tuple( fieldCi.xSize(), fieldCi.ySize(), fieldCi.zSize(),
                                        fieldPtr->fSize() );

   if ( fieldPtr->layout() == field::fzyx )
   {
      const uint_t dstAllocSizeZ = fieldPtr->zAllocSize();
      const uint_t srcAllocSizeZ = fieldCi.zSize();
      copyHostToDevFZYX( fieldPtr->pitchedPtr(), copyBufferPtr, dstOffset, srcOffset,
                         dstAllocSizeZ, srcAllocSizeZ, sizeof(FieldType),
                         intervalSize, unpackStream );
   }
   else
   {
      const uint_t dstAllocSizeY = fieldPtr->yAllocSize();
      const uint_t srcAllocSizeY = fieldCi.ySize();
      copyHostToDevZYXF( fieldPtr->pitchedPtr(), copyBufferPtr, dstOffset, srcOffset,
                         dstAllocSizeY, srcAllocSizeY, sizeof(FieldType),
                         intervalSize, unpackStream );
   }

   if ( copyAsync_ )
   {
      WALBERLA_GPU_CHECK( gpuStreamSynchronize( unpackStream ) );
   }
}


template<typename GPUField_T>
void GPUPackInfo<GPUField_T>::communicateLocal(const IBlock * sender, IBlock * receiver, stencil::Direction dir)
{
   const GPUField_T * sf = sender  ->getData< GPUField_T >( bdId_ );
         GPUField_T * rf = receiver->getData< GPUField_T >( bdId_ );

   WALBERLA_ASSERT_NOT_NULLPTR( sf )
   WALBERLA_ASSERT_NOT_NULLPTR( rf )

   WALBERLA_ASSERT_EQUAL(sf->xSize(), rf->xSize())
   WALBERLA_ASSERT_EQUAL(sf->ySize(), rf->ySize())
   WALBERLA_ASSERT_EQUAL(sf->zSize(), rf->zSize())
   WALBERLA_ASSERT_EQUAL(sf->fSize(), rf->fSize())

   WALBERLA_CHECK( sf->layout() == rf->layout(), "GPUPackInfo::communicateLocal: fields must have the same layout!" );

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( sf ) );

   CellInterval sCi = field::getSliceBeforeGhostLayer( *sf, dir, nrOfGhostLayers, false );
   CellInterval rCi = field::getGhostRegion( *rf, stencil::inverseDir[dir], nrOfGhostLayers, false );

   gpuStream_t & commStream = communicationStream_;

   auto dstOffset = std::make_tuple( uint_c(rCi.xMin() + nrOfGhostLayers),
                                     uint_c(rCi.yMin() + nrOfGhostLayers),
                                     uint_c(rCi.zMin() + nrOfGhostLayers),
                                     uint_c(0) );

   auto srcOffset = std::make_tuple( uint_c(sCi.xMin() + nrOfGhostLayers),
                                     uint_c(sCi.yMin() + nrOfGhostLayers),
                                     uint_c(sCi.zMin() + nrOfGhostLayers),
                                     uint_c(0) );

   auto intervalSize = std::make_tuple( rCi.xSize(), rCi.ySize(), rCi.zSize(), sf->fSize() );

   if ( sf->layout() == field::fzyx )
   {
      const uint_t dstAllocSizeZ = rf->zAllocSize();
      const uint_t srcAllocSizeZ = sf->zAllocSize();

      copyDevToDevFZYX( rf->pitchedPtr(), sf->pitchedPtr(), dstOffset, srcOffset,
                        dstAllocSizeZ, srcAllocSizeZ, sizeof(FieldType),
                        intervalSize, commStream );
   }
   else
   {
      const uint_t dstAllocSizeY = rf->yAllocSize();
      const uint_t srcAllocSizeY = sf->yAllocSize();

      copyDevToDevZYXF( rf->pitchedPtr(), sf->pitchedPtr(), dstOffset, srcOffset,
                        dstAllocSizeY, srcAllocSizeY, sizeof(FieldType),
                        intervalSize, commStream );
   }

   if ( copyAsync_ )
   {
      WALBERLA_GPU_CHECK( gpuStreamSynchronize( commStream ) )
   }
}


template<typename GPUField_T>
void GPUPackInfo<GPUField_T>::packDataImpl(const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer) const
{
   const GPUField_T * fieldPtr = sender->getData< GPUField_T >( bdId_ );
   WALBERLA_ASSERT_NOT_NULLPTR(fieldPtr)

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( fieldPtr ) );

   CellInterval fieldCi = field::getSliceBeforeGhostLayer( *fieldPtr, dir, nrOfGhostLayers, false );

   const std::size_t nrOfBytesToPack = fieldCi.numCells() * fieldPtr->fSize() * sizeof(FieldType);

   unsigned char * outBufferPtr = outBuffer.forward( nrOfBytesToPack );

   const gpuStream_t & packStream = communicationStream_;

   unsigned char * copyBufferPtr = outBufferPtr;
   if ( copyAsync_ )
   {
      PinnedMemoryBuffer & pinnedBuffer = pinnedSendBuffers_[dir];
      pinnedBuffer.clear();
      copyBufferPtr = pinnedBuffer.advance( nrOfBytesToPack );
   }

   auto dstOffset = std::make_tuple( uint_c(0), uint_c(0), uint_c(0), uint_c(0) );
   auto srcOffset = std::make_tuple( uint_c(fieldCi.xMin() + nrOfGhostLayers),
                                     uint_c(fieldCi.yMin() + nrOfGhostLayers),
                                     uint_c(fieldCi.zMin() + nrOfGhostLayers),
                                     uint_c(0) );

   auto intervalSize = std::make_tuple( fieldCi.xSize(), fieldCi.ySize(), fieldCi.zSize(),
                                        fieldPtr->fSize() );

   if ( fieldPtr->layout() == field::fzyx )
   {
      const uint_t dstAllocSizeZ = fieldCi.zSize();
      const uint_t srcAllocSizeZ = fieldPtr->zAllocSize();
      copyDevToHostFZYX( copyBufferPtr, fieldPtr->pitchedPtr(), dstOffset, srcOffset,
                         dstAllocSizeZ, srcAllocSizeZ, sizeof(FieldType),
                         intervalSize, packStream );
   }
   else
   {
      const uint_t dstAllocSizeZ = fieldCi.ySize();
      const uint_t srcAllocSizeZ = fieldPtr->yAllocSize();
      copyDevToHostZYXF( copyBufferPtr, fieldPtr->pitchedPtr(), dstOffset, srcOffset,
                         dstAllocSizeZ, srcAllocSizeZ, sizeof(FieldType),
                         intervalSize, packStream );
   }

   if ( copyAsync_ )
   {
      WALBERLA_GPU_CHECK( gpuStreamSynchronize( packStream ) )

      std::copy( copyBufferPtr, static_cast<unsigned char *>( copyBufferPtr + nrOfBytesToPack ), outBufferPtr );
   }
}


template<typename GPUField_T>
uint_t GPUPackInfo<GPUField_T>::numberOfGhostLayersToCommunicate( const GPUField_T * const field ) const
{
   if( communicateAllGhostLayers_ )
   {
      return field->nrOfGhostLayers();
   }
   else
   {
      WALBERLA_ASSERT_LESS_EQUAL( numberOfGhostLayers_, field->nrOfGhostLayers() )
      return numberOfGhostLayers_;
   }
}

} // namespace walberla::gpu::communication
