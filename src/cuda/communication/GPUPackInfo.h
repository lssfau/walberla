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
//! \ingroup cuda
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//======================================================================================================================

#pragma once

#include "cuda/GPUCopy.h"
#include "communication/UniformPackInfo.h"
#include "core/debug/Debug.h"
#include "field/GhostRegions.h"
#include "stencil/Directions.h"


namespace walberla {
namespace cuda {
namespace communication {


/**
 * Data packing/unpacking for ghost layer based communication of a cuda::GPUField
 * \ingroup cuda
 * Template Parameters:
 *    - GPUField_T   A fully qualified GPUField.
 */
template<typename GPUField_T>
class GPUPackInfo : public walberla::communication::UniformPackInfo
{
public:
   typedef typename GPUField_T::value_type T;

   GPUPackInfo( const BlockDataID & bdId ) : bdId_( bdId ), communicateAllGhostLayers_( true ),
            numberOfGhostLayers_( 0 ) {}

   GPUPackInfo( const BlockDataID & bdId, const uint_t numberOfGHostLayers ) : bdId_( bdId ),
            communicateAllGhostLayers_( false ), numberOfGhostLayers_(  numberOfGHostLayers ) {}

   virtual ~GPUPackInfo() {}

   bool constantDataExchange() const { return mpi::BufferSizeTrait<T>::constantSize; }
   bool threadsafeReceiving()  const { return true; }

   void unpackData(IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer);

   void communicateLocal(const IBlock * sender, IBlock * receiver, stencil::Direction dir);

protected:
   void packDataImpl(const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer) const;

   uint_t numberOfGhostLayersToCommunicate( const GPUField_T * const field ) const;

   const BlockDataID bdId_;
   bool   communicateAllGhostLayers_;
   uint_t numberOfGhostLayers_;
};


template<typename GPUField_T>
void GPUPackInfo<GPUField_T>::unpackData(IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer)
{
   GPUField_T * f = receiver->getData< GPUField_T >( bdId_ );
   WALBERLA_ASSERT_NOT_NULLPTR(f);

   if ( f->layout() != fzyx )
   {
      WALBERLA_ABORT( "GPUPackInfo currently only supports fzyx layout" );
   }

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( f ) );

   CellInterval ci = field::getGhostRegion( *f, dir, nrOfGhostLayers, false );

   uint_t nrOfBytesToRead = ci.xSize() * ci.ySize() * ci.zSize() * f->fSize() * sizeof(T);

   unsigned char * buf = buffer.skip(nrOfBytesToRead);

   copyHostToDevFZYX( f->pitchedPtr(), buf,
                      sizeof(T), f->zAllocSize(), ci.zSize(),
                      uint_c(ci.xMin() + nrOfGhostLayers),
                      uint_c(ci.yMin() + nrOfGhostLayers),
                      uint_c(ci.zMin() + nrOfGhostLayers), 0,
                      0, 0, 0, 0,
                      ci.xSize(), ci.ySize(), ci.zSize(), f->fSize() );
}


template<typename GPUField_T>
void GPUPackInfo<GPUField_T>::communicateLocal(const IBlock * sender, IBlock * receiver, stencil::Direction dir)
{
   const GPUField_T * sf = sender  ->getData< GPUField_T >( bdId_ );
         GPUField_T * rf = receiver->getData< GPUField_T >( bdId_ );

   if ( sf->layout() != fzyx || rf->layout() != fzyx )
   {
      WALBERLA_ABORT( "GPUPackInfo currently only supports fzyx layout" );
   }

   WALBERLA_ASSERT_EQUAL(sf->xSize(), rf->xSize());
   WALBERLA_ASSERT_EQUAL(sf->ySize(), rf->ySize());
   WALBERLA_ASSERT_EQUAL(sf->zSize(), rf->zSize());
   WALBERLA_ASSERT_EQUAL(sf->fSize(), rf->fSize());

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( sf ) );

   CellInterval sCi = field::getSliceBeforeGhostLayer( *sf, dir, nrOfGhostLayers, false );
   CellInterval rCi = field::getGhostRegion( *rf, stencil::inverseDir[dir], nrOfGhostLayers, false );

   copyDevToDevFZYX( rf->pitchedPtr(), sf->pitchedPtr(),
                     sizeof(T), rf->zAllocSize(), sf->zAllocSize(),
                     uint_c(rCi.xMin() + nrOfGhostLayers), uint_c(rCi.yMin() + nrOfGhostLayers), uint_c(rCi.zMin() + nrOfGhostLayers), 0,
                     uint_c(sCi.xMin() + nrOfGhostLayers), uint_c(sCi.yMin() + nrOfGhostLayers), uint_c(sCi.zMin() + nrOfGhostLayers), 0,
                     rCi.xSize(), rCi.ySize(), rCi.zSize(), sf->fSize() );
}


template<typename GPUField_T>
void GPUPackInfo<GPUField_T>::packDataImpl(const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer) const
{
   const GPUField_T * f = sender->getData< GPUField_T >( bdId_ );
   WALBERLA_ASSERT_NOT_NULLPTR(f);

   if ( f->layout() != fzyx )
   {
      WALBERLA_ABORT( "GPUPackInfo currently only supports fzyx layout" );
   }

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( f ) );

   CellInterval ci = field::getSliceBeforeGhostLayer( *f, dir, nrOfGhostLayers, false );

   uint_t nrOfBytesToPack = ci.xSize() * ci.ySize() * ci.zSize() * f->fSize() * sizeof(T);

   unsigned char * buf = outBuffer.forward( nrOfBytesToPack );

   copyDevToHostFZYX( buf, f->pitchedPtr(),
                      sizeof(T), ci.zSize(), f->zAllocSize(),
                      0, 0, 0, 0,
                      uint_c(ci.xMin() + nrOfGhostLayers),
                      uint_c(ci.yMin() + nrOfGhostLayers),
                      uint_c(ci.zMin() + nrOfGhostLayers), 0,
                      ci.xSize(), ci.ySize(), ci.zSize(), f->fSize() );
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
      WALBERLA_ASSERT_LESS_EQUAL( numberOfGhostLayers_, field->nrOfGhostLayers() );
      return numberOfGhostLayers_;
   }
}



} // namespace communication
} // namespace cuda
} // namespace walberla
