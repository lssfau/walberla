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
//! \file FieldIndexing.impl.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "FieldIndexing.h"
#include "GPUField.h"

#include "core/cell/CellInterval.h"
#include "core/debug/Debug.h"
#include "field/Layout.h"

#include <limits>
#include <cmath>

namespace walberla {
namespace gpu
{

template< typename T>
FieldIndexing<T>::FieldIndexing ( const GPUField<T> & field,
                                  dim3 _blockDim, dim3 _gridDim,
                                  const FieldAccessor<T> _gpuAccess )
   : field_( field ),
     blockDim_( _blockDim ),
     gridDim_( _gridDim ),
     gpuAccess_( _gpuAccess )
{
   WALBERLA_DEBUG_SECTION()
   {
      WALBERLA_DEVICE_SECTION()
      {
         gpuDeviceProp prop;
         int count;
         WALBERLA_GPU_CHECK(gpuGetDeviceCount(&count));
         int threadsPerBlock = std::numeric_limits< int >::max();
         for (int i = 0; i < count; i++)
         {
            WALBERLA_GPU_CHECK(gpuGetDeviceProperties(&prop, i));
            threadsPerBlock = std::min(prop.maxThreadsPerBlock, threadsPerBlock);
         }
         WALBERLA_ASSERT_LESS(int_c(blockDim_.x), threadsPerBlock,
                              "InnerCoordThreadIndexing works only for fields where each dimension x,y,z is smaller "
                                 << "than the maximal thread count per GPU block.")
      }
   }
}


template< typename T>
void shiftCoordinatesWhileFastestCoordHasSizeOne( typename FieldAccessor<T>::IndexingScheme & indexing, dim3 & gridDim, dim3 & blockDim )
{
   bool runLoop = true;
   while( blockDim.x == 1 && runLoop )
   {
      blockDim.x = gridDim.x;
      gridDim.x  = gridDim.y;
      gridDim.y  = gridDim.z;
      gridDim.z  = 1;


      switch ( indexing ) {
         case FieldAccessor<T>::FZYX: indexing = FieldAccessor<T>::FZY; break;
         case FieldAccessor<T>::FZY : indexing = FieldAccessor<T>::FZ ; break;
         case FieldAccessor<T>::FZ  : indexing = FieldAccessor<T>::F  ; break;

         case FieldAccessor<T>::ZYXF: indexing = FieldAccessor<T>::ZYX ; break;
         case FieldAccessor<T>::ZYX : indexing = FieldAccessor<T>::ZY  ; break;
         case FieldAccessor<T>::ZY  : indexing = FieldAccessor<T>::Z   ; break;

         // iteration goes only over a single element - stop the loop
         case FieldAccessor<T>::Z: runLoop = false; break;
         case FieldAccessor<T>::F: runLoop = false; break;
      }
   }
}

template< typename T>
FieldIndexing<T> FieldIndexing<T>::interval ( const GPUField<T> & f, const CellInterval & ci, int fBegin, int fEnd )
{
   uint_t xOffset;
   uint_t yOffset;
   uint_t zOffset;
   uint_t fOffset;

   if ( f.layout() == field::zyxf )
   {
      fOffset = sizeof(T);
      xOffset = f.pitchedPtr().pitch;
      yOffset = xOffset * f.xAllocSize();
      zOffset = yOffset * f.yAllocSize();
   }
   else
   {
      xOffset = sizeof(T);
      yOffset = f.pitchedPtr().pitch;
      zOffset = yOffset * f.yAllocSize();
      fOffset = zOffset * f.zAllocSize();
   }
   char * data = (char*)f.pitchedPtr().ptr;

   // Jump over ghost cells to first inner cell
   cell_idx_t gl = cell_idx_c( f.nrOfGhostLayers() );
   data += ( ci.xMin() + gl )* cell_idx_c(xOffset) +
           ( ci.yMin() + gl )* cell_idx_c(yOffset) +
           ( ci.zMin() + gl )* cell_idx_c(zOffset);


   dim3 gridDim;
   dim3 blockDim;
   typename FieldAccessor<T>::IndexingScheme firstCoord;
   if ( f.layout() == fzyx )
   {
      firstCoord = FieldAccessor<T>::FZYX;
      blockDim = dim3 ( (unsigned int)ci.xSize(), 1, 1 );
      gridDim  = dim3 ( (unsigned int)ci.ySize(), (unsigned int)ci.zSize(), (unsigned int)( fEnd - fBegin) );
   }
   else
   {
      firstCoord = FieldAccessor<T>::ZYXF;
      blockDim = dim3 ( (unsigned int)(fEnd - fBegin), 1, 1 );
      gridDim  = dim3 ( (unsigned int)ci.xSize(), (unsigned int)ci.ySize(), (unsigned int)ci.zSize() );
   }

   shiftCoordinatesWhileFastestCoordHasSizeOne<T>( firstCoord, gridDim, blockDim );

   return FieldIndexing<T> ( f, blockDim, gridDim,
                             FieldAccessor<T>( data, xOffset, yOffset, zOffset, fOffset, firstCoord ) );
}


template< typename T>
FieldIndexing<T> FieldIndexing<T>::xyz ( const GPUField<T> & f )
{
   CellInterval ci ( 0,0,0,
                    cell_idx_c( f.xSize() ) - 1,
                    cell_idx_c( f.ySize() ) - 1,
                    cell_idx_c( f.zSize() ) - 1 );

   return interval( f, ci, 0,1 );
}


template< typename T>
FieldIndexing<T> FieldIndexing<T>::withGhostLayerXYZ( const GPUField<T> & f, uint_t numGhostLayers )
{
   cell_idx_t gl = std::min( cell_idx_c(numGhostLayers), cell_idx_c( f.nrOfGhostLayers() ) );
   CellInterval ci ( -gl,-gl,-gl,
                    cell_idx_c( f.xSize() ) + gl - 1,
                    cell_idx_c( f.ySize() ) + gl - 1,
                    cell_idx_c( f.zSize() ) + gl - 1 );

   return interval( f, ci, 0, 1 );
}

template< typename T>
FieldIndexing<T> FieldIndexing<T>::ghostLayerOnlyXYZ( const GPUField<T> & f, uint_t thickness,
                                                      stencil::Direction dir, bool fullSlice )
{
   CellInterval ci;
   f.getGhostRegion( dir, ci, cell_idx_c(thickness), fullSlice );
   return interval( f, ci, 0, 1 );
}

template< typename T>
FieldIndexing<T> FieldIndexing<T>::sliceBeforeGhostLayerXYZ( const GPUField<T> & f, uint_t thickness,
                                                             stencil::Direction dir, bool fullSlice )
{
   CellInterval ci;
   f.getSliceBeforeGhostLayer( dir, ci, cell_idx_c(thickness), fullSlice );
   return interval( f, ci, 0, 1 );
}

template< typename T>
FieldIndexing<T> FieldIndexing<T>::sliceXYZ( const GPUField<T> & f, cell_idx_t distance, uint_t thickness,
                                             stencil::Direction dir, bool fullSlice )
{
   CellInterval ci;
   f.getSlice( dir, ci, distance, cell_idx_c(thickness), fullSlice );
   return interval( f, ci );
}

template< typename T>
FieldIndexing<T> FieldIndexing<T>::allInner ( const GPUField<T> & f )
{
   CellInterval ci ( 0,0,0,
                    cell_idx_c( f.xSize() ) - 1,
                    cell_idx_c( f.ySize() ) - 1,
                    cell_idx_c( f.zSize() ) - 1 );

   return interval( f, ci, 0, cell_idx_c( f.fSize() ) );
}

template< typename T>
FieldIndexing<T> FieldIndexing<T>::allWithGhostLayer ( const GPUField<T> & f )
{
   cell_idx_t gl = cell_idx_c( f.nrOfGhostLayers() );
   CellInterval ci ( -gl,-gl,-gl,
                    cell_idx_c( f.xSize() ) + gl - 1,
                    cell_idx_c( f.ySize() ) + gl - 1,
                    cell_idx_c( f.zSize() ) + gl - 1 );

   return interval( f, ci, 0, cell_idx_c( f.fSize() ) );

}

template< typename T>
FieldIndexing<T> FieldIndexing<T>::all ( const GPUField<T> & f, const cell::CellInterval & ci )
{
   return interval( f, ci, 0, cell_idx_c( f.fSize() ) );
}



} // namespace gpu
} // namespace walberla


