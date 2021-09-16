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
//! \file FieldIndexing3D.impl.h
//! \ingroup cuda
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//
//======================================================================================================================

#include "FieldIndexing3D.h"
#include "GPUField.h"

#include "core/cell/CellInterval.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "field/Layout.h"

#include <cuda_runtime.h>

#include <limits>
#include <cmath>

namespace walberla {
namespace cuda {

// Returns ( a % b != 0 ) ? ( a / b + 1 ) : ( a / b )
inline unsigned int iDivUp( unsigned int a, unsigned int b ) { return ( a + b - 1 ) / b; }


dim3 FieldIndexing3DBase::preferredBlockDim_( 32, 2, 2 );


template< typename T>
FieldIndexing3D<T>::FieldIndexing3D( const GPUField<T> & field,
                                     const dim3 & _blockDim, const dim3 & _gridDim,
                                     const FieldAccessor3D<T> & _gpuAccess )
   : field_( field ),
     blockDim_( _blockDim ),
     gridDim_( _gridDim ),
     gpuAccess_( _gpuAccess )
{
   WALBERLA_DEBUG_SECTION()
   {
      cudaDeviceProp prop;
      int count;
      cudaGetDeviceCount(&count);
      int threadsPerBlock = std::numeric_limits<int>::max();
      for (int i = 0; i < count; i++) {
         cudaGetDeviceProperties(&prop, i);
         threadsPerBlock = std::min( prop.maxThreadsPerBlock, threadsPerBlock );
      }
      WALBERLA_ASSERT_LESS( int_c( blockDim_.x ), threadsPerBlock,
                            "InnerCoordThreadIndexing works only for fields where each dimension x,y,z is smaller " <<
                            "than the maximal thread count per CUDA block." );
   }
}

template< typename T>
FieldIndexing3D<T> FieldIndexing3D<T>::interval( const GPUField<T> & f, const CellInterval & ci )
{
   uint_t xOffset, yOffset, zOffset, fOffset;

   WALBERLA_ASSERT( f.layout() == field::fzyx );

   xOffset = sizeof(T);
   yOffset = f.pitchedPtr().pitch;
   zOffset = yOffset * f.yAllocSize();
   fOffset = zOffset * f.zAllocSize();

   char * data = (char*)f.pitchedPtr().ptr;

   // position data according to ci
   cell_idx_t gl = cell_idx_c( f.nrOfGhostLayers() );
   data += ( ci.xMin() + gl ) * cell_idx_c(xOffset) +
           ( ci.yMin() + gl ) * cell_idx_c(yOffset) +
           ( ci.zMin() + gl ) * cell_idx_c(zOffset);


   dim3 idxDim( (unsigned int)ci.xSize(), (unsigned int)ci.ySize(), (unsigned int)ci.zSize() );
   unsigned int bx = std::min( preferredBlockDim_.x, idxDim.x );
   unsigned int by = std::min( preferredBlockDim_.y, idxDim.y );
   unsigned int bz = std::min( preferredBlockDim_.z, idxDim.z );
   dim3 gridDim( iDivUp( idxDim.x, bx ),
                 iDivUp( idxDim.y, by ),
                 iDivUp( idxDim.z, bz ) );
   dim3 blockDim( bx, by, bz );

   return FieldIndexing3D<T>( f, blockDim, gridDim,
                              FieldAccessor3D<T>( data, xOffset, yOffset, zOffset, fOffset,
                                                  idxDim, blockDim ) );
}

template< typename T>
FieldIndexing3D<T> FieldIndexing3D<T>::xyz ( const GPUField<T> & f )
{
   CellInterval ci( 0,0,0,
                    cell_idx_c( f.xSize() ) - 1,
                    cell_idx_c( f.ySize() ) - 1,
                    cell_idx_c( f.zSize() ) - 1 );

   return interval( f, ci );
}

template< typename T>
FieldIndexing3D<T> FieldIndexing3D<T>::withGhostLayerXYZ( const GPUField<T> & f, uint_t numGhostLayers )
{
   cell_idx_t gl = std::min( cell_idx_c( numGhostLayers ), cell_idx_c( f.nrOfGhostLayers() ) );
   CellInterval ci( -gl,-gl,-gl,
                    cell_idx_c( f.xSize() ) + gl - 1,
                    cell_idx_c( f.ySize() ) + gl - 1,
                    cell_idx_c( f.zSize() ) + gl - 1 );

   return interval( f, ci );
}

template< typename T>
FieldIndexing3D<T> FieldIndexing3D<T>::ghostLayerOnlyXYZ( const GPUField<T> & f, uint_t thickness,
                                                       stencil::Direction dir, bool fullSlice )
{
   CellInterval ci;
   f.getGhostRegion( dir, ci, cell_idx_c(thickness), fullSlice );
   return interval( f, ci );
}

template< typename T>
FieldIndexing3D<T> FieldIndexing3D<T>::sliceBeforeGhostLayerXYZ( const GPUField<T> & f, uint_t thickness,
                                                              stencil::Direction dir, bool fullSlice )
{
   CellInterval ci;
   f.getSliceBeforeGhostLayer( dir, ci, cell_idx_c(thickness), fullSlice );
   return interval( f, ci );
}

template< typename T>
FieldIndexing3D<T> FieldIndexing3D<T>::sliceXYZ( const GPUField<T> & f, cell_idx_t distance, uint_t thickness,
                                              stencil::Direction dir, bool fullSlice )
{
   CellInterval ci;
   f.getSlice( dir, ci, distance, cell_idx_c(thickness), fullSlice );
   return interval( f, ci );
}

template< typename T>
FieldIndexing3D<T> FieldIndexing3D<T>::intervalXYZ( const GPUField<T> & f, const cell::CellInterval & ci )
{
   return interval( f, ci );
}



} // namespace cuda
} // namespace walberla


