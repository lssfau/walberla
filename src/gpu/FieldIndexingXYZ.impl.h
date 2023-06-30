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
//! \file FieldIndexingXYZ.impl.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "FieldIndexingXYZ.h"
#include "GPUField.h"

#include "core/cell/CellInterval.h"


namespace walberla {
namespace gpu
{


template< typename T>
FieldIndexingXYZ<T>::FieldIndexingXYZ ( const GPUField<T> & field,
                                        dim3 _blockDim, dim3 _gridDim,
                                        const FieldAccessorXYZ<T> _gpuAccess )
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
         gpuGetDeviceCount(&count);
         int threadsPerBlock = std::numeric_limits< int >::max();
         for (int i = 0; i < count; i++)
         {
            gpuGetDeviceProperties(&prop, i);
            threadsPerBlock = std::min(prop.maxThreadsPerBlock, threadsPerBlock);
         }
         WALBERLA_ASSERT_LESS(int_c(blockDim_.x), threadsPerBlock,
                              "InnerCoordThreadIndexing works only for fields where each dimension x,y,z is smaller "
                                 << "than the maximal thread count per GPU block.")
      }
   }
}

template< typename T>
FieldIndexingXYZ<T> FieldIndexingXYZ<T>::interval ( const GPUField<T> & f, const CellInterval & ci )
{
   size_t xOffset;
   size_t yOffset;
   size_t zOffset;
   size_t fOffset;

   if ( f.layout() == field::zyxf )
   {
      fOffset = sizeof(T);
      xOffset = f.pitchedPtr().pitch;
      yOffset = xOffset * f.xSize();
      zOffset = yOffset * f.ySize();
   }
   else
   {
      xOffset = sizeof(T);
      yOffset = f.pitchedPtr().pitch;
      zOffset = yOffset * f.ySize();
      fOffset = zOffset * f.zSize();
   }
   char * data = (char*)f.pitchedPtr().ptr;

   // Jump over ghost cells to first inner cell
   cell_idx_t gl = cell_idx_c( f.nrOfGhostLayers() );
   data += ( ci.xMin() + gl )* int_c(xOffset) +
           ( ci.yMin() + gl )* int_c(yOffset) +
           ( ci.zMin() + gl )* int_c(zOffset);

   dim3 gridDim ( (unsigned int)ci.xSize(), 1, 1 );
   dim3 blockDim( (unsigned int)ci.ySize(), (unsigned int)ci.zSize(), 1 );
   return FieldIndexingXYZ<T> ( f, gridDim, blockDim,
                                FieldAccessorXYZ<T>( data, xOffset, yOffset, zOffset, fOffset ) );
}


template< typename T>
FieldIndexingXYZ<T> FieldIndexingXYZ<T>::xyz ( const GPUField<T> & f )
{
   CellInterval ci ( 0,0,0,
                    cell_idx_c( f.xSize() ) - 1,
                    cell_idx_c( f.ySize() ) - 1,
                    cell_idx_c( f.zSize() ) - 1 );

   return interval( f, ci );
}


template< typename T>
FieldIndexingXYZ<T> FieldIndexingXYZ<T>::withGhostLayerXYZ( const GPUField<T> & f, uint_t numGhostLayers )
{
   cell_idx_t gl = std::min( cell_idx_c(numGhostLayers), cell_idx_c( f.nrOfGhostLayers() ) );
   CellInterval ci ( -gl,-gl,-gl,
                    cell_idx_c( f.xSize() ) + gl - 1,
                    cell_idx_c( f.ySize() ) + gl - 1,
                    cell_idx_c( f.zSize() ) + gl - 1 );

   return interval( f, ci );
}



} // namespace gpu
} // namespace walberla


