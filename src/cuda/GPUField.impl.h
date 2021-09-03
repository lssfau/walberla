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
//! \file GPUField.impl.h
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "GPUField.h"
#include "ErrorChecking.h"
#include "AlignedAllocation.h"
#include "core/logging/Logging.h"

namespace walberla {
namespace cuda {


template<typename T>
GPUField<T>::GPUField( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t _fSize,
                       uint_t _nrOfGhostLayers, const Layout & _layout, bool usePitchedMem )
   : nrOfGhostLayers_( _nrOfGhostLayers ),
     xSize_( _xSize), ySize_( _ySize ), zSize_( _zSize ), fSize_( _fSize ),
     layout_( _layout ), usePitchedMem_( usePitchedMem )
{
   cudaExtent extent;
   if ( layout_ == zyxf )
   {
      extent.width  = _fSize * sizeof(T);
      extent.height = (_xSize + 2 * _nrOfGhostLayers );
      extent.depth  = (_ySize + 2 * _nrOfGhostLayers ) * ( _zSize + 2 * _nrOfGhostLayers );
   }
   else
   {
      extent.width  = (_xSize + 2 * _nrOfGhostLayers ) * sizeof(T);
      extent.height = (_ySize + 2 * _nrOfGhostLayers );
      extent.depth  = (_zSize + 2 * _nrOfGhostLayers ) * _fSize;
   }

   if ( usePitchedMem_ )
   {
      size_t pitch;
      const size_t alignment = 256;
      void * mem = allocate_pitched_with_offset( pitch, extent.width, extent.height * extent.depth, alignment,
                                                 sizeof(T) * nrOfGhostLayers_ );
      WALBERLA_ASSERT_EQUAL( size_t((char*)(mem) + sizeof(T) * nrOfGhostLayers_ ) % alignment, 0 );
      pitchedPtr_ = make_cudaPitchedPtr( mem, pitch, extent.width, extent.height );
   }
   else
   {
      pitchedPtr_ = make_cudaPitchedPtr( NULL, extent.width, extent.width, extent.height );
      WALBERLA_CUDA_CHECK ( cudaMalloc( &pitchedPtr_.ptr, extent.width * extent.height * extent.depth ) );
   }

   // allocation size is stored in pitched pointer
   // pitched pointer stores the amount of padded region in bytes
   // but we keep track of the size in #elements
   WALBERLA_ASSERT_EQUAL( pitchedPtr_.pitch % sizeof(T), 0 );
   if ( layout_ == field::fzyx )
   {
      xAllocSize_ = pitchedPtr_.pitch / sizeof(T);
      fAllocSize_ = fSize_;
   }
   else
   {
      fAllocSize_ = pitchedPtr_.pitch / sizeof(T);
      xAllocSize_ = xSize_ + 2 * nrOfGhostLayers_;
   }
}


template<typename T>
GPUField<T>::~GPUField()
{
   if( usePitchedMem_ )
      free_aligned_with_offset(pitchedPtr_.ptr );
   else
   {
      WALBERLA_CUDA_CHECK( cudaFree( pitchedPtr_.ptr ) );
   }
}


template<typename T>
T * GPUField<T>::dataAt(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f)
{
   auto offset = (x + cell_idx_c(nrOfGhostLayers_)) * xStride() +
                 (y + cell_idx_c(nrOfGhostLayers_)) * yStride() +
                 (z + cell_idx_c(nrOfGhostLayers_)) * zStride() +
                 f * fStride();
   return static_cast<T*>(pitchedPtr_.ptr) + offset;
}

template<typename T>
const T * GPUField<T>::dataAt(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f) const
{
   auto offset = (x + cell_idx_c(nrOfGhostLayers_)) * xStride() +
                 (y + cell_idx_c(nrOfGhostLayers_)) * yStride() +
                 (z + cell_idx_c(nrOfGhostLayers_)) * zStride() +
                 f * fStride();
   return static_cast<T*>(pitchedPtr_.ptr) + offset;
}


template<typename T>
void GPUField<T>::getGhostRegion(stencil::Direction d, CellInterval & ci,
                                 cell_idx_t thickness, bool fullSlice ) const
{
   const cell_idx_t sizeArr [] = { cell_idx_c( xSize() ),
                                   cell_idx_c( ySize() ),
                                   cell_idx_c( zSize() )};

   WALBERLA_ASSERT_GREATER( thickness, 0 );
   WALBERLA_ASSERT_LESS_EQUAL( uint_c(thickness), nrOfGhostLayers() );
   const cell_idx_t ghosts = cell_idx_c ( thickness );

   cell_idx_t fullSliceInc = fullSlice ? cell_idx_c( nrOfGhostLayers() ) : 0;

   for(unsigned int dim = 0; dim< 3; ++dim)
      switch ( stencil::c[dim][d] )
      {
         case -1: ci.min()[dim] =     -ghosts;     ci.max()[dim] =         0                   - 1; break;
         case  0: ci.min()[dim] = -fullSliceInc;   ci.max()[dim] =  sizeArr[dim]+fullSliceInc  - 1; break;
         case  1: ci.min()[dim] =   sizeArr[dim];  ci.max()[dim] =  sizeArr[dim]+ghosts        - 1; break;
      }
}


template<typename T>
inline CellInterval GPUField<T>::xyzSize() const
{
   return CellInterval (0,0,0,
                        cell_idx_c( xSize() )-1,
                        cell_idx_c( ySize() )-1,
                        cell_idx_c( zSize() )-1 );
}

template<typename T>
inline CellInterval GPUField<T>::xyzSizeWithGhostLayer() const
{
   CellInterval ci = GPUField<T>::xyzSize();
   for( uint_t i=0; i < 3; ++i ) {
      ci.min()[i] -= cell_idx_c( nrOfGhostLayers_ );
      ci.max()[i] += cell_idx_c( nrOfGhostLayers_ );
   }
   return ci;
}

template<typename T>
void GPUField<T>::getSlice(stencil::Direction d, CellInterval & ci,
                           cell_idx_t distance, cell_idx_t thickness, bool fullSlice ) const
{
   WALBERLA_ASSERT_GREATER( thickness, 0 );

   const cell_idx_t sizeArr [] = { cell_idx_c( xSize() ),
                                   cell_idx_c( ySize() ),
                                   cell_idx_c( zSize() ) };

   cell_idx_t fullSliceInc = fullSlice ? cell_idx_c(  nrOfGhostLayers() ) : 0;
   for(unsigned int dim = 0; dim< 3; ++dim)
   {
      switch ( stencil::c[dim][d] )
      {
         case -1:
            ci.min()[dim] = distance;
            ci.max()[dim] = distance + thickness - 1;
            break;
         case  0:
            ci.min()[dim] = -fullSliceInc;
            ci.max()[dim] =  sizeArr[dim] + fullSliceInc - 1;
            break;
         case  1:
            ci.min()[dim] = sizeArr[dim] - distance - thickness;
            ci.max()[dim] = sizeArr[dim] - distance - 1;
            break;
      }
   }
}

template<typename T>
inline uint_t GPUField<T>::size( uint_t coord ) const
{
   switch (coord) {
      case 0: return this->xSize();
      case 1: return this->ySize();
      case 2: return this->zSize();
      case 3: return this->fSize();
      default: WALBERLA_ASSERT(false); return 0;
   }
}

//*******************************************************************************************************************
/*! True if sizes of all dimensions match
 *******************************************************************************************************************/
template<typename T>
bool GPUField<T>::hasSameSize( const GPUField<T> & other ) const
{
   return xSize() == other.xSize() &&
          ySize() == other.ySize() &&
          zSize() == other.zSize();
}

//*******************************************************************************************************************
/*! True if allocation sizes of all dimensions match
 *******************************************************************************************************************/
template<typename T>
bool GPUField<T>::hasSameAllocSize( const GPUField<T> & other ) const
{
   return xAllocSize() == other.xAllocSize() &&
          yAllocSize() == other.yAllocSize() &&
          zAllocSize() == other.zAllocSize() &&
          fAllocSize() == other.fAllocSize();
}

//*******************************************************************************************************************
/*! Creates a new GPUField that has equal size, layout and memory type as this field but has uninitialized memory.
 *
 * \return a new FPUField that has to be freed by caller.
 *******************************************************************************************************************/
template <typename T>
GPUField<T> * GPUField<T>::cloneUninitialized() const
{
   GPUField<T> * res = new GPUField<T>( xSize(), ySize(), zSize(), fSize(),
                                        nrOfGhostLayers(), layout(), isPitchedMem() );

   WALBERLA_ASSERT( hasSameAllocSize( *res ) );
   WALBERLA_ASSERT( hasSameSize( *res ) );
   WALBERLA_ASSERT( layout() == res->layout() );
   WALBERLA_ASSERT( isPitchedMem() == res->isPitchedMem() );
   return res;
}

//*******************************************************************************************************************
/*! Implementation of equality operator, GPUField are equal if identical
*
* operator== is required in order to store GPUFields as block data.
* Implementing a correct equality check would require to call a kernel,
* which should be done manually.
* Therefore only identical GPUFields are considered equal.
*/
//*******************************************************************************************************************
template<typename T>
bool GPUField<T>::operator==( const GPUField & o ) const
{
   return pitchedPtr_.ptr  == o.pitchedPtr_.ptr  &&
          nrOfGhostLayers_ == o.nrOfGhostLayers_ &&
          xSize_           == o.xSize_           &&
          ySize_           == o.ySize_           &&
          zSize_           == o.zSize_           &&
          fSize_           == o.fSize_           &&
          layout_          == o.layout_          ;
}



template<typename T>
uint_t  GPUField<T>::xAllocSize() const
{
   return xAllocSize_;
}

template<typename T>
uint_t  GPUField<T>::yAllocSize() const
{
   return ySize_ + 2 * nrOfGhostLayers_;
}

template<typename T>
uint_t  GPUField<T>::zAllocSize() const
{
   return zSize_ + 2 * nrOfGhostLayers_;
}

template<typename T>
uint_t GPUField<T>::fAllocSize() const
{
   return fAllocSize_;
}

template<typename T>
void GPUField<T>::swapDataPointers( GPUField<T> & other )
{
   WALBERLA_ASSERT( hasSameAllocSize( other ) );
   WALBERLA_ASSERT( hasSameSize( other ) );
   WALBERLA_ASSERT( layout() == other.layout() );
   WALBERLA_ASSERT( isPitchedMem() == other.isPitchedMem() );
   std::swap( pitchedPtr_, other.pitchedPtr_ );
}



} // namespace cuda
} // namespace walberla


