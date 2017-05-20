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
//! \file GPUField.cpp
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "GPUField.h"
#include "ErrorChecking.h"
#include "GPUTypesExplicitInstantiation.h"

#include "core/logging/Logging.h"

namespace walberla {
namespace cuda {


template<typename T>
GPUField<T>::GPUField( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t _fSize,
                       uint_t _nrOfGhostLayers, const Layout & _layout )
   : nrOfGhostLayers_( _nrOfGhostLayers ),
     xSize_( _xSize), ySize_( _ySize ), zSize_( _zSize ), fSize_( _fSize ),
     layout_( _layout )
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

   WALBERLA_CUDA_CHECK ( cudaMalloc3D ( &pitchedPtr_, extent ) );
}


template<typename T>
GPUField<T>::~GPUField()
{
   cudaFree( pitchedPtr_.ptr );
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
void GPUField<T>::getSliceBeforeGhostLayer(stencil::Direction d, CellInterval & ci,
                                            cell_idx_t thickness, bool fullSlice ) const
{
   WALBERLA_ASSERT_GREATER( thickness, 0 );

   const cell_idx_t sizeArr [] = { cell_idx_c( xSize() ),
                                   cell_idx_c( ySize() ),
                                   cell_idx_c( zSize() )};

   cell_idx_t fullSliceInc = fullSlice ? cell_idx_c(  nrOfGhostLayers() ) : 0;
   for(unsigned int dim = 0; dim< 3; ++dim)
      switch ( stencil::c[dim][d] )
      {
         case -1: ci.min()[dim] =                      0;  ci.max()[dim] =     thickness              - 1; break;
         case  0: ci.min()[dim] =          -fullSliceInc;  ci.max()[dim] =  sizeArr[dim] +fullSliceInc- 1; break;
         case  1: ci.min()[dim] = sizeArr[dim]-thickness;  ci.max()[dim] =  sizeArr[dim]              - 1; break;
      }
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
   if ( layout_ == field::fzyx )
   {
      // allocation size is stored in pitched pointer
      // pitched pointer stores the amount of padded region in bytes
      // but we have to return the size in #elements
      WALBERLA_ASSERT_EQUAL( pitchedPtr_.pitch % sizeof(T), 0 );
      return pitchedPtr_.pitch / sizeof(T);
   }
   return xSize_ + 2 * nrOfGhostLayers_;
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
   if ( layout_ == field::zyxf )
   {
      WALBERLA_ASSERT_EQUAL( pitchedPtr_.pitch % sizeof(T), 0 );
      return pitchedPtr_.pitch / sizeof(T);
   }
   return fSize_;
}

template<typename T>
void GPUField<T>::swapDataPointers( GPUField<T> & other )
{
   WALBERLA_ASSERT_EQUAL( xAllocSize(), other.xAllocSize() );
   WALBERLA_ASSERT_EQUAL( yAllocSize(), other.yAllocSize() );
   WALBERLA_ASSERT_EQUAL( zAllocSize(), other.zAllocSize() );
   WALBERLA_ASSERT_EQUAL( fAllocSize(), other.fAllocSize() );
   WALBERLA_ASSERT_EQUAL( layout(),     other.layout()     );
   std::swap( pitchedPtr_, other.pitchedPtr_ );
}



GPU_CLASS_TEMPLATE_INSTANTIATION( GPUField )


} // namespace cuda
} // namespace walberla


