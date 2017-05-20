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
//! \file GPUField.h
//! \ingroup moduleName
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "field/Layout.h"
#include "stencil/Directions.h"

#include <cuda_runtime.h>



namespace walberla {
namespace cuda {

   using field::Layout;
   using field::fzyx;
   using field::zyxf;


   //*******************************************************************************************************************
   /*! GhostLayerField stored on a CUDA GPU
   *
   *  Basically a wrapper around a CUDA device pointer together with size information about the field
   *  i.e. sizes in x,y,z,f directions and number of ghost layers.
   *
   *  Internally represented by a cudaPitchedPtr which is allocated with cudaMalloc3D to take padding of the
   *  innermost coordinate into account.
   *
   *  Supports Array-of-Structures (AoS,zyxf) layout and Structure-of-Arrays (SoA, fzyx) layout, in a similiar way
   *  to field::Field
   *
   *  To work with the GPUField look at the cuda::fieldCpy functions to transfer a field::Field to a cuda::GPUField
   *  and vice versa.
   *  When writing CUDA kernels for GPUFields have a look at the FieldIndexing and FieldAccessor concepts.
   *  These simplify the "iteration" i.e. indexing of cells in GPUFields.
   */
   //*******************************************************************************************************************
   template<typename T>
   class GPUField
   {
   public:
      typedef T value_type;

      GPUField( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t _fSize,
                uint_t _nrOfGhostLayers, const Layout & _layout = zyxf );

      ~GPUField();

      Layout layout() const { return layout_; }

      cudaPitchedPtr pitchedPtr() const { return pitchedPtr_; }


      inline uint_t  xSize() const  { return xSize_; }
      inline uint_t  ySize() const  { return ySize_; }
      inline uint_t  zSize() const  { return zSize_; }
      inline uint_t  fSize() const  { return fSize_; }

      cell_idx_t xOff() const { return cell_idx_c( nrOfGhostLayers_ ); }
      cell_idx_t yOff() const { return cell_idx_c( nrOfGhostLayers_ ); }
      cell_idx_t zOff() const { return cell_idx_c( nrOfGhostLayers_ ); }


      uint_t xAllocSize() const;
      uint_t yAllocSize() const;
      uint_t zAllocSize() const;
      uint_t fAllocSize() const;


      void swapDataPointers( GPUField<T> & other );
      void swapDataPointers( GPUField<T> * other ) { swapDataPointers( *other); }


      inline uint_t  nrOfGhostLayers() const { return nrOfGhostLayers_; }

      bool operator==( const GPUField & other ) const;

      void getGhostRegion( stencil::Direction d, CellInterval & ci,
                           cell_idx_t thickness, bool fullSlice ) const;
      void getSliceBeforeGhostLayer(stencil::Direction d, CellInterval & ci,
                                    cell_idx_t thickness, bool fullSlice ) const;


            void * data()            { return pitchedPtr_.ptr; }
      const void * data() const      { return pitchedPtr_.ptr; }

   protected:
      cudaPitchedPtr pitchedPtr_;
      uint_t         nrOfGhostLayers_;
      uint_t         xSize_;
      uint_t         ySize_;
      uint_t         zSize_;
      uint_t         fSize_;
      Layout         layout_;
   };


} // namespace cuda
} // namespace walberla


