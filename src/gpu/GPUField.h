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
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"

#include "field/Layout.h"

#include "stencil/Directions.h"

#include "gpu/DeviceWrapper.h"

namespace walberla {
namespace gpu
{

   using field::Layout;
   using field::fzyx;
   using field::zyxf;


   //*******************************************************************************************************************
   /*! GhostLayerField stored on a CUDA/HIP GPU
   *
   *  Basically a wrapper around a CUDA/HIP device pointer together with size information about the field
   *  i.e. sizes in x,y,z,f directions and number of ghost layers.
   *
   *  Internally represented by a \c gpuPitchedPtr which is allocated with extra padding for the
   *  innermost coordinate.
   *  Pitched memory is a type of non-linear memory where padding is introduced
   *  to optimize data alignment and thus reduce data access latency,
   *  for example by avoiding shared memory bank conflicts.
   *
   *  Supports Array-of-Structures (AoS,zyxf) layout and Structure-of-Arrays (SoA, fzyx) layout, in a similar way
   *  to \ref field::Field
   *
   *  To work with the \ref gpu::GPUField look at the \ref gpu::fieldCpy functions to transfer a \ref field::Field to a \ref gpu::GPUField
   *  and vice versa.
   *
   *  When writing device kernels for a \ref GPUField, have a look at the \ref FieldIndexing and \ref FieldAccessor concepts.
   *  These simplify the "iteration" i.e. indexing of cells in a \ref GPUField.
   */
   //*******************************************************************************************************************
   template<typename T>
   class GPUField
   {
   public:
      typedef T value_type;

      GPUField( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t _fSize,
                uint_t _nrOfGhostLayers, const Layout & _layout = fzyx, bool usePitchedMem = true );

      ~GPUField();

      Layout layout() const { return layout_; }

      bool isPitchedMem() const { return usePitchedMem_; }

      gpuPitchedPtr pitchedPtr() const { return pitchedPtr_; }


      inline uint_t  xSize() const  { return xSize_; }
      inline uint_t  ySize() const  { return ySize_; }
      inline uint_t  zSize() const  { return zSize_; }
      inline uint_t  fSize() const  { return fSize_; }
      inline uint_t  size()  const  { return fSize() * xSize() * ySize() * zSize(); }
      inline uint_t  size( uint_t coord )  const;

      inline uint_t       xSizeWithGhostLayer()        const  { return xSize() + uint_c(2)*nrOfGhostLayers_; }
      inline uint_t       ySizeWithGhostLayer()        const  { return ySize() + uint_c(2)*nrOfGhostLayers_; }
      inline uint_t       zSizeWithGhostLayer()        const  { return zSize() + uint_c(2)*nrOfGhostLayers_; }
      inline uint_t       sizeWithGhostLayer(uint_t i) const  { return i==3 ? fSize_ :
                                                                              size(i) + uint_c(2)*nrOfGhostLayers_; }

      cell_idx_t xOff() const { return cell_idx_c( nrOfGhostLayers_ ); }
      cell_idx_t yOff() const { return cell_idx_c( nrOfGhostLayers_ ); }
      cell_idx_t zOff() const { return cell_idx_c( nrOfGhostLayers_ ); }

      cell_idx_t xStride() const { return (layout_ == fzyx) ? cell_idx_t(1) :
                                                              cell_idx_c(fAllocSize()); }
      cell_idx_t yStride() const { return (layout_ == fzyx) ? cell_idx_t(xAllocSize()) :
                                                              cell_idx_c(fAllocSize() * xAllocSize()); }
      cell_idx_t zStride() const { return (layout_ == fzyx) ? cell_idx_t(xAllocSize() * yAllocSize()) :
                                                              cell_idx_c(fAllocSize() * xAllocSize() * yAllocSize()); }
      cell_idx_t fStride() const { return (layout_ == fzyx) ? cell_idx_t(xAllocSize() * yAllocSize() * zAllocSize()) :
                                                              cell_idx_c(1); }


      uint_t xAllocSize() const;
      uint_t yAllocSize() const;
      uint_t zAllocSize() const;
      uint_t fAllocSize() const;
      inline uint_t allocSize() const { return fAllocSize() * xAllocSize() * yAllocSize() * zAllocSize(); }

      bool hasSameAllocSize( const GPUField<T> & other ) const;
      bool hasSameSize( const GPUField<T> & other ) const;

      GPUField<T> * cloneUninitialized() const;

      void swapDataPointers( GPUField<T> & other );
      void swapDataPointers( GPUField<T> * other ) { swapDataPointers( *other ); }


      inline uint_t  nrOfGhostLayers() const { return nrOfGhostLayers_; }

      inline CellInterval xyzSize()               const;
      inline CellInterval xyzSizeWithGhostLayer() const;

      bool operator==( const GPUField & other ) const;

      void getGhostRegion( stencil::Direction d, CellInterval & ci,
                           cell_idx_t thickness, bool fullSlice = false ) const;
      void getSliceBeforeGhostLayer(stencil::Direction d, CellInterval & ci,
                                    cell_idx_t thickness, bool fullSlice = false ) const
      {
         getSlice( d, ci, 0, thickness, fullSlice );
      }
      void getSlice(stencil::Direction d, CellInterval & ci,
                    cell_idx_t distance, cell_idx_t thickness, bool fullSlice ) const;

            void * data()            { return pitchedPtr_.ptr; }
      const void * data() const      { return pitchedPtr_.ptr; }

      T       * dataAt(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f);
      const T * dataAt(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f) const;

      //** TimestepInformation *****************************************************************************************
      /*! \name TimestepCounter */
      //@{
      inline uint8_t advanceTimestep()
      {
         timestepCounter_ = (timestepCounter_ + 1) & 1;
         return timestepCounter_;
      }
      inline uint8_t getTimestep() const { return timestepCounter_; }
      inline uint8_t getTimestepPlusOne() const { return (timestepCounter_ + 1) & 1; }
      inline bool isEvenTimeStep() const {return (((timestepCounter_) &1) ^ 1); }
      //@}
      //****************************************************************************************************************

   protected:
      gpuPitchedPtr  pitchedPtr_;
      uint_t         nrOfGhostLayers_;
      uint_t         xSize_;
      uint_t         ySize_;
      uint_t         zSize_;
      uint_t         fSize_;

      uint_t         xAllocSize_;
      uint_t         fAllocSize_;
      Layout         layout_;
      bool           usePitchedMem_;
      uint8_t        timestepCounter_;
   };


} // namespace gpu
} // namespace walberla


#include "GPUField.impl.h"