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
//! \file FieldIndexing3D.h
//! \ingroup cuda
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \brief Indexing Scheme that executes all elements of inner coordinate within on thread block
//
//======================================================================================================================

#pragma once

#include "FieldAccessor3D.h"

#include "stencil/Directions.h"
#include <cuda_runtime.h>

namespace walberla { namespace cell {  class CellInterval;  } }

namespace walberla {
namespace cuda {

   // Forward Declarations
   template< typename T> class GPUField;


   class FieldIndexing3DBase
   {
   public:

      static void setPreferredBlockDim( dim3 blockDim ) { preferredBlockDim_ = blockDim; }
      static void setPreferredBlockDim( unsigned int x, unsigned int y, unsigned int z ) { preferredBlockDim_ = dim3( x, y, z ); }

   protected:

      static dim3 preferredBlockDim_;
   };


   template<typename T>
   class FieldIndexing3D: public FieldIndexing3DBase
   {
   public:

      //** Kernel call        ******************************************************************************************
      /*! \name Kernel call  */
      //@{
      dim3 blockDim() const                        { return blockDim_; }
      dim3 gridDim () const                        { return gridDim_;  }

      const FieldAccessor3D<T> & gpuAccess() const { return gpuAccess_; }
      //@}
      //****************************************************************************************************************




      //** Creation        *********************************************************************************************
      /*! \name Creation  */
      //@{
      static FieldIndexing3D<T> interval ( const GPUField<T> & f,
                                           const cell::CellInterval & ci );


      static FieldIndexing3D<T> xyz                      ( const GPUField<T> & f );
      static FieldIndexing3D<T> withGhostLayerXYZ        ( const GPUField<T> & f, uint_t numGhostLayers );
      static FieldIndexing3D<T> ghostLayerOnlyXYZ        ( const GPUField<T> & f, uint_t thickness,
                                                           stencil::Direction dir, bool fullSlice = false );
      static FieldIndexing3D<T> sliceBeforeGhostLayerXYZ ( const GPUField<T> & f, uint_t thickness,
                                                           stencil::Direction dir, bool fullSlice = false );
      static FieldIndexing3D<T> sliceXYZ                 ( const GPUField<T> & f, cell_idx_t distance, uint_t thickness,
                                                        stencil::Direction dir, bool fullSlice = false );

      static FieldIndexing3D<T> intervalXYZ              ( const GPUField<T> & f, const cell::CellInterval & ci );
      //@}
      //****************************************************************************************************************

   protected:
      FieldIndexing3D ( const GPUField<T> & field,
                        const dim3 & _blockDim, const dim3 & _gridDim,
                        const FieldAccessor3D<T> & _gpuAccess );

      const GPUField<T> &  field_;
      dim3 blockDim_;
      dim3 gridDim_;
      FieldAccessor3D<T> gpuAccess_;
   };


} // namespace cuda
} // namespace walberla


#include "FieldIndexing3D.impl.h"