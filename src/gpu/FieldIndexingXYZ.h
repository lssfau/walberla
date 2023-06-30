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
//! \file FieldIndexingXYZ.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "DeviceWrapper.h"
#include "FieldAccessorXYZ.h"

namespace walberla { namespace cell {  class CellInterval;  } }


namespace walberla {
namespace gpu
{

// Forward Declarations
template< typename T> class GPUField;


   template<typename T>
   class FieldIndexingXYZ
   {
   public:

      //** Kernel call        ******************************************************************************************
      /*! \name Kernel call  */
      //@{
      dim3 blockDim() const                        { return blockDim_; }
      dim3 gridDim () const                        { return gridDim_;  }

      const FieldAccessorXYZ<T> & gpuAccess() const { return gpuAccess_; }
      //@}
      //****************************************************************************************************************


      //** Creation        *********************************************************************************************
      /*! \name Creation  */
      //@{

      static FieldIndexingXYZ<T> interval ( const GPUField<T> & f, const cell::CellInterval & ci  );


      static FieldIndexingXYZ<T> xyz ( const GPUField<T> & f );
      static FieldIndexingXYZ<T> withGhostLayerXYZ       ( const GPUField<T> & f, uint_t numGhostLayers );
      //@}
      //****************************************************************************************************************

   protected:
      FieldIndexingXYZ<T> ( const GPUField<T> & field, dim3 _blockDim, dim3 _gridDim, const FieldAccessorXYZ<T> _gpuAccess );

      const GPUField<T> &  field_;
      dim3 blockDim_;
      dim3 gridDim_;
      FieldAccessorXYZ<T> gpuAccess_;
   };


} // namespace gpu
} // namespace walberla


#include "FieldIndexingXYZ.impl.h"