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
//! \file FieldAccessor.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "gpu/GPUWrapper.h"
#include "gpu/DeviceWrapper.h"

namespace walberla {
namespace gpu
{



   /**
    * \brief Handle to the underlying device data of a \ref GPUField.
    *
    * Encapsulate the device memory pointer and offsets necessary
    * to calculate the address of a cell from a GPU kernel's thread
    * coordinates in the thread block.
    */
   template<typename T>
   class FieldAccessor
   {
   public:
      enum IndexingScheme { FZYX, FZY, FZ, F,
                            ZYXF, ZYX, ZY, Z
                            };

      FieldAccessor( char * ptr,
                     uint_t xOffset,
                     uint_t yOffset,
                     uint_t zOffset,
                     uint_t fOffset,
                     IndexingScheme indexingScheme )
            : ptr_(ptr), xOffset_(xOffset), yOffset_(yOffset), zOffset_(zOffset),
              fOffset_(fOffset), indexingScheme_(indexingScheme )
      {}

      __device__ void set( uint3 _blockIdx, uint3 _threadIdx )
      {
         switch ( indexingScheme_)
         {
            case FZYX: ptr_ += _blockIdx.z * fOffset_ + _blockIdx.y * zOffset_ + _blockIdx.x * yOffset_ + _threadIdx.x * xOffset_; break;
            case FZY : ptr_ +=                          _blockIdx.y * fOffset_ + _blockIdx.x * zOffset_ + _threadIdx.x * yOffset_; break;
            case FZ  : ptr_ +=                                                   _blockIdx.x * fOffset_ + _threadIdx.x * zOffset_; break;
            case F   : ptr_ +=                                                                            _threadIdx.x * fOffset_; break;

            case ZYXF: ptr_ += _blockIdx.z * zOffset_ + _blockIdx.y * yOffset_ + _blockIdx.x * xOffset_ + _threadIdx.x * fOffset_; break;
            case ZYX : ptr_ +=                          _blockIdx.y * zOffset_ + _blockIdx.x * yOffset_ + _threadIdx.x * xOffset_; break;
            case ZY  : ptr_ +=                                                   _blockIdx.x * zOffset_ + _threadIdx.x * yOffset_; break;
            case Z   : ptr_ +=                                                                            _threadIdx.x * zOffset_; break;
         }
      }


      __device__ uint_t getLinearIndex( uint3 _blockIdx, uint3 _threadIdx, uint3 _gridDim, uint3 _blockDim )
      {
         return _threadIdx.x                                +
                _blockIdx.x * _blockDim.x                   +
                _blockIdx.y * _blockDim.x * _gridDim.x      +
                _blockIdx.z * _blockDim.x * _gridDim.x * _gridDim.y ;
      }

      // This is always true for this specific field indexing class.
      __device__ __forceinline__ bool isValidPosition()  { return true; }

      __device__ T & get()       { return * (T*)(ptr_);                }
      __device__ T & get( uint_t f) { return * (T*)(ptr_ + f * fOffset_); }


      __device__ T & getNeighbor( int cx, int cy, int cz ) const
      {
         return * (T*)( ptr_ + cx * xOffset_ +
                               cy * yOffset_ +
                               cz * zOffset_ );
      }

      __device__ T & getNeighbor( int cx, int cy, int cz, uint_t cf )
      {
         return * (T*)( ptr_ + cx * xOffset_ +
                               cy * yOffset_ +
                               cz * zOffset_ +
                               cf * fOffset_ );
      }


   protected:
      char * ptr_;
      uint_t xOffset_;
      uint_t yOffset_;
      uint_t zOffset_;
      uint_t fOffset_;
      IndexingScheme indexingScheme_;
   };




} // namespace gpu
} // namespace walberla


