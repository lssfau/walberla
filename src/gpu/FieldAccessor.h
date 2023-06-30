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

      __device__ void set( uint3 blockIdx, uint3 threadIdx )
      {
         switch ( indexingScheme_)
         {
            case FZYX: ptr_ += blockIdx.z * fOffset_ + blockIdx.y * zOffset_ + blockIdx.x * yOffset_ + threadIdx.x * xOffset_; break;
            case FZY : ptr_ +=                         blockIdx.y * fOffset_ + blockIdx.x * zOffset_ + threadIdx.x * yOffset_; break;
            case FZ  : ptr_ +=                                                 blockIdx.x * fOffset_ + threadIdx.x * zOffset_; break;
            case F   : ptr_ +=                                                                         threadIdx.x * fOffset_; break;

            case ZYXF: ptr_ += blockIdx.z * zOffset_ + blockIdx.y * yOffset_ + blockIdx.x * xOffset_ + threadIdx.x * fOffset_; break;
            case ZYX : ptr_ +=                         blockIdx.y * zOffset_ + blockIdx.x * yOffset_ + threadIdx.x * xOffset_; break;
            case ZY  : ptr_ +=                                                 blockIdx.x * zOffset_ + threadIdx.x * yOffset_; break;
            case Z   : ptr_ +=                                                                         threadIdx.x * zOffset_; break;
         }
      }


      __device__ uint_t getLinearIndex( uint3 blockIdx, uint3 threadIdx, uint3 gridDim, uint3 blockDim )
      {
         return threadIdx.x                              +
                blockIdx.x * blockDim.x                  +
                blockIdx.y * blockDim.x * gridDim.x      +
                blockIdx.z * blockDim.x * gridDim.x * gridDim.y ;
      }

      // This is always true for this specific field indexing class.
      __device__ __forceinline__ bool isValidPosition()  { return true; }

      __device__ T & get()       { return * (T*)(ptr_);                }
      __device__ T & get( int f) { return * (T*)(ptr_ + f * fOffset_); }


      __device__ T & getNeighbor( int cx, int cy, int cz ) const
      {
         return * (T*)( ptr_ + cx * xOffset_ +
                               cy * yOffset_ +
                               cz * zOffset_ );
      }

      __device__ T & getNeighbor( int cx, int cy, int cz, int cf )
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


