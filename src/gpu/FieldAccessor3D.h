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
//! \file FieldAccessor3D.h
//! \ingroup gpu
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

namespace walberla {
namespace gpu
{


   template<typename T>
   class FieldAccessor3D
   {
   public:

      FieldAccessor3D( char * ptr,
                       uint_t xOffset,
                       uint_t yOffset,
                       uint_t zOffset,
                       uint_t fOffset,
                       const dim3 & _idxDim,
                       const dim3 & _blockDim )
            : ptr_( ptr ), xOffset_( xOffset ), yOffset_( yOffset ), zOffset_( zOffset ), fOffset_( fOffset ),
              idxDim_( _idxDim ), blockDim_( _blockDim )
      {}

      __device__ __forceinline__ void set( const uint3& _blockIdx, const uint3& _threadIdx )
      {
         uint_t x = _blockIdx.x * blockDim_.x + _threadIdx.x;
         uint_t y = _blockIdx.y * blockDim_.y + _threadIdx.y;
         uint_t z = _blockIdx.z * blockDim_.z + _threadIdx.z;

         if ( x < idxDim_.x && y < idxDim_.y && z < idxDim_.z )
         {
            ptr_ += x * xOffset_ + y * yOffset_ + z * zOffset_;
            isValidPosition_ = true;
         }
      }

      __device__ __forceinline__ bool isValidPosition()  { return isValidPosition_; }

      __device__ __forceinline__ T & get()               { return * (T*)(ptr_);                }
      __device__ __forceinline__ T & get( int f )        { return * (T*)(ptr_ + f * fOffset_); }


      __device__ __forceinline__ T & getNeighbor( int cx, int cy, int cz ) const
      {
         return * (T*)( ptr_ + cx * xOffset_ +
                               cy * yOffset_ +
                               cz * zOffset_ );
      }

      __device__ __forceinline__ T & getNeighbor( int cx, int cy, int cz, int cf )
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
      dim3 idxDim_;
      dim3 blockDim_;
      bool isValidPosition_{ false };
   };



} // namespace gpu
} // namespace walberla


