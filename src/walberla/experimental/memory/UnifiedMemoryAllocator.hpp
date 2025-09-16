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
//! \file UnifiedMemoryAllocator.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"

#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)

#include "gpu/GPUWrapper.h"
#include "gpu/ErrorChecking.h"

namespace walberla::experimental::memory
{

/**
 * @brief Allocator for managed GPU memory.
 * 
 * This allocator implementes the Allocator named requirement and
 * manages allocation and disposal of memory blocks with unified
 * memory semantics, using the CUDA or HIP runtime (depending on which is active).
 */
template< typename T >
class UnifiedMemoryAllocator
{
 public:
   using value_type = T;
   using pointer_type = T *;


   pointer_type allocate(size_t n) {
      pointer_type ptr;
      size_t bytes { n * sizeof(value_type) };
      WALBERLA_GPU_CHECK(
         gpuMallocManaged(&ptr, bytes)
      );
      return ptr;
   }

   void deallocate(pointer_type ptr, size_t /*n*/) {
      WALBERLA_GPU_CHECK(
         gpuFree(ptr)
      );
   }

   bool operator== (const UnifiedMemoryAllocator< value_type > &) {
      return true;
   }
};

} // namespace walberla::experimental::memory

#endif
