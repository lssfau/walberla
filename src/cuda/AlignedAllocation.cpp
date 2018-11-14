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
//! \file AlignedAllocation.cpp
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "AlignedAllocation.h"
#include "cuda/ErrorChecking.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

#include <map>

namespace walberla {
namespace cuda {

   static std::map<void *, void*> freePointers_;

   void *allocate_aligned_with_offset( uint_t size, uint_t alignment, uint_t offset )
   {
      // With 0 alignment this function makes no sense
      // use normal malloc instead
      WALBERLA_ASSERT_GREATER( alignment, 0 );
      // Tests if alignment is power of two (assuming alignment>0)
      WALBERLA_ASSERT( !(alignment & (alignment - 1)) );

      WALBERLA_ASSERT_LESS( offset, alignment );

      if( offset == 0 )
      {
         void * result = nullptr;
         WALBERLA_CUDA_CHECK( cudaMalloc( &result, size ) );
         freePointers_[result] = result;
         return result;
      }

      void *pa;  // pointer to allocated memory
      void *ptr; // pointer to usable aligned memory

      WALBERLA_CUDA_CHECK( cudaMalloc( &pa, size + alignment ));
      WALBERLA_CHECK_EQUAL(size_t(pa) % alignment, 0 , "CUDA malloc did not return memory with requested alignment");
      ptr = (void *) ((char *) (pa) + alignment - offset);
      freePointers_[ptr] = pa;

      WALBERLA_ASSERT_EQUAL(((size_t) ptr + offset) % alignment, 0 );
      return ptr;
   }


   void free_aligned_with_offset( void *ptr )
   {
      // assume that pointer to real allocated chunk is stored just before
      // chunk that was given to user
      WALBERLA_CUDA_CHECK( cudaFree( freePointers_[ptr] ));
      freePointers_.erase(ptr);
   }


   void *allocate_pitched_with_offset( size_t &pitchOut, size_t width, size_t height,
                                              size_t alignment, size_t alignmentOffset )
   {
      if( width % alignment == 0)
         pitchOut = width;
      else
         pitchOut = ((width + alignment) / alignment ) * alignment;

      WALBERLA_ASSERT_GREATER_EQUAL( pitchOut, width );
      WALBERLA_ASSERT_EQUAL( pitchOut  % alignment, 0 );

      return allocate_aligned_with_offset( pitchOut * height, alignment, alignmentOffset );
   }


} // namespace cuda
} // namespace walberla

