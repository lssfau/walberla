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
//! \file AlignedMalloc.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Implementation of aligned memory allocation
//
//======================================================================================================================

#include "AlignedMalloc.h"

#include "core/debug/Debug.h"


namespace walberla {
namespace field {


   void *aligned_malloc( uint_t size, uint_t alignment )
   {
      // With 0 alignment this function makes no sense
      // use normal malloc instead
      WALBERLA_ASSERT_GREATER( alignment, 0 );
      // Tests if alignment is power of two (assuming alignment>0)
      WALBERLA_ASSERT( !(alignment & (alignment - 1)) );


      void *pa;  // pointer to allocated memory
      void *ptr; // pointer to usable aligned memory

      pa = std::malloc((size + alignment - 1) + sizeof( void * ));
      if(!pa)
         return nullptr;

      // Find next aligned position, starting at pa+sizeof(void*)-1
      ptr=(void*)( ((size_t)pa+sizeof(void *)+alignment-1) & ~(alignment-1));

      // Store pointer to real allocated chunk just before usable chunk
      *((void **)ptr-1)=pa;

      WALBERLA_ASSERT_EQUAL( ((size_t)ptr) % alignment, 0 );

      return ptr;
   }

   void *aligned_malloc_with_offset( uint_t size, uint_t alignment, uint_t offset )
   {
      // With 0 alignment this function makes no sense
      // use normal malloc instead
      WALBERLA_ASSERT_GREATER( alignment, 0 );
      // Tests if alignment is power of two (assuming alignment>0)
      WALBERLA_ASSERT( !(alignment & (alignment - 1)) );

      WALBERLA_ASSERT_LESS( offset, alignment );

      if( offset == 0 )
         return aligned_malloc( size, alignment );


      void *pa;  // pointer to allocated memory
      void *ptr; // pointer to usable aligned memory

      pa=std::malloc( (size+2*alignment-1 )+sizeof(void *));
      if(!pa)
         return nullptr;

      // Find next aligned position, starting at pa+sizeof(void*)-1
      ptr=(void*)( ((size_t)pa+sizeof(void *)+alignment-1) & ~(alignment-1));
      ptr=(void*) ( (char*)(ptr) + alignment - offset);

      // Store pointer to real allocated chunk just before usable chunk
      *((void **)ptr-1)=pa;

      WALBERLA_ASSERT_EQUAL( ((size_t)ptr+offset) % alignment, 0 );

      return ptr;
   }


   void aligned_free( void *ptr )
   {
      // assume that pointer to real allocated chunk is stored just before
      // chunk that was given to user
      if(ptr)
         std::free(*((void **)ptr-1));
   }

} // namespace field
} // namespace walberla


