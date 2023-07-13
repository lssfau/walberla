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
//! \file CustomMemoryBuffer.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Basic Buffer supporting different memory spaces
//
//======================================================================================================================

#pragma once

#include <algorithm>
#include <cstring>

#include "gpu/ErrorChecking.h"
#include "gpu/GPUWrapper.h"

namespace walberla {
namespace gpu
{
namespace communication {


   struct HostMemoryAllocator;
   struct DeviceMemoryAllocator;

   //*******************************************************************************************************************
   /*!
    * Simple buffer class that supports memory allocators, e.g. for pinned host memory or GPU memory
    *
    * \ingroup gpu
    *
    * In contrast to core::mpi::Buffer this class does not support stream operators "<<" and ">>" because these
    * operators imply serial (un)packing which is not feasible on the GPU.
    * The allocator template has to provide:
    *   - static void *allocate( size_t size )
    *   - void deallocate( void *ptr )
    *   - void memcpy( void *dst, void *src, size_t count )
    *
    * The buffer has a beginning, a current position and an end position. Here is an overview of the most important
    * operations:
    *   - clear: reset current position to begin, does not change size
    *   - advance: moves current position number of bytes forward and returns pointer to the old current position
    *              two versions are available, one that automatically resizes and reallocates the buffer, and one that
    *              fails if not enough space is available
    */
   //*******************************************************************************************************************
   template<typename Allocator>
   class CustomMemoryBuffer
   {
   public:
      using ElementType = uint8_t;

      explicit CustomMemoryBuffer();
      explicit CustomMemoryBuffer( std::size_t initSize );
      explicit CustomMemoryBuffer( const CustomMemoryBuffer &pb );
      ~CustomMemoryBuffer();
      CustomMemoryBuffer &operator=( const CustomMemoryBuffer &pb );

      void resize( std::size_t newSize );
      inline std::size_t allocSize() const { return std::size_t(end_ - begin_); }
      inline std::size_t size() const { return std::size_t(cur_ - begin_); }
      inline std::size_t remainingSize() const { return std::size_t(end_ - cur_); }
      ElementType *ptr() const { return begin_; }
      ElementType *cur() const { return cur_; }

      inline void clear() { cur_ = begin_; }

      ElementType *advance( std::size_t bytes );
      ElementType *advanceNoResize( std::size_t bytes );

      template<typename T>
      T *advance( std::size_t bytes ) { return reinterpret_cast<T *>( advance( bytes * sizeof( T ))); }
      template<typename T>
      T *advanceNoResize( std::size_t bytes ) { return reinterpret_cast<T *>( advanceNoResize( bytes * sizeof( T ))); }

   private:
      ElementType *begin_ = nullptr;
      ElementType *cur_ = nullptr;
      ElementType *end_ = nullptr;
   };


   using PinnedMemoryBuffer = CustomMemoryBuffer<HostMemoryAllocator>;
   using GPUMemoryBuffer    = CustomMemoryBuffer<DeviceMemoryAllocator>;


   struct HostMemoryAllocator
   {
      static void *allocate( size_t size )
      {
         void *p;
         WALBERLA_GPU_CHECK( gpuMallocHost( &p, size ))
         return p;
      }

      static void deallocate( void *ptr )
      {
         WALBERLA_GPU_CHECK( gpuFreeHost( ptr ))
      }

      static void memcpy( void *dst, void *src, size_t count )
      {
         std::memcpy( dst, src, count );
      }
   };

   struct DeviceMemoryAllocator
   {
      static void *allocate( size_t size )
      {
         void *p;
         WALBERLA_GPU_CHECK( gpuMalloc( &p, size ))
         return p;
      }

      static void deallocate( void *ptr )
      {
         WALBERLA_GPU_CHECK( gpuFree( ptr ))
      }

      static void memcpy( void *dst, void *src, size_t count )
      {
         WALBERLA_GPU_CHECK( gpuMemcpy( dst, src, count, gpuMemcpyDeviceToHost ) )
      }
   };


} // namespace communication
} // namespace gpu
} // namespace walberla

#include "CustomMemoryBuffer.impl.h"
