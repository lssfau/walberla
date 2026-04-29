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
//! \file MemoryTraits.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include <memory>

#include "./MemoryTags.hpp"
#include "./Allocators.hpp"

namespace walberla::v8::memory
{

constexpr size_t hostAlignment() {
#if defined(__ARM_FEATURE_SVE) && defined(__ARM_FEATURE_SVE_BITS) && __ARM_FEATURE_SVE_BITS > 0
   return static_cast<size_t>(__ARM_FEATURE_SVE_BITS / 8);
#elif defined(__ARM_FEATURE_SVE)
   return static_cast<size_t>(64);
#elif defined(__ARM_NEON)
   return static_cast<size_t>(16);
#elif defined(__AVX512F__)
   return static_cast<size_t>(64);
#elif defined(__AVX__)
   return static_cast<size_t>(32);
#elif defined(__SSE__)
   return static_cast<size_t>(16);
#elif defined(__BIGGEST_ALIGNMENT__)
   return static_cast<size_t>(__BIGGEST_ALIGNMENT__);
#else
   return static_cast<size_t>(64);
#endif
}

constexpr size_t deviceAlignment() {
#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
   return static_cast<size_t>(256);
#else
   return hostAlignment();
#endif
}


/**
 * Memory traits for a given memory tag.
 *
 * This struct template must be specialized for every instantiable memory tag.
 * Specializations must define
 *  - `AllocatorType`: A type that satisfies the standard Allocator named requirements for type `T`
 *  - `alignment()`: Static member function returning `size_t` compile-time constant holding the expected memory
 * alignment
 *  - `getAllocator()`: Static member function returning an instance of `AllocatorType`
 */
template< MemTag mem_tag_t, typename T >
struct MemoryTraits
{
#if !defined(__NVCC__)
   static_assert(false && "MemoryTraits not specialized for this memory tag");
#endif

   using AllocatorType = void;
   static constexpr size_t alignment();
   static AllocatorType getAllocator(size_t /* offset */);
};

// ------------ stdmem ------------

template< typename T >
struct MemoryTraits< memtag::stdmem, T >
{
   using AllocatorType = std::allocator< T >;
   static constexpr size_t alignment() { return alignof(T); };
   static AllocatorType getAllocator(size_t /* offset */) { return {}; }   
};

// ------------ host ------------

template< typename T >
struct MemoryTraits< memtag::host, T >
{
   using AllocatorType = HostAllocator< T >;
   static constexpr size_t alignment() { return hostAlignment(); };
   static AllocatorType getAllocator(size_t offset) { return {offset, alignment()}; }
};

// ------------ unified ------------

template< typename T >
struct MemoryTraits< memtag::unified, T >
{
#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
   using AllocatorType = UnifiedAllocator< T >;
   static constexpr size_t alignment() { return deviceAlignment(); };
#else
   using AllocatorType = HostAllocator< T >;
   static constexpr size_t alignment() { return hostAlignment(); };
#endif
   static AllocatorType getAllocator(size_t offset) { return {offset, alignment()}; }
};

// ------------ device ------------

template< typename T >
struct MemoryTraits< memtag::device, T >
{
#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
   using AllocatorType = DeviceAllocator< T >;
   static constexpr size_t alignment() { return deviceAlignment(); };
   static AllocatorType getAllocator(size_t offset) { return {offset, alignment()}; }
#else
   static_assert(never_true<T>::value, "memtag::device requires CUDA or HIP");
#endif
};

// ------------ pinned ------------

template< typename T >
struct MemoryTraits< memtag::pinned, T >
{
#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
   using AllocatorType = PinnedAllocator< T >;
   static constexpr size_t alignment() { return deviceAlignment(); };
#else
   using AllocatorType = HostAllocator< T >;
   static constexpr size_t alignment() { return hostAlignment(); };
#endif
   static AllocatorType getAllocator(size_t offset) { return {offset, alignment()}; }
};


} // namespace walberla::v8::memory