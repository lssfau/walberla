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
//! \file MemoryTags.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"

#include <concepts>
#include <memory>

#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)

#include "./UnifiedMemoryAllocator.hpp"

#endif


namespace walberla::experimental::memory
{

namespace memtag
{

/**
 * @brief Common base class for all memory tags.
 */
struct _mem_tag
{};

/**
 * @brief Regular host memory tag.
 */
struct host : public _mem_tag
{};
inline host host_v;

/**
 * @brief Memory tag indicating managed GPU memory with unified memory semantics.
 */
struct unified : public _mem_tag
{};
inline unified unified_v;

/**
 * @brief Memory tag indicating GPU device memory.
 */
struct device : public _mem_tag
{};
inline device device_v;

} // namespace memtag

template< typename T >
concept MemTag = std::derived_from< T, memtag::_mem_tag >;

namespace detail
{
template< MemTag memtag_t >
struct SelectStandardAllocator;

template<>
struct SelectStandardAllocator< memtag::host >
{
   template< typename T >
   using allocator = std::allocator< T >;
};

#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)

/* If CUDA or HIP are active, use the managed memory allocator for unified memory */

template<>
struct SelectStandardAllocator< memtag::unified >
{
   template< typename T >
   using allocator = UnifiedMemoryAllocator< T >;
};

#else

/* If no GPU support is active, use std::allocator */

template<>
struct SelectStandardAllocator< memtag::unified >
{
   template< typename T >
   using allocator = std::allocator< T >;
};

#endif

} // namespace detail

/**
 * @brief Standard allocator implementation for a given memory tag and value type.
 *
 * This typedef alias references a type implementing the Allocator named requirement for type T
 * which implements the memory semantics required by the given memory tag.
 */
template< MemTag memtag_t, typename T >
using StandardAllocator_T = detail::SelectStandardAllocator< memtag_t >::template allocator< T >;

} // namespace walberla::experimental::memory
