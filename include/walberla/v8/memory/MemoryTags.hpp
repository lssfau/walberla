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

#include <concepts>
#include <memory>

#include "waLBerlaDefinitions.h"

/**
 * @defgroup MemoryTags Memory Tags
 * @ingroup v8core-memory
 * @brief Control allocation and placement of buffers between host and accelerator memory spaces.
 *
 * The *memory tag* protocol is used by all memory-managing simulation data structures of the V8 core library.
 * It is a compile-time protocol where the desired memory semantics are communicated to a class by passing a *memory
 * tag* as a template parameter.
 * 
 * Each memory tag specifies whether the allocated memory is accessible from the host and/or the device 
 * via two compile-time boolean flags: `isHostAccessible` and `isDeviceAccessible`.
 * 
 * Depending on waLBerla's build settings, different memory tags are available.
 */

namespace walberla::v8
{
namespace memory
{

namespace memtag
{

/**
 * @brief Common base class for all memory tags.
 */
struct _mem_tag
{};

/**
 * @brief Standard pageable host memory.
 * @ingroup MemoryTags
 *
 * Represents CPU memory without alignment requirements, as provided by `std::allocator`.
 *
 * @note `stdmem` is mostly meant for debugging and testing purposes. For actual CPU simulation,
 *       prefer aligned host memory via `memtag::host`.
 */
struct stdmem : public _mem_tag
{
    static constexpr bool isHostAccessible = true;
    static constexpr bool isDeviceAccessible = false;
};
inline stdmem stdmem_v;

/**
 * @brief Aligned pageable host memory.
 * @ingroup MemoryTags
 *
 * Represents CPU memory aligned according to the specific alignment requirements of the target architecture.
 */
struct host : public _mem_tag
{
    static constexpr bool isHostAccessible = true;
    static constexpr bool isDeviceAccessible = false;
};
inline host host_v;

/**
 * @brief Managed GPU memory with unified memory semantics.
 * @ingroup MemoryTags
 *
 */
struct unified : public _mem_tag
{
    static constexpr bool isHostAccessible = true;
    static constexpr bool isDeviceAccessible = true;
};
inline unified unified_v;

/**
 * @brief Memory tag indicating GPU device memory.
 * @ingroup MemoryTags
 */
struct device : public _mem_tag
{
    static constexpr bool isHostAccessible = false;
    static constexpr bool isDeviceAccessible = true;
};
inline device device_v;

/**
 * @brief Memory tag indicating pinned (page-locked) host memory.
 * @ingroup MemoryTags
 */
struct pinned : public _mem_tag
{
    static constexpr bool isHostAccessible = true;
    static constexpr bool isDeviceAccessible = false;
};
inline pinned pinned_v;

/**
 * @brief Automatically select a memory tag according to the current build configuration
 * @ingroup MemoryTags
 * 
 * Select a memory tag according to the build configuration:
 *  - If CUDA or HIP are enabled, use `unified` memory
 *  - If no GPU language is enabled, use `host` memory
 */
using automatic =
#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
    unified;
#else 
    host;
#endif

} // namespace memtag

/**
 * @brief Concept identifying memory tag types.
 * @ingroup MemoryTags
 */
template< typename T >
concept MemTag = std::derived_from< T, memtag::_mem_tag > &&
    requires {
        { T::isHostAccessible   } -> std::same_as<const bool&>;
        { T::isDeviceAccessible } -> std::same_as<const bool&>;
    };

} // namespace memory

namespace memtag = memory::memtag;

} // namespace walberla::v8
