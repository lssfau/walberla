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
//! \file Allocators.hpp
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "gpu/ErrorChecking.h"

#include <new>
#include <map>

namespace walberla::v8::memory {

/**
 * @defgroup v8core-allocators Allocators
 * @ingroup v8core-memory
 * @brief Aligned memory allocation with offset, with support for differnet backends.
 * 
 * The `Allocator` class handles alignment and pointer bookkeeping, 
 * while the backends are responsible only for raw memory allocation and deallocation. 
 * 
 */

/**
 * @brief Concept for allocator backends.
 * @ingroup v8core-allocators
 *
 * A type satisfies `AllocatorBackend` if it provides:
 *  - `allocMemory(void*& ptr, size_t bytes)`   : allocates `bytes` of memory and stores the pointer in `ptr`. May throw.
 *  - `freeMemory(void* ptr)`                   : frees the memory pointed to by `ptr`. Must not throw.
 */
template< typename T >
concept AllocatorBackend = requires(T backend, void*& ptr, size_t bytes) {
    {backend.allocMemory(ptr, bytes)} -> std::same_as<void>;
    {backend.freeMemory(ptr)} noexcept -> std::same_as<void>;     
};

/**
 * @brief Concept for types that can be allocated by `Allocator`.
 * @ingroup v8core-allocators
 *
 * A type satisfies `AllocatableType` if it is a cv-unqualified object type.
 */
template<typename T>
concept AllocatableType =
    std::is_object_v<T> &&
    !std::is_const_v<T> &&
    !std::is_volatile_v<T>;

/**
 * @brief Aligned memory allocator with offset.
 * @ingroup v8core-allocators
 *
 * @tparam T       Element type to allocate.
 * @tparam Backend Allocator backend satisfying the `AllocatorBackend` concept.
 */
template< AllocatableType T, AllocatorBackend Backend >
class Allocator
{
public:
    Allocator() noexcept = default; 

    /**
     * @brief Construct an allocator with a given offset and alignment.
     * @param offset    Number of bytes into the data block which sepcifies the aligned region.
     * @param alignment Required byte alignment.
     * @throws std::runtime_error if `alignment` is less than `alignof(T)` or not a power of two.
     */
    Allocator(size_t offset, size_t alignment) :  offset_(offset), align_(alignment){
        if(align_ < alignof(T))     throw std::runtime_error("Alignment must be at least as large as alignof(T)");
        if(align_ & (align_ - 1))   throw std::runtime_error("Alignment must be a power of two");
    }

    using value_type = T;

    /**
     * @brief Allocate memory for `n` elements of type `T`.
     *
     * The returned (outer) pointer satisfies the offset and alignment requirement such that:
     * @code
     * (reinterpret_cast<uintptr_t>(outer) + offset) % alignment == 0
     * @endcode
     *  
     * The original (base) pointer is stored internally and recovered during deallocation.
     * Memory layout:
     * @code
     *                       <-------------------------- (data) n*sizeof(T) -------------------------->
     * <------ unused ------><--------- offset --------->
     * [base]================[outer]=====================[inner (aligned)]=============================
     * @endcode
     * 
     * @param n Number of elements to allocate.
     * @return  Pointer to the outer memory region.
     * @throws std::bad_array_new_length if the allocation size overflows.
     * @throws std::bad_alloc if the offset is larger than the data itself.
     */
    [[nodiscard]] T* allocate(size_t n) {
        if(n > std::numeric_limits<size_t>::max() / sizeof(T))  throw std::bad_array_new_length();
        if(offset_ > n*sizeof(T))   throw std::bad_alloc();

        // We need at most (alignment - 1) extra bytes to shift the 
        // inner pointer to the next aligned address (only align up)
        size_t totalSize =  n*sizeof(T) + align_ - 1;
        void* basePtr = nullptr;
        backend_.allocMemory(basePtr, totalSize);
        uintptr_t base = reinterpret_cast<uintptr_t>(basePtr);
        uintptr_t inner = base + offset_;                    // find inital (unaligned) inner address
        uintptr_t innerAligned = (inner + align_ - 1) & ~(align_ -1); // align up inner address
        uintptr_t outer = innerAligned - offset_;            // find outer for the aligned inner
        
        // store the basePtr corresponding to the returned (outer) pointer in the map
        void* outerPtr = reinterpret_cast<void*>(outer);
        outerToBase[outerPtr] = basePtr;
        return reinterpret_cast<T*>(outer);
    }

    /**
     * @brief Deallocate memory previously allocated by `allocate`.
     *
     * Recovers the original base pointer from the internal map and frees it. 
     *
     * @param ptr Pointer previously returned by `allocate`.
     */
    void deallocate(T* ptr, size_t /* n */) noexcept {
        void* outerPtr = reinterpret_cast<void*>(ptr);
        auto it = outerToBase.find(outerPtr);

        if(it == outerToBase.end()) 
            WALBERLA_ABORT("Deacllocation failed: Unkown pointer")

        void* basePtr = it->second;
        backend_.freeMemory(basePtr);
        outerToBase.erase(it);
    }

    size_t offset() const noexcept { return offset_; }
    size_t alignment() const noexcept { return align_; }

    // Allocators of the same type compare equal since they share the same outerToBase map
    // meaning memory allocated by one can be deallocated by the other.
    bool operator==(const Allocator&) const noexcept { return true; }
    bool operator!=(const Allocator&) const noexcept { return false; }

    template <typename U>
    Allocator(const Allocator<U, Backend>& other) noexcept
        : Allocator(other.offset(), other.alignment()) {}

    template <typename U>
    struct rebind { using other = Allocator<U, Backend>; };

private:
    size_t offset_ = 0;
    size_t align_ = alignof(T); 
    Backend backend_;
    static inline std::map<void*, void*> outerToBase;
};


/**
 * @defgroup v8core-allocator-backends Allocator Backends
 * @ingroup v8core-allocators
 * @brief Raw memory allocation/deallocation for different backends.
 *
 * Each backend satisfies the `AllocatorBackend` concept and is responsible only
 * for raw allocation and deallocation. They are not meant to be used directly;
 * use the provided allocator type aliases instead (e.g. `HostAllocator`, `DeviceAllocator`).
 */
namespace allocator_backends {

/// @brief Allocator backend for pageable host memory. @ingroup v8core-allocator-backends
struct HostMemory{
    void allocMemory(void*& ptr, size_t bytes) {
        ptr = ::operator new(bytes);
    }

    void freeMemory(void* ptr) noexcept {
        ::operator delete(ptr);
    }
};

#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)

/// @brief Allocator backend for unified (managed) GPU memory. @ingroup v8core-allocator-backends
struct UnifiedMemory{
    void allocMemory(void*& ptr, size_t bytes) noexcept {
        WALBERLA_GPU_CHECK( gpuMallocManaged( &ptr, bytes ) )
    }

    void freeMemory(void* ptr) noexcept {
        WALBERLA_GPU_CHECK( gpuFree( ptr ) )
    }
};

/// @brief Allocator backend for GPU device memory. @ingroup v8core-allocator-backends
struct DeviceMemory {
    void allocMemory(void*& ptr, size_t bytes) noexcept {
        WALBERLA_GPU_CHECK( gpuMalloc( &ptr, bytes ) )
    }

    void freeMemory(void* ptr) noexcept {
        WALBERLA_GPU_CHECK( gpuFree( ptr ) )
    }
};

/// @brief Allocator backend for pinned (page-locked) host memory. @ingroup v8core-allocator-backends
struct PinnedMemory{
    void allocMemory(void*& ptr, size_t bytes) noexcept {
        WALBERLA_GPU_CHECK( gpuMallocHost( &ptr, bytes ) )
    }

    void freeMemory(void* ptr) noexcept {
        WALBERLA_GPU_CHECK( gpuFreeHost( ptr ) )
    }
};

#endif
} // namespace allocator_backends

/// @brief Aligned allocator with offset for pageable host memory. @ingroup v8core-allocators
template<AllocatableType T>
using HostAllocator = Allocator<T, allocator_backends::HostMemory>;

#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
/// @brief Aligned allocator with offset for unified (managed) GPU memory. @ingroup v8core-allocators
template<AllocatableType T>
using UnifiedAllocator = Allocator<T, allocator_backends::UnifiedMemory>;

/// @brief Aligned allocator with offset for GPU device memory. @ingroup v8core-allocators
template<AllocatableType T>
using DeviceAllocator = Allocator<T, allocator_backends::DeviceMemory>;

/// @brief Aligned allocator with offset for pinned (page-locked) host memory. @ingroup v8core-allocators
template<AllocatableType T>
using PinnedAllocator = Allocator<T, allocator_backends::PinnedMemory>;

#endif

} // namespace walberla::v8::memory
