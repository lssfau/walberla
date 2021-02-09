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
//! \file FieldAllocator.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Allocation Strategies for Fields
//
//======================================================================================================================

#pragma once

#include "AlignedMalloc.h"
#include "core/debug/Debug.h"
#include "field/CMakeDefs.h"
#include "core/VectorTrait.h"

#include <map>
#include <new>


namespace walberla {
namespace field {


   template <typename T>
   class FieldAllocatorBase
   {
   protected:
      static std::map<T*, uint_t> referenceCounts_;
   };

   //*******************************************************************************************************************
   /*! Allocation Strategy base class for fields
   *
   * \ingroup field
   *
   * The base class handles reference counting, and has purely virtual member functions
   * that have to be implemented by concrete strategy.
   */
   //*******************************************************************************************************************
   template<typename T>
   class FieldAllocator : FieldAllocatorBase<typename std::conditional<VectorTrait<T>::F_SIZE!=0, typename VectorTrait<T>::OutputType, T>::type>
   {
      public:
         using BaseType = typename std::conditional<VectorTrait<T>::F_SIZE!=0, typename VectorTrait<T>::OutputType, T>::type;

         virtual ~FieldAllocator() = default;

         /**
          * \brief Allocate memory for a field of given sizes and initializes reference counter with one
          *
          * \param size*       The size of the field. size0*size1*size2*size*3 is the minimum amount of
          *                    memory that has to be allocated. size0 is the size of the outermost (slowest)
          *                    coordinate, size3 the size of the fastest.
          *
          * \param allocSize*  Output parameters. The actual size of the allocation may be bigger. Here
          *                    allocSize can be bigger than size if alignment >0. There is no allocSize0,
          *                    since extra allocations (padding) for the outer coordinate do not make sense.
          */
         T * allocate (  uint_t size0, uint_t size1, uint_t size2, uint_t size3,
                         uint_t & allocSize1, uint_t & allocSize2, uint_t & allocSize3)
         {
            T * mem = allocateMemory(size0,size1,size2,size3,allocSize1,allocSize2,allocSize3);
            BaseType * bmem = reinterpret_cast<BaseType *>(mem);
            #ifdef WALBERLA_THREAD_SAFE_FIELD_ALLOCATION
            #ifdef _OPENMP
            #pragma omp critical( walberla_field_allocator_refcount )
            #endif
            #endif
            {
               WALBERLA_ASSERT( referenceCounts_.find(bmem) == referenceCounts_.end() || referenceCounts_[bmem] == 0 );
               referenceCounts_[bmem] = 1;
            }
            return mem;
         }


         virtual void setInnerGhostLayerSize( uint_t /*innerGhostLayerSize*/ ) {}

         /**
          * \brief Allocate memory of the given size
          *
          * This is used by the field when a memory portion of the same size as an already allocated
          * region (for example in clone() )
          */
         T * allocate ( uint_t allocSize )
         {
            T * mem = allocateMemory( allocSize );
            BaseType * bmem = reinterpret_cast<BaseType *>(mem);
            
            #ifdef WALBERLA_THREAD_SAFE_FIELD_ALLOCATION
            #ifdef _OPENMP
            #pragma omp critical( walberla_field_allocator_refcount )
            #endif
            #endif
            {
               WALBERLA_ASSERT( referenceCounts_.find(bmem) == referenceCounts_.end() || referenceCounts_[bmem] == 0 );
               referenceCounts_[bmem] = 1;
            }
            return mem;
         }

         /**
          * \brief Increments the reference count for the given pointer
          *
          * \param allocatedPointer a pointer that was returned by the allocate function
          *
          * The reference counter is the number of field objects that currently use the allocated
          * memory region
          */
         void incrementReferenceCount( T * mem )
         {
            BaseType * bmem = reinterpret_cast<BaseType *>(mem);

            #ifdef WALBERLA_THREAD_SAFE_FIELD_ALLOCATION
            #ifdef _OPENMP
            #pragma omp critical( walberla_field_allocator_refcount )
            #endif
            #endif
            {
               WALBERLA_ASSERT( referenceCounts_.find(bmem) != referenceCounts_.end() );
               WALBERLA_ASSERT_GREATER( referenceCounts_[bmem], 0 );
               referenceCounts_[bmem]++;
            }
         }


         /**
          * \brief Decrements the reference count for the given pointer
          *
          * \param allocatedPointer a pointer that was returned by the allocate function
          *
          * A memory region is freed only when the reference count is zero.
          * \return true if memory was freed,
          *         false if some other field still uses the memory
          */
         bool decrementReferenceCount( T * mem )
         {
            BaseType * bmem = reinterpret_cast<BaseType *>(mem);
            bool memoryFreed = false;
            
            uint_t refCount = 0;

            #ifdef WALBERLA_THREAD_SAFE_FIELD_ALLOCATION
            #ifdef _OPENMP
            #pragma omp critical( walberla_field_allocator_refcount )
            #endif
            #endif
            {
               WALBERLA_ASSERT( referenceCounts_.find(bmem) != referenceCounts_.end() );
               WALBERLA_ASSERT_GREATER( referenceCounts_[bmem], 0 );

               refCount = --referenceCounts_[bmem];
            }

            if( refCount == 0 ) {
               deallocate( mem );
               memoryFreed = true;
            }
            
            return memoryFreed;
         }


         uint_t referenceCount ( T * mem ) const
         {
            BaseType * bmem = reinterpret_cast<BaseType *>(mem);
            uint_t refCount = 0;
            #ifdef WALBERLA_THREAD_SAFE_FIELD_ALLOCATION
            #ifdef _OPENMP
            #pragma omp critical( walberla_field_allocator_refcount )
            #endif
            #endif
            {
               WALBERLA_ASSERT( referenceCounts_.find(bmem) != referenceCounts_.end()  );
               refCount = referenceCounts_[bmem];
            }
            
            return refCount;
         }

      protected:
         /**
          * \brief Same as allocated, without handling of reference counter
          */
         virtual T * allocateMemory (  uint_t size0, uint_t size1, uint_t size2, uint_t size3,
                                       uint_t & allocSize1, uint_t & allocSize2, uint_t & allocSize3 ) = 0;


         /**
          * \brief Same as allocated, without handling of reference counter
          */
         virtual T * allocateMemory ( uint_t size ) = 0;

         /**
          * \brief Deallocates memory that has been allocated using allocate()
          *
          * \param values      Return value of allocate function()
          */
         virtual void deallocate( T *& values ) = 0;

      private:
         using FieldAllocatorBase<BaseType>::referenceCounts_;
   };

   template<typename T>
   std::map<T*, uint_t> FieldAllocatorBase<T>::referenceCounts_ = std::map<T*,uint_t>();



   /****************************************************************************************************************//**
   * Aligned Allocation strategy for Fields
   *
   * \ingroup field
   *
   * Template parameters:
   *  - T          type that is stored in field
   *  - alignment  the beginning of each row of the field is placed such that the memory
   *               address 'a' of each row fulfills: a % alignment == 0
   *               alignment has to be a power of 2
   ********************************************************************************************************************/
   template <typename T, uint_t alignment>
   class AllocateAligned : public FieldAllocator<T>
   {

      protected:

         T * allocateMemory (  uint_t size0, uint_t size1, uint_t size2, uint_t size3,
                                       uint_t & allocSize1, uint_t & allocSize2, uint_t & allocSize3) override
         {
            allocSize1=size1;
            allocSize2=size2;
            allocSize3=size3;
            uint_t lineLength = size3 * static_cast<uint_t>( sizeof(T) );
            if(lineLength % alignment !=0 )
               allocSize3 = ((lineLength + alignment) / alignment ) * (alignment / sizeof(T));

            WALBERLA_ASSERT_GREATER_EQUAL( allocSize3, size3 );
            WALBERLA_ASSERT_EQUAL( (allocSize3 * sizeof(T)) % alignment, 0 );

            return allocateMemory ( size0 * allocSize1 * allocSize2 * allocSize3 );
         }

         T * allocateMemory (  uint_t size ) override
         {
            void * ptr = aligned_malloc_with_offset(size * sizeof(T) + alignment, alignment, offset_ % alignment );
            if(!ptr)
               throw std::bad_alloc();

            // placement new
            new (ptr) T[ size ];

            T * ret = reinterpret_cast<T*>( ptr );

            #ifdef _OPENMP
            #pragma omp critical( walberla_field_aligned_allocator_nrOfElements )
            #endif
            {
               nrOfElements_[ret] = size;
            }
            return ret;
         }

         void setInnerGhostLayerSize( uint_t innerGhostLayerSize ) override {
            offset_ = sizeof(T) * innerGhostLayerSize;
         }

         void deallocate(T *& values ) override
         {
            WALBERLA_ASSERT ( nrOfElements_.find(values) != nrOfElements_.end() );

            size_t nrOfValues = 0;

            #ifdef _OPENMP
            #pragma omp critical( walberla_field_aligned_allocator_nrOfElements )
            #endif
            {
               nrOfValues = nrOfElements_[values];
            }

            for( uint_t i = 0; i < nrOfValues; ++i )
               values[i].~T();

            #ifdef _OPENMP
            #pragma omp critical( walberla_field_aligned_allocator_nrOfElements )
            #endif
            {
               nrOfElements_.erase( values );
            }

            aligned_free(values);
         }

         static_assert(alignment > 0, "Use StdFieldAlloc");
         static_assert(!(alignment & (alignment - 1)) , "Alignment has to be power of 2");

      private:
         /// Nr of elements per allocated pointer has to be stored to call the destructor on each element
         static std::map<T*, uint_t> nrOfElements_;

         uint_t offset_;
   };
   template <typename T, uint_t alignment>
   std::map<T*,uint_t> AllocateAligned<T,alignment>::nrOfElements_ = std::map<T*,uint_t>();



   /****************************************************************************************************************//**
   *  Allocator without alignment using new and delete[]
   *
   * \ingroup field
   *
   ********************************************************************************************************************/
   template <typename T>
   class StdFieldAlloc : public FieldAllocator<T>
   {
      public:
         T * allocateMemory (  uint_t size0, uint_t size1, uint_t size2, uint_t size3,
                                       uint_t & allocSize1, uint_t & allocSize2, uint_t & allocSize3) override
         {
            allocSize1=size1;
            allocSize2=size2;
            allocSize3=size3;
            return new T[size0*size1*size2*size3];
         }

         T * allocateMemory ( uint_t size ) override
         {
            return new T[size];
         }

         void deallocate(T *& values) override {
            delete[] values;
            values = nullptr;
         }
   };


} // namespace field
} // namespace walberla


