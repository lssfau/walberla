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
//! \file HostFieldAllocator.h
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Allocator that allocates a CPU! field using cudaHostAlloc
//
//======================================================================================================================

#pragma once

#include "ErrorChecking.h"
#include "field/allocation/FieldAllocator.h"

#include <cuda_runtime.h>


namespace walberla {
namespace cuda {


   //*******************************************************************************************************************
   /*!
   * Allocator that allocates a CPU! field using cudaHostAlloc without padding
   *
   * Uses cudaHostAlloc for the allocation - which allocates page-locked memory that is faster to transfer to the GPU
   * This allocator should be used for CPU fields that are often transfered to GPU and back
   *
   * \ingroup cuda
   *
   */
   //*******************************************************************************************************************
   template<typename T, unsigned int cudaHostAllocFlags = cudaHostAllocDefault>
   class HostFieldAllocator : public field::FieldAllocator<T>
   {
   public:
      virtual ~HostFieldAllocator() {}

      virtual T * allocateMemory (  uint_t size0, uint_t size1, uint_t size2, uint_t size3,
                                    uint_t & allocSize1, uint_t & allocSize2, uint_t & allocSize3 )
      {
         allocSize1=size1;
         allocSize2=size2;
         allocSize3=size3;
         void * result;
         WALBERLA_CUDA_CHECK( cudaHostAlloc( &result, size0*size1*size2*size3*sizeof(T), cudaHostAllocFlags ) );
         return (T*)(result);
      }

      virtual T * allocateMemory ( uint_t size )
      {
         T* result;
         cudaHostAlloc( &result, size*sizeof(T), cudaHostAllocFlags );
         return result;
      }

      virtual void deallocate(T *& values) {
         WALBERLA_CUDA_CHECK( cudaFreeHost( values ) );
      }
   };




} // namespace cuda
} // namespace walberla


