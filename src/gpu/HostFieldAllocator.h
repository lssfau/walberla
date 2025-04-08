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
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Allocator that allocates a CPU! field using gpuHostAlloc
//
//======================================================================================================================

#pragma once

#include "gpu/ErrorChecking.h"
#include "gpu/DeviceWrapper.h"
#include "field/allocation/FieldAllocator.h"


namespace walberla {
namespace gpu
{


   //*******************************************************************************************************************
   /*!
   * Allocator that allocates a CPU! field using gpuHostAlloc without padding.
   *
   * Uses gpuHostAlloc for the allocation - which allocates page-locked memory that is faster to transfer to the GPU.
   * This allocator should be used for CPU fields that are often transferred to GPU and back.
   *
   * \ingroup gpu
   *
   */
   //*******************************************************************************************************************
   template<typename T, unsigned int HostAllocFlags = gpuHostAllocDefault>
   class HostFieldAllocator : public field::FieldAllocator<T>
   {
   public:
      ~HostFieldAllocator() override = default;

      T * allocateMemory (  uint_t size0, uint_t size1, uint_t size2, uint_t size3,
                            uint_t & allocSize1, uint_t & allocSize2, uint_t & allocSize3 ) override
      {
         WALBERLA_NON_DEVICE_SECTION()
         {
            WALBERLA_ABORT(__FUNCTION__ << "Using GPU method without WALBERLA_BUILD_WITH_GPU_SUPPORT being enabled.")
         }

         allocSize1=size1;
         allocSize2=size2;
         allocSize3=size3;
         void * result = nullptr;
         WALBERLA_GPU_CHECK(gpuHostAlloc(&result, size0 * size1 * size2 * size3 * sizeof(T), HostAllocFlags))
         return (T*)(result);
      }

      T * allocateMemory ( uint_t size ) override
      {
         WALBERLA_NON_DEVICE_SECTION()
         {
            WALBERLA_ABORT(__FUNCTION__ << "Using GPU method without WALBERLA_BUILD_WITH_GPU_SUPPORT being enabled.")
         }

         void * result = nullptr;
         WALBERLA_GPU_CHECK(gpuHostAlloc(&result, size*sizeof(T), HostAllocFlags))
         return (T*)(result);
      }

      void deallocate(T *& values) override {
         WALBERLA_NON_DEVICE_SECTION() {
            WALBERLA_ABORT(__FUNCTION__ << "Using GPU method without WALBERLA_BUILD_WITH_GPU_SUPPORT being enabled.")
         }
         WALBERLA_GPU_CHECK(gpuFreeHost(values))
      }
   };




} // namespace gpu
} // namespace walberla


