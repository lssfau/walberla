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
//! \file Initializers.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include <algorithm>

#include "./BufferSystem.hpp"
#include "./MemoryTags.hpp"

namespace walberla::v8::memory::initializers
{

template< MemTag TMemTag, typename T >
struct BufferInitializer
{
#if !defined(__NVCC__)
   static_assert(false && "BufferInitializer not specialized for this memory tag");
#endif
   
   using type = void;
};

template< typename T >
class SerialBufferInitializer
{
 public:
   SerialBufferInitializer(T initVal) : initVal_{ initVal } {}

   template< IBufferView TView >
      requires(std::same_as< typename TView::element_type, T >)
   void operator()(const TView& view)
   {
      std::ranges::fill(view.linearAllocView(), initVal_);
   }

 private:
   T initVal_;
};


#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)

namespace detail{
#if defined(__CUDACC__) || defined(__HIPCC__)
   template< typename T >
   __global__ void initKernel(stdlib::span< T > linView,  T val){
      size_t i = blockIdx.x * blockDim.x + threadIdx.x;
      if (i < linView.size()) linView[i] = val;
   }

   template< typename T >
   void initKernelWrapper(size_t bpg, size_t tpb, stdlib::span< T > linView, T val){
      detail::initKernel<<<bpg, tpb>>>(linView, val);
      WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
   }
#else
   template< typename T >
   void initKernelWrapper(size_t /*bpg*/, size_t /*tpb*/, stdlib::span< T > /*linView*/, T /*val*/){
      static_assert(never_true<T>::value, 
         "Instantiation of 'initKernelWrapper' during host code compilation is invalid. "
         "Make sure your source file is compiled with nvcc or hipcc.");
   }
#endif
} // namespace detail

template< typename T >
class GPUBufferInitializer
{
 public:
   GPUBufferInitializer(T initVal) : initVal_{ initVal } {}

   template< IBufferView TView >
      requires(std::same_as< typename TView::element_type, T >)
   void operator()(const TView& view)
   {
      size_t allocSize{ view.linearAllocView().size() };
      size_t tpb{ 256 };
      size_t bpg{ (tpb + allocSize - 1) / tpb };
      detail::initKernelWrapper(bpg, tpb, view.linearAllocView(), initVal_);
   }

 private:
   T initVal_;
};

#endif

template< typename T >
struct BufferInitializer< memtag::stdmem, T >
{
   using type = SerialBufferInitializer< T >;
};

template< typename T >
struct BufferInitializer< memtag::host, T >
{
   using type = SerialBufferInitializer< T >;
};

template< typename T >
struct BufferInitializer< memtag::unified, T >
{
#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
   using type = GPUBufferInitializer< T >;
#else
   using type = SerialBufferInitializer< T >;
#endif

};

template< typename T >
struct BufferInitializer< memtag::device, T >
{
#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
   using type = GPUBufferInitializer< T >;
#else
   static_assert(never_true<T>::value, "BufferInitializer: memtag::device requires CUDA or HIP");
#endif
};

template< typename T >
struct BufferInitializer< memtag::pinned, T >
{
   using type = SerialBufferInitializer< T >;
};

template< MemTag TMemTag, typename T >
using BufferInitializer_T = BufferInitializer< TMemTag, T >::type;

} // namespace walberla::v8::memory::initializers