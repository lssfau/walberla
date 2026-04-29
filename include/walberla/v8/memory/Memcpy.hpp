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
//! \file Memcpy.hpp
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "./BufferSystem.hpp"
#include "./Field.hpp"

namespace walberla::v8::memory
{

enum class MemcpyKind {
   HostToHost,
   HostToDevice,
   DeviceToHost,
   DeviceToDevice,
   Default
};

namespace detail {

template< MemTag DstTag, MemTag SrcTag >
consteval MemcpyKind deduceMemcpyKind()
{
   constexpr bool srcHostOnly = !SrcTag::isDeviceAccessible &&  SrcTag::isHostAccessible;
   constexpr bool srcDevOnly =   SrcTag::isDeviceAccessible && !SrcTag::isHostAccessible;
   constexpr bool dstHostOnly = !DstTag::isDeviceAccessible &&  DstTag::isHostAccessible;
   constexpr bool dstDevOnly =   DstTag::isDeviceAccessible && !DstTag::isHostAccessible;

   if constexpr (srcHostOnly && dstHostOnly) return MemcpyKind::HostToHost;
   if constexpr (srcHostOnly && dstDevOnly ) return MemcpyKind::HostToDevice;
   if constexpr (srcDevOnly  && dstHostOnly) return MemcpyKind::DeviceToHost;
   if constexpr (srcDevOnly  && dstDevOnly ) return MemcpyKind::DeviceToDevice;

   return MemcpyKind::Default;   // Default only if either src or dst is unified
}

template< MemcpyKind kind >
inline void memcpyLinear(void* dst, const void* src, size_t bytes)
{
   if      constexpr (kind == MemcpyKind::HostToHost)     std::memcpy(dst, src, bytes);
   else if constexpr (kind == MemcpyKind::HostToDevice)   WALBERLA_GPU_CHECK(gpuMemcpy(dst, src, bytes, gpuMemcpyHostToDevice))
   else if constexpr (kind == MemcpyKind::DeviceToHost)   WALBERLA_GPU_CHECK(gpuMemcpy(dst, src, bytes, gpuMemcpyDeviceToHost))
   else if constexpr (kind == MemcpyKind::DeviceToDevice) WALBERLA_GPU_CHECK(gpuMemcpy(dst, src, bytes, gpuMemcpyDeviceToDevice))
   else if constexpr (kind == MemcpyKind::Default)        WALBERLA_GPU_CHECK(gpuMemcpy(dst, src, bytes, gpuMemcpyDefault))
}

} // detail

template< IBufferView TBViewDst, IBufferView TBViewSrc >
   requires(!TBViewDst::IS_CONST && std::is_same_v< typename TBViewDst::value_type, typename TBViewSrc::value_type>)
void bufferCpy(const TBViewDst& bViewDst, const TBViewSrc& bViewSrc)
{
   if (bViewDst.indexing() != bViewSrc.indexing())
      throw std::runtime_error("bufferCpy: BufferIndexing of source and destination must be identical");
   
   detail::memcpyLinear< detail::deduceMemcpyKind< typename TBViewDst::memory_tag, typename TBViewSrc::memory_tag >() >(
      static_cast< void* >(bViewDst.allocData()), 
      static_cast< const void* >(bViewSrc.allocData()), 
      bViewDst.indexing().linearBufferAllocSize() * sizeof(typename TBViewDst::value_type));
}

template< IFieldView TFViewDst, IFieldView TFViewSrc >
void fieldCpy(const TFViewDst& fViewDst, const TFViewSrc& fViewSrc){
   bufferCpy(fViewDst.bufferView(), fViewSrc.bufferView());
}


} // namespace walberla::v8::memory