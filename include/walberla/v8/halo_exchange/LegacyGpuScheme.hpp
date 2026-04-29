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
//! \file LegacyGpuScheme.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Stdlib.hpp"
#include WALBERLA_STDLIB(span)

#include "blockforest/BlockID.h"
#include "blockforest/StructuredBlockForest.h"

#include "domain_decomposition/IBlock.h"

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
#include "gpu/communication/GeneratedGPUPackInfo.h"
#include "gpu/communication/UniformGPUScheme.h"
#endif

#include "stencil/Directions.h"

#include "./AbstractCommScheme.hpp"
#include "./IPackInfo.hpp"
#include "walberla/v8/memory/Allocators.hpp"
#include "walberla/v8/sweep/ExecutionTags.hpp"

namespace walberla::v8::halo_exchange::legacy_gpu
{

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)

template< typename TPackInfo >
concept IDevicePackInfo = packinfo_interface::IPackInfo< TPackInfo, exectag::GPU >;

template< IDevicePackInfo TImpl >
class LegacyGpuPackInfoWrapper : public gpu::GeneratedGPUPackInfo
{
 public:
   using value_type   = typename TImpl::value_type;
   using pointer_type = value_type*;

   explicit LegacyGpuPackInfoWrapper(TImpl impl) : impl_{ std::move(impl) } {}

   void pack(stencil::Direction dir, unsigned char* buffer, IBlock* sender, gpuStream_t stream) override
   {
      void* const rawBuffer = buffer;

      WALBERLA_DEBUG_SECTION()
      {
         const size_t expectedAlignment{ alignof(value_type) };
         if (std::ptrdiff_t(rawBuffer) % expectedAlignment != 0)
         {
            throw std::runtime_error{ std::format(
               "Buffer alignment did not match requirements of message data type.\n"
               "Expected alignment: {}, Actual address: {}\n"
               "This is a bug in the legacy UniformGpuScheme - we apologize for the inconvenience.",
               expectedAlignment, rawBuffer) };
         }
      }

      const size_t numValues{ packinfo_interface::sendPacketSize(impl_, *sender, dir) };
      const pointer_type typedBuffer = static_cast< pointer_type >(rawBuffer);
      const stdlib::span< value_type > bufferWindow{ typedBuffer, numValues };
      impl_.pack(*sender, dir, bufferWindow, exectag::GPU{ stream });
   }

   void unpack(stencil::Direction dirToSender, unsigned char* buffer, IBlock* receiver, gpuStream_t stream) override
   {
      void* const rawBuffer = buffer;

      WALBERLA_DEBUG_SECTION()
      {
         const size_t expectedAlignment{ alignof(value_type) };
         if (std::ptrdiff_t(rawBuffer) % expectedAlignment != 0)
         {
            throw std::runtime_error{ std::format(
               "Buffer alignment did not match requirements of message data type.\n"
               "Expected alignment: {}, Actual address: {}\n"
               "This is a bug in the legacy UniformGpuScheme - we apologize for the inconvenience.",
               expectedAlignment, rawBuffer) };
         }
      }

      const stencil::Direction dir{ stencil::inverseDir[dirToSender] };
      const size_t numValues{ packinfo_interface::receivePacketSize(impl_, *receiver, dir) };
      const pointer_type typedBuffer = static_cast< pointer_type >(rawBuffer);
      const stdlib::span< const value_type > bufferWindow{ typedBuffer, numValues };
      impl_.unpack(*receiver, dir, bufferWindow, exectag::GPU{ stream });
   }

   void communicateLocal(stencil::Direction dir, const IBlock* sender, IBlock* receiver, gpuStream_t stream) override
   {
      if constexpr (packinfo_interface::HasLocalCopy< TImpl, exectag::GPU >)
      {
         impl_.localCopy(*sender, dir, *receiver, exectag::GPU{ stream });
      }
      else
      {
         const size_t numValues{ packinfo_interface::sendPacketSize(impl_, *sender, dir) };

         memory::DeviceAllocator< value_type > alloc;
         pointer_type buffer = alloc.allocate(numValues);

         exectag::GPU gex{ stream };

         stdlib::span< value_type > packView{ buffer, numValues };
         impl_.pack(*sender, dir, packView, gex);

         stdlib::span< const value_type > unpackView{ buffer, numValues };
         impl_.unpack(*receiver, dir, unpackView, gex);

         gex.sync();

         alloc.deallocate(buffer, numValues);
      }
   }

   uint_t size(stencil::Direction dir, IBlock* sender) override
   {
      return sizeof(value_type) * packinfo_interface::sendPacketSize(impl_, *sender, dir);
   }

 private:
   TImpl impl_;
};

template< typename TCommStencil >
class LegacyGpuScheme : public AbstractCommScheme, public gpu::communication::UniformGPUScheme< TCommStencil >
{
 public:
   using BaseCommScheme = gpu::communication::UniformGPUScheme< TCommStencil >;
   using BaseCommScheme::BaseCommScheme;

   void startCommunication() override { static_cast< BaseCommScheme& >(*this).startCommunication(); }
   void wait() override { static_cast< BaseCommScheme& >(*this).wait(); }
};

template< typename TCommStencil >
struct LegacyGpuSchemeTraits
{
   using CommStencil     = TCommStencil;
   using CommScheme      = LegacyGpuScheme< TCommStencil >;
   using PackInfoWrapper = gpu::GeneratedGPUPackInfo;

   template< typename T >
   static consteval bool isPackInfoType()
   {
      return IDevicePackInfo< T >;
   }

   template< IDevicePackInfo TPackInfo >
   static std::shared_ptr< PackInfoWrapper > wrapPackInfo(TPackInfo pInfo)
   {
      return std::make_shared< LegacyGpuPackInfoWrapper< TPackInfo > >(std::move(pInfo));
   }

   static std::unique_ptr< AbstractCommScheme >
      createCommScheme(const std::shared_ptr< StructuredBlockForest >& blocks,
                       std::vector< std::shared_ptr< PackInfoWrapper > > pInfos)
   {
      auto scheme = std::make_unique< CommScheme >(blocks);
      for (auto& pInfo : pInfos)
      {
         scheme->addPackInfo(pInfo);
      }

      return std::move(scheme);
   }
};

#endif

} // namespace walberla::v8::halo_exchange::legacy_gpu