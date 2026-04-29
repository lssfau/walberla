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
//! \file LegacyCpuScheme.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Stdlib.hpp"
#include WALBERLA_STDLIB(span)

#include "blockforest/BlockID.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "communication/UniformPackInfo.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include "./AbstractCommScheme.hpp"
#include "./IPackInfo.hpp"
#include "walberla/v8/sweep/ExecutionTags.hpp"

namespace walberla::v8::halo_exchange::legacy_cpu
{

template< typename TPackInfo >
concept IHostPackInfo = packinfo_interface::IPackInfo< TPackInfo, exectag::Serial >;

template< IHostPackInfo TImpl >
class LegacyHostPackInfoWrapper : public communication::UniformPackInfo
{
 public:
   using value_type   = typename TImpl::value_type;
   using pointer_type = value_type*;

   explicit LegacyHostPackInfoWrapper(TImpl impl) : impl_{ std::move(impl) } {}

   bool constantDataExchange() const override { return true; }
   bool threadsafeReceiving() const override { return true; }

   void unpackData(IBlock* receiver, stencil::Direction dirToSender, mpi::RecvBuffer& buffer) override
   {
      const stencil::Direction dir{ stencil::inverseDir[dirToSender] };
      const size_t numValues{ packinfo_interface::receivePacketSize(impl_, *receiver, dir) };
      size_t reservedBytes{ (numValues + 1) * sizeof(value_type) };
      void* bytesPtr                   = static_cast< void* >(buffer.skip(reservedBytes));
      const pointer_type bufferPointer = static_cast< pointer_type >(
         std::align(alignof(value_type), numValues * sizeof(value_type), bytesPtr, reservedBytes));
      const stdlib::span< const value_type > bufferWindow{ bufferPointer, numValues };

      impl_.unpack(*receiver, dir, bufferWindow, exectag::Serial{});
   }

   void communicateLocal(const IBlock* sender, IBlock* receiver, stencil::Direction dir) override
   {
      if constexpr (packinfo_interface::HasLocalCopy< TImpl, exectag::Serial >)
      {
         impl_.localCopy(*sender, dir, *receiver, exectag::Serial{});
      }
      else
      {
         const size_t numValues{ packinfo_interface::sendPacketSize(impl_, *sender, dir) };
         std::vector< value_type > buffer;
         buffer.resize(numValues);
         impl_.pack(*sender, dir, stdlib::span< value_type >(buffer), exectag::Serial{});
         impl_.unpack(*receiver, dir, stdlib::span< const value_type >(buffer), exectag::Serial{});
      }
   }

 protected:
   void packDataImpl(const IBlock* sender, stencil::Direction dir, mpi::SendBuffer& buffer) const override
   {
      const size_t numValues{ packinfo_interface::sendPacketSize(impl_, *sender, dir) };
      size_t reservedBytes{ (numValues + 1) * sizeof(value_type) };
      void* bytesPtr                   = static_cast< void* >(buffer.forward(reservedBytes));
      const pointer_type bufferPointer = static_cast< pointer_type >(
         std::align(alignof(value_type), numValues * sizeof(value_type), bytesPtr, reservedBytes));
      const stdlib::span< value_type > bufferWindow{ bufferPointer, numValues };

      impl_.pack(*sender, dir, bufferWindow, exectag::Serial{});
   }

 private:
   TImpl impl_;
};

template< typename TCommStencil >
class LegacyCpuScheme : public AbstractCommScheme,
                        public blockforest::communication::UniformBufferedScheme< TCommStencil >
{
 public:
   using BaseCommScheme = blockforest::communication::UniformBufferedScheme< TCommStencil >;
   using BaseCommScheme::BaseCommScheme;

   void startCommunication() override { static_cast< BaseCommScheme& >(*this).startCommunication(); }
   void wait() override { static_cast< BaseCommScheme& >(*this).wait(); }
};

template< typename TCommStencil >
struct LegacyCpuSchemeTraits
{
   using CommStencil     = TCommStencil;
   using CommScheme      = LegacyCpuScheme< TCommStencil >;
   using PackInfoWrapper = communication::UniformPackInfo;

   template< typename T >
   static consteval bool isPackInfoType()
   {
      return IHostPackInfo< T >;
   }

   template< IHostPackInfo TPackInfo >
   static std::shared_ptr< PackInfoWrapper > wrapPackInfo(TPackInfo pInfo)
   {
      return std::make_shared< LegacyHostPackInfoWrapper< TPackInfo > >(std::move(pInfo));
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

} // namespace walberla::v8::halo_exchange::legacy_cpu