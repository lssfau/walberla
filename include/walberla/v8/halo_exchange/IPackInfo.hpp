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
//! \file IPackInfo.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Stdlib.hpp"
#include WALBERLA_STDLIB(span)

#include "blockforest/BlockID.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include <concepts>
#include <span>

#include "walberla/v8/sweep/ExecutionTags.hpp"

namespace walberla::v8::halo_exchange
{

namespace packinfo_interface
{

template< typename TPackInfo, typename ExecTag >
concept HasPack = requires(const TPackInfo& pInfo, const IBlock& srcBlock, stencil::Direction commDir,
                           stdlib::span< typename TPackInfo::value_type > buffer, ExecTag xtag) {
   typename TPackInfo::value_type;
   { pInfo.pack(srcBlock, commDir, buffer, xtag) } -> std::same_as< void >;
};

template< typename TPackInfo, typename ExecTag >
concept HasUnpack = requires(const TPackInfo& pInfo, IBlock& dstBlock, stencil::Direction commDir,
                             stdlib::span< const typename TPackInfo::value_type > buffer, ExecTag xtag) {
   { pInfo.unpack(dstBlock, commDir, buffer, xtag) } -> std::same_as< void >;
};

template< typename TPackInfo, typename ExecTag >
concept HasLocalCopy = requires(const TPackInfo& pInfo, const IBlock& srcBlock, IBlock& dstBlock,
                                stencil::Direction commDir, ExecTag xtag) {
   { pInfo.localCopy(srcBlock, commDir, dstBlock, xtag) } -> std::same_as< void >;
};

template< typename TPackInfo >
concept HasPacketSize = requires(const TPackInfo& pInfo, stencil::Direction commDir) {
   { pInfo.packetSize(commDir) } -> std::same_as< size_t >;
};

template< typename TPackInfo >
concept HasSendAndReceivePacketSize =
   requires(const TPackInfo& pInfo, const IBlock& block, stencil::Direction commDir) {
      { pInfo.sendPacketSize(block, commDir) } -> std::same_as< size_t >;
      { pInfo.receivePacketSize(block, commDir) } -> std::same_as< size_t >;
   };

template< typename TPackInfo, typename ExecTag >
concept IPackInfo = (                                                          //
   std::semiregular< typename TPackInfo::value_type >                          //
   && HasPack< TPackInfo, ExecTag >                                            //
   && HasUnpack< TPackInfo, ExecTag >                                          //
   && (HasPacketSize< TPackInfo > || HasSendAndReceivePacketSize< TPackInfo >) //
   &&!(HasPacketSize< TPackInfo > && HasSendAndReceivePacketSize< TPackInfo >) //
);

template< HasPacketSize TPackInfo >
inline size_t sendPacketSize(const TPackInfo& packInfo, const IBlock& /*sender*/, stencil::Direction dir)
{
   return packInfo.packetSize(dir);
}

template< HasPacketSize TPackInfo >
inline size_t receivePacketSize(const TPackInfo& packInfo, const IBlock& /*receiver*/, stencil::Direction dir)
{
   return packInfo.packetSize(dir);
}

template< HasSendAndReceivePacketSize TPackInfo >
inline size_t sendPacketSize(const TPackInfo& packInfo, const IBlock& sender, stencil::Direction dir)
{
   return packInfo.sendPacketSize(sender, dir);
}

template< HasSendAndReceivePacketSize TPackInfo >
inline size_t receivePacketSize(const TPackInfo& packInfo, const IBlock& receiver, stencil::Direction dir)
{
   return packInfo.receivePacketSize(receiver, dir);
}

} // namespace packinfo_interface

} // namespace walberla::v8::halo_exchange