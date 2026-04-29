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
//! \file CommSchemeTraits.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "./AbstractCommScheme.hpp"
#include "./LegacyCpuScheme.hpp"
#include "./LegacyGpuScheme.hpp"
#include "walberla/v8/memory/MemoryTags.hpp"

namespace walberla::v8::halo_exchange
{

template< typename TCommStencil, memory::MemTag TMemTag >
struct CommSchemeTraits
{
#if !defined(__NVCC__)
   static_assert(false && "CommSchemeTraits not specialized for this memory tag");
#endif

   using CommStencil     = void;
   using CommScheme      = void;
   using PackInfoWrapper = void;

   template< typename T >
   static consteval bool isPackInfoType();

   template< typename TPackInfo >
   static std::shared_ptr< PackInfoWrapper > wrapPackInfo(TPackInfo);

   static std::unique_ptr< AbstractCommScheme > createCommScheme(const std::shared_ptr< StructuredBlockForest >&,
                                                                 std::vector< std::shared_ptr< PackInfoWrapper > >);
};

template< typename TCommStencil >
struct CommSchemeTraits< TCommStencil, memtag::stdmem > : public legacy_cpu::LegacyCpuSchemeTraits< TCommStencil >
{};

template< typename TCommStencil >
struct CommSchemeTraits< TCommStencil, memtag::host > : public legacy_cpu::LegacyCpuSchemeTraits< TCommStencil >
{};

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
template< typename TCommStencil >
struct CommSchemeTraits< TCommStencil, memtag::device > : public legacy_gpu::LegacyGpuSchemeTraits< TCommStencil >
{};

template< typename TCommStencil >
struct CommSchemeTraits< TCommStencil, memtag::unified > : public legacy_gpu::LegacyGpuSchemeTraits< TCommStencil >
{};
#endif

} // namespace walberla::v8::halo_exchange