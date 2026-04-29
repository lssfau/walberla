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
//! \file PackInfoSelection.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "./GenericFieldPackInfos.hpp"

namespace walberla::v8::halo_exchange
{

template< memory::IField TField >
struct HaloExchangeTraits< TField, memtag::stdmem >
{
   using PackInfoType = fields::GenericFieldPackInfo< TField, exectag::Serial >;

   static PackInfoType createPackInfo(const TField& field) { return PackInfoType(field); }
};

template< memory::IField TField >
struct HaloExchangeTraits< TField, memtag::host > : public HaloExchangeTraits< TField, memtag::stdmem >
{};

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
template< memory::IField TField >
struct HaloExchangeTraits< TField, memtag::device >
{
   using PackInfoType = fields::GenericFieldPackInfo< TField, exectag::GPU >;

   static PackInfoType createPackInfo(const TField& field) { return PackInfoType(field); }
};
#endif

template< memory::IField TField >
struct HaloExchangeTraits< TField, memtag::unified > : public HaloExchangeTraits< TField, memtag::device >
{};

namespace detail
{

template< typename TStencil >
struct _StreamPullSyncFactory
{
   template< IField TField >
      requires(std::same_as< typename TField::memory_tag, memtag::stdmem > ||
               std::same_as< typename TField::memory_tag, memtag::host >)
   auto operator()(const TField& field)
   {
      return fields::GenericStreamPullPackInfo< TField, TStencil, exectag::Serial >{ field };
   }

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   template< IField TField >
      requires(std::same_as< typename TField::memory_tag, memtag::device > ||
               std::same_as< typename TField::memory_tag, memtag::unified >)
   auto operator()(const TField& field)
   {
      return fields::GenericStreamPullPackInfo< TField, TStencil, exectag::GPU >{ field };
   }
#endif

};
} // namespace detail

template< typename TStencil >
inline detail::_StreamPullSyncFactory< TStencil > streamPullSync;

} // namespace walberla::v8::halo_exchange