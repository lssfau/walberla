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
//! \file HaloExchange.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include <memory>

#include "./AbstractCommScheme.hpp"
#include "./CommSchemeTraits.hpp"

#include "walberla/v8/memory/MemoryTags.hpp"

namespace walberla::v8
{

namespace halo_exchange
{

class HaloExchange;

template< typename TEntity, memory::MemTag TMemTag >
struct HaloExchangeTraits;

template< typename TEntity, typename TMemTag >
concept HasHaloExchangeTraits = requires(const TEntity& entity) {
   {
      HaloExchangeTraits< TEntity, TMemTag >::createPackInfo(entity)
   } -> std::same_as< typename HaloExchangeTraits< TEntity, TMemTag >::PackInfoType >;
};

/**
 * @brief Builder for the `HaloExchange` class
 * @ingroup v8core-haloexchange
 */
template< typename TCommStencil, memory::MemTag TMemTag >
class HaloExchangeBuilder
{
 public:
   using CsTraits        = CommSchemeTraits< TCommStencil, TMemTag >;
   using PackInfoWrapper = typename CsTraits::PackInfoWrapper;

   using Self = HaloExchangeBuilder< TCommStencil, TMemTag >;

   using CommStencil = TCommStencil;
   using memory_tag  = TMemTag;

   explicit HaloExchangeBuilder(const std::shared_ptr< StructuredBlockForest >& blocks) : blocks_{ blocks } {}

   template< typename TPackInfo >
      requires(CsTraits::template isPackInfoType< TPackInfo >())
   Self& sync(TPackInfo pi)
   {
      packInfos_.push_back(CsTraits::wrapPackInfo(std::move(pi)));
      return *this;
   }

   template< HasHaloExchangeTraits< memory_tag > TEntity >
      requires(CsTraits::template isPackInfoType< typename HaloExchangeTraits< TEntity, TMemTag >::PackInfoType >())
   Self& sync(TEntity entity)
   {
      packInfos_.push_back(CsTraits::wrapPackInfo(HaloExchangeTraits< TEntity, TMemTag >::createPackInfo(entity)));
      return *this;
   }

   HaloExchange build();

   std::shared_ptr< HaloExchange > makeShared() {
      return std::make_shared< HaloExchange >(build());
   }

 private:
   std::shared_ptr< StructuredBlockForest > blocks_;
   std::vector< std::shared_ptr< PackInfoWrapper > > packInfos_;
};

/**
 * @brief Primary facilitator for exchange of halo-layer data between adjacent blocks
 * @ingroup v8core-haloexchange
 */
class HaloExchange
{
 public:

   /**
    * @brief Create a new halo-exchange object through its builder
    * 
    * @tparam TCommStencil Stencil defining the communication neighborhood;
    *         i.e. the set of neighbors with whom each block should exchange data
    * @tparam TMemTag Memory tag that controls the placement of the communication buffers
    */
   template< typename TCommStencil, memory::MemTag TMemTag >
   static HaloExchangeBuilder< TCommStencil, TMemTag > create(const std::shared_ptr< StructuredBlockForest >& blocks)
   {
      return HaloExchangeBuilder< TCommStencil, TMemTag >{ blocks };
   }

   explicit HaloExchange(std::unique_ptr< AbstractCommScheme > scheme) : scheme_{ std::move(scheme) } {}

   void startCommunication() { scheme_->startCommunication(); }
   void wait() { scheme_->wait(); }

   void communicate()
   {
      startCommunication();
      wait();
   }

   void operator()() { communicate(); }

 private:
   std::unique_ptr< AbstractCommScheme > scheme_;
};

template< typename TCommStencil, memory::MemTag TMemTag >
inline HaloExchange HaloExchangeBuilder< TCommStencil, TMemTag >::build()
{
   return HaloExchange(CsTraits::createCommScheme(blocks_, std::move(packInfos_)));
}

} // namespace halo_exchange

using halo_exchange::HaloExchange;

} // namespace walberla::v8
