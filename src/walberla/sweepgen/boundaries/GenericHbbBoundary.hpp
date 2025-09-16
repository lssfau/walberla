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
//! \file GenericHbbBoundary.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"
#include "core/cell/Cell.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"

#include "stencil/Directions.h"

#include <memory>
#include <tuple>

#include "walberla/experimental/sweep/SparseIndexList.hpp"
#include "walberla/experimental/sweep/Sweeper.hpp"

namespace walberla::sweepgen
{

struct HbbLink
{
   int32_t x;
   int32_t y;
   int32_t z;
   int32_t dir;

   /**
    * Create a halfway bounce-back boundary link from a fluid cell `c` to a boundary cell in direction `d`.
    */
   HbbLink(const Cell& c, const stencil::Direction d)
      : x{ int32_c(c.x()) }, y{ int32_c(c.y()) }, z{ int32_c(c.z()) }, dir{ int32_t(d) }
   {}

   bool operator==(const HbbLink& other) const
   {
      return std::tie(x, y, z, dir) == std::tie(other.x, other.y, other.z, other.dir);
   }
};

template< typename Impl >
struct GenericHbbFactoryImplTraits
{
   using Stencil_T         = typename Impl::Stencil;
   using IdxStruct_T       = HbbLink;
   using IndexListMemTag_T = typename Impl::memtag_t;
   using IndexList_T       = walberla::experimental::sweep::SparseIndexList< IdxStruct_T, IndexListMemTag_T >;
};

/**
 * @brief Represents a potential boundary link.
 */
struct PotentialBoundaryLink
{
   const Block & block;
   Cell fluidCell;
   stencil::Direction dir;
   Cell wallCell;
};

template< typename F >
concept LinkPredicate = requires(F f, PotentialBoundaryLink link) {
   { f(link) } -> std::same_as< bool >;
};

template< typename Impl >
class GenericHbbFactory
{
   using ImplTraits = GenericHbbFactoryImplTraits< Impl >;

 public:
   GenericHbbFactory(const shared_ptr< StructuredBlockForest > blocks)
      : blocks_{ blocks }
   {}

   template< typename FlagField_T >
   auto fromFlagField(BlockDataID flagFieldID, field::FlagUID boundaryFlagUID, field::FlagUID fluidFlagUID)
   {
      using flag_t      = typename FlagField_T::flag_t;
      using IndexList_T = typename ImplTraits::IndexList_T;
      using Stencil     = typename ImplTraits::Stencil_T;

      IndexList_T indexList{ *blocks_ };

      for (auto& block : *blocks_)
      {
         FlagField_T* flagField = block.template getData< FlagField_T >(flagFieldID);
         flag_t boundaryFlag{ flagField->getFlag(boundaryFlagUID) };
         flag_t fluidFlag{ flagField->getFlag(fluidFlagUID) };

         if (!(flagField->flagExists(boundaryFlagUID) && flagField->flagExists(fluidFlagUID))) continue;

         auto& idxVector = indexList.getVector(block);

         for (auto it = flagField->beginXYZ(); it != flagField->end(); ++it)
         {
            Cell fluidCell{ it.cell() };
            if (!isFlagSet(it, fluidFlag)) { continue; }

            for (auto dIt = Stencil::beginNoCenter(); dIt != Stencil::end(); ++dIt)
            {
               stencil::Direction dir{ *dIt };
               Cell neighbor{ fluidCell + dir };
               if (flagField->isFlagSet(neighbor, boundaryFlag)) { idxVector.emplace_back(fluidCell, dir); }
            }
         }
      }

      return impl().irregularFromIndexVector(indexList);
   }

   /**
    * @brief Construct a boundary condition sweep by selecting boundary links via a predicate
    * 
    * Construct a boundary condition sweep by iterating over all lattice links
    * and setting the boundary condition on all links where the given predicate `pred`
    * evaluates to true.
    * 
    * @param pred Predicate to select boundary links.
    */
   template< LinkPredicate P >
   auto selectLinks(P pred)
   {
      using IndexList_T = typename ImplTraits::IndexList_T;
      using Stencil     = typename ImplTraits::Stencil_T;

      IndexList_T indexList{ *blocks_ };

      walberla::experimental::sweep::SerialSweeper sweeper{ blocks_ };

      for (auto& ib : *blocks_)
      {
         Block & block = dynamic_cast< Block & >(ib);
         auto& idxVector = indexList.getVector(block);

         sweeper.forAllCells([&](Cell fluidCell) {
            for (auto dIt = Stencil::beginNoCenter(); dIt != Stencil::end(); ++dIt)
            {
               stencil::Direction dir{ *dIt };
               Cell wallCell{ fluidCell + dir };
               PotentialBoundaryLink link { .block=block, .fluidCell=fluidCell, .dir=dir, .wallCell=wallCell };
               if (pred(link)) { idxVector.emplace_back(fluidCell, dir); }
            }
         });
      }

      return impl().irregularFromIndexVector(indexList);
   }

 protected:
   shared_ptr< StructuredBlockForest > blocks_;

 private:
   Impl& impl() { return *static_cast< Impl* >(this); }

   const Impl& impl() const { return *static_cast< const Impl* >(this); }
};

} // namespace walberla::sweepgen
