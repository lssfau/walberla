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
//! \file IrregularFreeSlip.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/Cell.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"

#include "field/FlagField.h"
#include "field/GhostLayerField.h"

#include "stencil/Directions.h"

#include <memory>
#include <tuple>
#include <type_traits>

#include "walberla/experimental/Memory.hpp"
#include "walberla/experimental/sweep/SparseIndexList.hpp"

namespace walberla::sweepgen
{
struct IrregularFreeSlipLinkInfo
{
   int32_t x;
   int32_t y;
   int32_t z;
   int32_t dir;
   int32_t source_offset_x;
   int32_t source_offset_y;
   int32_t source_offset_z;
   int32_t source_dir;
   IrregularFreeSlipLinkInfo(walberla::Cell fluidCell, walberla::stencil::Direction linkDir, walberla::Cell sourceCell,
                             walberla::stencil::Direction sourceDir)
      : x{ fluidCell.x() }, y{ fluidCell.y() }, z{ fluidCell.z() }, dir{ int32_t(linkDir) },
        source_offset_x{ sourceCell.x() - fluidCell.x() }, source_offset_y{ sourceCell.y() - fluidCell.y() },
        source_offset_z{ sourceCell.z() - fluidCell.z() }, source_dir{ int32_t(sourceDir) }
   {}

   bool operator==(const IrregularFreeSlipLinkInfo& other) const
   {
      return std::tie(x, y, z, dir, source_offset_x, source_offset_y, source_offset_z, source_dir) ==
             std::tie(other.x, other.y, other.z, other.dir, other.source_offset_x, other.source_offset_y,
                      other.source_offset_z, other.source_dir);
   }
};

namespace detail
{
template< typename Stencil, typename FlagField_T, typename IndexList_T >
class FreeSlipLinksFromFlagField
{
 public:
   using flag_t    = typename FlagField_T::flag_t;

   FreeSlipLinksFromFlagField(StructuredBlockForest& sbfs, BlockDataID flagFieldID, field::FlagUID boundaryFlagUID,
                              field::FlagUID fluidFlagUID)
      : sbfs_{ sbfs }, flagFieldID_(flagFieldID), boundaryFlagUID_(boundaryFlagUID), fluidFlagUID_(fluidFlagUID)
   {}

   IndexList_T collectLinks()
   {
      IndexList_T indexList{ sbfs_ };

      for (auto& block : sbfs_)
      {
         FlagField_T* flagField = block.template getData< FlagField_T >(flagFieldID_);
         flag_t boundaryFlag{ flagField->getFlag(boundaryFlagUID_) };
         flag_t fluidFlag{ flagField->getFlag(fluidFlagUID_) };

         auto& idxVector = indexList.getVector(block);

         for (auto it = flagField->beginXYZ(); it != flagField->end(); ++it)
         {
            Cell c{ it.cell() };
            if (!isFlagSet(it, fluidFlag)) { continue; }

            for (auto dIt = Stencil::beginNoCenter(); dIt != Stencil::end(); ++dIt)
            {
               stencil::Direction dir{ *dIt };
               Cell neighbor{ c + dir };
               if (flagField->isFlagSet(neighbor, boundaryFlag))
               {
                  idxVector.emplace_back(createLink(flagField, c, dir));
               }
            }
         }
      }

      return indexList;
   }

 private:
   IrregularFreeSlipLinkInfo createLink(FlagField_T* flagField, const Cell& fluidCell, const stencil::Direction dir)
   {
      const flag_t fluidFlag{ flagField->getFlag(fluidFlagUID_) };
      const Cell wallCell{ fluidCell + dir };

      // inverse direction of 'dir' as lattice vector

      const cell_idx_t ix = stencil::cx[stencil::inverseDir[dir]];
      const cell_idx_t iy = stencil::cy[stencil::inverseDir[dir]];
      const cell_idx_t iz = stencil::cz[stencil::inverseDir[dir]];

      stencil::Direction sourceDir =
         stencil::inverseDir[dir]; // compute reflected (mirrored) of inverse direction of 'dir'

      cell_idx_t wnx = 0; // compute "normal" vector of free slip wall
      cell_idx_t wny = 0;
      cell_idx_t wnz = 0;

      if (flagField->isFlagSet(wallCell.x() + ix, wallCell.y(), wallCell.z(), fluidFlag))
      {
         wnx       = ix;
         sourceDir = stencil::mirrorX[sourceDir];
      }
      if (flagField->isFlagSet(wallCell.x(), wallCell.y() + iy, wallCell.z(), fluidFlag))
      {
         wny       = iy;
         sourceDir = stencil::mirrorY[sourceDir];
      }
      if (flagField->isFlagSet(wallCell.x(), wallCell.y(), wallCell.z() + iz, fluidFlag))
      {
         wnz       = iz;
         sourceDir = stencil::mirrorZ[sourceDir];
      }

      // concave corner (neighbors are non-fluid)
      if (wnx == 0 && wny == 0 && wnz == 0)
      {
         wnx       = ix;
         wny       = iy;
         wnz       = iz;
         sourceDir = dir;
      }

      const Cell sourceCell{ wallCell.x() + wnx, wallCell.y() + wny, wallCell.z() + wnz };

      return { fluidCell, dir, sourceCell, sourceDir };
   }

   StructuredBlockForest& sbfs_;
   const BlockDataID flagFieldID_;
   const field::FlagUID boundaryFlagUID_;
   const field::FlagUID fluidFlagUID_;
};
} // namespace detail

template< typename Impl >
struct IrregularFreeSlipFactoryImplTraits
{
   using Stencil_T   = typename Impl::Stencil;
   using IdxStruct_T = IrregularFreeSlipLinkInfo;
   using IndexListMemTag_T = typename Impl::memtag_t;
   using IndexList_T = walberla::experimental::sweep::SparseIndexList< IrregularFreeSlipLinkInfo, IndexListMemTag_T >;
};

template< typename Impl >
class IrregularFreeSlipFactory
{
   using ImplTraits = IrregularFreeSlipFactoryImplTraits< Impl >;

 public:
   IrregularFreeSlipFactory(const shared_ptr< StructuredBlockForest > blocks, BlockDataID pdfFieldID)
      : blocks_{ blocks }, pdfFieldID_{ pdfFieldID }
   {}

   template< typename FlagField_T >
   auto fromFlagField(BlockDataID flagFieldID, field::FlagUID boundaryFlagUID, field::FlagUID fluidFlagUID)
   {
      detail::FreeSlipLinksFromFlagField< typename ImplTraits::Stencil_T, FlagField_T, typename ImplTraits::IndexList_T > linksFromFlagField{
         *blocks_, flagFieldID, boundaryFlagUID, fluidFlagUID
      };
      auto indexVector = linksFromFlagField.collectLinks();
      return impl().irregularFromIndexVector(indexVector);
   }

 protected:
   shared_ptr< StructuredBlockForest > blocks_;
   BlockDataID pdfFieldID_;

 private:
   Impl& impl() { return *static_cast< Impl* >(this); }

   const Impl& impl() const { return *static_cast< const Impl* >(this); }
};
} // namespace walberla::experimental::lbm
