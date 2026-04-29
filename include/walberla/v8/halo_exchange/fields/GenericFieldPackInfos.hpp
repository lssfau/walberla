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
//! \file GenericFieldPackInfo.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Stdlib.hpp"
#include WALBERLA_STDLIB(span)

#include "stencil/Directions.h"

#include "../HaloExchange.hpp"
#include "walberla/v8/StencilRanges.hpp"
#include "walberla/v8/Sweep.hpp"
#include "walberla/v8/memory/Field.hpp"
#include "walberla/v8/memory/MemoryTags.hpp"

namespace walberla::v8::halo_exchange
{

namespace fields
{

template< memory::IField TField, sweep::ExecTag XTag >
   requires(!TField::IS_CONST)
class GenericFieldPackInfo
{
 public:
   using FieldType                = TField;
   using value_type               = typename FieldType::value_type;
   static constexpr size_t F_SIZE = FieldType::F_SIZE;

   /**
    * @param field The field whose ghost layers to exchange
    * @param numGhostLayers Number of ghost layers to exchange. Must not exceed the field's total number of ghost
    * layers.
    */
   GenericFieldPackInfo(const FieldType& field, size_t numGhostLayers)
      : field_{ field }, numGhostLayers_{ numGhostLayers }
   {
      if (numGhostLayers_ == 0 || numGhostLayers_ > field.numGhostLayers())
      {
         throw std::invalid_argument{ "Invalid number of ghost layers" };
      }
   }

   explicit GenericFieldPackInfo(const FieldType& field) : GenericFieldPackInfo(field, field.numGhostLayers()) {}

   void pack(const IBlock& srcBlock, stencil::Direction commDir, stdlib::span< value_type > buffer, XTag xtag) const
   {
      const ConstFieldView srcView{ field_, srcBlock };
      const CellInterval srcRegion{ srcView.slices().borderSlice(commDir, numGhostLayers_) };

      constexpr size_t valuesPerCell{ F_SIZE };

      // TODO place mdspan over buffer for md-indexing
      // TODO Honor memory layout of field (fzyx or zyxf)

      sweep::forAllCells(xtag, srcRegion, [=] WALBERLA_HOST_DEVICE (Cell cell) {
         const Cell rc{ cell - srcRegion.min() };
         const size_t segmentStart(
            valuesPerCell * (uint_c(rc.x()) + srcRegion.xSize() * (uint_c(rc.y()) + srcRegion.ySize() * uint_c(rc.z())))
         );

         for (cell_idx_t q = 0; q < cell_idx_t(F_SIZE); ++q)
         {
            buffer[segmentStart + size_t(q)] = srcView(cell, q);
         }
      });
   }

   void unpack(IBlock& dstBlock, stencil::Direction commDir, stdlib::span< const value_type > buffer, XTag xtag) const
   {
      const FieldView dstView{ field_, dstBlock };
      const CellInterval dstRegion{ dstView.slices().ghostSlice(stencil::inverseDir[commDir], numGhostLayers_) };

      constexpr size_t valuesPerCell{ F_SIZE };

      sweep::forAllCells(xtag, dstRegion, [=] WALBERLA_HOST_DEVICE (Cell cell) {
         const Cell rc{ cell - dstRegion.min() };
         const size_t segmentStart(
            valuesPerCell * (uint_c(rc.x()) + dstRegion.xSize() * (uint_c(rc.y()) + dstRegion.ySize() * uint_c(rc.z())))
         );

         for (cell_idx_t q = 0; q < cell_idx_t(F_SIZE); ++q)
         {
            dstView(cell, q) = buffer[segmentStart + size_t(q)];
         }
      });
   }

   void localCopy(const IBlock& srcBlock, stencil::Direction commDir, IBlock& dstBlock, XTag xtag) const
   {
      const ConstFieldView srcView{ field_, srcBlock };
      const CellInterval srcRegion{ srcView.slices().borderSlice(commDir, numGhostLayers_) };

      const FieldView dstView{ field_, dstBlock };
      const CellInterval dstRegion{ dstView.slices().ghostSlice(stencil::inverseDir[commDir], numGhostLayers_) };

      sweep::forAllCellPairs(xtag, srcRegion, dstRegion, [=] WALBERLA_HOST_DEVICE (Cell srcCell, Cell dstCell) {
         for (cell_idx_t q = 0; q < cell_idx_t(F_SIZE); ++q)
         {
            dstView(dstCell, q) = srcView(srcCell, q);
         }
      });
   }

   size_t packetSize(stencil::Direction commDir) const
   {
      size_t numCells{ field_.slices().borderSlice(commDir, numGhostLayers_).numCells() };
      return numCells * F_SIZE;
   }

 private:
   FieldType field_;
   size_t numGhostLayers_;
};

template< IField TField, typename TStencil, sweep::ExecTag XTag >
   requires(TField::F_SIZE == TStencil::Q)
class GenericStreamPullPackInfo
{
 public:
   using FieldType  = TField;
   using Stencil    = TStencil;
   using value_type = typename FieldType::value_type;

   /**
    * @param field The field whose ghost layers to exchange
    * @param numGhostLayers Number of ghost layers to exchange. Must not exceed the field's total number of ghost
    * layers.
    */
   explicit GenericStreamPullPackInfo(const FieldType& field) : field_{ field }
   {
      if (field.numGhostLayers() == 0) { throw std::invalid_argument{ "Field must have at least one ghost layer" }; }
   }

   void pack(const IBlock& srcBlock, stencil::Direction commDir, stdlib::span< value_type > buffer, XTag xtag) const
   {
      const ConstFieldView srcView{ field_, srcBlock };
      const CellInterval srcRegion{ srcView.slices().borderSlice(commDir, 1) };

      sweep::forAllCells(xtag, srcRegion, [=] WALBERLA_HOST_DEVICE (Cell cell) {
         const Cell rc{ cell - srcRegion.min() };

         const auto subdirs = stencil_ranges::subdirections< Stencil >(commDir);
         const size_t valuesPerCell{ subdirs.size() };
         const size_t segmentStart{
            valuesPerCell * size_t(rc.x() + cell_idx_c(srcRegion.xSize()) * (rc.y() + cell_idx_c(srcRegion.ySize()) * rc.z()))
         };

         for (size_t i = 0; i < subdirs.size(); ++i)
         {
            const auto q                     = subdirs[i];
            buffer[segmentStart + i] = srcView(cell, q);
         }
      });
   }

   void unpack(IBlock& dstBlock, stencil::Direction commDir, stdlib::span< const value_type > buffer, XTag xtag) const
   {
      const FieldView dstView{ field_, dstBlock };
      const CellInterval dstRegion{ dstView.slices().ghostSlice(stencil::inverseDir[commDir], 1) };

      sweep::forAllCells(xtag, dstRegion, [=] WALBERLA_HOST_DEVICE (Cell cell) {
         const Cell rc{ cell - dstRegion.min() };

         const auto subdirs = stencil_ranges::subdirections< Stencil >(commDir);
         const size_t valuesPerCell{ subdirs.size() };
         const size_t segmentStart{
            valuesPerCell * size_t(rc.x() + cell_idx_c(dstRegion.xSize()) * (rc.y() + cell_idx_c(dstRegion.ySize()) * rc.z()))
         };

         for (size_t i = 0; i < subdirs.size(); ++i)
         {
            const auto q     = subdirs[i];
            dstView(cell, q) = buffer[segmentStart + i];
         }
      });
   }

   void localCopy(const IBlock& srcBlock, stencil::Direction commDir, IBlock& dstBlock, XTag xtag) const
   {
      const ConstFieldView srcView{ field_, srcBlock };
      const CellInterval srcRegion{ srcView.slices().borderSlice(commDir, 1) };

      const FieldView dstView{ field_, dstBlock };
      const CellInterval dstRegion{ dstView.slices().ghostSlice(stencil::inverseDir[commDir], 1) };

      sweep::forAllCellPairs(xtag, srcRegion, dstRegion, [=] WALBERLA_HOST_DEVICE (Cell srcCell, Cell dstCell) {
         for (cell_idx_t q : stencil_ranges::subdirections< Stencil >(commDir))
         {
            dstView(dstCell, q) = srcView(srcCell, q);
         }
      });
   }

   size_t packetSize(stencil::Direction commDir) const
   {
      size_t numCells{ field_.slices().borderSlice(commDir, 1).numCells() };
      return numCells * stencil_ranges::subdirections< Stencil >(commDir).size();
   }

 private:
   FieldType field_;
};

} // namespace fields

} // namespace walberla::v8::halo_exchange