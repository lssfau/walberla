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
//! \file GridGeometry.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/cell/Cell.h"
#include "core/cell/CellInterval.h"
#include "core/math/AABB.h"

#include "domain_decomposition/IBlock.h"

namespace walberla::v8::domain
{

/**
 * @brief Geometric properties of a cell grid
 * @ingroup v8core-device
 *
 * Represents a grid of cells with a given bounding box (`CellGrid::aabb()`), covered by a given interval of cells
 * (`CellGrid::cells()`) and exposes functions to compute the bounding boxes and center points of cells that are part of
 * this grid.
 *
 * @note This class is device-copyable. All its public member functions can be called from device code.
 */
class CellGrid
{
 public:
   WALBERLA_HOST_DEVICE CellGrid(AABB aabb, CellInterval cellBb, Vector3< real_t > cellExtents)
      : aabb_{ aabb }, cellBb_{ cellBb }, cellExtents_{ cellExtents }
   {}

   WALBERLA_HOST_DEVICE const AABB& aabb() const { return aabb_; }

   WALBERLA_HOST_DEVICE const CellInterval& cells() const { return cellBb_; }

   WALBERLA_HOST_DEVICE const Vector3< real_t >& cellExtents() const { return cellExtents_; }

   WALBERLA_HOST_DEVICE Vector3< real_t > cellCenter(const Cell& c) const
   {
      const auto cc = Vector3< real_t >(cellExtents_[0] * (real_c(c.x()) + real_t(0.5)),
                                        cellExtents_[1] * (real_c(c.y()) + real_t(0.5)),
                                        cellExtents_[2] * (real_c(c.z()) + real_t(0.5)));
      return aabb_.min() + cc;
   }

 private:
   AABB aabb_;
   CellInterval cellBb_;
   Vector3< real_t > cellExtents_;
};

/**
 * @brief Container for geometric properties of a simulation block.
 * @ingroup v8core-device
 *
 * @note This class is device-copyable. All its public member functions (except constructors) can be called from device
 * code.
 */
class BlockGeometry
{
 public:
   BlockGeometry(const StructuredBlockForest& blocks, const IBlock& b)
      : aabb_{ blocks.getAABB(b.getId()) },   //
        cellBb_{ blocks.getBlockCellBB(b) },  //
        refinementLevel_(blocks.getLevel(b)), //
        cellExtents_{ blocks.dx(refinementLevel_), blocks.dy(refinementLevel_), blocks.dz(refinementLevel_) }
   {}

   WALBERLA_HOST_DEVICE const AABB& aabb() const { return aabb_; }

   WALBERLA_HOST_DEVICE const CellInterval& cellBoundingBox() const { return cellBb_; }

   WALBERLA_HOST_DEVICE size_t refinementLevel() const { return refinementLevel_; }

   WALBERLA_HOST_DEVICE const Vector3< real_t >& cellExtents() const { return cellExtents_; }

   WALBERLA_HOST_DEVICE CellGrid localGrid() const
   {
      return {
         aabb_,                                                                                 //
         CellInterval{ {}, { cellBb_.xSize() - 1, cellBb_.ySize() - 1, cellBb_.zSize() - 1 } }, //
         cellExtents_                                                                           //
      };
   }

 private:
   AABB aabb_;
   CellInterval cellBb_;
   size_t refinementLevel_;
   Vector3< real_t > cellExtents_;
};

} // namespace walberla::v8::domain