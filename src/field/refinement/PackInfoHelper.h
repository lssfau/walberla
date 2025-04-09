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
//! \file PackInfoHelper.h
//! \ingroup field
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/cell/CellInterval.h"
#include "stencil/Directions.h"
#include "blockforest/BlockNeighborhoodSection.h"

namespace walberla {
namespace field {
namespace refinement {

///////////////////////////////////////////////////////////////////////
// Helper functions for determining packing/unpacking cell intervals //
///////////////////////////////////////////////////////////////////////

inline bool divisibleByTwo(const CellInterval& cellBB)
{
   return ((cellBB.xSize() & uint_t(1)) == uint_t(0)) && ((cellBB.ySize() & uint_t(1)) == uint_t(0)) &&
          ((cellBB.zSize() & uint_t(1)) == uint_t(0));
}

//////////////////////////////
// General helper functions //
//////////////////////////////

inline bool isFaceDirection(stencil::Direction dir)
{
   return (dir == stencil::N) || (dir == stencil::S) || (dir == stencil::W) || (dir == stencil::E) ||
          (dir == stencil::T) || (dir == stencil::B);
}

inline bool isEdgeDirection(stencil::Direction dir)
{
   return (dir == stencil::NW) || (dir == stencil::NE) || (dir == stencil::SW) || (dir == stencil::SE) ||
          (dir == stencil::TN) || (dir == stencil::TS) || (dir == stencil::TW) || (dir == stencil::TE) ||
          (dir == stencil::BN) || (dir == stencil::BS) || (dir == stencil::BW) || (dir == stencil::BE);
}

inline bool isCornerDirection(stencil::Direction dir)
{
   return (dir == stencil::TNE) || (dir == stencil::TNW) || (dir == stencil::TSE) || (dir == stencil::TSW) ||
          (dir == stencil::BNE) || (dir == stencil::BNW) || (dir == stencil::BSE) || (dir == stencil::BSW);
}

inline bool blocksConnectedByFaces(const Block* block, const BlockID& neighbor)
{
   const std::array<uint_t, 6> face = { uint_t(4), uint_t(10), uint_t(12), uint_t(13), uint_t(15), uint_t(21) };
   for (unsigned long i : face)
   {
      for (uint_t n = 0; n != block->getNeighborhoodSectionSize(i); ++n)
         if (block->getNeighborId(i, n) == neighbor) return true;
   }
   return false;
}

inline bool blocksConnectedByEdges(const Block* block, const BlockID& neighbor)
{
   const std::array<uint_t, 12> face = { uint_t(1),  uint_t(3),  uint_t(5),  uint_t(7),  uint_t(9),  uint_t(11),
                           uint_t(14), uint_t(16), uint_t(18), uint_t(20), uint_t(22), uint_t(24) };

   for (unsigned long i : face)
   {
      for (uint_t n = 0; n != block->getNeighborhoodSectionSize(i); ++n)
         if (block->getNeighborId(i, n) == neighbor) return true;
   }
   return false;
}

inline bool coarserNeighborExists(const Block* block, stencil::Direction dir)
{
   if (block->getLevel() == uint_t(0)) return false;

   Vector3< int > min(-1);
   Vector3< int > max(1);

   for (uint_t i = 0; i != 3; ++i)
   {
      if (stencil::c[i][dir] == -1) max[i] = 0;
      if (stencil::c[i][dir] == 1) min[i] = 0;
   }

   for (int z = min[2]; z <= max[2]; ++z)
   {
      for (int y = min[1]; y <= max[1]; ++y)
      {
         for (int x = min[0]; x <= max[0]; ++x)
         {
            if (x == 0 && y == 0 && z == 0) continue;
            if (block->neighborhoodSectionHasLargerBlock(blockforest::getBlockNeighborhoodSectionIndex(x, y, z)))
               return true;
         }
      }
   }

   return false;
}

inline Vector3< cell_idx_t > getNeighborShift(const BlockID& smallBlock,
                                              stencil::Direction dir) // dir: direction from big to small block
{
   Vector3< cell_idx_t > shift;

   uint_t branchId = smallBlock.getBranchId();

   shift[0] = (stencil::cx[dir] == 0) ? (((branchId & uint_t(1)) == uint_t(0)) ? cell_idx_t(-1) : cell_idx_t(1)) :
                                        cell_idx_t(0);
   shift[1] = (stencil::cy[dir] == 0) ? (((branchId & uint_t(2)) == uint_t(0)) ? cell_idx_t(-1) : cell_idx_t(1)) :
                                        cell_idx_t(0);
   shift[2] = (stencil::cz[dir] == 0) ? (((branchId & uint_t(4)) == uint_t(0)) ? cell_idx_t(-1) : cell_idx_t(1)) :
                                        cell_idx_t(0);

   return shift;
}

inline CellInterval equalLevelUnpackInterval(stencil::Direction dir, const CellInterval& cellBB, const uint_t numberOfLayers)
{
   CellInterval interval(cellBB);
   interval.expand(cell_idx_c(numberOfLayers));

   for (uint_t i = 0; i != 3; ++i)
   {
      const auto c = cell_idx_c(stencil::c[i][dir]);

      if (c == -1)
         interval.max()[i] = interval.min()[i] + cell_idx_c(numberOfLayers - 1);
      else if (c == 1)
         interval.min()[i] = interval.max()[i] - cell_idx_c(numberOfLayers - 1);
      else
      {
         WALBERLA_ASSERT_EQUAL(c, cell_idx_t(0));
         interval.min()[i] += cell_idx_c(numberOfLayers);
         interval.max()[i] -= cell_idx_c(numberOfLayers);
      }
   }

   return interval;
}

inline CellInterval equalLevelPackInterval(stencil::Direction dir, const CellInterval& cellBB,
                                           const uint_t numberOfLayers)
{
   CellInterval interval = equalLevelUnpackInterval(dir, cellBB, numberOfLayers);

   for (uint_t i = 0; i != 3; ++i)
   {
      const auto offset = cell_idx_c(stencil::c[i][dir]) * cell_idx_c(numberOfLayers);
      interval.min()[i] -= offset;
      interval.max()[i] -= offset;
   }

   return interval;
}

inline CellInterval coarseToFinePackInterval(stencil::Direction dir, const CellInterval& cellBB, const BlockID& smallBlock)
{
   CellInterval interval       = equalLevelPackInterval(dir, cellBB, uint_t(1));
   Vector3< cell_idx_t > shift = getNeighborShift(smallBlock, dir);

   WALBERLA_ASSERT(divisibleByTwo(cellBB));

   for (uint_t i = 0; i != 3; ++i)
   {
      if (shift[i] == cell_idx_t(-1)) interval.max()[i] = interval.min()[i] + cell_idx_c(cellBB.size(i) / uint_t(2));
      if (shift[i] == cell_idx_t(1)) interval.min()[i] = interval.max()[i] - cell_idx_c(cellBB.size(i) / uint_t(2));
   }

   WALBERLA_ASSERT(cellBB.contains(interval));

   return interval;
}

inline CellInterval coarseToFineUnpackInterval(stencil::Direction dir, const CellInterval& cellBB, const BlockID& smallBlock)
{
   CellInterval interval       = equalLevelUnpackInterval(dir, cellBB, uint_t(2));
   Vector3< cell_idx_t > shift = getNeighborShift(smallBlock, dir);

   for (uint_t i = 0; i != 3; ++i)
   {
      if (shift[i] == cell_idx_t(-1)) interval.max()[i] += cell_idx_t(2);
      if (shift[i] == cell_idx_t(1)) interval.min()[i] -= cell_idx_t(2);
   }

#ifndef NDEBUG
   CellInterval expandedCellBB(cellBB);
   expandedCellBB.expand(cell_idx_t(2));
   WALBERLA_ASSERT(expandedCellBB.contains(interval));
#endif

   return interval;
}

inline CellInterval fineToCoarsePackInterval(stencil::Direction dir, const CellInterval& cellBB)
{
   return equalLevelPackInterval(dir, cellBB, uint_t(2));
}

inline CellInterval fineToCoarseUnpackInterval(stencil::Direction dir, const CellInterval& cellBB, const BlockID& smallBlock)
{
   CellInterval interval       = equalLevelUnpackInterval(dir, cellBB, uint_t(1));
   Vector3< cell_idx_t > shift = getNeighborShift(smallBlock, dir);

   WALBERLA_ASSERT(divisibleByTwo(cellBB));

   for (uint_t i = 0; i != 3; ++i)
   {
      if (shift[i] == cell_idx_t(-1))
         interval.max()[i] = interval.min()[i] + cell_idx_c(cellBB.size(i) / uint_t(2)) - cell_idx_t(1);
      if (shift[i] == cell_idx_t(1))
         interval.min()[i] = interval.max()[i] - cell_idx_c(cellBB.size(i) / uint_t(2)) + cell_idx_t(1);
   }

#ifndef NDEBUG
   CellInterval expandedCellBB(cellBB);
   expandedCellBB.expand(cell_idx_t(1));
   WALBERLA_ASSERT(expandedCellBB.contains(interval));
#endif

   return interval;
}

}
}
} //namespace walberla