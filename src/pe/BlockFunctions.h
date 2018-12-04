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
//! \file BlockFunctions.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/rigidbody/Owner.h"

#include "core/math/Vector3.h"
#include "core/DataTypes.h"

namespace walberla {
namespace pe {

/**
 * Checks if one of the neighbor blocks is owned by process with rank \a rank.
 */
template <class BlockT>
bool hasNeighborOwner(const BlockT& block, const Owner& owner)
{
   for( uint_t i = uint_t(0); i != block.getNeighborhoodSize(); ++i )
   {
      if ((owner.rank_ == int_c( block.getNeighborProcess(i)) && (owner.blockID_ == block.getNeighborId(i)) )) return true;
   }
   return false;
}

/**
 * Looks through all neighboring blocks to find the one whose AABB contains \a pt.
 * Also checks if \a pt is located in the block itself.
 * Returns -1 if no valid block is found otherwise the process rank of the containing block is returned.
 *
 * \attention If periodic boundaries are used you have to make sure the point is mapped to the domain before calling this function!
 */
template <class BlockT>
Owner findContainingProcess(const BlockT& block, const math::Vector3<real_t> pt)
{
   WALBERLA_DEBUG_SECTION()
   {
      auto pt2 = pt;
      block.getBlockStorage().mapToPeriodicDomain(pt2);
      WALBERLA_ASSERT_EQUAL(pt, pt2);
   }

   if (block.getAABB().contains(pt)) return Owner(int_c(block.getProcess()), block.getId().getID());
   for( uint_t i = uint_t(0); i != block.getNeighborhoodSize(); ++i )
   {
      if (block.getNeighborAABB(i).contains(pt)) return Owner(int_c(block.getNeighborProcess(i)), block.getNeighborId(i).getID());
   }
   return Owner();
}

}  // namespace pe
}  // namespace walberla
