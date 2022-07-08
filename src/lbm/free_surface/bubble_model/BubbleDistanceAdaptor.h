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
//! \file BubbleDistanceAdaptor.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Get the distance to a certain bubble ID.
//
//======================================================================================================================

#pragma once

#include "field/adaptors/GhostLayerFieldAdaptor.h"

#include "DisjoiningPressureBubbleModel.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
 * Get the distance to a certain bubble ID. The search region is limited by maxDistance.
 **********************************************************************************************************************/
class BubbleDistanceAdaptionFunction
{
 public:
   using basefield_t        = DisjoiningPressureBubbleModel::DistanceField;
   using basefield_iterator = basefield_t::const_base_iterator;
   using value_type         = real_t;

   static const uint_t F_SIZE = 1u;

   BubbleDistanceAdaptionFunction(BubbleID ownBubbleID, real_t maxDistance)
      : bubbleID_(ownBubbleID), maxDistance_(maxDistance)
   {
      WALBERLA_ASSERT_GREATER(maxDistance_, real_c(1.0));
   }

   value_type operator()(const basefield_t& baseField, cell_idx_t x, cell_idx_t y, cell_idx_t z,
                         cell_idx_t /*f*/ = 0) const
   {
      WALBERLA_ASSERT_GREATER(maxDistance_, real_c(1.0));
      return baseField.get(x, y, z).getDistanceToNearestBubble(bubbleID_, maxDistance_);
   }

   value_type operator()(const basefield_iterator& it) const
   {
      WALBERLA_ASSERT_GREATER(maxDistance_, real_c(1.0));
      return it->getDistanceToNearestBubble(bubbleID_, maxDistance_);
   }

 private:
   BubbleID bubbleID_;
   real_t maxDistance_;
}; // class BubbleDistanceAdaptionFunction

using BubbleDistanceAdaptor = field::GhostLayerFieldAdaptor< BubbleDistanceAdaptionFunction, 0 >;

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
