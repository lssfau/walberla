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
//! \file DisjoiningPressureBubbleModel.h
//! \ingroup bubble_model
//! \author Daniela Anderl
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Bubble Model with additional pressure term (disjoining pressure).
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"

#include <cmath>

#include "BubbleModel.h"
#include "DistanceInfo.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
 * Bubble Model with additional pressure term (disjoining pressure) which depends on distance to the nearest bubble. The
 * disjoining pressure is computed as in the paper from Koerner et al., 2005, where variable
 * - "disjPressConst_" is called "c_pi"
 * - "maxDistance_" is called "d_range"
 * - "dist" is called "d_int".
 *
 * Be aware that the value of the first two variables are phenomenological, i.e., they have to be chosen according to
 * comparisons with experiments. The current default values are more or less randomly chosen.
 **********************************************************************************************************************/
template< typename Stencil_T >
class DisjoiningPressureBubbleModel : public BubbleModel< Stencil_T >
{
 public:
   using DistanceField_T = GhostLayerField< DistanceInfo, 1 >;

   explicit DisjoiningPressureBubbleModel(const std::shared_ptr< StructuredBlockForest >& blockStorage,
                                          const real_t maxDistance    = real_c(7),
                                          const real_t disjPressConst = real_c(0.05), bool enableBubbleSplits = true,
                                          uint_t distanceFieldUpdateInterval = uint_c(1));

   ~DisjoiningPressureBubbleModel() override = default;

   // compute the gas density with disjoining pressure according to the paper from Koerner et al., 2005
   real_t getDensity(IBlock* block, const Cell& cell) const override
   {
      // get bubble field and bubble ID
      const BubbleField_T* bf = BubbleModel< Stencil_T >::getBubbleField(block);
      const BubbleID id       = bf->get(cell);

      WALBERLA_ASSERT_UNEQUAL(id, INVALID_BUBBLE_ID, "Cell " << cell << " does not contain a bubble.");

      const real_t gasDensity = BubbleModel< Stencil_T >::bubbles_[id].getDensity();

      // cache a pointer to the block and to the distance field since block->getData<DistanceField>() is expensive
      static IBlock* blockCache         = nullptr;
      static DistanceField_T* distField = nullptr;

      // update cache
      if (block != blockCache)
      {
         blockCache = block;
         distField  = block->getData< DistanceField_T >(distanceFieldSrcID_);
      }

      const DistanceInfo& distanceInfo = distField->get(cell);

      // get the distance to the nearest neighboring interface cell from a different bubble
      const real_t dist = distanceInfo.getDistanceToNearestBubble(id, maxDistance_);

      real_t disjoiningPressure = real_c(0.0);

      // computation of disjoining pressure according to the paper from Koerner et al., 2005 where variable
      // - "disjPressConst_" is called "c_pi"
      // - "maxDistance_" is called "d_range"
      // - "dist" is called "d_int"
      if (dist > maxDistance_)
      {
         // disjoining pressure is zero for bubbles outside of range "maxDistance"
         disjoiningPressure = real_c(0);
      }
      else
      {
         // disjoining pressure as Koerner et al., 2005
         disjoiningPressure = disjPressConst_ * real_c(std::fabs(maxDistance_ - dist));
      }

      return gasDensity - disjoiningPressure;
   }

   // update the bubble model
   void update() override
   {
      static uint_t step = uint_c(0);

      // update the distance field
      if (step % distanceFieldUpdateInterval_ == 0) updateDistanceField();

      // update regular bubble model
      BubbleModel< Stencil_T >::update();

      ++step;
   }

   ConstBlockDataID getDistanceFieldID() { return distanceFieldSrcID_; }

 protected:
   void updateDistanceField();

   real_t maxDistance_;    // d_range in the paper from Koerner et al., 2005
   real_t disjPressConst_; // c_pi in the paper from Koerner et al., 2005

   BlockDataID distanceFieldSrcID_;
   BlockDataID distanceFieldDstID_;

   uint_t distanceFieldUpdateInterval_;
}; // class DisjoiningPressureBubbleModel

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

#include "DisjoiningPressureBubbleModel.impl.h"
