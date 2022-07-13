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
//! \file DisjoiningPressureBubbleModel.impl.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Bubble Model with additional pressure term (disjoining pressure).
//
//======================================================================================================================

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "DisjoiningPressureBubbleModel.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
template< typename Stencil_T >
DisjoiningPressureBubbleModel< Stencil_T >::DisjoiningPressureBubbleModel(
   const std::shared_ptr< StructuredBlockForest >& blockStorage, real_t maxDistance, real_t disjPressConst,
   bool enableBubbleSplits, uint_t distanceFieldUpdateInterval)
   : BubbleModel< Stencil_T >(blockStorage, enableBubbleSplits), maxDistance_(maxDistance),
     disjPressConst_(disjPressConst), distanceFieldSrcID_(field::addToStorage< DistanceField_T >(
                                         blockStorage, "BubbleDistance Src", DistanceInfo(), field::fzyx, uint_c(1))),
     distanceFieldDstID_(field::addToStorage< DistanceField_T >(blockStorage, "BubbleDistance Dst", DistanceInfo(),
                                                                field::fzyx, uint_c(1))),
     distanceFieldUpdateInterval_(distanceFieldUpdateInterval)
{
   BubbleModel< Stencil_T >::bubbleFieldCommunication_.addPackInfo(
      std::make_shared< field::communication::PackInfo< DistanceField_T > >(distanceFieldSrcID_));
}

template< typename Stencil_T >
void DisjoiningPressureBubbleModel< Stencil_T >::updateDistanceField()
{
   for (auto blockIt = BubbleModel< Stencil_T >::blockStorage_->begin();
        blockIt != BubbleModel< Stencil_T >::blockStorage_->end(); ++blockIt)
   {
      // get fields
      const BubbleField_T* const bubbleField =
         blockIt->template getData< const BubbleField_T >(BubbleModel< Stencil_T >::getBubbleFieldID());
      DistanceField_T* const distSrcField = blockIt->template getData< DistanceField_T >(distanceFieldSrcID_);
      DistanceField_T* const distDstField = blockIt->template getData< DistanceField_T >(distanceFieldDstID_);

      WALBERLA_FOR_ALL_CELLS(distDstIt, distDstField, distSrcIt, distSrcField, bubbleIt, bubbleField, {
         distDstIt->clear();

         // if current cell is part of a bubble, set the distance to zero
         if (*bubbleIt != INVALID_BUBBLE_ID) { distDstIt->setToZero(*bubbleIt); }

         // loop over src field neighborhood and combine them into this cell
         for (auto d = Stencil_T::beginNoCenter(); d != Stencil_T::end(); ++d)
         {
            // sum distance with already known distances from neighborhood
            distDstIt->combine(distSrcIt.neighbor(*d), d.length(), maxDistance_);
         }
      }); // WALBERLA_FOR_ALL_CELLS

      // swap pointers to distance fields
      distSrcField->swapDataPointers(*distDstField);
   }
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
