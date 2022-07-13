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
//! \file BubbleBodyMover.impl.h
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#include "geometry/bodies/AABBBody.h"
#include "geometry/bodies/Ellipsoid.h"
#include "geometry/bodies/Sphere.h"

#include "lbm/free_surface/bubble_model/Geometry.h"

#include "BubbleBodyMover.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
template< typename Body_T, typename Stencil_T >
BubbleBodyMover< Body_T, Stencil_T >::BubbleBodyMover(const std::shared_ptr< StructuredBlockStorage >& blockStorage,
                                                      const std::shared_ptr< BubbleModel< Stencil_T > >& bubbleModel)
   : BubbleModelTester< Stencil_T >(blockStorage, bubbleModel),
     srcFillInitializer_(*blockStorage, BubbleModelTester< Stencil_T >::srcFillLevelFieldID_),
     dstFillInitializer_(*blockStorage, BubbleModelTester< Stencil_T >::dstFillLevelFieldID_)
{}

template< typename Body_T, typename Stencil_T >
void BubbleBodyMover< Body_T, Stencil_T >::initAddedBodies()
{
   // initialize fill levels
   for (auto body = bodies_.begin(); body != bodies_.end(); ++body)
   {
      srcFillInitializer_.init(*body, false);
      dstFillInitializer_.init(*body, false);
   }

   // set flags
   BubbleModelTester< Stencil_T >::srcFlagsFromSrcFills();
   BubbleModelTester< Stencil_T >::dstFlagsFromDstFills();

   // initialize bubble model
   BubbleModelTester< Stencil_T >::bubbleModel_->initFromFillLevelField(
      BubbleModelTester< Stencil_T >::dstFillLevelFieldID_);
}

template< typename Body_T, typename Stencil_T >
void BubbleBodyMover< Body_T, Stencil_T >::updateDestinationFields()
{
   using BubbleModelTester_T = BubbleModelTester< Stencil_T >;

   // move bodies according to moveFunctions_
   for (uint_t i = uint_c(0); i < bodies_.size(); ++i)
   {
      moveFunctions_[i](bodies_[i]);
   }

   // clear destination field
   for (auto blockIt = this->blockStorage_->begin(); blockIt != this->blockStorage_->end(); ++blockIt)
   {
      typename BubbleModelTester_T::FlagField_T* const flagField =
         blockIt->template getData< typename BubbleModelTester_T::FlagField_T >(this->dstFlagFieldID_);
      typename BubbleModelTester_T::ScalarField_T* const fillField =
         blockIt->template getData< typename BubbleModelTester_T::ScalarField_T >(this->dstFillLevelFieldID_);

      flagField->setWithGhostLayer(typename BubbleModelTester_T::flag_t(0));
      fillField->setWithGhostLayer(real_c(1.0));
   }

   // update fill level field with new body position
   for (auto body = bodies_.begin(); body != bodies_.end(); ++body)
   {
      dstFillInitializer_.init(*body, false);
   }

   // update flags
   this->dstFlagsFromDstFills();
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
