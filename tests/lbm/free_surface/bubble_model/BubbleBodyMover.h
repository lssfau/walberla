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
//! \file BubbleBodyMover.h
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Helper class for bubble model test cases.
//!
//! Add bubbles that are represented by geometric bodies. Move the bubbles by updating the fill level and flag field.
//
//======================================================================================================================

#pragma once

#include "geometry/initializer/OverlapFieldFromBody.h"

#include "BubbleModelTester.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
template< typename Body_T, typename Stencil_T >
class BubbleBodyMover : public BubbleModelTester< Stencil_T >
{
 public:
   BubbleBodyMover(const std::shared_ptr< StructuredBlockStorage >& blockStorage,
                   const std::shared_ptr< BubbleModel< Stencil_T > >& bubbleModel);

   inline const std::vector< Body_T >& getBodies() const { return bodies_; }

   void addBody(const Body_T& b, std::function< void(Body_T&) > moveFunction)
   {
      bodies_.push_back(b);
      moveFunctions_.push_back(moveFunction);
   }

   // initialize the added bodies in the fill level and flag fields, and in the bubble model
   void initAddedBodies();

 protected:
   void updateDestinationFields() override;

   std::vector< Body_T > bodies_;

   std::vector< std::function< void(Body_T&) > > moveFunctions_;

   geometry::initializer::OverlapFieldFromBody srcFillInitializer_;
   geometry::initializer::OverlapFieldFromBody dstFillInitializer_;
};

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

#include "BubbleBodyMover.impl.h"