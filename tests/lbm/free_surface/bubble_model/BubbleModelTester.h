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
//! \file BubbleModelTester.h
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Helper class for bubble model test cases.
//!
//! Provides a source and destination field for fill levels and flags. A derived class is supposed to work on the
//! destination fields. This (base) class then updates the bubble model and converts cells with respect to any changes
//! in the destination fields.
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "field/FlagField.h"
#include "field/GhostLayerField.h"

#include "lbm/free_surface/bubble_model/BubbleModel.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
template< typename Stencil_T >
class BubbleModelTester
{
 public:
   // define fields
   using flag_t        = uint32_t;
   using FlagField_T   = FlagField< flag_t >;
   using ScalarField_T = GhostLayerField< real_t, 1 >;

   BubbleModelTester(const std::shared_ptr< StructuredBlockStorage >& blockStorage,
                     const std::shared_ptr< BubbleModel< Stencil_T > >& bubbleModel);
   virtual ~BubbleModelTester() = default;

   // swap source and destination fill level fields;
   // swap source and destination flag fields;
   // update bubble model;
   // update destination fields;
   virtual void operator()();

   inline ConstBlockDataID getFlagFieldID() const { return srcFlagFieldID_; }
   inline ConstBlockDataID getFillLevelFieldID() const { return srcFillLevelFieldID_; }

 private:
   // report cell conversions from liquid to interface;
   // report changes in the fill level (destination fill level - source fill level);
   // report cell conversions from interface to liquid;
   // check interface layer for correctness
   void updateBubbleModel();

 protected:
   virtual void updateDestinationFields() = 0;

   // initialize flag field from fill level
   void srcFlagsFromSrcFills();
   void dstFlagsFromDstFills();

   std::shared_ptr< StructuredBlockStorage > blockStorage_;
   std::shared_ptr< BubbleModel< Stencil_T > > bubbleModel_;

   uint32_t liquidFlag_;
   uint32_t gasFlag_;
   uint32_t interfaceFlag_;

   BlockDataID srcFlagFieldID_;
   BlockDataID dstFlagFieldID_;
   BlockDataID srcFillLevelFieldID_;
   BlockDataID dstFillLevelFieldID_;
};

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

#include "BubbleModelTester.impl.h"