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
//! \file BubbleModelTester.impl.h
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#include "field/AddToStorage.h"

#include "lbm/free_surface/bubble_model/Geometry.h"

#include "BubbleModelTester.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
template< typename Stencil_T >
BubbleModelTester< Stencil_T >::BubbleModelTester(const std::shared_ptr< StructuredBlockStorage >& blockStorage,
                                                  const std::shared_ptr< BubbleModel< Stencil_T > >& bubbleModel)
   : blockStorage_(blockStorage), bubbleModel_(bubbleModel)
{
   // add flag fields
   srcFlagFieldID_ = field::addFlagFieldToStorage< FlagField_T >(blockStorage_, "FlagFieldSrc");
   dstFlagFieldID_ = field::addFlagFieldToStorage< FlagField_T >(blockStorage_, "FlagFieldDst");

   // add fill level fields
   srcFillLevelFieldID_ =
      field::addToStorage< ScalarField_T >(blockStorage_, "FillLevelSrc", real_c(1.0), field::fzyx, uint_c(1));
   dstFillLevelFieldID_ =
      field::addToStorage< ScalarField_T >(blockStorage_, "FillLevelDst", real_c(1.0), field::fzyx, uint_c(1));

   // register flags
   for (auto blockIt = blockStorage->begin(); blockIt != blockStorage->end(); ++blockIt)
   {
      FlagField_T* const srcFlagField = blockIt->getData< FlagField_T >(srcFlagFieldID_);
      FlagField_T* const dstFlagField = blockIt->getData< FlagField_T >(dstFlagFieldID_);

      liquidFlag_    = srcFlagField->registerFlag("liquid", uint_c(1));
      gasFlag_       = srcFlagField->registerFlag("gas", uint_c(2));
      interfaceFlag_ = srcFlagField->registerFlag("interface", uint_c(3));

      dstFlagField->registerFlag("liquid", uint_c(1));
      dstFlagField->registerFlag("gas", uint_c(2));
      dstFlagField->registerFlag("interface", uint_c(3));
   }
}

template< typename Stencil_T >
void BubbleModelTester< Stencil_T >::srcFlagsFromSrcFills()
{
   // initialize flag field from fill level
   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      FlagField_T* const srcFlagField         = blockIt->getData< FlagField_T >(srcFlagFieldID_);
      const ScalarField_T* const srcFillLevel = blockIt->getData< const ScalarField_T >(srcFillLevelFieldID_);

      bubble_model::setFlagFieldFromFillLevels(srcFlagField, srcFillLevel, "liquid", "gas", "interface");
   }
}

template< typename Stencil_T >
void BubbleModelTester< Stencil_T >::dstFlagsFromDstFills()
{
   // initialize flag field from fill level
   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      FlagField_T* const dstFlagField         = blockIt->getData< FlagField_T >(dstFlagFieldID_);
      const ScalarField_T* const dstFillLevel = blockIt->getData< const ScalarField_T >(dstFillLevelFieldID_);

      bubble_model::setFlagFieldFromFillLevels(dstFlagField, dstFillLevel, "liquid", "gas", "interface");
   }
}

template< typename Stencil_T >
void BubbleModelTester< Stencil_T >::operator()()
{
   // swap source and destination fill level fields, and flag fields
   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      ScalarField_T* const srcFillLevel = blockIt->getData< ScalarField_T >(srcFillLevelFieldID_);
      ScalarField_T* const dstFillLevel = blockIt->getData< ScalarField_T >(dstFillLevelFieldID_);

      FlagField_T* const srcFlagField = blockIt->getData< FlagField_T >(srcFlagFieldID_);
      FlagField_T* const dstFlagField = blockIt->getData< FlagField_T >(dstFlagFieldID_);

      srcFlagField->swapDataPointers(*dstFlagField);
      srcFillLevel->swapDataPointers(*dstFillLevel);
   }

   updateDestinationFields();
   updateBubbleModel();
}

template< typename Stencil_T >
void BubbleModelTester< Stencil_T >::updateBubbleModel()
{
   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      IBlock* thisBlock = &(*blockIt);

      ScalarField_T* const srcFillLevel = blockIt->getData< ScalarField_T >(srcFillLevelFieldID_);
      ScalarField_T* const dstFillLevel = blockIt->getData< ScalarField_T >(dstFillLevelFieldID_);
      WALBERLA_ASSERT(srcFillLevel->hasSameSize(*dstFillLevel));
      WALBERLA_ASSERT_EQUAL_2(srcFillLevel->xyzSize(), dstFillLevel->xyzSize());

      FlagField_T* const srcFlagField = blockIt->getData< FlagField_T >(srcFlagFieldID_);
      FlagField_T* const dstFlagField = blockIt->getData< FlagField_T >(dstFlagFieldID_);
      WALBERLA_ASSERT_EQUAL_2(srcFlagField->xyzSize(), dstFlagField->xyzSize());

      WALBERLA_ASSERT_EQUAL_2(srcFlagField->xyzSize(), srcFillLevel->xyzSize());

      // report cell conversion from liquid to interface; explicitly avoid OpenMP when setting bubble IDs
      WALBERLA_FOR_ALL_CELLS_OMP(srcFlagFieldIt, srcFlagField, dstFlagFieldIt, dstFlagField, omp critical, {
         if (isFlagSet(srcFlagFieldIt, liquidFlag_) && isFlagSet(dstFlagFieldIt, interfaceFlag_))
         {
            bubbleModel_->reportLiquidToInterfaceConversion(thisBlock, srcFlagFieldIt.cell());
         }
      }) // WALBERLA_FOR_ALL_CELLS_OMP

      // report changes in the fill level
      WALBERLA_FOR_ALL_CELLS_OMP(
         srcFlagFieldIt, srcFlagField, dstFlagFieldIt, dstFlagField, srcFillIt, srcFillLevel, dstFillIt, dstFillLevel,
         omp critical, {
            if (isFlagSet(srcFlagFieldIt, interfaceFlag_) || isFlagSet(dstFlagFieldIt, interfaceFlag_))
            {
               bubbleModel_->reportFillLevelChange(thisBlock, srcFillIt.cell(), *dstFillIt - *srcFillIt);
            }
         }) // WALBERLA_FOR_ALL_CELLS_OMP

      // report cell conversion from interface to liquid
      WALBERLA_FOR_ALL_CELLS_OMP(srcFlagFieldIt, srcFlagField, dstFlagFieldIt, dstFlagField, omp critical, {
         if (isFlagSet(srcFlagFieldIt, interfaceFlag_) && isFlagSet(dstFlagFieldIt, liquidFlag_))
         {
            bubbleModel_->reportInterfaceToLiquidConversion(thisBlock, srcFlagFieldIt.cell());
         }
      }) // WALBERLA_FOR_ALL_CELLS_OMP
   }
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
