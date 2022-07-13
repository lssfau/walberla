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
//! \file BlockStateDetectorSweep.h
//! \ingroup free_surface
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Set block states according to their content, i.e., free surface flags.
//
//======================================================================================================================

#pragma once

#include "core/uid/SUID.h"

#include "FlagInfo.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Set block states according to free surface flag field:
 *  - blocks that contain at least one interface cell are marked as "fullFreeSurface" (computationally most expensive
 *    type of blocks)
 *  - cells without interface cell and with at least one liquid cell are marked "onlyLBM"
 *  - all other blocks are marked as "onlyGasAndBoundary"
 **********************************************************************************************************************/
template< typename FlagField_T >
class BlockStateDetectorSweep
{
 public:
   static const SUID fullFreeSurface;
   static const SUID onlyGasAndBoundary;
   static const SUID onlyLBM;

   BlockStateDetectorSweep(const std::weak_ptr< StructuredBlockForest >& blockForestPtr,
                           FlagInfo< FlagField_T > flagInfo, ConstBlockDataID flagFieldID)
      : flagFieldID_(flagFieldID), flagInfo_(flagInfo)
   {
      const auto blockForest = blockForestPtr.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         (*this)(&(*blockIt)); // initialize block states
      }
   }

   void operator()(IBlock* const block)
   {
      bool liquidFound    = false;
      bool interfaceFound = false;

      const FlagField_T* const flagField = block->getData< const FlagField_T >(flagFieldID_);

      // search the flag field for interface and liquid cells
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flagField, uint_c(1), {
         const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);

         // at inflow boundaries, interface cells are generated; therefore, the block must be fullFreeSurface even if it
         // does not yet contain an interface cell
         if (flagInfo_.isInterface(flagFieldPtr) || flagInfo_.isInflow(flagFieldPtr))
         {
            interfaceFound = true;
            break; // stop search, as this block belongs already to the computationally most expensive block type
         }

         if (flagInfo_.isLiquid(flagFieldPtr)) { liquidFound = true; }
      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

      block->clearState();
      if (interfaceFound)
      {
         block->addState(fullFreeSurface);
         return;
      }
      if (liquidFound)
      {
         block->addState(onlyLBM);
         return;
      }
      block->addState(onlyGasAndBoundary);
   }

 protected:
   ConstBlockDataID flagFieldID_;
   FlagInfo< FlagField_T > flagInfo_;
}; // class BlockStateDetectorSweep

template< typename FlagField_T >
const SUID BlockStateDetectorSweep< FlagField_T >::fullFreeSurface = SUID("fullFreeSurface");
template< typename FlagField_T >
const SUID BlockStateDetectorSweep< FlagField_T >::onlyGasAndBoundary = SUID("onlyGasAndBoundary");
template< typename FlagField_T >
const SUID BlockStateDetectorSweep< FlagField_T >::onlyLBM = SUID("onlyLBM");

} // namespace free_surface
} // namespace walberla
