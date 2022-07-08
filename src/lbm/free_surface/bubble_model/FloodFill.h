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
//! \file FloodFill.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Flood fill algorithm to identify connected gas volumes as bubbles.
//
//======================================================================================================================

#pragma once

#include "core/cell/Cell.h"

#include <queue>

#include "BubbleDefinitions.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
 * Use the flood fill (also called seed fill) algorithm to identify connected gas volumes as bubbles and mark the whole
 * region belonging to the bubble with a common bubble ID in bubbleField.
 **********************************************************************************************************************/
class FloodFillInterface
{
 public:
   virtual ~FloodFillInterface() = default;

   /********************************************************************************************************************
    * Marks the region of a bubble in the BubbleField
    *
    * bubbleField   Field where the bubble is marked. BubbleField must not contain entries equal to newBubbleID.
    * startCell     Cell that belongs to bubble. It is used as starting point for flood fill algorithm.
    * newBubbleID   An integer, greater than zero. This value is used as the bubble marker (ID) in the bubble field.
    * volume        Volume of the gas phase, i.e., the sum of (1-fillLevel) for all cells inside the bubble.
    * nrOfCells     The number of cells that belong to this bubble.
    ********************************************************************************************************************/
   virtual void run(IBlock& block, BlockDataID bubbleIDField, const Cell& startCell, BubbleID newBubbleID,
                    real_t& volume, uint_t& nrOfCells) = 0;
}; // class FloodFillInterface

/***********************************************************************************************************************
 * Use only fill level to identify bubbles.
 * Problem: isInterfaceFromFillLevel() has to be called to identify interface cells; this is expensive and requires
 * information from neighboring cells.
 **********************************************************************************************************************/
template< typename Stencil_T >
class FloodFillUsingFillLevel : public FloodFillInterface
{
 public:
   FloodFillUsingFillLevel(ConstBlockDataID fillFieldID) : fillFieldID_(fillFieldID) {}
   ~FloodFillUsingFillLevel() override = default;

   void run(IBlock& block, BlockDataID bubbleFieldID, const Cell& startCell, BubbleID newBubbleID, real_t& volume,
            uint_t& nrOfCells) override;

 private:
   ConstBlockDataID fillFieldID_;
}; // FloodFillUsingFillLevel

/***********************************************************************************************************************
 * Use flag field to identify interface cells which should be faster than using only the fill level.
 **********************************************************************************************************************/
template< typename FlagField_T >
class FloodFillUsingFlagField : public FloodFillInterface
{
 public:
   FloodFillUsingFlagField(ConstBlockDataID fillFieldID, ConstBlockDataID flagFieldID)
      : fillFieldID_(fillFieldID), flagFieldID_(flagFieldID)
   {}

   ~FloodFillUsingFlagField() override = default;

   void run(IBlock& block, BlockDataID bubbleFieldID, const Cell& startCell, BubbleID newBubbleID, real_t& volume,
            uint_t& nrOfCells) override;

 private:
   ConstBlockDataID fillFieldID_;
   ConstBlockDataID flagFieldID_;
}; // FloodFillUsingFlagField

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

#include "FloodFill.impl.h"