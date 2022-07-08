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
//! \file FloodFill.impl.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Flood fill algorithm to identify connected gas volumes as bubbles.
//
//======================================================================================================================

#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/InterfaceFromFillLevel.h"

#include "FloodFill.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
// Sets newBubbleID in all cells that belong to the same bubble as startCell. Only uses the fillField but needs to
// call isInterfaceFromFillLevel() to identify interface cells. This is expensive and requires information from
// neighboring cells.
template< typename Stencil_T >
void FloodFillUsingFillLevel< Stencil_T >::run(IBlock& block, BlockDataID bubbleFieldID, const Cell& startCell,
                                               BubbleID newBubbleID, real_t& volume, uint_t& nrOfCells)
{
   // get fields
   BubbleField_T* const bubbleField     = block.getData< BubbleField_T >(bubbleFieldID);
   const ScalarField_T* const fillField = block.getData< const ScalarField_T >(fillFieldID_);

   WALBERLA_ASSERT_EQUAL(fillField->xyzSize(), bubbleField->xyzSize());

   volume    = real_c(0);
   nrOfCells = uint_c(0);

   std::queue< Cell > cellQueue;
   cellQueue.push(startCell);

   std::vector< bool > addedLast;

   using namespace stencil;
   const int dirs[4]    = { N, S, T, B };
   const uint_t numDirs = uint_c(4);

   CellInterval fieldSizeInterval = fillField->xyzSize();

   while (!cellQueue.empty())
   {
      // process first cell in queue
      Cell cell = cellQueue.front();

      // go to the beginning of the x-line, i.e., find the minimum x-coordinate that still belongs to this bubble
      while (cell.x() > cell_idx_c(0) &&
             bubbleField->get(cell.x() - cell_idx_c(1), cell.y(), cell.z()) != newBubbleID &&
             (fillField->get(cell.x() - cell_idx_c(1), cell.y(), cell.z()) < real_c(1.0) ||
              isInterfaceFromFillLevel< Stencil_T >(*fillField, cell.x() - cell_idx_c(1), cell.y(), cell.z())))
      {
         --cell.x();
      }

      // vector stores whether a cell in a direction from dirs was added already to make sure that for the whole x-line,
      // only one neighboring cell per direction is added to the cellQueue
      addedLast.assign(numDirs, false);

      // loop to the end of the x-line and mark any cell that belongs to this bubble
      while (cell.x() < cell_idx_c(fillField->xSize()) && bubbleField->get(cell) != newBubbleID &&
             (fillField->get(cell) < real_c(1) || isInterfaceFromFillLevel< Stencil_T >(*fillField, cell)))
      {
         // set bubble ID to this cell
         bubbleField->get(cell) = newBubbleID;
         volume += real_c(1.0) - fillField->get(cell);
         nrOfCells++;

         // iterate over all directions in dirs to check which cells in y- and z-direction should be processed next,
         // i.e., which cells should be added to cellQueue
         for (uint_t i = uint_c(0); i < numDirs; ++i)
         {
            Cell neighborCell(cell.x(), cell.y() + cy[dirs[i]], cell.z() + cz[dirs[i]]);

            // neighboring cell is not inside the field
            if (!fieldSizeInterval.contains(neighborCell)) { continue; }

            // neighboring cell is part of the bubble
            if ((fillField->get(neighborCell) < real_c(1) ||
                 isInterfaceFromFillLevel< Stencil_T >(*fillField, neighborCell)) &&
                bubbleField->get(neighborCell) != newBubbleID)
            {
               // make sure that for the whole x-line, only one neighboring cell per direction is added to the cellQueue
               if (!addedLast[i])
               {
                  addedLast[i] = true;

                  // add neighboring cell to queue such that it serves as starting cell in one of the next iterations
                  cellQueue.push(neighborCell);
               }
            }
            else { addedLast[i] = false; }
         }

         ++cell.x();
      }

      // remove the processed cell from cellQueue
      cellQueue.pop();
   }
}

// Sets newBubbleID in all cells that belong to the same bubble as startCell. Uses the fillField and the flagField
// and is therefore less expensive than the approach that only uses the fillField.
template< typename FlagField_T >
void FloodFillUsingFlagField< FlagField_T >::run(IBlock& block, BlockDataID bubbleFieldID, const Cell& startCell,
                                                 BubbleID newBubbleID, real_t& volume, uint_t& nrOfCells)
{
   // get fields
   BubbleField_T* const bubbleField     = block.getData< BubbleField_T >(bubbleFieldID);
   const ScalarField_T* const fillField = block.getData< const ScalarField_T >(fillFieldID_);
   const FlagField_T* const flagField   = block.getData< const FlagField_T >(flagFieldID_);

   WALBERLA_ASSERT_EQUAL(fillField->xyzSize(), bubbleField->xyzSize());
   WALBERLA_ASSERT_EQUAL(fillField->xyzSize(), flagField->xyzSize());

   using flag_t = typename FlagField_T::flag_t;
   // gasInterfaceFlagMask masks interface and gas cells (using bitwise OR)
   flag_t gasInterfaceFlagMask = flag_t(0);
   gasInterfaceFlagMask        = flag_t(gasInterfaceFlagMask | flagField->getFlag(flagIDs::interfaceFlagID));
   gasInterfaceFlagMask        = flag_t(gasInterfaceFlagMask | flagField->getFlag(flagIDs::gasFlagID));

   volume    = real_c(0);
   nrOfCells = uint_c(0);

   std::queue< Cell > cellQueue;
   cellQueue.push(startCell);

   std::vector< bool > addedLast;

   using namespace stencil;
   const int dirs[4]    = { N, S, T, B };
   const uint_t numDirs = uint_c(4);

   CellInterval fieldSizeInterval = flagField->xyzSize();

   while (!cellQueue.empty())
   {
      // process first cell in queue
      Cell& cell = cellQueue.front();

      // go to the beginning of the x-line, i.e., find the minimum x-coordinate that still belongs to this bubble
      while (isPartOfMaskSet(flagField->get(cell.x() - cell_idx_c(1), cell.y(), cell.z()), gasInterfaceFlagMask) &&
             cell.x() > cell_idx_c(0) && bubbleField->get(cell.x() - cell_idx_c(1), cell.y(), cell.z()) != newBubbleID)
      {
         --cell.x();
      }

      // vector stores whether a cell in a direction from dirs was added already to make sure that for the whole x-line,
      // only one neighboring cell per direction is added to the cellQueue
      addedLast.assign(numDirs, false);

      // loop to the end of the x-line and mark any cell that belongs to this bubble
      while (isPartOfMaskSet(flagField->get(cell), gasInterfaceFlagMask) && bubbleField->get(cell) != newBubbleID &&
             cell.x() < cell_idx_c(flagField->xSize()))
      {
         // set bubble ID to this cell
         bubbleField->get(cell) = newBubbleID;
         volume += real_c(1.0) - fillField->get(cell);
         nrOfCells++;

         // iterate over all directions in dirs to check which cells in y- and z-direction should be processed next,
         // i.e., which cells should be added to cellQueue
         for (uint_t i = uint_c(0); i < numDirs; ++i)
         {
            Cell neighborCell(cell.x(), cell.y() + cy[dirs[i]], cell.z() + cz[dirs[i]]);

            // neighboring cell is not inside the field
            if (!fieldSizeInterval.contains(neighborCell)) { continue; }

            // neighboring cell is part of the bubble
            if (!isPartOfMaskSet(flagField->get(neighborCell), gasInterfaceFlagMask) &&
                bubbleField->get(neighborCell) != newBubbleID)
            {
               // make sure that for the whole x-line, only one neighboring cell per direction is added to the cellQueue
               if (!addedLast[i])
               {
                  addedLast[i] = true;

                  // add neighboring cell to queue such that it serves as starting cell in one of the next iterations
                  cellQueue.push(neighborCell);
               }
            }
            else { addedLast[i] = false; }
         }
         ++cell.x();
      }

      // remove the processed cell from cellQueue
      cellQueue.pop();
   }
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
