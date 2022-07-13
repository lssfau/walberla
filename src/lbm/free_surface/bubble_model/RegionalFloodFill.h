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
//! \file RegionalFloodFill.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Flood fill algorithm in a constrained neighborhood around a cell.
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
 *
 * Flood fill algorithm in a constrained neighborhood around a cell.
 *
 * This algorithm is used to check if a bubble has split up. A split detection is triggered
 * when an interface cell becomes a fluid cell. The normal case is to check in a 1-neighborhood if
 * the bubble has split up: ("-" means no bubble, numbers indicate bubble IDs)
 * Example:  0 0 -
 *           - x -
 *           0 0 -
 * In this case, the bubble has split up.
 *
 * This check in a 1-neighborhood is performed by the function BubbleModel::mapNeighborhood().
 *
 * To look only in 1-neighborhoods has the drawback that we trigger splits when no split has occurred. This is OK, since
 * later it is detected that the bubble has not really split. However, this process requires communication and is
 * therefore expensive.
 *
 * The RegionalFloodFill has the task to rule out splits before they are triggered by looking at larger neighborhoods.
 *
 * Example: (same as above but with larger neighborhood)
 *          0 0 0 - -
 *          0 0 0 - -
 *          0 - x - -
 *          0 0 0 - -
 *          0 0 0 - -
 * By looking at the 2-neighborhood we see that the bubble has not really split up.
 * The RegionalFloodFill needs as input the current cell (marked with x), a direction to start, e.g., N, and the size of
 * the neighborhood to look at, e.g., 2. In this case, it creates a small 5x5x5 field where a flood fill is executed.
 * The field might be smaller if the cell is near the boundary, as the algorithm is capable of detecting boundaries.
 *
 **********************************************************************************************************************/
template< typename T, typename Stencil_T >
class RegionalFloodFill
{
 public:
   RegionalFloodFill(const GhostLayerField< T, 1 >* externField, const Cell& startCell,
                     stencil::Direction startDirection, const T& searchValue, cell_idx_t neighborhood = 2);

   RegionalFloodFill(const GhostLayerField< T, 1 >* externField, const Cell& startCell, const T& emptyValue,
                     cell_idx_t neighborhood = 2);

   ~RegionalFloodFill() { delete workingField_; }

   inline bool connected(stencil::Direction d);
   inline bool connected(cell_idx_t xOff, cell_idx_t yOff, cell_idx_t zOff);

   const Field< bool, 1 >& workingField() const { return *workingField_; }

 protected:
   inline bool cellInsideAndNotMarked(const Cell& cell);

   void runFloodFill(const Cell& startCell, stencil::Direction startDirection, const T& searchValue,
                     cell_idx_t neighborhood);

   const GhostLayerField< T, 1 >* externField_;
   Field< bool, 1 >* workingField_;

   // the size of the working field is 2*neighborhood+1 for storing startCell (+1) and neighborhood in each direction;
   // boundaries in externField_ reduce the size of workingField_ accordingly
   CellInterval workingFieldSize_;

   T searchValue_;
   Cell offset_; // offset of workingField_, allows coordinate transformations:
                 // externFieldCoordinates = workingFieldCoordinates + offset
   Cell startCellInWorkingField_;
}; // class RegionalFloodFill

// iterate over all neighboring cells and run flood fill with starting direction of first neighboring cell whose value
// is not equal to emptyValue
template< typename T, typename Stencil_T >
RegionalFloodFill< T, Stencil_T >::RegionalFloodFill(const GhostLayerField< T, 1 >* externField, const Cell& startCell,
                                                     const T& emptyValue, cell_idx_t neighborhood)
   : externField_(externField), workingField_(nullptr)
{
   for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
   {
      const T& neighborValue = externField_->getNeighbor(startCell, *d);
      if (neighborValue != emptyValue)
      {
         // run flood fill with starting direction of first non-empty neighboring cell
         runFloodFill(startCell, *d, neighborValue, neighborhood);
         return;
      }
   }

   WALBERLA_ASSERT(false); // should not happen, only empty values in neighborhood
}

/***********************************************************************************************************************
 * Create a small overlay field in given neighborhood and execute a flood fill in this small field.
 *
 * The flood fill starts in startCell and startDirection. The neighboring cell of startCell in startDirection has to be
 * set to searchValue. The size of the neighborhood is specified by neighborhood. In case of nearby boundary in a
 * certain direction, the size of the neighborhood is automatically reduced in this direction.
 **********************************************************************************************************************/
template< typename T, typename Stencil_T >
RegionalFloodFill< T, Stencil_T >::RegionalFloodFill(const GhostLayerField< T, 1 >* externField, const Cell& startCell,
                                                     stencil::Direction startDirection, const T& searchValue,
                                                     cell_idx_t neighborhood)
   : externField_(externField), workingField_(nullptr)
{
   runFloodFill(startCell, startDirection, searchValue, neighborhood);
}

template< typename T, typename Stencil_T >
void RegionalFloodFill< T, Stencil_T >::runFloodFill(const Cell& startCell, stencil::Direction startDirection,
                                                     const T& searchValue, cell_idx_t neighborhood)
{
   searchValue_ = searchValue;

   const cell_idx_t externFieldGhostLayers = cell_idx_c(externField_->nrOfGhostLayers());

   // size of workingField_ is defined by the number of cells around startCell
   cell_idx_t extendLower[3]; // number of cells in negative x-, y-, and z-direction of startCell
   cell_idx_t extendUpper[3]; // number of cells in negative x-, y-, and z-direction of startCell

   // determine the size of the working field with respect to boundaries of externField_
   for (uint_t i = uint_c(0); i < uint_c(3); ++i)
   {
      const cell_idx_t minCoord = -externFieldGhostLayers;
      const cell_idx_t maxCoord = cell_idx_c(externField_->size(i)) - cell_idx_c(1) + externFieldGhostLayers;

      // in case of nearby boundaries in externField_, the size of workingField_ is adjusted such that these boundaries
      // are respected
      extendLower[i] = std::min(neighborhood, startCell[i] - minCoord);
      extendUpper[i] = std::min(neighborhood, maxCoord - startCell[i]);
   }

   // offset_ can be used to transform coordinates from externField_ to coordinates of workingField_
   offset_ = Cell(startCell[0] - extendLower[0], startCell[1] - extendLower[1], startCell[2] - extendLower[2]);

   startCellInWorkingField_ = startCell - offset_;

   // create workingField_
   workingField_ = new Field< bool, 1 >(uint_c(1) + uint_c(extendLower[0] + extendUpper[0]),
                                        uint_c(1) + uint_c(extendLower[1] + extendUpper[1]),
                                        uint_c(1) + uint_c(extendLower[2] + extendUpper[2]), false);

   workingFieldSize_ = workingField_->xyzSize();

   // mark the startCell such that flood fill can not take a path across startCell
   workingField_->get(startCellInWorkingField_) = true;

   // use stack to store cells that still need to be searched
   std::vector< Cell > stack;
   stack.reserve(Stencil_T::Size);
   stack.emplace_back(startCellInWorkingField_[0] + stencil::cx[startDirection],
                      startCellInWorkingField_[1] + stencil::cy[startDirection],
                      startCellInWorkingField_[2] + stencil::cz[startDirection]);

   while (!stack.empty())
   {
      // next search cell is the last entry in stack
      Cell currentCell = stack.back();

      // remove the searched cell from stack
      stack.pop_back();

      WALBERLA_ASSERT_EQUAL(externField_->get(currentCell + offset_), searchValue_);

      // mark the searched cell to be connected in workingField_
      workingField_->get(currentCell) = true;

      // iterate over all existing neighbors that are not yet marked and push them onto the stack
      for (auto d = Stencil_T::beginNoCenter(); d != Stencil_T::end(); ++d)
      {
         Cell neighbor(currentCell[0] + d.cx(), currentCell[1] + d.cy(), currentCell[2] + d.cz());

         // check if neighboring cell is:
         // - inside workingField_
         // - not yet marked as connected
         // - connected to this cell
         if (cellInsideAndNotMarked(neighbor))
         {
            // add neighbor to stack such that it gets marked as connected and used as start cell in the next iteration;
            // cells that are not connected to startCell will never get marked as connected
            stack.push_back(neighbor);
         }
      }
   }
}

// check if startCell is connected to the neighboring cell in direction d
template< typename T, typename Stencil_T >
bool RegionalFloodFill< T, Stencil_T >::connected(stencil::Direction d)
{
   using namespace stencil;
   return workingField_->get(startCellInWorkingField_[0] + cx[d], startCellInWorkingField_[1] + cy[d],
                             startCellInWorkingField_[2] + cz[d]);
}

// check if startCell is connected to the cell with some offset from startCell
template< typename T, typename Stencil_T >
bool RegionalFloodFill< T, Stencil_T >::connected(cell_idx_t xOff, cell_idx_t yOff, cell_idx_t zOff)
{
   return workingField_->get(startCellInWorkingField_[0] + xOff, startCellInWorkingField_[1] + yOff,
                             startCellInWorkingField_[2] + zOff);
}

// check if cell is:
// - inside workingField_
// - not yet marked as connected
// - connected to this cell
template< typename T, typename Stencil_T >
bool RegionalFloodFill< T, Stencil_T >::cellInsideAndNotMarked(const Cell& cell)
{
   if (!workingFieldSize_.contains(cell)) { return false; }

   // check if cell is already marked as connected
   if (!workingField_->get(cell))
   {
      // test if cell is connected
      return (externField_->get(cell + offset_) == searchValue_);
   }
   else { return false; }
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
