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
//! \file Geometry.impl.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Geometrical helper functions for the bubble model.
//
//======================================================================================================================

#include "core/logging/Logging.h"

#include "field/FlagField.h"
#include "field/GhostLayerField.h"

#include <limits>

#include "Geometry.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
template< typename Body_T >
void addBodyToFillLevelField(StructuredBlockStorage& blockStorage, BlockDataID fillFieldID, const Body_T& body,
                             bool isGas)
{
   const real_t dx = blockStorage.dx();

   for (auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt)
   {
      // get block
      IBlock* block = &(*blockIt);

      // get fill level field
      GhostLayerField< real_t, 1 >* const fillField = block->getData< GhostLayerField< real_t, 1 > >(fillFieldID);

      // get the block's bounding box
      AABB blockBB = block->getAABB();

      // extend the bounding box with the ghost layer
      blockBB.extend(dx * real_c(fillField->nrOfGhostLayers()));

      // skip blocks that do not intersect the body
      if (geometry::fastOverlapCheck(body, blockBB) == geometry::COMPLETELY_OUTSIDE) { continue; }

      // get the global coordinates (via the bounding box) of the block's first cell
      AABB firstCellBB;
      blockStorage.getBlockLocalCellAABB(*block, fillField->beginWithGhostLayer().cell(), firstCellBB);

      // get the global coordinates of the midpoint of the block's first cell
      Vector3< real_t > firstCellMidpoint;
      for (uint_t i = uint_c(0); i < uint_c(3); ++i)
      {
         firstCellMidpoint[i] = firstCellBB.min(i) + real_c(0.5) * firstCellBB.size(i);
      }

      // get the number of ghost layers
      const uint_t numGl     = fillField->nrOfGhostLayers();
      cell_idx_t glCellIndex = cell_idx_c(numGl);

      // starting from a block's first cell, iterate over all cells and determine the overlap of the body with each cell
      // to set the cell's fill level accordingly
      Vector3< real_t > currentMidpoint;
      currentMidpoint[2] = firstCellMidpoint[2];
      for (cell_idx_t z = -glCellIndex; z < cell_idx_c(fillField->zSize() + numGl); ++z, currentMidpoint[2] += dx)
      {
         currentMidpoint[1] = firstCellMidpoint[1];
         for (cell_idx_t y = -glCellIndex; y < cell_idx_c(fillField->ySize() + numGl); ++y, currentMidpoint[1] += dx)
         {
            currentMidpoint[0] = firstCellMidpoint[0];
            for (cell_idx_t x = -glCellIndex; x < cell_idx_c(fillField->xSize() + numGl); ++x, currentMidpoint[0] += dx)
            {
               // get the current cell's overlap with the body (full overlap=1, no overlap=0, else in between)
               real_t overlapFraction = geometry::overlapFraction(body, currentMidpoint, dx);

               // get the current fill level
               real_t& fillLevel = fillField->get(x, y, z);
               WALBERLA_ASSERT(fillLevel >= 0 && fillLevel <= 1);

               // initialize a gas bubble (fill level=0)
               if (isGas)
               {
                  // subtract the fraction of the body's overlap from the fill level
                  fillLevel -= overlapFraction;

                  // limit the fill level such that it does not become negative during initialization
                  fillLevel = std::max(fillLevel, real_c(0.0));
               }
               else // initialize a liquid drop (fill level=1)
               {
                  // add the fraction of the body's overlap from the fill level
                  fillLevel += overlapFraction;

                  // limit the fill level such that it does not become too large during initialization
                  fillLevel = std::min(fillLevel, real_c(1.0));
               }
            }
         }
      }
   }
}

template< typename flag_t >
void setFlagFieldFromFillLevels(FlagField< flag_t >* flagField, const GhostLayerField< real_t, 1 >* fillField,
                                const FlagUID& liquid, const FlagUID& gas, const FlagUID& interFace)
{
   WALBERLA_ASSERT(flagField->xyzSize() == fillField->xyzSize());

   // get flags from flag field
   flag_t liquidFlag    = flagField->getFlag(liquid);
   flag_t gasFlag       = flagField->getFlag(gas);
   flag_t interfaceFlag = flagField->getFlag(interFace);

   flag_t allMask = liquidFlag | gasFlag | interfaceFlag;

   using FillField_T = GhostLayerField< real_t, 1 >;

   // iterate over all cells (including ghost layer) in fill level field and flag field
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(fillField, {
      const typename FillField_T::ConstPtr fillFieldPtr(*fillField, x, y, z);
      const typename FlagField< flag_t >::Ptr flagFieldPtr(*flagField, x, y, z);

      // clear liquid, gas, and interface flags in flag field
      removeMask(flagFieldPtr, allMask);

      if (*fillFieldPtr <= real_c(0)) { addFlag(flagFieldPtr, gasFlag); }
      else
      {
         if (*fillFieldPtr >= real_c(1)) { addFlag(flagFieldPtr, liquidFlag); }
         else
         {
            // add interface flag for fill levels between 0 and 1
            addFlag(flagFieldPtr, interfaceFlag);
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}

template< typename flag_t, typename Stencil_T >
bool checkForValidInterfaceLayer(const FlagField< flag_t >* flagField, const FlagUID& liquid, const FlagUID& gas,
                                 const FlagUID& interFace, bool printWarnings)
{
   // variable only used when printWarnings is true; avoids premature termination of the search such that any non-valid
   // cell can be print
   bool valid = true;

   // get flags
   flag_t liquidFlag    = flagField->getFlag(liquid);
   flag_t gasFlag       = flagField->getFlag(gas);
   flag_t interfaceFlag = flagField->getFlag(interFace);

   // iterate flag field
   for (auto flagFieldIt = flagField->begin(); flagFieldIt != flagField->end(); ++flagFieldIt)
   {
      // check that not more than one flag is set for a cell
      if (isMaskSet(flagFieldIt, flag_t(liquidFlag | gasFlag)) ||
          isMaskSet(flagFieldIt, flag_t(liquidFlag | interfaceFlag)) ||
          isMaskSet(flagFieldIt, flag_t(gasFlag | interfaceFlag)))
      {
         if (!printWarnings) { return false; }
         valid = false;
         WALBERLA_LOG_WARNING("More than one free surface flag is set in cell " << flagFieldIt.cell() << ".");
      }

      // if current cell is a gas cell it must not have a fluid cell in its neighborhood
      if (isFlagSet(flagFieldIt, liquidFlag))
      {
         for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
         {
            if (isFlagSet(flagFieldIt.neighbor(*dir), gasFlag))
            {
               if (!printWarnings) { return false; }
               valid = false;
               WALBERLA_LOG_WARNING("Fluid cell " << flagFieldIt.cell()
                                                  << " has a gas cell in its direct neighborhood.");
            }
         }
      }
      // if current cell is a fluid cell it must not have a gas cell in its neighborhood
      else
      {
         if (isFlagSet(flagFieldIt, gasFlag))
         {
            for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
            {
               if (isFlagSet(flagFieldIt.neighbor(*dir), liquidFlag))
               {
                  if (!printWarnings) { return false; }
                  valid = false;
                  WALBERLA_LOG_WARNING("Gas cell " << flagFieldIt.cell()
                                                   << " has a fluid cell in its direct neighborhood.");
               }
            }
         }
      }
   }

   return valid;
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla