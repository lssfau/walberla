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
//! \file BubbleModel.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief System for tracking pressure/density in gas volumes.
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/math/Vector3.h"

#include "field/GhostLayerField.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q27.h"

#include "Bubble.h"
#include "BubbleDefinitions.h"
#include "FloodFill.h"
#include "MergeInformation.h"
#include "NewBubbleCommunication.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
class BubbleModelBase
{
 public:
   virtual ~BubbleModelBase() = default;

   virtual real_t getDensity(IBlock* block, const Cell& cell) const       = 0;
   virtual void setDensity(IBlock* block, const Cell& cell, real_t value) = 0;
   virtual void setDensityOfAllBubbles(real_t val)                        = 0;

   // call updateVolumeDiff() with fillLevelDifference
   virtual void reportFillLevelChange(IBlock* block, const Cell& cell, real_t fillLevelDifference) = 0;
   virtual void reportLiquidToInterfaceConversion(IBlock* block, const Cell& cell)                 = 0;
   virtual void reportInterfaceToLiquidConversion(IBlock* block, const Cell& cell)                 = 0;

   virtual void update() = 0;
}; // class BubbleModelBase

/***********************************************************************************************************************
 * Implementation for setups in which no bubble model is required. Reports always a constant pressure
 **********************************************************************************************************************/
class BubbleModelConstantPressure : public BubbleModelBase
{
 public:
   BubbleModelConstantPressure(real_t constantLatticeDensity) : constantLatticeDensity_(constantLatticeDensity) {}
   ~BubbleModelConstantPressure() override = default;

   real_t getDensity(IBlock*, const Cell&) const override { return constantLatticeDensity_; }
   void setDensity(IBlock*, const Cell&, real_t) override {}
   void setDensityOfAllBubbles(real_t val) override { constantLatticeDensity_ = val; }

   void reportFillLevelChange(IBlock*, const Cell&, real_t) override{};
   void reportLiquidToInterfaceConversion(IBlock*, const Cell&) override{};
   void reportInterfaceToLiquidConversion(IBlock*, const Cell&) override{};

   void update() override {}

 private:
   real_t constantLatticeDensity_;
}; // class BubbleModelConstantPressure

/***********************************************************************************************************************
 * System for tracking pressure/density in gas volumes.
 *
 * The pure volume of fluid code calculates how mass/fluid moves across the domain. As input it needs the gas density or
 * pressure. The density is equal in the same bubble. To track the pressure, individual gas volumes have to be tracked.
 * The bubbles can split or merge, can possibly range across multiple blocks and across multiple processes. The handling
 * of this is implemented in this class.
 **********************************************************************************************************************/
template< typename Stencil_T >
class BubbleModel : public BubbleModelBase
{
 public:
   BubbleModel(const std::shared_ptr< StructuredBlockForest >& blockStorage, bool enableBubbleSplits);
   ~BubbleModel() override = default;

   // initialize bubble model from fill level field; bubble model is cleared and bubbles are created in cells with fill
   // level less than 1
   void initFromFillLevelField(const ConstBlockDataID& fillField);

   void setDensityOfAllBubbles(real_t rho) override;

   // mark the specified (gas) cell for belonging to the atmosphere bubble with constant pressure; the atmosphere's
   // bubble ID is set to the highest ID that was found on any of the processes at cells that belong to the atmosphere;
   // WARNING: This function must be called on all processes for the same cells, even if the cells are not located on
   // the current block.
   void setAtmosphere(const Cell& cellInGlobalCoordinates, real_t constantRho = real_c(1.0));

   // accessing
   real_t getDensity(IBlock* block, const Cell& cell) const override
   {
      // get the bubble containing cell
      const Bubble* bubble = getBubble(block, cell);
      WALBERLA_ASSERT_NOT_NULLPTR(bubble, "Cell " << cell << " does not belong to a bubble.");

      return bubble->getDensity();
   }

   void setDensity(IBlock* block, const Cell& cell, real_t value) override
   {
      // get the bubble containing cell
      Bubble* bubble = getBubble(block, cell);
      WALBERLA_ASSERT_NOT_NULLPTR(bubble, "Cell " << cell << " does not belong to a bubble.");

      bubble->setDensity(value);
   }

   const BubbleID& getBubbleID(IBlock* block, const Cell& cell) const;
   BubbleID& getBubbleID(IBlock* block, const Cell& cell);

   // cell and fill level update
   void reportFillLevelChange(IBlock* block, const Cell& cell, real_t fillLevelDifference) override;

   // assign a bubble ID (from the first found neighboring cell) to the newly created interface cell; if multiple bubble
   // IDs are found in neighborhood, register a merge
   void reportLiquidToInterfaceConversion(IBlock* block, const Cell& cell) override;

   // invalidate the bubble ID of the converted liquid cell and check for bubble splits
   void reportInterfaceToLiquidConversion(IBlock* block, const Cell& cell) override;

   ConstBlockDataID getBubbleFieldID() const { return bubbleFieldID_; }

   // combine information about a bubble
   struct BubbleInfo
   {
      BubbleInfo() : nrOfCells(uint_c(0)) {}
      Vector3< real_t > centerOfMass;
      uint_t nrOfCells;
      Bubble* bubble;
   };

   // compute bubbleInfo for each bubble and (MPI) reduce it on root
   std::vector< BubbleInfo > computeBubbleStats();

   // write bubbleInfo to terminal on root process
   void logBubbleStatsOnRoot();

   // update the bubble model:
   // - communicate bubble ID field
   // - merge and split bubbles (involves global MPI communications)
   void update() override;

 protected:
   const Bubble* getBubble(IBlock* block, const Cell& cell) const;
   Bubble* getBubble(IBlock* block, const Cell& cell);

   const std::vector< Bubble >& getBubbles() const { return bubbles_; };

   // check a 3x3x3 neighborhood whether a bubble could have split; calling is function is relatively inexpensive and
   // can be used as indicator whether the more expensive extendedSplitCheck() makes sense
   static bool checkForSplit(BubbleField_T* bf, const Cell& cell, BubbleID prevBubbleID);

   // check "neighborhood" cells in each direction around "cell" to ensure that a bubble has really split
   static bool extendedSplitCheck(BubbleField_T* bf, const Cell& cell, BubbleID oldBubbleID,
                                  cell_idx_t neighborhood = 2);

   using StencilForSplit_T =
      typename std::conditional< Stencil_T::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;

   // split bubbles, assign new global bubble IDs, and handle potential bubble merges resulting from splitting (involves
   // global MPI communication)
   void handleSplits();

   // delete potentially splitting bubbles and create new ones from them; new bubbles are added to newBubbleCom
   void markAndCreateSplittedBubbles(NewBubbleCommunication& newBubbleComm, const std::vector< bool >& splitIndicator);

   // in a 3x3x3 neighborhood, find all directions that are connected to startDir (by having same bubbleID) using a
   // flood fill algorithm; flood fill is more efficient than simply iterating over all neighbors since the latter would
   // require lots of extra logic
   static uint32_t mapNeighborhood(BubbleField_T* bf, stencil::Direction startDir, const Cell& cell, BubbleID bubbleID);

   // block storage to access bubbleField and fillField
   std::shared_ptr< StructuredBlockStorage > blockStorage_;

   // field that stores every cell's bubble ID; if a cell does not belong to any bubble, its ID is set to
   // INVALID_BUBBLE_ID; this field is managed by the BubbleModel and should not be passed outside
   BlockDataID bubbleFieldID_;

   // vector that stores all bubbles; it is kept synchronized across all processes
   std::vector< Bubble > bubbles_;

   // helper class to manage bubble merges
   MergeInformation mergeInformation_;

   // communication scheme for the bubble field
   blockforest::communication::UniformBufferedScheme< StencilForSplit_T > bubbleFieldCommunication_;

   // store split information, i.e., hints for splitting; store only hints since merges have to be processed first
   struct SplitHint
   {
      SplitHint(IBlock* _block, const Cell& _cell) : block(_block), cell(_cell) {}
      IBlock* block;
      Cell cell;
   };

   // vector with outstanding splits that (still) need to be processed
   std::vector< SplitHint > splitsToProcess_;

   std::shared_ptr< FloodFillInterface > floodFill_;

   // disable splits to decrease computational costs
   bool enableBubbleSplits_;

   inline BubbleField_T* getBubbleField(IBlock* block) const { return block->getData< BubbleField_T >(bubbleFieldID_); }

}; // class BubbleModel

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

using walberla::free_surface::bubble_model::BubbleModelBase;

#include "BubbleModel.impl.h"
