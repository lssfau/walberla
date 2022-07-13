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
//! \file BubbleModel.impl.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief System for tracking pressure/density in gas volumes.
//
//======================================================================================================================

#include "field/AddToStorage.h"

#include "lbm/free_surface/InterfaceFromFillLevel.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q27.h"

#include "BubbleIDFieldPackInfo.h"
#include "BubbleModel.h"
#include "RegionalFloodFill.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
template< typename Stencil_T >
BubbleModel< Stencil_T >::BubbleModel(const std::shared_ptr< StructuredBlockForest >& blockStorage,
                                      bool enableBubbleSplits)
   : blockStorage_(blockStorage), bubbleFieldID_(field::addToStorage< BubbleField_T >(
                                     blockStorage, "BubbleIDs", BubbleID(INVALID_BUBBLE_ID), field::fzyx, uint_c(1))),
     bubbleFieldCommunication_(blockStorage), enableBubbleSplits_(enableBubbleSplits)
{
   bubbleFieldCommunication_.addPackInfo(
      std::make_shared< BubbleIDFieldPackInfo< Stencil_T > >(bubbleFieldID_, &mergeInformation_));
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::initFromFillLevelField(const ConstBlockDataID& fillFieldID)
{
   // mark regions belonging to the same bubble
   floodFill_ = std::make_shared< FloodFillUsingFillLevel< Stencil_T > >(fillFieldID);

   bubbles_.clear();

   NewBubbleCommunication newBubbleComm;

   // start numbering the bubbles from 0
   BubbleID firstNewBubbleID = 0;
   BubbleID nextID           = firstNewBubbleID;

   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      // get fields
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(fillFieldID);
      BubbleField_T* const bubbleField     = blockIt->getData< BubbleField_T >(bubbleFieldID_);

      // initialize bubbleField with invalid IDs
      bubbleField->set(INVALID_BUBBLE_ID);

      auto bubbleIt    = bubbleField->begin();
      auto fillFieldIt = fillField->begin();

      while (bubbleIt != bubbleField->end())
      {
         // only consider cells
         // - that are either gas or interface
         // - for which no bubble ID is set yet (each call to floodFill_->run() sets new bubbleIDs to cells)
         if ((*fillFieldIt < real_c(1) || isInterfaceFromFillLevel< Stencil_T >(*fillField, bubbleIt.cell())) &&
             *bubbleIt == INVALID_BUBBLE_ID)
         {
            real_t volume;
            uint_t nrOfCells;

            // set (block local, preliminary) bubble IDs in the bubble field
            floodFill_->run(*blockIt, bubbleFieldID_, bubbleIt.cell(), nextID++, volume, nrOfCells);

            // create new bubble with preliminary bubbleIDs; the final global bubbleIDs are set in communicateAndApply()
            newBubbleComm.createBubble(Bubble(volume));
         }

         ++bubbleIt;
         ++fillFieldIt;
      }
   }

   // set global bubble IDs
   newBubbleComm.communicateAndApply(bubbles_, *blockStorage_, bubbleFieldID_);

   // clear merge information, i.e., make sure that no merge is already registered
   mergeInformation_.resizeAndClear(bubbles_.size());

   // communicate bubble field IDs
   bubbleFieldCommunication_();

   // communicate bubble merges
   mergeInformation_.communicateMerges();

   if (mergeInformation_.hasMerges())
   {
      // merge bubbles (add bubble volumes and delete/rename bubble IDs)
      mergeInformation_.mergeAndReorderBubbleVector(bubbles_);

      // rename bubble IDs on all blocks
      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
         mergeInformation_.renameOnBubbleField(blockIt->getData< BubbleField_T >(bubbleFieldID_));
   }

   // clear merge information after bubbles have been merged
   mergeInformation_.resizeAndClear(bubbles_.size());
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::setDensityOfAllBubbles(real_t rho)
{
   for (auto it = bubbles_.begin(); it != bubbles_.end(); ++it)
   {
      if (!it->hasConstantDensity()) { it->setDensity(rho); }
   }
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::setAtmosphere(const Cell& cellInGlobalCoordinates, real_t rho)
{
   // temporarily set the atmosphere's bubble ID to invalid
   BubbleID atmosphereBubbbleID = INVALID_BUBBLE_ID;

   // get the cell's block (or nullptr if the block is not on this process)
   IBlock* blockWithAtmosphereBubble = blockStorage_->getBlock(cellInGlobalCoordinates);

   // set atmosphere bubble ID to this cell's bubble ID
   if (blockWithAtmosphereBubble) // else: block does not exist locally
   {
      Cell localCell;
      blockStorage_->transformGlobalToBlockLocalCell(localCell, *blockWithAtmosphereBubble, cellInGlobalCoordinates);

      const BubbleField_T* const bf = blockWithAtmosphereBubble->getData< const BubbleField_T >(bubbleFieldID_);

      // get this cell's bubble ID
      atmosphereBubbbleID = bf->get(localCell);

      // cell must be a gas cell; therefore, a valid bubble ID must be set
      WALBERLA_ASSERT_UNEQUAL(atmosphereBubbbleID, INVALID_BUBBLE_ID);
   }

   // variable for (MPI) reducing the bubble ID; value of -1 is set in order to ignore this bubble in maximum reduction
   int reducedBubbleID = atmosphereBubbbleID != INVALID_BUBBLE_ID ? int_c(atmosphereBubbbleID) : -1;

   // get the highest of all processes' bubble IDs
   WALBERLA_MPI_SECTION()
   {
      MPI_Allreduce(MPI_IN_PLACE, &reducedBubbleID, 1, MPITrait< int >::type(), MPI_MAX, MPI_COMM_WORLD);
   }

   if (reducedBubbleID < 0)
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Could not set atmosphere in the non-gas cell " << cellInGlobalCoordinates);
   }
   else
   {
      // set atmosphere bubble ID to highest of all processes' bubble IDs
      atmosphereBubbbleID = BubbleID(reducedBubbleID);
      bubbles_[atmosphereBubbbleID].setConstantDensity(rho);
   }
}

template< typename Stencil_T >
const Bubble* BubbleModel< Stencil_T >::getBubble(IBlock* blockIt, const Cell& cell) const
{
   const BubbleField_T* bf = getBubbleField(blockIt);
   const BubbleID id       = bf->get(cell);

   return (id == INVALID_BUBBLE_ID) ? nullptr : &bubbles_[id];
}

template< typename Stencil_T >
Bubble* BubbleModel< Stencil_T >::getBubble(IBlock* blockIt, const Cell& cell)
{
   const BubbleField_T* bf = getBubbleField(blockIt);
   const BubbleID id       = bf->get(cell);

   return (id == INVALID_BUBBLE_ID) ? nullptr : &bubbles_[id];
}

template< typename Stencil_T >
const BubbleID& BubbleModel< Stencil_T >::getBubbleID(IBlock* blockIt, const Cell& cell) const
{
   const BubbleField_T* bf = getBubbleField(blockIt);
   return bf->get(cell);
}

template< typename Stencil_T >
BubbleID& BubbleModel< Stencil_T >::getBubbleID(IBlock* blockIt, const Cell& cell)
{
   BubbleField_T* bf = getBubbleField(blockIt);
   return bf->get(cell);
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::reportFillLevelChange(IBlock* blockIt, const Cell& cell, real_t fillLevelDifference)
{
   Bubble* b = getBubble(blockIt, cell);
   WALBERLA_ASSERT_NOT_NULLPTR(b,
                               "Reporting fill level change in cell " << cell << " where no bubble ID is registered.");

   // update the bubble volume change; fillLevelDifference is negated because variable is:
   // - positive if fill level increased => bubble volume has to decrease
   // - negative if fill level decreased => bubble volume has to increase
   if (b) { b->updateVolumeDiff(-fillLevelDifference); }
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::reportLiquidToInterfaceConversion(IBlock* blockIt, const Cell& cell)
{
   // get bubble field
   BubbleField_T* bf = getBubbleField(blockIt);

   // get this cell's bubble ID
   BubbleID& thisCellID = bf->get(cell);

   // this cell is converted from liquid to interface; liquid cells have no bubble ID such that this cell can not have
   // a bubble ID, yet
   WALBERLA_ASSERT_EQUAL(thisCellID, INVALID_BUBBLE_ID);

   // iterate neighborhood and assign the first found bubble ID to this new interface cell
   using SearchStencil_T = typename std::conditional< Stencil_T::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;
   for (auto d = SearchStencil_T::beginNoCenter(); d != SearchStencil_T::end(); ++d)
   {
      // get bubble ID of neighboring cell
      BubbleID neighborID = bf->get(cell[0] + d.cx(), cell[1] + d.cy(), cell[2] + d.cz());
      if (neighborID != INVALID_BUBBLE_ID)
      {
         // assign the first found neighbor's bubble ID to this cell
         if (thisCellID == INVALID_BUBBLE_ID) { thisCellID = neighborID; }
         else
         {
            // if multiple different bubble IDs are in neighborhood, trigger merging
            if (thisCellID != neighborID) { mergeInformation_.registerMerge(thisCellID, neighborID); }
         }
      }
   }
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::reportInterfaceToLiquidConversion(IBlock* blockIt, const Cell& cell)
{
   // get bubble field
   BubbleField_T* bf = getBubbleField(blockIt);

   // get this cell's bubble ID
   BubbleID oldBubbleID = bf->get(cell);

   // this cell is converted from interface to liquid; the interface cell must already have a valid bubble ID
   WALBERLA_ASSERT_UNEQUAL(oldBubbleID, INVALID_BUBBLE_ID);

   // invalidate the converted cell's bubble ID (liquid cells must not have a bubble ID)
   bf->get(cell) = INVALID_BUBBLE_ID;

   if (enableBubbleSplits_)
   {
      // check a 3x3x3 neighborhood whether a bubble could have split
      if (checkForSplit(bf, cell, oldBubbleID))
      {
         WALBERLA_LOG_INFO("Possible bubble split detected due to conversion in cell " << cell << ".");

         // check a larger neighborhood to ensure that the bubble has really split
         if (extendedSplitCheck(bf, cell, oldBubbleID, cell_idx_c(3)))
         {
            WALBERLA_LOG_INFO("Extended split check confirmed split.");
            // register this bubble split
            splitsToProcess_.emplace_back(blockIt, cell);
         }
         else { WALBERLA_LOG_INFO("Extended split check ruled out split."); }
      }
   }
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::update()
{
   // communicate field with bubble IDs
   bubbleFieldCommunication_();

   // vector for (MPI) reducing each bubble's volume change and indicators for whether merges and splits occurred
   static std::vector< real_t > reduceVector;

   // indicators that are (MPI) reduced to identify merges and splits
   real_t mergeIndicator = mergeInformation_.hasMerges() ? real_c(1) : real_c(0);
   real_t splitIndicator = splitsToProcess_.empty() ? real_c(0) : real_c(1);

   // extend the vector for storing the merge and split indicator
   reduceVector.resize(bubbles_.size() + 2);

   uint_t i = uint_c(0);
   for (; i < bubbles_.size(); ++i)
   {
      // get each bubble's volume change
      reduceVector[i] = bubbles_[i].getAndResetVolumeDiff();
   }

   // append the indicators at the end of reduceVector
   reduceVector[i++] = mergeIndicator;
   reduceVector[i++] = splitIndicator;
   WALBERLA_ASSERT_EQUAL(i, reduceVector.size()); // make sure that indexing is correct

   WALBERLA_MPI_SECTION()
   {
      // globally (MPI) reduce each bubble's volume change, the number of merges, and the number of splits
      MPI_Allreduce(MPI_IN_PLACE, &reduceVector[0], int_c(reduceVector.size()), MPITrait< real_t >::type(), MPI_SUM,
                    MPI_COMM_WORLD);
   }

   uint_t j = uint_c(0);
   for (; j < bubbles_.size(); ++j)
   {
      // update each bubble's volume and density
      bubbles_[j].applyVolumeDiff(reduceVector[j]);
   }

   // check for merges and splits
   bool mergeHappened = (reduceVector[j++] > real_c(0));
   bool splitHappened = (reduceVector[j++] > real_c(0));
   WALBERLA_ASSERT_EQUAL(j, reduceVector.size()); // make sure that indexing is correct

   // treat bubble merges
   if (mergeHappened)
   {
      WALBERLA_ROOT_SECTION()
      {
         // std::stringstream ss;
         // mergeInformation_.print(ss);
         // WALBERLA_LOG_INFO("Merge detected, " << ss.str());
         WALBERLA_LOG_INFO("Merge detected");
      }

      // globally communicate bubble merges and rename bubble IDs accordingly
      mergeInformation_.communicateMerges();

      // merge bubbles
      mergeInformation_.mergeAndReorderBubbleVector(bubbles_);

      // update, i.e., rename bubble IDs in the bubble ID field
      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
      {
         mergeInformation_.renameOnBubbleField(getBubbleField(&(*blockIt)));
      }
   }

   // treat bubble splits
   if (splitHappened) { handleSplits(); }

   mergeInformation_.resizeAndClear(bubbles_.size());
   splitsToProcess_.clear();
}

template< typename Stencil_T >
bool BubbleModel< Stencil_T >::checkForSplit(BubbleField_T* bf, const Cell& cell, BubbleID previousBubbleID)
{
   using namespace stencil;

   // Variable neighborHoodInfo has bits set in each connected neighboring direction where the cell belongs to bubble
   // previousBubbleID (=ID of the bubble that got "lost" by converting cell from interface to liquid). The
   // neighborHoodInfo is built via flood fill exactly once from a starting direction (e.g., first direction in which
   // previousBubbleID was found). Thus, the connected neighborhood is built for only one direction and cells that
   // are not connected to this starting direction are not part of the connected neighborhood. A split occurs, i.e., is
   // detected when previousBubbleID is found in an unconnected region.
   // Example:
   // previousBubbleID is 1, the cells in the center were converted from bubble ID=1 to liquid (ID=f)
   //   1 1 1
   //   f f f
   //   1 1 1
   // The first direction in which previousBubbleID is set shall be N (north). The connected neighborhood is then built
   // (c=connected, n=not connected):
   //   c c c
   //   n n n
   //   n n n
   // At neighbor N, neighbors NW and NE are found: no split is detected since we are in the connected neighborhood.
   // Then, neighbor S is found. Since S is not in the (first) connected neighborhood, a split is detected.
   uint32_t neighborHoodInfo = uint32_c(0);

   for (auto d = StencilForSplit_T::beginNoCenter(); d != StencilForSplit_T::end(); ++d)
   {
      // get the bubble ID of the neighboring cell
      BubbleID neighborID = bf->getNeighbor(cell, *d);

      // neighboring bubble is a different one (or no bubble) than this cell's bubble
      if (neighborID != previousBubbleID) { continue; }
      // => from here: bubbles are the same, i.e., neighborID == previousBubbleID

      if (neighborHoodInfo > uint32_c(0)) // the neighborhood map has already been created
      {
         // "connected" bit is set in this direction, i.e., the neighbor is connected
         if (neighborHoodInfo & dirToBinary[*d]) { continue; }
         else // "connected" bit is not set in this direction, i.e., the neighbor is not connected
         {
            // since neighborID == previousBubbleID but neighbor is not connected, a split had to occur
            return true;
         }
      }
      else // the neighborhood map has not been created, yet
      {
         // create connected neighborhood starting from direction d
         neighborHoodInfo = mapNeighborhood(bf, *d, cell, neighborID);
      }
   }
   return false;
}

template< typename Stencil_T >
bool BubbleModel< Stencil_T >::extendedSplitCheck(BubbleField_T* bf, const Cell& cell, BubbleID previousBubbleID,
                                                  cell_idx_t neighborhood)
{
   // RegionalFloodFill is used to find connected regions in a larger (>3x3x3) neighborhood
   RegionalFloodFill< BubbleID, StencilForSplit_T >* neighborHoodInfo = nullptr;

   for (auto d = StencilForSplit_T::beginNoCenter(); d != StencilForSplit_T::end(); ++d)
   {
      // get the bubble ID of the neighboring cell
      BubbleID neighborID = bf->getNeighbor(cell, *d);

      // neighboring bubble is a different one (or no bubble) than this cell's bubble
      if (neighborID != previousBubbleID) { continue; }
      // => from here: bubbles are the same, i.e., neighborID == previousBubbleID

      if (neighborHoodInfo) // the neighborhood map has already been created
      {
         if (neighborHoodInfo->connected(*d)) { continue; }
         else // bubble is not connected in direction d and a split occurred
         {
            delete neighborHoodInfo;
            return true;
         }
      }
      else // the neighborhood map has not been created, yet
      {
         // create connected neighborhood starting from direction d
         neighborHoodInfo =
            new RegionalFloodFill< BubbleID, StencilForSplit_T >(bf, cell, *d, neighborID, neighborhood);
      }
   }

   delete neighborHoodInfo;
   return false;
}

template< typename Stencil_T >
uint32_t BubbleModel< Stencil_T >::mapNeighborhood(BubbleField_T* bf, stencil::Direction startDir, const Cell& cell,
                                                   BubbleID bubbleID)
{
   using namespace stencil;

   uint32_t result = uint32_c(0);

   // use stack to store directions that still need to be searched
   std::vector< Direction > stack;
   stack.push_back(startDir);

   while (!stack.empty())
   {
      // next search direction is the last entry in stack
      Direction d = stack.back();

      // remove the current search direction from stack
      stack.pop_back();

      WALBERLA_ASSERT(d != C); // do not search in center direction
      WALBERLA_ASSERT(bf->get(cell[0] + cx[d], cell[1] + cy[d], cell[2] + cz[d]) ==
                      bubbleID); // cell must belong to the same bubble

      // add this direction to result, i.e., to the "connected neighborhood" using bitwise OR
      result |= dirToBinary[d];

      // in direction d, iterate over d's neighboring directions i, i.e., iterate over a (at maximum) 3x3x3 neighborhood
      // from the viewpoint of cell
      for (uint_t i = uint_c(1); i < StencilForSplit_T::dir_neighbors_length[d]; ++i)
      {
         // transform direction d's neighbor in direction i to a direction from the viewpoint of cell, e.g., d=N, i=W =>
         // nDir=NW
         Direction nDir = StencilForSplit_T::dir_neighbors[d][i];

         // cell in direction nDir belongs to the bubble and is not already in result
         if (bf->get(cell[0] + cx[nDir], cell[1] + cy[nDir], cell[2] + cz[nDir]) == bubbleID &&
             (result & dirToBinary[nDir]) == 0)
         {
            // add nDir to stack such that it gets added to result and used as start cell in the next iteration;
            // the connected bit will never be set for directions that are not connected to startDir
            stack.push_back(nDir);
         }
      }
   }

   return result;
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::handleSplits()
{
   // splitIndicator is set to true for bubbles for which a split was registered; this can not be done earlier since
   // bubble IDs may have changed during merging
   std::vector< bool > splitIndicator(bubbles_.size());

   // process all registered splits
   for (auto i = splitsToProcess_.begin(); i != splitsToProcess_.end(); ++i)
   {
      BubbleField_T* bubbleField = getBubbleField(i->block);

      const Cell& c = i->cell;

      // cell c was transformed from interface to liquid and has no valid bubble ID anymore; mark remaining gas cells
      // with valid bubble IDs for being split
      for (auto d = StencilForSplit_T::begin(); d != StencilForSplit_T::end(); ++d)
      {
         BubbleID id = bubbleField->get(c.x() + d.cx(), c.y() + d.cy(), c.z() + d.cz());

         // mark bubble's cell for splitting in splitIndicator
         if (id != INVALID_BUBBLE_ID) { splitIndicator[id] = true; }
      }
   }

   // communicate (MPI reduce) splitIndicator among all processes
   allReduceInplace(splitIndicator, mpi::BITWISE_OR);

   // create communication for new bubbles
   NewBubbleCommunication newBubbleComm(bubbles_.size());

   // treat split and create new (splitted) bubbles with new IDs
   markAndCreateSplittedBubbles(newBubbleComm, splitIndicator);

   // communicate all bubbles and assign new global bubble IDs
   newBubbleComm.communicateAndApply(bubbles_, splitIndicator, *blockStorage_, bubbleFieldID_);

   // clear merge information
   mergeInformation_.resizeAndClear(bubbles_.size());

   // communicate bubble field
   bubbleFieldCommunication_();

   // communicate merges
   mergeInformation_.communicateMerges();

   // merge new splitted bubbles (merges can occur at the block border after the bubbleField was communicated;
   // bubbleFieldCommunication is linked to mergeInformation_ to report the merges there)
   if (mergeInformation_.hasMerges())
   {
      // merge bubbles
      mergeInformation_.mergeAndReorderBubbleVector(bubbles_);
      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
      {
         // assign new bubble IDs
         mergeInformation_.renameOnBubbleField(blockIt->getData< BubbleField_T >(bubbleFieldID_));
      }
   }

   // clear merge information
   mergeInformation_.resizeAndClear(bubbles_.size());
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::markAndCreateSplittedBubbles(NewBubbleCommunication& newBubbleComm,
                                                            const std::vector< bool >& splitIndicator)
{
   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      BubbleField_T* bubbleField = getBubbleField(&(*blockIt));

      // loop over the whole domain and start a flood fill at positions where bubbles are marked in splitIndicator
      for (auto bubbleIt = bubbleField->begin(); bubbleIt != bubbleField->end(); ++bubbleIt)
      {
         // skip cells that
         // - do not belong to a bubble
         // - belong to a bubble that might have been created during this function call (i.e., this bubble ID is not yet
         // known to the vector spliIndicator)
         // - belong to a bubble for which no split was detected
         if (*bubbleIt == INVALID_BUBBLE_ID || *bubbleIt >= splitIndicator.size() || !splitIndicator[*bubbleIt])
         {
            continue;
         }

         const Cell& curCell       = bubbleIt.cell();
         real_t densityOfOldBubble = getDensity(&(*blockIt), curCell);
         real_t volume;
         uint_t nrOfCells;

         // mark the whole region of the bubble (to which this cells belongs) in bubbleField
         floodFill_->run(*blockIt, bubbleFieldID_, curCell, newBubbleComm.nextFreeBubbleID(), volume, nrOfCells);

         // create new bubble
         newBubbleComm.createBubble(Bubble(volume, densityOfOldBubble));
      }
   }
}

template< typename Stencil_T >
std::vector< typename BubbleModel< Stencil_T >::BubbleInfo > BubbleModel< Stencil_T >::computeBubbleStats()
{
   std::vector< BubbleInfo > bubbleStats;
   bubbleStats.assign(bubbles_.size(), BubbleInfo());

   // iterate all bubbles on each block
   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      const BubbleField_T* bubbleField = getBubbleField(&(*blockIt));
      WALBERLA_FOR_ALL_CELLS(bubbleFieldIt, bubbleField, {
         if (*bubbleFieldIt == INVALID_BUBBLE_ID) { continue; }

         Vector3< real_t > cellCenter;
         Cell globalCell;
         blockStorage_->transformBlockLocalToGlobalCell(globalCell, *blockIt, bubbleFieldIt.cell());
         blockStorage_->getCellCenter(cellCenter[0], cellCenter[1], cellCenter[2], globalCell);

         // center of mass of this bubble on this block
         bubbleStats[*bubbleFieldIt].centerOfMass += cellCenter;

         // bubble's number of cells
         bubbleStats[*bubbleFieldIt].nrOfCells++;
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // store bubble information in reduceVec; reducing a vector of real_t is significantly easier than reducing
   // bubbleStats (vector of bubbleInfo)
   std::vector< real_t > reduceVec;
   for (auto statsIter = bubbleStats.begin(); statsIter != bubbleStats.end(); ++statsIter)
   {
      reduceVec.push_back(statsIter->centerOfMass[0]);
      reduceVec.push_back(statsIter->centerOfMass[1]);
      reduceVec.push_back(statsIter->centerOfMass[2]);
      reduceVec.push_back(real_c(statsIter->nrOfCells));
   }

   // (MPI) reduce bubble information on root
   mpi::reduceInplace(reduceVec, mpi::SUM, 0, MPI_COMM_WORLD);

   WALBERLA_ROOT_SECTION()
   {
      uint_t idx = uint_c(0);
      uint_t i   = uint_c(0);
      for (auto statsIter = bubbleStats.begin(); statsIter != bubbleStats.end(); ++statsIter)
      {
         statsIter->centerOfMass[0] = reduceVec[idx++];
         statsIter->centerOfMass[1] = reduceVec[idx++];
         statsIter->centerOfMass[2] = reduceVec[idx++];
         statsIter->nrOfCells       = uint_c(reduceVec[idx++]);

         // compute the bubble's global center of mass
         statsIter->centerOfMass /= real_c(statsIter->nrOfCells);

         // store the bubble ID
         statsIter->bubble = &(bubbles_[i]);
         ++i;
      }
      WALBERLA_ASSERT_EQUAL(idx, reduceVec.size());

      return bubbleStats;
   }
   else // else belongs to macro WALBERLA_ROOT_SECTION()
   {
      return std::vector< BubbleInfo >();
   }
}

template< typename Stencil_T >
void BubbleModel< Stencil_T >::logBubbleStatsOnRoot()
{
   std::vector< BubbleInfo > bubbleStats = computeBubbleStats();

   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_LOG_RESULT("Bubble Status:");
      for (auto it = bubbleStats.begin(); it != bubbleStats.end(); ++it)
      {
         WALBERLA_LOG_RESULT("\tPosition:" << it->centerOfMass << "  #Cells: " << it->nrOfCells);
      }
   }
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
