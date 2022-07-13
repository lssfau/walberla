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
//! \file MergeInformation.cpp
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Manage merging of bubbles.
//
//======================================================================================================================

#include "MergeInformation.h"

#include "core/mpi/MPIManager.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
MergeInformation::MergeInformation(uint_t numberOfBubbles)
{
   hasMerges_ = false;
   resize(numberOfBubbles);
}

/***********************************************************************************************************************
 * Implementation Note (see also MergeInformationTest.cpp):
 * ---------------------
 *    renameVec_ = { 0, 1, 2, 3, 4, 5 }
 *
 *    position 5 is renamed to the ID at position 3:
 *    (1)   registerMerge( 5, 3 );
 *    renameVec_ = { 0, 1, 2, 3, 4, 3 }
 *
 *    position 3 is renamed to the ID at position 1:
 *    (2)   registerMerge(3, 1);
 *    renameVec_ = { 0, 1, 2, 1, 4, 3 }
 *
 *    since position 3 was already renamed to the ID at position 1, registerMerge(1, 0) is called internally and
 *    position 1 is renamed to the ID at position 0; position 3 is renamed to the ID at position 0
 *    (3)   registerMerge( 3, 0 );
 *          -> leads to recursive registerMerge(1, 0);
 *    renameVec_ = { 0, 0, 2, 0, 4, 3 }
 *    since position 5 was already renamed to the ID of position 3, registerMerge(3, 2) is called internally; since
 *    position 3 was already renamed to the ID of position 0, registerMerge(2, 0) is called internally; position 2 is
 *    renamed to the ID at position 0, position 5 is renamed to the ID at position 0
 *    (4)   registerMerge( 5, 2)
 *          -> leads to recursive registerMerge( 3, 2)
 *                                registerMerge( 0, 2)
 *    renameVec_ = { 0, 0, 0, 0, 4, 0 }
 *
 *  Recursive call in Step(3) is necessary because otherwise the mapping from 3->1 would be "forgotten".
 *  However, not all transitive renames are resolved by this function as seen in step (3):  5->3->1
 *  Thus, the additional function resolveTransitiveRenames() is required.
 **********************************************************************************************************************/
void MergeInformation::registerMerge(BubbleID b0, BubbleID b1)
{
   WALBERLA_ASSERT_LESS(b0, renameVec_.size());
   WALBERLA_ASSERT_LESS(b1, renameVec_.size());

   // identical bubbles can not be merged
   if (b0 == b1) { return; }

   // register the merge using hasMerges_
   hasMerges_ = true;

   // ensure that b0 < b1 (increasing ID ordering is required for some functions later on)
   if (b1 < b0) { std::swap(b0, b1); }

   WALBERLA_ASSERT_LESS(b0, b1);

   // if bubble b1 is also marked for merging with another bubble, e.g., bubble b2
   if (isRenamed(b1))
   {
      // mark bubble b2 for merging with b0 (b2 = renameVec_[b1])
      registerMerge(b0, renameVec_[b1]);

      // mark bubble b1 for merging with bubble b0 or bubble b2 (depending on the ID order found in
      // registerMerge(b0,b2) above)
      renameVec_[b1] = renameVec_[renameVec_[b1]];
   }
   else
   {
      // mark bubble b1 for merging with bubble b0
      renameVec_[b1] = b0;
   }
}

void MergeInformation::mergeAndReorderBubbleVector(std::vector< Bubble >& bubbles)
{
   // no bubble merge possible with less than 2 bubbles
   if (bubbles.size() < uint_c(2)) { return; }

   // rename merged bubbles (always keep smallest bubble ID)
   resolveTransitiveRenames();

#ifndef NDEBUG
   uint_t numberOfMerges = countNumberOfRenames();
#endif

   WALBERLA_ASSERT_EQUAL(bubbles.size(), renameVec_.size());
   WALBERLA_ASSERT_EQUAL(renameVec_[0], 0);

   for (size_t i = uint_c(0); i < bubbles.size(); ++i)
   {
      // any entry that does not point to itself needs to be merged
      if (renameVec_[i] != i)
      {
         WALBERLA_ASSERT_LESS(renameVec_[i], i);

         // merge bubbles
         bubbles[renameVec_[i]].merge(bubbles[i]);
      }
   }

   // create temporary vector with "index = bubble ID" to store the exchanged bubble IDs
   std::vector< BubbleID > exchangeVector(bubbles.size());
   for (size_t i = uint_c(0); i < exchangeVector.size(); ++i)
   {
      exchangeVector[i] = BubbleID(i);
   }

   // gapPointer searches for the renamed bubbles in increasing order; vector entries found by gapPointer are "gaps" and
   // can be overwritten
   uint_t gapPointer = uint_c(1);
   while (gapPointer < bubbles.size() && !isRenamed(gapPointer)) // find first renamed bubble
   {
      ++gapPointer;
   }

   // lastPointer searches for not-renamed bubbles in decreasing order; vector entries found by lastPointer remain in
   // the vector and will be copied to the gaps
   uint_t lastPointer = uint_c(renameVec_.size() - 1);
   while (isRenamed(lastPointer) && lastPointer > uint_c(0)) // find last non-renamed bubble
   {
      --lastPointer;
   }

   while (lastPointer > gapPointer)
   {
      // exchange the last valid (non-renamed) bubble ID with the first non-valid (renamed) bubble ID; anything
      // gapPointer points to will be deleted later; this reorders the bubble vector
      std::swap(bubbles[gapPointer], bubbles[lastPointer]);

      // store the above exchange
      exchangeVector[lastPointer] = BubbleID(gapPointer);
      exchangeVector[gapPointer]  = BubbleID(lastPointer);

      // update lastPointer, i.e., find next non-renamed bubble
      do
      {
         --lastPointer;
         // important condition: "lastPointer > gapPointer" since gapPointer is valid now
      } while (isRenamed(lastPointer) && lastPointer > gapPointer);

      // update gapPointer, i.e., find next renamed bubble
      do
      {
         ++gapPointer;
      } while (gapPointer < bubbles.size() && !isRenamed(gapPointer));
   }

   // shrink bubble vector (any element after lastPointer is not valid and can be removed)
   uint_t newSize = lastPointer + uint_c(1);
   WALBERLA_ASSERT_EQUAL(newSize, bubbles.size() - numberOfMerges);
   WALBERLA_ASSERT_LESS_EQUAL(newSize, bubbles.size());
   bubbles.resize(newSize);

   // update renameVec_ with exchanged bubble IDs with the highest unnecessary bubble IDs being dropped; this ensures
   // that bubble IDs are always numbered continuously from 0 upwards
   for (size_t i = 0; i < renameVec_.size(); ++i)
   {
      renameVec_[i] = exchangeVector[renameVec_[i]];
      WALBERLA_ASSERT_LESS(renameVec_[i], newSize);
   }
}

void MergeInformation::communicateMerges()
{
   // merges can only be communicated if they occurred and are registered in renameVec_
   if (renameVec_.empty()) { return; }

   // rename process local bubble IDs
   resolveTransitiveRenames();

   // globally combine all rename vectors
   int numProcesses = MPIManager::instance()->numProcesses();
   std::vector< BubbleID > allRenameVectors(renameVec_.size() * uint_c(numProcesses));
   WALBERLA_MPI_SECTION()
   {
      MPI_Allgather(&renameVec_[0], int_c(renameVec_.size()), MPITrait< BubbleID >::type(), &allRenameVectors[0],
                    int_c(renameVec_.size()), MPITrait< BubbleID >::type(), MPI_COMM_WORLD);
   }

   // check for inter-process bubble merges
   for (size_t i = renameVec_.size() - 1; i > 0; --i)
      for (int process = 0; process < numProcesses; ++process)
      {
         if (process == MPIManager::instance()->rank()) { continue; } // local merges have already been treated

         size_t idx = uint_c(process) * renameVec_.size() + i;

         // register inter-process bubble merge (this updated renameVec_)
         if (allRenameVectors[idx] != i) { registerMerge(allRenameVectors[idx], BubbleID(i)); }
      }

   // rename global bubble IDs
   resolveTransitiveRenames();
}

void MergeInformation::renameOnBubbleField(BubbleField_T* bubbleField) const
{
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(bubbleField, {
      const typename BubbleField_T::Ptr bubblePtr(*bubbleField, x, y, z);
      if (*bubblePtr != INVALID_BUBBLE_ID) { *bubblePtr = renameVec_[*bubblePtr]; }
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}

void MergeInformation::resizeAndClear(uint_t numberOfBubbles)
{
   renameVec_.resize(numberOfBubbles);
   for (uint_t i = uint_c(0); i < numberOfBubbles; ++i)
   {
      renameVec_[i] = BubbleID(i);
   }

   hasMerges_ = false;
}

void MergeInformation::resize(uint_t numberOfBubbles)
{
   if (numberOfBubbles > renameVec_.size())
   {
      renameVec_.reserve(numberOfBubbles);
      for (size_t i = renameVec_.size(); i < numberOfBubbles; ++i)
      {
         renameVec_.push_back(BubbleID(i));
      }
   }
}

void MergeInformation::resolveTransitiveRenames()
{
   WALBERLA_ASSERT_GREATER(renameVec_.size(), 0);
   WALBERLA_ASSERT_EQUAL(renameVec_[0], 0);

   for (size_t i = renameVec_.size() - 1; i > 0; --i)
   {
      // create new bubble ID for each entry in renameVec_
      BubbleID& newBubbleID = renameVec_[i];

      // example 1: "renameVec_[4] = 2" means that bubble 4 has merged with bubble 2
      // => bubble 4 is renamed to bubble 2 (always the smaller ID is kept)
      // example 2: "renameVec_[4] = 2" and "renameVec_[2] = 1" means above and that bubble 2 has merged with bubble 1
      // => bubble 4 should finally be renamed to bubble 1

      // this loop ensures that renaming is correct even if multiple bubbles have merged (as in example 2 from above)
      while (renameVec_[newBubbleID] != newBubbleID)
      {
         newBubbleID = renameVec_[newBubbleID];
      }
   }

   WALBERLA_ASSERT(transitiveRenamesResolved()); // ensure that renaming was resolved correctly
}

bool MergeInformation::transitiveRenamesResolved() const
{
   for (size_t i = 0; i < renameVec_.size(); ++i)
   {
      // A = renameVec_[i]
      // B = renameVec_[A]
      // A must be equal to B after bubble renaming
      if (renameVec_[i] != renameVec_[renameVec_[i]]) { return false; }
   }

   // ensure that no entry points to an element that has been renamed
   for (size_t i = 0; i < renameVec_.size(); ++i)
   {
      if (renameVec_[i] != i) // bubble at entry i has been renamed
      {
         for (size_t j = 0; j < renameVec_.size(); ++j)
         {
            // renameVec_[j] points to entry i although i has been renamed
            if (i != j && renameVec_[j] == i) { return false; }
         }
      }
   }

   return true;
}

uint_t MergeInformation::countNumberOfRenames() const
{
   uint_t renames = uint_c(0);
   for (size_t i = 0; i < renameVec_.size(); ++i)
   {
      // only bubbles with "bubble ID != index" have been renamed
      if (renameVec_[i] != i) { ++renames; }
   }

   return renames;
}

void MergeInformation::print(std::ostream& os) const
{
   os << "Merge Information: ";

   for (size_t i = 0; i < renameVec_.size(); ++i)
      os << "( " << i << " -> " << renameVec_[i] << " ), ";

   os << std::endl;
}

std::ostream& operator<<(std::ostream& os, const MergeInformation& mi)
{
   mi.print(os);
   return os;
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
