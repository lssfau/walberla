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
//! \file MergeInformation.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Manage merging of bubbles.
//
//======================================================================================================================

#pragma once

#include <iostream>

#include "Bubble.h"
#include "BubbleDefinitions.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
 * Class for managing bubble merging.
 *
 * - assumption: each process holds a synchronized, global vector of bubbles
 * - each process creates a new MergeInformation class, which internally holds a mapping from BubbleID -> BubbleID
 * - then during the time step, several merges can be registered using registerMerge() member function
 * - this information is then communicated between all processes using communicateMerges()
 * - Then the merges are executed on each process locally. This means that bubble volumes have to be added,
 *   and bubbles have to be deleted. To still have a continuous array, the last elements are copied to the place where
 *   the deleted bubbles have been -> renaming of bubbles.
 * - the renaming has to be done also in the bubble field using renameOnBubbleField()
 **********************************************************************************************************************/
class MergeInformation
{
 public:
   MergeInformation(uint_t numberOfBubbles = uint_c(0));

   void registerMerge(BubbleID b0, BubbleID b1);

   // globally communicate bubble merges and rename bubble IDs accordingly in renameVec_
   void communicateMerges();

   // merge bubbles according to renameVec_, and reorder and shrink the bubble vector such that bubble IDs are always
   // numbered continuously from 0 upwards
   void mergeAndReorderBubbleVector(std::vector< Bubble >& bubbles);

   // rename all bubbles in the bubble field according to renameVec_
   void renameOnBubbleField(BubbleField_T* bf) const;

   void resizeAndClear(uint_t numberOfBubbles);

   bool hasMerges() const { return hasMerges_; }

   void print(std::ostream& os) const;

   // vector for tracking bubble renaming:
   // - "renameVec_[4] = 2": bubble 4 has merged with bubble 2, i.e., any 4 can be replaced by 2 in the bubble field
   // - "renameVec_[i] <= i" for all i: when two bubbles merge, the smaller bubble ID is chosen for the resulting bubble
   std::vector< BubbleID > renameVec_;

 private:
   // adapt size of permutation vector, and init new elements with "identity"
   void resize(uint_t numberOfBubbles);

   inline bool isRenamed(size_t i) const { return renameVec_[i] != i; }

   // ensure that renaming of bubble IDs is correct even if multiple bubbles merge during the same time step
   void resolveTransitiveRenames();

   // check whether renaming has been done correctly, i.e., check the correctness of renameVec_
   bool transitiveRenamesResolved() const;

   uint_t countNumberOfRenames() const;

   // signals that merges have to be treated in this time step, true whenever renameVec_ is not the identity mapping
   bool hasMerges_;

   // friend function used in unit tests (MergeInformationTest)
   friend void checkRenameVec(const MergeInformation& mi, const std::vector< BubbleID >& vecCompare);
}; // MergeInformation

std::ostream& operator<<(std::ostream& os, const MergeInformation& mi);

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
