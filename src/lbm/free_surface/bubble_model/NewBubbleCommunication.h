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
//! \file NewBubbleCommunication.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Communication for the creation of new bubbles.
//
//======================================================================================================================

#pragma once

#include "core/logging/Logging.h"

#include <vector>

#include "Bubble.h"
#include "BubbleDefinitions.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
// forward declaration of internal data structure
class CommunicatedNewBubbles;

/***********************************************************************************************************************
 * Communication for the creation of new bubbles.
 *
 *  - Manages the distribution of newly created bubbles.
 *  - Bubbles must have the same BubbleID everywhere, since BubbleID is an index in the global bubble array. Thus, a
 *    problem arises when two (or more) processes add a new bubble in the same time step.
 *  - Every bubble is first created via this class and assigned a temporary BubbleID that is not used, yet.
 *  - These data are then synchronized in a communication step.
 *  - New bubbles are sorted by process rank such that there is then a sorted array of new bubbles that is identical on
 *    each process.
 *  - This vector of new bubbles can then be appended/inserted into the global bubble array.
 *  - Finally, temporary bubble IDs in the bubble field are renamed to global IDs.
 *
 *  Usage example:
 *  \code
     NewBubbleCommunication newBubbleComm;

     //the following lines are usually different on different processes:
     BubbleID newID     = newBubbleComm.createBubble( Bubble ( ... ) );
     BubbleID anotherID = newBubbleComm.createBubble( Bubble ( ... ) );
     // use the returned temporary bubble IDs in the BubbleField

    // The following call, updates the global bubble vector and the bubble field.
    // Here the new bubbles are appended to the vector and the temporary bubble IDs in the bubble field are
    // mapped to global IDs.
    newBubbleComm.communicateAndApply( globalBubbleVector, blockStorage, bubbleFieldID );

 *  \endcode
 *
 *
 **********************************************************************************************************************/
class NewBubbleCommunication
{
 public:
   explicit NewBubbleCommunication(uint_t nrOfExistingBubbles = uint_c(0))
   {
      // initialize the temporary bubble ID with the currently registered number of bubbles
      nextFreeBubbleID_ = numeric_cast< BubbleID >(nrOfExistingBubbles);
   }

   // add a new bubble which gets assigned a temporary bubble ID that should be stored in bubbleField; temporary IDs are
   // converted to valid global IDs in communicateAndApply() later
   BubbleID createBubble(const Bubble& newBubble)
   {
      bubblesToCreate_.push_back(newBubble);

      // get a temporary bubble ID
      return nextFreeBubbleID_++;
   }

   // get the temporary bubble ID that is returned at the next call to createBubble()
   BubbleID nextFreeBubbleID() const { return nextFreeBubbleID_; }

   /********************************************************************************************************************
    * Communicate all locally added bubbles and append them to a global bubble array. Convert the preliminary bubble IDs
    * in the bubble field to global IDs.
    *******************************************************************************************************************/
   void communicateAndApply(std::vector< Bubble >& vectorToAddBubbles, StructuredBlockStorage& blockStorage,
                            BlockDataID bubbleFieldID);

   /********************************************************************************************************************
    * Communicate all locally added bubbles and delete bubbles marked inside a boolean array. Convert the preliminary
    * bubble IDs in the bubble field to global IDs. Instead of appending all new bubbles, first all "holes" in the
    * globalVector are filled
    *
    * \Important: - the bubblesToOverwrite vector must be the same on all processes, i.e., it has to be communicated/
    *               reduced before calling this function
    *             - bubblesToCreate_.size() > number of "true"s in bubblesToOverwrite, i.e., more new bubbles have to be
    *               created than deleted
    *******************************************************************************************************************/
   void communicateAndApply(std::vector< Bubble >& vectorToAddBubbles, const std::vector< bool >& bubblesToOverwrite,
                            StructuredBlockStorage& blockStorage, BlockDataID bubbleFieldID);

 private:
   // communicate all new local bubbles and store the (MPI) gathered new bubbles in the returned object;
   // bubblesToCreate_ is cleared
   std::shared_ptr< CommunicatedNewBubbles > communicate(size_t nrOfBubblesBefore);

   BubbleID nextFreeBubbleID_; // temporary bubble ID
   std::vector< Bubble > bubblesToCreate_;
}; // class NewBubbleCommunication

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
