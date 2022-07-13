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
//! \file NewBubbleCommunication.cpp
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Communication for the creation of new bubbles.
//
//======================================================================================================================

#include "NewBubbleCommunication.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
class CommunicatedNewBubbles
{
 public:
   CommunicatedNewBubbles(size_t nrOfBubblesBefore, size_t locallyCreatedBubbles, size_t localOffset)
      : nextBubbleCtr_(0)
   {
      temporalIDToNewIDMap_.resize(locallyCreatedBubbles);
      nrOfBubblesBefore_         = numeric_cast< BubbleID >(nrOfBubblesBefore);
      nrOfLocallyCreatedBubbles_ = numeric_cast< BubbleID >(locallyCreatedBubbles);
      localOffset_               = numeric_cast< BubbleID >(localOffset);
   }

   void storeNextBubble(BubbleID newBubbleID, Bubble& bubbleToStoreTo)
   {
      // store newBubbleID in map
      if (wasNextBubbleCreatedOnThisProcess()) { temporalIDToNewIDMap_[nextBubbleCtr_ - localOffset_] = newBubbleID; }

      // store bubble
      recvBuffer_ >> bubbleToStoreTo;

      ++nextBubbleCtr_;
   }

   // return whether there are new bubbles that have to be processed
   bool hasMoreBubbles() const { return !recvBuffer_.isEmpty(); }

   // map the preliminary bubble IDs to global new bubble IDs
   void mapTemporalToNewBubbleID(BubbleID& id) const
   {
      if (id >= nrOfBubblesBefore_) { id = temporalIDToNewIDMap_[id - nrOfBubblesBefore_]; }
      // else: bubble ID is already mapped correctly
   }

   mpi::RecvBuffer& recvBuffer() { return recvBuffer_; }

 private:
   // return whether the next bubble was created on this process
   bool wasNextBubbleCreatedOnThisProcess() const
   {
      // return whether counter of current bubble is between offset and offset+numberOfLocallyCreatedBubbles
      return nextBubbleCtr_ >= localOffset_ && nextBubbleCtr_ < localOffset_ + nrOfLocallyCreatedBubbles_;
   }

   mpi::RecvBuffer recvBuffer_;

   BubbleID nrOfBubblesBefore_;         // size of global bubble array without new/deleted bubbles
   BubbleID nrOfLocallyCreatedBubbles_; // number of bubbles that were added on this process
   BubbleID localOffset_;               // the offset of the locally created bubbles in recvBuffer
   BubbleID nextBubbleCtr_;             // number of bubbles that have already been stored with storeNextBubble()

   std::vector< BubbleID > temporalIDToNewIDMap_; // maps the temporal bubble IDs to new global bubble IDs
};                                                // class CommunicatedNewBubbles

std::shared_ptr< CommunicatedNewBubbles > NewBubbleCommunication::communicate(size_t nrOfBubblesBefore)
{
   std::shared_ptr< MPIManager > mpiManager = MPIManager::instance();

   // cast number of every process' new bubbles to int (int is used in MPI_Allgather)
   int localNewBubbles = int_c(bubblesToCreate_.size());

   std::vector< int > numNewBubblesPerProcess(uint_t(mpiManager->numProcesses()), 0);

   // communicate, i.e., gather each process' number of new bubbles
   WALBERLA_MPI_SECTION()
   {
      MPI_Allgather(&localNewBubbles, 1, MPITrait< int >::type(), &numNewBubblesPerProcess[0], 1,
                    MPITrait< int >::type(), MPI_COMM_WORLD);
   }
   WALBERLA_NON_MPI_SECTION() { numNewBubblesPerProcess[0] = localNewBubbles; }

   // offsetVector[i] is the number of new bubbles created on processes with rank smaller than i
   std::vector< int > offsetVector;
   offsetVector.push_back(0);
   for (size_t i = 0; i < numNewBubblesPerProcess.size() - 1; ++i)
   {
      offsetVector.push_back(offsetVector.back() + numNewBubblesPerProcess[i]);
   }

   WALBERLA_ASSERT_EQUAL(offsetVector.size(), numNewBubblesPerProcess.size());

   // get each process' individual offset
   size_t offset = uint_c(offsetVector[uint_c(mpiManager->worldRank())]);

   // get the total number of new bubbles
   size_t numNewBubbles = uint_c(offsetVector.back() + numNewBubblesPerProcess.back());

   mpi::SendBuffer sendBuffer;

   size_t bytesPerBubble = mpi::BufferSizeTrait< Bubble >::size;

   // reserve space for a bubble's size in bytes (clean way to create send buffer)
   sendBuffer.reserve(bytesPerBubble);

   // pack local new bubble data into sendBuffer
   for (auto i = bubblesToCreate_.begin(); i != bubblesToCreate_.end(); ++i)
   {
      sendBuffer << *i;
   }

   // compute the byte size of numNewBubblesPerProcess[i] and offsetVector[i]
   for (uint_t i = uint_c(0); i < uint_c(mpiManager->numProcesses()); ++i)
   {
      numNewBubblesPerProcess[i] *= int_c(bytesPerBubble);
      offsetVector[i] *= int_c(bytesPerBubble);
   }

   // create new CommunicatedNewBubbles object to store the gathered new bubbles (see below)
   auto result = std::make_shared< CommunicatedNewBubbles >(nrOfBubblesBefore, bubblesToCreate_.size(), offset);

   // resize recvBuffer for storing the total number of new bubbles
   result->recvBuffer().resize(bytesPerBubble * numNewBubbles);

   WALBERLA_MPI_SECTION()
   {
      // communicate, i.e., gather all locally added bubbles (and store them in the new CommunicatedNewBubbles object)
      MPI_Allgatherv(sendBuffer.ptr(), int_c(sendBuffer.size()), MPI_UNSIGNED_CHAR, result->recvBuffer().ptr(),
                     &numNewBubblesPerProcess[0], &offsetVector[0], MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
   }
   WALBERLA_NON_MPI_SECTION() { result->recvBuffer() = sendBuffer; }

   // all bubbles have been "created" and vector can be cleared
   bubblesToCreate_.clear();

   return result;
}

void NewBubbleCommunication::communicateAndApply(std::vector< Bubble >& vectorToAddBubbles,
                                                 StructuredBlockStorage& blockStorage, BlockDataID bubbleFieldID)
{
   // communicate
   std::shared_ptr< CommunicatedNewBubbles > newBubbleContainer = communicate(vectorToAddBubbles.size());

   // append new bubbles to vectorToAddBubbles
   BubbleID nextBubbleID = numeric_cast< BubbleID >(vectorToAddBubbles.size());
   while (newBubbleContainer->hasMoreBubbles())
   {
      // create and append a new empty bubble
      vectorToAddBubbles.emplace_back();

      // update the empty bubble with the received bubble
      newBubbleContainer->storeNextBubble(nextBubbleID, vectorToAddBubbles.back());

      ++nextBubbleID;
   }

   // rename bubble IDs
   for (auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt)
   {
      BubbleField_T* const bubbleField = blockIt->getData< BubbleField_T >(bubbleFieldID);
      for (auto fieldIt = bubbleField->begin(); fieldIt != bubbleField->end(); ++fieldIt)
      {
         // rename an existing bubble ID in bubbleField
         if (*fieldIt != INVALID_BUBBLE_ID) { newBubbleContainer->mapTemporalToNewBubbleID(*fieldIt); }
      }
   }
}

void NewBubbleCommunication::communicateAndApply(std::vector< Bubble >& vectorToAddBubbles,
                                                 const std::vector< bool >& bubblesToOverwrite,
                                                 StructuredBlockStorage& blockStorage, BlockDataID bubbleFieldID)
{
   WALBERLA_ASSERT_EQUAL(bubblesToOverwrite.size(), vectorToAddBubbles.size());

   // communicate
   std::shared_ptr< CommunicatedNewBubbles > newBubbleContainer = communicate(vectorToAddBubbles.size());

   // overwrite existing bubbles with new bubbles; there must be at least the same number of new bubbles than in
   // bubblesToOverwrite
   for (size_t i = 0; i < bubblesToOverwrite.size(); ++i)
   {
      if (bubblesToOverwrite[i]) { newBubbleContainer->storeNextBubble(BubbleID(i), vectorToAddBubbles[i]); }
   }

   // append remaining new bubbles to vectorToAddBubbles
   BubbleID nextBubbleID = numeric_cast< BubbleID >(vectorToAddBubbles.size());
   while (newBubbleContainer->hasMoreBubbles())
   {
      // create and append a new empty bubble
      vectorToAddBubbles.emplace_back();

      // update the empty bubble with the received bubble
      newBubbleContainer->storeNextBubble(nextBubbleID, vectorToAddBubbles.back());

      ++nextBubbleID;
   }

   // rename bubble IDs
   for (auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt)
   {
      BubbleField_T* const bubbleField = blockIt->getData< BubbleField_T >(bubbleFieldID);
      for (auto fieldIt = bubbleField->begin(); fieldIt != bubbleField->end(); ++fieldIt)
      {
         // rename an existing bubble ID in bubbleField
         if (*fieldIt != INVALID_BUBBLE_ID) { newBubbleContainer->mapTemporalToNewBubbleID(*fieldIt); }
      }
   }
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
