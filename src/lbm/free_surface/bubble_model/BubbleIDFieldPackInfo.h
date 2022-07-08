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
//! \file BubbleIDFieldPackInfo.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Pack/unpack information for a field containing bubble IDs.
//
//======================================================================================================================

#pragma once

#include "field/communication/PackInfo.h"

#include "stencil/D3Q27.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
 * Pack/unpack information for a field containing bubble IDs.
 *
 * The bubble ID field requires a special pack info, since whenever the ghost layers are updated, the bubble model has
 * to look for possible bubble merges.
 *
 * This could also be implemented with a regular FieldPackInfo and an immediate loop over the ghost layer only. However,
 * it is more efficient by directly looking for bubble merges while unpacking the ghost layer, since the elements
 * do not have to be loaded twice.
 ***********************************************************************************************************************/
template< typename Stencil_T >
class BubbleIDFieldPackInfo : public field::communication::PackInfo< GhostLayerField< BubbleID, 1 > >
{
 public:
   using CommunicationStencil_T =
      typename std::conditional< Stencil_T::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;

   using field_t = GhostLayerField< BubbleID, 1 >;

   BubbleIDFieldPackInfo(const BlockDataID& bdId, MergeInformation* mergeInfo)
      : field::communication::PackInfo< field_t >(bdId), mergeInfo_(mergeInfo)
   {}

   bool threadsafeReceiving() const override { return false; }

   void unpackData(IBlock* receiver, stencil::Direction dir, mpi::RecvBuffer& buffer) override
   {
      field_t* const field = receiver->getData< field_t >(this->bdId_);
      WALBERLA_ASSERT_NOT_NULLPTR(field);

#ifndef NDEBUG
      uint_t xSize;
      uint_t ySize;
      uint_t zSize;
      buffer >> xSize >> ySize >> zSize;
      WALBERLA_ASSERT_EQUAL(xSize, field->xSize());
      WALBERLA_ASSERT_EQUAL(ySize, field->ySize());
      WALBERLA_ASSERT_EQUAL(zSize, field->zSize());
#endif

      for (auto fieldIt = field->beginGhostLayerOnly(dir); fieldIt != field->end(); ++fieldIt)
      {
         // update ghost layer with received values
         buffer >> *fieldIt;

         // look for bubble merges with bubbles from other blocks, i.e., analyze the just received ghost layer
         lookForMerges(fieldIt, dir, field);
      }
   }

   // communicate bubble IDs locally (between blocks on the same process) and immediately check for bubble merges
   void communicateLocal(const IBlock* sender, IBlock* receiver, stencil::Direction dir) override
   {
      // get sender and receiver fields
      const field_t* const senderField = sender->getData< const field_t >(this->bdId_);
      field_t* const receiverField     = receiver->getData< field_t >(this->bdId_);

      WALBERLA_ASSERT_EQUAL(senderField->xSize(), receiverField->xSize());
      WALBERLA_ASSERT_EQUAL(senderField->ySize(), receiverField->ySize());
      WALBERLA_ASSERT_EQUAL(senderField->zSize(), receiverField->zSize());

      auto srcIt = senderField->beginSliceBeforeGhostLayer(dir); // iterates only over last slice before ghost layer
      auto dstIt = receiverField->beginGhostLayerOnly(stencil::inverseDir[dir]); // iterates only over ghost layer

      while (srcIt != senderField->end())
      {
         // fill receiver's ghost layer with values from the sender's outermost inner layer
         *dstIt = *srcIt;

         // look for bubble merges with bubbles from other blocks, i.e., analyze the just received ghost layer
         lookForMerges(dstIt, stencil::inverseDir[dir], receiverField);

         ++srcIt;
         ++dstIt;
      }

      WALBERLA_ASSERT(srcIt == senderField->end() && dstIt == receiverField->end());
   }

 protected:
   // looks for bubble merges with bubbles from other blocks from the view of a ghost layer; argument "iterator i"
   // should iterate ghost layer only
   void lookForMerges(const field_t::iterator& i, stencil::Direction dir, const field_t* field)
   {
      using namespace stencil;

      // only iterate the relevant neighborhood, for example:
      // in ghost layer at "W", check all neighbors in directions containing "E"
      // in ghost layer at "NE", check all neighbors in directions containing "SW" (=> SW, TSW and BSW)
      // in ghost layer at "TNE", only check the neighbor in direction "BSW"
      for (uint_t d = uint_c(0); d < uint_c(CommunicationStencil_T::d_per_d_length[inverseDir[dir]]); ++d)
      {
         const Direction neighborDirection = CommunicationStencil_T::d_per_d[inverseDir[dir]][d];
         auto neighborVal                  = i.neighbor(neighborDirection);

         // merge only occurs if bubble IDs from both blocks are valid and different
         if (neighborVal != INVALID_BUBBLE_ID && *i != INVALID_BUBBLE_ID && neighborVal != *i)
         {
            Cell neighborCell = i.cell() + Cell(cx[neighborDirection], cy[neighborDirection], cz[neighborDirection]);

            // make sure that the neighboring cell is not in a different ghost layer but part of the domain
            if (field->isInInnerPart(neighborCell)) { mergeInfo_->registerMerge(*i, neighborVal); }
         }
      }
   }

   MergeInformation* mergeInfo_;
}; // class BubbleIDFieldPackInfo

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
