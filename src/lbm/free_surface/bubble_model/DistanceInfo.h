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
//! \file DistanceInfo.h
//! \ingroup bubble_model
//! \author Daniela Anderl
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute the distances between bubbles.
//
//======================================================================================================================

#pragma once

#include "core/mpi/BufferDataTypeExtensions.h"

#include "field/GhostLayerField.h"

#include <map>

#include "BubbleDefinitions.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
 * Compute the distances from other bubbles in a diffusive manner.
 *
 * Each bubble spreads its distance until a defined maximum distance is reached. Cells that are within this
 * maximum distance from one ore more bubbles know the distance from itself to each of these bubbles (or the single
 * bubble).
 *
 **********************************************************************************************************************/
class DistanceInfo
{
 public:
   DistanceInfo() = default;

   // get the distance to the bubble that is closest to bubble "ownID" within the range of maxDistance
   real_t getDistanceToNearestBubble(const BubbleID& ownID, real_t maxDistance) const;

   // combine "bubbleInfos" (std::map) of the current cell with "bubbleInfos" of another cell; the distance between
   // these two cells is linkDistance:
   // - if the other cell's bubble is not yet registered, insert it and calculate the distance
   // - if the other cell's bubble is already registered, update the distance
   void combine(const DistanceInfo& other, real_t linkDistance, real_t maxDistance);

   // set the distance to bubble "id" to zero; create an entry for this bubble, if non yet exists
   void setToZero(bubble_model::BubbleID id);

   // remove any entry in "bubbleInfos" with distance <= 0
   void removeZeroDistanceFromLiquid();

   bool operator==(DistanceInfo& other);

   void clear() { bubbleInfos_.clear(); }

   // print the content of bubbleInfos
   void print(std::ostream& os) const;

 protected:
   friend mpi::SendBuffer& operator<<(mpi::SendBuffer&, const DistanceInfo&);
   friend mpi::RecvBuffer& operator>>(mpi::RecvBuffer&, DistanceInfo&);

   std::map< bubble_model::BubbleID, real_t > bubbleInfos_;
}; // class DistanceInfo

mpi::SendBuffer& operator<<(mpi::SendBuffer& buf, const DistanceInfo& di);
mpi::RecvBuffer& operator>>(mpi::RecvBuffer& buf, DistanceInfo& di);

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

namespace walberla
{
namespace mpi
{
template<>
struct BufferSizeTrait< free_surface::bubble_model::DistanceInfo >
{
   static const bool constantSize = false;
};
} // namespace mpi
} // namespace walberla
