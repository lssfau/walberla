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
//! \file DistanceInfo.cpp
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute the distances between bubbles.
//
//======================================================================================================================

#include "DistanceInfo.h"

#include "core/mpi/BufferDataTypeExtensions.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
mpi::SendBuffer& operator<<(mpi::SendBuffer& buf, const DistanceInfo& di)
{
   buf << di.bubbleInfos_;
   return buf;
}

mpi::RecvBuffer& operator>>(mpi::RecvBuffer& buf, DistanceInfo& di)
{
   di.clear();
   buf >> di.bubbleInfos_;
   return buf;
}

real_t DistanceInfo::getDistanceToNearestBubble(const BubbleID& ownID, real_t maxDistance) const
{
   // limit search region to maxDistance
   real_t minDistance = maxDistance;

   // find the minimum distance
   for (auto it = bubbleInfos_.begin(); it != bubbleInfos_.end(); it++)
   {
      if (minDistance > it->second && it->first != ownID) { minDistance = it->second; }
   }

   return minDistance;
}

void DistanceInfo::combine(const DistanceInfo& other, real_t linkDistance, real_t maxDistance)
{
   // sum of another cell's distance to the nearest bubble and linkDistance
   real_t sumDistance = real_c(0);

   // loop over the "bubbleInfos" (std::map) of another cell
   for (auto it = other.bubbleInfos_.begin(); it != other.bubbleInfos_.end(); it++)
   {
      // the other cell's bubble is not in contained in this cell's bubbleInfos, yet
      if (bubbleInfos_.find(it->first) == bubbleInfos_.end())
      {
         sumDistance = it->second + linkDistance;

         // add an entry for the other cell's nearest bubble, if the summed distance is below maxDistance
         if (sumDistance < maxDistance) { bubbleInfos_[it->first] = sumDistance; }
      }
      else // the other cell's bubble is contained in this cell's bubbleInfos
      {
         sumDistance   = it->second + linkDistance;
         real_t& entry = bubbleInfos_[it->first];

         // update the distance for this bubble if it is less than the currently known distance
         if (entry > sumDistance) { entry = sumDistance; }
      }
   }
}

void DistanceInfo::setToZero(bubble_model::BubbleID id)
{
   // "bubbleInfos" is of type std::map, if no entry with key "id" exists, a new entry is created
   bubbleInfos_[id] = real_c(0.0);
}

void DistanceInfo::removeZeroDistanceFromLiquid()
{
   for (auto it = bubbleInfos_.begin(); it != bubbleInfos_.end(); it++)
      if (it->second <= real_c(0.0)) { bubbleInfos_.erase(it); }
}

bool DistanceInfo::operator==(DistanceInfo& other) { return (this->bubbleInfos_ == other.bubbleInfos_); }

void DistanceInfo::print(std::ostream& os) const
{
   for (auto it = bubbleInfos_.begin(); it != bubbleInfos_.end(); it++)
   {
      os << "( " << it->first << "," << it->second << ")\n";
   }
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
