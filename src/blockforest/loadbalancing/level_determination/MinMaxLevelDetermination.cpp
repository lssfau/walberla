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
//! \file MinMaxLevelDetermination.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "MinMaxLevelDetermination.h"

#include "core/logging/Logging.h"

namespace walberla {
namespace blockforest {

void MinMaxLevelDetermination::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                           std::vector< const Block * > &,
                                           const BlockForest & /*forest*/ )
{
   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      const auto infoIt = ic_->find(it->first->getId());
      WALBERLA_ASSERT_UNEQUAL( infoIt, ic_->end() );

      it->second = it->first->getLevel(); //keep everything as it is

      //check for refinement
      if (infoIt->second.computationalWeight > maxBodies_)
      {
         it->second = it->first->getLevel() + uint_t(1);
         continue;
      }

      //check for coarsening
      if ((it->first->getLevel() > 0) && (infoIt->second.computationalWeight < minBodies_))
      {
         if (getOrCreateCoarseInfo(it->first->getId())->second.computationalWeight < maxBodies_)
         {
            it->second = it->first->getLevel() - uint_t(1);
         }
         continue;
      }
   }
}

blockforest::InfoCollection::const_iterator MinMaxLevelDetermination::getOrCreateCoarseInfo( const blockforest::BlockID& id )
{
   auto fatherId = id.getFatherId();
   auto infoIt   = ic_->find( fatherId );
   if (infoIt != ic_->end()) return infoIt;

   blockforest::BlockInfo newWeight( 0, 0);
   for (uint_t child = 0; child < 8; ++child)
   {
      blockforest::BlockID childId(fatherId, child);
      auto childIt = ic_->find( childId );
      //might be not available if not all blocks are on the same level
      //return giant number to prevent coarsening
      if (childIt == ic_->end())
      {
         newWeight = blockforest::BlockInfo( std::numeric_limits<uint_t>::max(),
                                             std::numeric_limits<uint_t>::max());
         break;
      } else
      {
         newWeight += childIt->second;
      }
   }
   WALBERLA_LOG_DETAIL("creating coarse weights (" << newWeight << ")");
   return ic_->insert( std::make_pair(fatherId, newWeight) ).first;
}

} // namespace blockforest
} // namespace walberla
