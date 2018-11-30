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
//! \file GlobalBodyPresenceLevelDetermination.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "GlobalBodyPresenceLevelDetermination.h"

#include "pe_coupling/geometry/PeOverlapFraction.h"

namespace walberla {
namespace pe_coupling {
namespace amr {

void GlobalBodyPresenceLevelDetermination::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                                       std::vector< const Block * > &, const BlockForest & /*forest*/ )
{
   for (auto &minTargetLevel : minTargetLevels) {
      uint_t currentLevelOfBlock = minTargetLevel.first->getLevel();

      auto blockExtendedAABB = minTargetLevel.first->getAABB().getExtended(blockExtensionLength_);
      bool blockPartiallyOverlapsWithGlobalBodies = checkForPartialOverlapWithGlobalBodies(blockExtendedAABB);

      uint_t targetLevelOfBlock = currentLevelOfBlock; //keep everything as it is
      if ( blockPartiallyOverlapsWithGlobalBodies )
      {
         // set block to finest level since overlap with at least one global body is present
         targetLevelOfBlock = finestLevel_;
      }
      else
      {
         // block could coarsen since there are no overlaps with global bodies
         if( currentLevelOfBlock > uint_t(0) )
            targetLevelOfBlock = currentLevelOfBlock - uint_t(1);
      }

      WALBERLA_ASSERT_LESS_EQUAL(std::abs(int_c(targetLevelOfBlock) - int_c(currentLevelOfBlock)), uint_t(1), "Only level difference of maximum 1 allowed!");
      minTargetLevel.second = targetLevelOfBlock;
   }
}

bool GlobalBodyPresenceLevelDetermination::checkForPartialOverlapWithGlobalBodies(const AABB& box)
{
   const Vector3<real_t> boxMidPoint( box.min() + real_t(0.5) * box.sizes());
   const Vector3<real_t> dxVec( box.sizes() );
   const auto maxDepthSuperSampling = uint_t(2);

   bool partialOverlapWithAllBodies = false;

   for( auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt )
   {
      auto bodyID = bodyIt.getBodyID();
      if( globalBodySelectorFct_(bodyID))
      {
         real_t overlapFraction = pe_coupling::overlapFractionPe(*bodyID, boxMidPoint, dxVec, maxDepthSuperSampling);

         // check for partial overlap
         if( overlapFraction > real_t(0) && overlapFraction < real_t(1))
         {
            partialOverlapWithAllBodies = true;
         }
         // if fully contained in at least one body, no partial overlap possible
         if( floatIsEqual(overlapFraction, real_t(1) ))
         {
            partialOverlapWithAllBodies = false;
            break;
         }
      }
   }
   return partialOverlapWithAllBodies;
}

} // namespace amr
} // namespace pe_coupling
} // namespace walberla
