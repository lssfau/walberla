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
//! \file BodyPresenceLevelDetermination.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "BodyPresenceLevelDetermination.h"

namespace walberla {
namespace pe_coupling {
namespace amr {


void BodyPresenceLevelDetermination::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                                 std::vector< const Block * > &, const BlockForest & /*forest*/ )
{
   for (auto &minTargetLevel : minTargetLevels) {
      uint_t currentLevelOfBlock = minTargetLevel.first->getLevel();

      const uint_t numberOfParticlesInDirectNeighborhood = getNumberOfLocalAndShadowBodiesInNeighborhood(minTargetLevel.first);

      uint_t targetLevelOfBlock = currentLevelOfBlock; //keep everything as it is
      if ( numberOfParticlesInDirectNeighborhood > uint_t(0) )
      {
         // set block to finest level if there are bodies nearby
         targetLevelOfBlock = finestLevel_;
      }
      else
      {
         // block could coarsen since there are no bodies nearby
         if( currentLevelOfBlock > uint_t(0) )
            targetLevelOfBlock = currentLevelOfBlock - uint_t(1);
      }

      WALBERLA_CHECK_LESS_EQUAL(std::abs(int_c(targetLevelOfBlock) - int_c(currentLevelOfBlock)), uint_t(1), "Only level difference of maximum 1 allowed!");
      minTargetLevel.second = targetLevelOfBlock;
   }
}

uint_t BodyPresenceLevelDetermination::getNumberOfLocalAndShadowBodiesInNeighborhood(const Block * block)
{
   auto numBodies = uint_t(0);

   // add bodies of current block
   const auto infoIt = infoCollection_->find(block->getId());
   WALBERLA_CHECK_UNEQUAL(infoIt, infoCollection_->end(), "Block with ID " << block->getId() << " not found in info collection!");

   numBodies += infoIt->second.numberOfLocalBodies;
   numBodies += infoIt->second.numberOfShadowBodies;

   // add bodies of all neighboring blocks
   for(uint_t i = 0; i < block->getNeighborhoodSize(); ++i)
   {
      const BlockID &neighborBlockID = block->getNeighborId(i);
      const auto infoItNeighbor = infoCollection_->find(neighborBlockID);
      WALBERLA_CHECK_UNEQUAL(infoItNeighbor, infoCollection_->end(), "Neighbor block with ID " << neighborBlockID << " not found in info collection!");

      numBodies += infoItNeighbor->second.numberOfLocalBodies;
      numBodies += infoItNeighbor->second.numberOfShadowBodies;
   }
   return numBodies;
}


} // namespace amr
} // namespace pe_coupling
} // namespace walberla
