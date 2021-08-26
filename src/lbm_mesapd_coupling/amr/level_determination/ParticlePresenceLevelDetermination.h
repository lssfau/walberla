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
//! \file ParticlePresenceLevelDetermination.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockForest.h"
#include "lbm_mesapd_coupling/amr/InfoCollection.h"

namespace walberla {
namespace lbm_mesapd_coupling {
namespace amr
{
/*
 * Class to determine the minimum level a block can be.
 * For coupled LBM-PE simulations the following rules apply:
 *  - a moving particle will always remain on the finest block
 *  - a moving particle is not allowed to extend into an area with a coarser block
 *  - if no moving particle is present, the level can be as coarse as possible (restricted by the 2:1 rule)
 * Therefore, if a particle, local or remote (due to particles that are larger than a block), is present on any of the
 * neighboring blocks of a certain block, this block's target level is the finest level.
 * This, together with a refinement checking frequency that depends on the maximum translational particle velocity,
 * ensures the above given requirements.
 */
class ParticlePresenceLevelDetermination
{
 public:
   ParticlePresenceLevelDetermination(const shared_ptr< InfoCollection > infoCollection,
                                      uint_t finestLevel)
      : infoCollection_(infoCollection), finestLevel_(finestLevel)
   {}

   void operator()(std::vector< std::pair< const Block*, uint_t > >& minTargetLevels, std::vector< const Block* >&,
                   const BlockForest& /*forest*/)
   {
      for (auto& minTargetLevel : minTargetLevels)
      {
         uint_t currentLevelOfBlock = minTargetLevel.first->getLevel();

         const uint_t numberOfParticlesInDirectNeighborhood =
            getNumberOfLocalAndShadowParticlesInNeighborhood(minTargetLevel.first);

         uint_t targetLevelOfBlock = currentLevelOfBlock; // keep everything as it is
         if (numberOfParticlesInDirectNeighborhood > uint_t(0))
         {
            // set block to finest level if there are particles nearby
            targetLevelOfBlock = finestLevel_;
         }
         else
         {
            // block could coarsen since there are no particles nearby
            if (currentLevelOfBlock > uint_t(0)) targetLevelOfBlock = currentLevelOfBlock - uint_t(1);
         }

         WALBERLA_CHECK_LESS_EQUAL(std::abs(int_c(targetLevelOfBlock) - int_c(currentLevelOfBlock)), uint_t(1),
                                   "Only level difference of maximum 1 allowed!");
         minTargetLevel.second = targetLevelOfBlock;
      }
   }

 private:
   uint_t getNumberOfLocalAndShadowParticlesInNeighborhood(const Block* block)
   {
      auto numParticles = uint_t(0);

      // add particles of current block
      const auto infoIt = infoCollection_->find(block->getId());
      WALBERLA_CHECK_UNEQUAL(infoIt, infoCollection_->end(),
                             "Block with ID " << block->getId() << " not found in info collection!");

      numParticles += infoIt->second.numberOfLocalParticles + infoIt->second.numberOfGhostParticles;

      // add particles of all neighboring blocks
      for (uint_t i = 0; i < block->getNeighborhoodSize(); ++i)
      {
         const BlockID& neighborBlockID = block->getNeighborId(i);
         const auto infoItNeighbor      = infoCollection_->find(neighborBlockID);
         WALBERLA_CHECK_UNEQUAL(infoItNeighbor, infoCollection_->end(),
                                "Neighbor block with ID " << neighborBlockID << " not found in info collection!");

         numParticles += infoItNeighbor->second.numberOfLocalParticles + infoItNeighbor->second.numberOfGhostParticles;
      }
      return numParticles;
   }

   shared_ptr< InfoCollection > infoCollection_;
   uint_t finestLevel_;
};

} // namespace amr
} // namespace lbm_mesapd_coupling
} // namespace walberla
