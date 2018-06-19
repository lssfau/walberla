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
//! \file BodyPresenceLevelDetermination.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockForest.h"
#include "pe_coupling/amr/InfoCollection.h"

namespace walberla {
namespace pe_coupling {
namespace amr {

/*
 * Class to determine the minimum level a block can be.
 * For coupled LBM-PE simulations the following rules apply:
 *  - a moving body will always remain on the finest block
 *  - a moving body is not allowed to extend into an area with a coarser block
 *  - if no moving body is present, the level can be as coarse as possible (restricted by the 2:1 rule)
 * Therefore, if a body, local or remote (due to bodies that are larger than a block), is present on any of the
 * neighboring blocks of a certain block, this block's target level is the finest level.
 * This, together with a refinement checking frequency that depends on the maximum translational body velocity,
 * ensures the above given requirements.
 */
class BodyPresenceLevelDetermination
{
public:

   BodyPresenceLevelDetermination( const shared_ptr<pe_coupling::InfoCollection> & infoCollection, uint_t finestLevel) :
         infoCollection_( infoCollection ), finestLevel_( finestLevel)
   {}

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > &, const BlockForest & /*forest*/ );

private:

   uint_t getNumberOfLocalAndShadowBodiesInNeighborhood(const Block * block);

   shared_ptr<pe_coupling::InfoCollection> infoCollection_;
   uint_t finestLevel_;
};

} // namespace amr
} // namespace pe_coupling
} // namespace walberla
