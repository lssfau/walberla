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
//! \file GlobalBodyPresenceLevelDetermination.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockForest.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe_coupling/utility/BodySelectorFunctions.h"

namespace walberla {
namespace pe_coupling {
namespace amr {

/*
 * Class to determine the minimum level a block can be based on presence of global bodies.
 * Only the global bodies that are given by the body selection function are considered.
 * If a global body partially overlaps with a block this is block is determined to be on the finest grid.
 * The block can be extended by the given extension length (e.g. the number of ghost layers * dx).
 * This ensures correctness of the body mapping across block borders.
 * Note: Blocks that are fully contained inside a global body can have any level.
 */
class GlobalBodyPresenceLevelDetermination
{
public:

   GlobalBodyPresenceLevelDetermination( const shared_ptr<pe::BodyStorage> & globalBodyStorage,
                                         uint_t finestLevel, real_t blockExtensionLength = real_t(0),
                                         const std::function<bool(pe::BodyID)> & globalBodySelectorFct = selectAllBodies) :
         globalBodyStorage_( globalBodyStorage ), finestLevel_( finestLevel ),
         blockExtensionLength_( blockExtensionLength ),
         globalBodySelectorFct_( globalBodySelectorFct )
   {}

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > &, const BlockForest & /*forest*/ );

private:

   bool checkForPartialOverlapWithGlobalBodies(const AABB& box);

   shared_ptr<pe::BodyStorage> globalBodyStorage_;
   uint_t finestLevel_;
   real_t blockExtensionLength_;
   std::function<bool(pe::BodyID)> globalBodySelectorFct_;
};

} // namespace amr
} // namespace pe_coupling
} // namespace walberla
