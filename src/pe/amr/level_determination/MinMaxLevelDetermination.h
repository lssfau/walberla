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
//! \file MinMaxLevelDetermination.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <pe/Types.h>
#include <pe/amr/InfoCollection.h>

#include <blockforest/Block.h>
#include <blockforest/BlockForest.h>
#include <core/logging/Logging.h>
#include <domain_decomposition/BlockDataID.h>

namespace walberla {
namespace pe {
namespace amr {

class MinMaxLevelDetermination
{
public:

   MinMaxLevelDetermination( const shared_ptr<blockforest::InfoCollection>& ic,
                             const size_t minBodies,
                             const size_t maxBodies) :
      ic_( ic ), minBodies_(minBodies), maxBodies_(maxBodies)
   {}

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > &,
                    const BlockForest & forest );

public:
   const shared_ptr<blockforest::InfoCollection> ic_;
   size_t      minBodies_;
   size_t      maxBodies_;

   blockforest::InfoCollection::const_iterator getOrCreateCoarseInfo( const blockforest::BlockID& id );
};

} // namespace amr
} // namespace pe
} // namespace walberla
