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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/Flags.h>

namespace walberla {
namespace mesa_pd {
namespace collision_detection {

/**
 * Checks if particle idx1 is within interaction distance with particle idx2.
 * Collisions of infinite particles are also filtered.
 */
template <typename Accessor>
bool isInInteractionDistance(const size_t idx1, const size_t idx2, Accessor& ac)
{
   using namespace data::particle_flags;
   if (isSet(ac.getFlags(idx1), INFINITE)) return true;
   if (isSet(ac.getFlags(idx2), INFINITE)) return true;
   auto separationDist = ac.getInteractionRadius(idx1) + ac.getInteractionRadius(idx2);
   auto realDist2 = (ac.getPosition(idx1) - ac.getPosition(idx2)).sqrLength();
   return realDist2 < separationDist*separationDist;
}

} //namespace collision_detection
} //namespace mesa_pd
} //namespace walberla
