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
//! \file AABBConversion.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/math/AABB.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/IAccessor.h>

#include <limits>

namespace walberla {
namespace mesa_pd {

math::AABB getAABBFromInteractionRadius(const Vector3<real_t> & pos, const real_t interactionRadius )
{
   return math::AABB( pos[0]-interactionRadius, pos[1]-interactionRadius, pos[2]-interactionRadius,
                      pos[0]+interactionRadius, pos[1]+interactionRadius, pos[2]+interactionRadius );
}


// requires: interactionRadius
// flags: infinite
template<typename ParticleAccessor_T>
math::AABB getParticleAABB(const size_t particleIdx, const ParticleAccessor_T& ac)
{
   static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

   if( mesa_pd::data::particle_flags::isSet( ac.getFlags(particleIdx), mesa_pd::data::particle_flags::INFINITE) )
   {
      auto inf = std::numeric_limits<real_t>::infinity();
      return math::AABB( -inf, -inf, -inf,
                          inf,  inf,  inf);
   }
   else
   {
      return getAABBFromInteractionRadius(ac.getPosition(particleIdx), ac.getInteractionRadius(particleIdx));
   }
}

} //namespace mesa_pd
} //namespace walberla
