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
//! \file AnalyticCollisionFunctions.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>

#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {
namespace collision_detection {
namespace analytic {

inline
bool detectSphereSphereCollision( const Vec3&   pos1,
                                  const real_t& radius1,
                                  const Vec3&   pos2,
                                  const real_t& radius2,
                                        Vec3&   contactPoint,
                                        Vec3&   contactNormal,
                                        real_t& penetrationDepth,
                                  const real_t& contactThreshold)
{
   contactNormal         = pos1 - pos2;
   auto contactNormalSq  = contactNormal.sqrLength();
   real_t separationDist = radius1 + radius2 + contactThreshold;

   if( contactNormalSq < separationDist * separationDist ) {
      penetrationDepth = std::sqrt(contactNormalSq) - radius1 - radius2;
      math::normalize(contactNormal);
      const real_t k( radius2 + real_c(0.5) * penetrationDepth );
      contactPoint = ( pos2 + contactNormal * k );
      return true;
   }
   return false;
}

inline
bool detectSphereHalfSpaceCollision( const Vec3&   pos1,
                                     const real_t& radius1,
                                     const Vec3&   pos2,
                                     const Vec3&   normal2,
                                           Vec3&   contactPoint,
                                           Vec3&   contactNormal,
                                           real_t& penetrationDepth,
                                     const real_t& contactThreshold)
{
   /**
    * Plane displacement from the origin.
    *
    * The displacement can be categorized in the following way:\n
    * - > 0: The global origin is inside the plane\n
    * - < 0: The global origin is outside the plane\n
    * - = 0: The global origin is on the surface of the plane
   **/
   const real_t d = math::dot(normal2, pos2);

   const real_t k = math::dot(normal2, pos1);
   penetrationDepth = ( k - radius1 - d );

   if( penetrationDepth < contactThreshold )
   {
      contactPoint = ( pos1 - ( radius1 + penetrationDepth ) * normal2 );
      contactNormal = normal2;
      return true;
   }
   return false;
}


} //namespace collision_detection
} //namespace analytic
} //namespace mesa_pd
} //namespace walberla
