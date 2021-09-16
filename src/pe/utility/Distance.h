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
//! \file Distance.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>

#include "core/DataTypes.h"
#include "core/math/Limits.h"

#include "pe/Types.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"


namespace walberla {
namespace pe {

//=================================================================================================
//
//  DISTANCE CALCULATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
//* These functions can be used to calculate the distance between two rigid bodies (surface to surface).
//*************************************************************************************************


//*************************************************************************************************
/*!\name Distance calculation functions */
//@{
inline real_t getSurfaceDistance( ConstSphereID s1  , ConstSphereID s2   );
inline real_t getSurfaceDistance( ConstSphereID s   , ConstPlaneID p     );
inline real_t getSurfaceDistance( ConstPlaneID p1   , ConstPlaneID p2    );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Distance calculation between two Sphere primitives.
 *
 * \param s1 The first sphere.
 * \param s2 The second sphere.
 * \return The minimum distance between the two spheres.
 *
 * This function returns the distance between the two spheres \a s1 and \a s2. Note that a
 * negative distance indicates that the two spheres are overlapping.
 */
inline real_t getSurfaceDistance( ConstSphereID s1, ConstSphereID s2 )
{
   const Vec3 normal( s1->getPosition() - s2->getPosition() );
   return normal.length() - s1->getRadius() - s2->getRadius();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Distance calculation between a Sphere and a Plane.
 *
 * \param s The sphere.
 * \param p The plane.
 * \return The minimum distance between the sphere and the plane.
 *
 * This function returns the distance between the sphere \a s and the plane \a p. Note that
 * a negative distance indicates that the two bodies are overlapping.
 */
inline real_t getSurfaceDistance( ConstSphereID s, ConstPlaneID p )
{
   return (  s->getPosition() * p->getNormal() ) - s->getRadius() - p->getDisplacement();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Distance calculation between two Plane primitives.
 *
 * \param p1 The first plane.
 * \param p2 The second plane.
 * \return The minimum distance between the two planes.
 *
 * This function returns the distance between the two planes \a p1 and \a p2. Note that in case
 * the two (infinite) planes overlap the returned distance is -math::Limits<real_t>::inf().
 */
inline real_t getSurfaceDistance( ConstPlaneID p1, ConstPlaneID p2 )
{
   if( (  p1->getNormal() * p2->getNormal() ) <= real_t(1E-12-1.0) )
      return std::fabs( p1->getDisplacement() - p2->getDisplacement() );
   else return -math::Limits<real_t>::inf();
}
//*************************************************************************************************

} // namespace pe
} // namespace walberla
