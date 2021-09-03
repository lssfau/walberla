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
//! \file IntersectionRatio.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "geometry/bodies/Sphere.h"

#include <exception>
#include <limits>
#include <string>

namespace walberla {
namespace lbm {



template <typename Body>
real_t intersectionRatioBisection( const Body & body,
                                   const Vector3<real_t> & fluidPoint,
                                   const Vector3<real_t> & direction,
                                   const real_t epsilon );



real_t intersectionRatioSphere( const geometry::Sphere & sphere,
                                const Vector3<real_t> & fluidPoint,
                                const Vector3<real_t> & direction );


//*******************************************************************************************************************
/*!
* \brief Computes the intersection of a ray segment with a body surface
*
* Let P_i be the intersection point
*
* \param body       surface
* \param fluidPoint (P_f), start point of ray
* \param direction  (d) of ray (length of d is the length of the ray segment)
* \param epsilon    abortion criterion for iterative methods. Epsilon relates to the distance of P_i to the surface
* \return Intersection ratio: | P_i - P_f | / | d |
*/
//*******************************************************************************************************************
template <typename Body>
real_t intersectionRatio( const Body & body,
                          const Vector3<real_t> & fluidPoint,
                          const Vector3<real_t> & direction,
                          const real_t epsilon )
{
   return intersectionRatioBisection( body, fluidPoint, direction, epsilon );
}



inline real_t intersectionRatio( const geometry::Sphere & sphere,
                                 const Vector3<real_t> & fluidPoint,
                                 const Vector3<real_t> & direction,
                                 const real_t /*epsilon*/ = real_t(0) )
{
   return intersectionRatioSphere( sphere, fluidPoint, direction );
}






} // namespace lbm
} // namespace walberla


#include "IntersectionRatio.impl.h"


