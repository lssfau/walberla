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
//! \file MapPointToPeriodicDomain.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"

#include <array>



namespace walberla {
namespace domain_decomposition {



//**********************************************************************************************************************
/*!
*   This function can be used to transform any point in 3D space into the periodic simulation space. For example, if the
*   simulation is periodic in x direction and the simulation domain spans from x = 0 to x = 10, then a point located at
*   x = 38 is mapped to x = 8, and a point located at x = -13 is mapped to x = 7.
*   The min points of the domain are included in the simulation space, the max points are excluded!
*/
//**********************************************************************************************************************
void mapPointToPeriodicDomain( const std::array< bool, 3 > & periodic, const AABB & domain, real_t & x, real_t & y, real_t & z );



/// see documentation of 'void mapPointToPeriodicDomain( const std::array< bool, 3 > & periodic, const AABB & domain, real_t & x, real_t & y, real_t & z )'
inline void mapPointToPeriodicDomain( const std::array< bool, 3 > & periodic, const AABB & domain, Vector3< real_t > & p )
{
   mapPointToPeriodicDomain( periodic, domain, p[0], p[1], p[2] );
}



/// see documentation of 'void mapPointToPeriodicDomain( const std::array< bool, 3 > & periodic, const AABB & domain, real_t & x, real_t & y, real_t & z )'
inline Vector3< real_t > mapPointToPeriodicDomain( const std::array< bool, 3 > & periodic, const AABB & domain, const Vector3< real_t > & p )
{
   Vector3< real_t > point( p );
   mapPointToPeriodicDomain( periodic, domain, point[0], point[1], point[2] );
   return point;
}



} // namespace domain_decomposition
} // namespace walberla
