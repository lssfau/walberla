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
//! \file AABBBody.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "BodyOverlapFunctions.h"
#include "core/math/AABB.h"


namespace walberla {
namespace geometry {


template<>
inline real_t overlapFraction ( const AABB & body, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx, uint_t )
{
   AABB box = AABB::createFromMinMaxCorner( cellMidpoint[0] - real_t(0.5)*dx[0], cellMidpoint[1] - real_t(0.5)*dx[1], cellMidpoint[2] - real_t(0.5)*dx[2],
                                            cellMidpoint[0] + real_t(0.5)*dx[0], cellMidpoint[1] + real_t(0.5)*dx[1], cellMidpoint[2] + real_t(0.5)*dx[2]);

   return body.intersectionVolume( box ) / ( dx[0] * dx[1] * dx[2] );
}


template<>
inline FastOverlapResult fastOverlapCheck ( const AABB & body, const AABB & cell )
{
   if ( body.contains( cell ) )
      return CONTAINED_INSIDE_BODY;
   if ( ! body.intersects( cell ) )
      return COMPLETELY_OUTSIDE;

   return DONT_KNOW;
}


template<>
inline bool contains ( const AABB & aabb, const Vector3<real_t> & point )
{
   return aabb.contains( point );
}


} // namespace geometry
} // namespace walberla


