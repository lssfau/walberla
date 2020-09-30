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
//! \file PeBodyOverlapFunctions.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Implementation of BodyOverlapFunctions for the pe bodies
//
//======================================================================================================================

#pragma once

#include "geometry/bodies/BodyOverlapFunctions.h"

#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Squirmer.h"

namespace walberla {
namespace geometry {

/////////////////////////////////////////////////
// specialization for pe::Sphere               //
// duplication from geometry/bodies/Sphere.cpp //
/////////////////////////////////////////////////

template<> inline FastOverlapResult fastOverlapCheck( const pe::Sphere & peSphere, const AABB & box )
{
   if ( ! peSphere.getAABB().intersects( box ) )
   {
      return COMPLETELY_OUTSIDE;
   }
   else
   {
      Vector3<real_t> midPoint = peSphere.getPosition();
      const real_t oneOverSqrt3 = real_t(1) / std::sqrt( real_t(3) );
      real_t halfBoxEdge = peSphere.getRadius() * oneOverSqrt3;
      AABB innerBox = AABB::createFromMinMaxCorner( midPoint[0] - halfBoxEdge, midPoint[1] - halfBoxEdge, midPoint[2] - halfBoxEdge,
                                                    midPoint[0] + halfBoxEdge, midPoint[1] + halfBoxEdge, midPoint[2] + halfBoxEdge );
      if( innerBox.contains(box) )
      {
         return CONTAINED_INSIDE_BODY;
      }
   }

   return DONT_KNOW;
}

template<> inline FastOverlapResult fastOverlapCheck( const pe::Sphere & peSphere, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx )
{
   const real_t sqrt3half = std::sqrt( real_t(3) ) / real_t(2);

   const real_t midPointDistSq = (peSphere.getPosition() - cellMidpoint).sqrLength();

   const real_t dxMax = dx.max();
   // Check against outer circle of box
   const real_t dist1 = peSphere.getRadius() + sqrt3half * dxMax;
   if ( midPointDistSq > dist1 * dist1 )
      return COMPLETELY_OUTSIDE;

   // Check against inner circle of box
   const real_t dist2 = peSphere.getRadius() - sqrt3half * dxMax;
   if ( midPointDistSq < dist2 * dist2 )
      return CONTAINED_INSIDE_BODY;

   return DONT_KNOW;
}

template<> inline bool contains( const pe::Sphere & peSphere, const Vector3<real_t> & point )
{
   return peSphere.containsPoint( point[0], point[1], point[2] );
}

template<> inline FastOverlapResult fastOverlapCheck( const pe::Squirmer & peSquirmer, const AABB & box )
{
   return fastOverlapCheck( static_cast<const pe::Sphere &>(peSquirmer), box );
}

template<> inline FastOverlapResult fastOverlapCheck( const pe::Squirmer & peSquirmer, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx )
{
   return fastOverlapCheck( static_cast<const pe::Sphere &>(peSquirmer), cellMidpoint, dx );
}

template<> inline bool contains( const pe::Squirmer & peSquirmer, const Vector3<real_t> & point )
{
   return contains( static_cast<const pe::Sphere &>(peSquirmer), point );
}

/////////////////////////////////////////////////
// specialization for pe::Plane                //
/////////////////////////////////////////////////

template<> inline FastOverlapResult fastOverlapCheck( const pe::Plane & pePlane, const AABB & box )
{
   if ( ! pePlane.getAABB().intersects( box ) )
   {
      // if axis aligned plane, its AABB is not inf
      return COMPLETELY_OUTSIDE;
   }

   uint_t numberOfContainedCorners( 0 );
   for( const Vector3<real_t> & aabbCorner : box.corners() )
   {
      if( pePlane.containsPoint(aabbCorner))
      {
         ++numberOfContainedCorners;
      }
   }

   if( numberOfContainedCorners == uint_t(0) )
   {
      return COMPLETELY_OUTSIDE;
   }
   if( numberOfContainedCorners == uint_t(8) )
   {
      return CONTAINED_INSIDE_BODY;
   }
   // else
   return PARTIAL_OVERLAP;

}

template<> inline FastOverlapResult fastOverlapCheck( const pe::Plane & pePlane, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx )
{
   AABB box = AABB::createFromMinMaxCorner( cellMidpoint[0] - real_t(0.5)*dx[0], cellMidpoint[1] - real_t(0.5)*dx[1], cellMidpoint[2] - real_t(0.5)*dx[2],
                                            cellMidpoint[0] + real_t(0.5)*dx[0], cellMidpoint[1] + real_t(0.5)*dx[1], cellMidpoint[2] + real_t(0.5)*dx[2]);
   return fastOverlapCheck(pePlane, box);
}

template<> inline bool contains( const pe::Plane & pePlane, const Vector3<real_t> & point )
{
   return pePlane.containsPoint( point[0], point[1], point[2] );
}

////////////////////////////////////
// general pe body implementation //
////////////////////////////////////

template<> inline FastOverlapResult fastOverlapCheck( const pe::RigidBody & peBody, const AABB & box )
{
   if ( ! peBody.getAABB().intersects( box ) )
   {
      return COMPLETELY_OUTSIDE;
   }
   else
   {
      return DONT_KNOW;
   }
}

template<> inline FastOverlapResult fastOverlapCheck( const pe::RigidBody & peBody, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx )
{

   AABB box = AABB::createFromMinMaxCorner( cellMidpoint[0] - real_t(0.5)*dx[0], cellMidpoint[1] - real_t(0.5)*dx[1], cellMidpoint[2] - real_t(0.5)*dx[2],
                                            cellMidpoint[0] + real_t(0.5)*dx[0], cellMidpoint[1] + real_t(0.5)*dx[1], cellMidpoint[2] + real_t(0.5)*dx[2]);

   if ( ! peBody.getAABB().intersects( box ) )
   {
      return COMPLETELY_OUTSIDE;
   }
   else
   {
      return DONT_KNOW;
   }
}

template<> inline bool contains( const pe::RigidBody & peBody, const Vector3<real_t> & point )
{
   return peBody.containsPoint( point[0], point[1], point[2] );
}

} // namespace geometry
} // namespace walberla
