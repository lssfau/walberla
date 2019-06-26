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
bool detectSphereBoxCollision( const Vec3&   pos1,
                               const real_t& radius1,
                               const Vec3&   pos2,
                               const Vec3&   edgeLength2,
                               const Rot3&   rot2,
                               Vec3&   contactPoint,
                               Vec3&   contactNormal,
                               real_t& penetrationDepth,
                               const real_t& contactThreshold)
{
   //WALBERLA_ABORT(s << b);
   const Vec3  l( real_t(0.5)*edgeLength2 );  // Box half side lengths
   const Mat3& R( rot2.getMatrix() );         // Rotation of the box

   bool outside( false );

   // Calculating the distance between the sphere and the box
   const Vec3 d( pos1 - pos2 );

   // Calculating the the sphere/box distance in the box's frame of reference
   Vec3 p( d[0]*R[0] + d[1]*R[3] + d[2]*R[6],
         d[0]*R[1] + d[1]*R[4] + d[2]*R[7],
         d[0]*R[2] + d[1]*R[5] + d[2]*R[8] );

   // Projection of the x-component
   if     ( p[0] < -l[0] ) { p[0] = -l[0]; outside = true; }
   else if( p[0] >  l[0] ) { p[0] =  l[0]; outside = true; }

   // Projection of the y-component
   if     ( p[1] < -l[1] ) { p[1] = -l[1]; outside = true; }
   else if( p[1] >  l[1] ) { p[1] =  l[1]; outside = true; }

   // Projection of the z-component
   if     ( p[2] < -l[2] ) { p[2] = -l[2]; outside = true; }
   else if( p[2] >  l[2] ) { p[2] =  l[2]; outside = true; }

   // Special treatment if the sphere's center of mass is inside the box
   // In this case, a contact at the sphere's center of mass with a normal away from
   // the closest face of the box is generated.
   if( !outside )
   {
      real_t dist( std::fabs(p[0]) - l[0] );
      size_t face( 0 );

      for( size_t i=1; i<3; ++i ) {
         const real_t tmp( std::fabs(p[i]) - l[i] );
         if( dist < tmp ) {
            dist = tmp;
            face = i;
         }
      }

      // Calculating the normal direction of the contact
      Vec3 n;
      n[face] = math::sign( p[face] );
      n = R * n;

      // Calculating the contact distance
      dist -= radius1;

      // Creating a single contact between the sphere and the box
      contactPoint     = pos1;
      contactNormal    = n;
      penetrationDepth = dist;
      return true;
   }

   const Vec3 q( R * p );  // Transformation from the projection to the global world frame
   const Vec3 n( d - q );  // Normal direction of the contact (pointing from the box to the sphere)

   const real_t dist( n.length() - radius1 );  // Distance between the sphere and the box

   // Creating a single contact between the sphere and the box
   if( dist < contactThreshold )
   {
      contactPoint     = pos2+q;
      contactNormal    = n.getNormalizedOrZero();
      if (contactNormal.sqrLength() < real_t(0.1))
      {
         contactNormal = Vec3(1,0,0);
      }
      penetrationDepth = dist;
      return true;
   }
   return false;
}

inline
bool detectSphereCylindricalBoundaryCollision( const Vec3&   pos1,
                                               const real_t& radius1,
                                               const Vec3&   pos2,
                                               const real_t& radius2,
                                               const Vec3&   axis2,
                                               Vec3&         contactPoint,
                                               Vec3&         contactNormal,
                                               real_t&       penetrationDepth,
                                               const real_t& contactThreshold)
{
   WALBERLA_ASSERT_FLOAT_EQUAL( axis2.sqrLength(), real_t(1));
   WALBERLA_ASSERT_GREATER( radius2, radius1 );

   const Vec3   dist      = (pos1 - pos2) - ((pos1 - pos2) * axis2) * axis2;
   const real_t effRadius = radius2 - radius1;
   if( effRadius * effRadius - dist.sqrLength() < contactThreshold ) {
      contactNormal    = -dist.getNormalized();
      penetrationDepth = effRadius - dist.length();
      contactPoint     = ( pos1 - ( radius1 + penetrationDepth ) * contactNormal );

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


} //namespace analytic
} //namespace collision_detection
} //namespace mesa_pd
} //namespace walberla
