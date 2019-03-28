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
//! \file Collide.h
//! \author Klaus Iglberger
//! \author Tobias Scharpff
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/Types.h"
#include "pe/contact/Contact.h"
#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/CylindricalBoundary.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Union.h"
#include "pe/utility/BodyCast.h"

#include "core/debug/Debug.h"
#include "core/math/RotationMatrix.h"
#include "core/math/Shims.h"
#include "geometry/GeometricalFunctions.h"

namespace walberla {
namespace pe {
namespace fcd {
namespace analytic {

template <typename Container>
inline
bool collide( GeomID bd1, GeomID bd2, Container& container );

template <typename Container>
inline
bool collide( SphereID s1, SphereID s2, Container& container );

template <typename Container>
inline
bool collide( SphereID s, PlaneID p, Container& container );
template <typename Container>
inline
bool collide( PlaneID p, SphereID s, Container& container );

template <typename Container>
inline
bool collide( SphereID s, CylindricalBoundaryID cb, Container& container );
template <typename Container>
inline
bool collide( CylindricalBoundaryID cb, SphereID s, Container& container );

template <typename Container>
inline
bool collide( SphereID s, BoxID b, Container& container );
template <typename Container>
inline
bool collide( BoxID b, SphereID s, Container& container );

template <typename Container>
inline
bool collide( BoxID b1, BoxID b2, Container& container );

template <typename Container>
inline
bool collide( BoxID b, PlaneID p, Container& container );
template <typename Container>
inline
bool collide( PlaneID p, BoxID b, Container& container );

template <typename Container>
inline
bool collide( CapsuleID c1, CapsuleID c2, Container& container );

template <typename Container>
inline
bool collide( CapsuleID c, PlaneID p, Container& container );
template <typename Container>
inline
bool collide( PlaneID p, CapsuleID c, Container& container );

template <typename Container>
inline
bool collide( SphereID s, CapsuleID c, Container& container );
template <typename Container>
inline
bool collide( CapsuleID c, SphereID s, Container& container );


template <typename Container>
inline
bool collide( BoxID b, CapsuleID c, Container& container );
template <typename Container>
inline
bool collide( CapsuleID c, BoxID b, Container& container );

template <typename... BodyTypes, typename BodyB, typename Container>
inline
bool collide( Union<BodyTypes...>* bd1, BodyB* bd2, Container& container );

template <typename BodyA, typename... BodyTypes, typename Container>
inline
bool collide( BodyA* bd1, Union<BodyTypes...>* bd2, Container& container );

template <typename... BodyTypesA, typename... BodyTypesB, typename Container>
inline
bool collide( Union<BodyTypesA...>* bd1, Union<BodyTypesB...>* bd2, Container& container );

} //namespace analytic

template <typename Container>
struct AnalyticCollideFunctor
{
   Container& contacts_;

   AnalyticCollideFunctor(Container& contacts) : contacts_(contacts) {}

   template< typename BodyType1, typename BodyType2 >
   bool operator()( BodyType1* bd1, BodyType2* bd2) { using namespace analytic; return collide( bd1, bd2, contacts_); }
};

template <typename BodyType1, typename Container>
struct AnalyticSingleCollideFunctor
{
   BodyType1* bd1_;
   Container& contacts_;

   AnalyticSingleCollideFunctor(BodyType1* bd1, Container& contacts) : bd1_(bd1), contacts_(contacts) {}

   template< typename BodyType2 >
   bool operator()( BodyType2* bd2) { using namespace analytic; return collide( bd1_, bd2, contacts_); }
};

namespace analytic {

//*************************************************************************************************
/*!\brief Contact generation between two colliding rigid bodies.
 *
 * \param s1 The first colliding rigid body.
 * \param s2 The second colliding rigid body.
 * \param contacts Contact container for the generated contacts.
 * \return true if contact is detected, false otherwise
 *
 * \attention If no contact was detected the value of penetrationDepth, contactPoint, contactNormal is undefined!
 */
template <typename Container>
inline
bool collide( GeomID /*bd1*/, GeomID /*bd2*/, Container& /*container*/ )
{
   WALBERLA_ABORT("UNSUPPORTED COLLISION!");
   return false;
}

template <typename Container>
inline
bool collide( SphereID s1, SphereID s2, Container& container )
{
   WALBERLA_ASSERT_UNEQUAL(s1->getSystemID(), s2->getSystemID(), "colliding with itself!\n" << s1 << "\n" << s2);

   Vec3 contactNormal = ( s1->getPosition() - s2->getPosition() );
   real_t penetrationDepth = ( contactNormal.length() - s1->getRadius() - s2->getRadius() );

   if( penetrationDepth < contactThreshold ) {
      normalize(contactNormal);
      const real_t k( s2->getRadius() + real_c(0.5) * penetrationDepth );
      Vec3 contactPoint = ( s2->getPosition() + contactNormal * k );

//      if ((-penetrationDepth > s1->getRadius()) || (-penetrationDepth > s2->getRadius()))
//      {
//         WALBERLA_ABORT("DEEP SPHERE-SPHERE COLLISION!!!\n" <<
//                              s1 << "\n" << s2 <<
//                              "\ncontactPoint: " << contactPoint <<
//                              "\ncontactNormal: " << contactNormal <<
//                              "\npenetrationDepth: " << penetrationDepth );
//      }

      WALBERLA_LOG_DETAIL( "      Contact created between sphere " << s1->getSystemID()
             << " and sphere " << s2->getSystemID() << " (dist=" << penetrationDepth << ")" );
      container.push_back( Contact(s1, s2, contactPoint, contactNormal, penetrationDepth) );
      return true;
   }
   return false;
}

template <typename Container>
inline
bool collide( SphereID s, PlaneID p, Container& container )
{
   const real_t k( p->getNormal() * s->getPosition() );
   real_t penetrationDepth = ( k - s->getRadius() - p->getDisplacement() );

   if( penetrationDepth < contactThreshold ) {
      Vec3 contactPoint = ( s->getPosition() - ( s->getRadius() + penetrationDepth ) * p->getNormal() );

      WALBERLA_LOG_DETAIL( "      Contact created between sphere " << s->getID()
             << " and plane " << p->getID() << " (dist=" << penetrationDepth << ")");

      Vec3 contactNormal = p->getNormal();

      container.push_back( Contact(s, p, contactPoint, contactNormal, penetrationDepth) );
      return true;
   }
   return false;
}

template <typename Container>
inline
bool collide( PlaneID p, SphereID s, Container& container )
{
   return collide(s, p, container);
}

template <typename Container>
inline
bool collide( SphereID s, CylindricalBoundaryID cb, Container& container )
{
   WALBERLA_ASSERT_GREATER( cb->getRadius(), s->getRadius() );

   const Vec3   dist      = (s->getPosition() - cb->getPosition()) - ((s->getPosition() - cb->getPosition()) * cb->getAxis()) * cb->getAxis();
   const real_t effRadius = cb->getRadius() - s->getRadius();
   if( effRadius * effRadius - dist.sqrLength() < contactThreshold ) {
      const Vec3   contactNormal    = -dist.getNormalized();
      const real_t penetrationDepth = effRadius - dist.length();
      const Vec3   contactPoint     = ( s->getPosition() - ( s->getRadius() + penetrationDepth ) * contactNormal );

      WALBERLA_LOG_DETAIL( "      Contact created between sphere " << s->getID()
             << " and cylindrical boundary " << cb->getID() << " (dist=" << penetrationDepth << ")");

      container.push_back( Contact(s, cb, contactPoint, contactNormal, penetrationDepth) );
      return true;
   }
   return false;
}

template <typename Container>
inline
bool collide( CylindricalBoundaryID cb, SphereID s, Container& container )
{
   return collide(s, cb, container);
}

//*************************************************************************************************
/*!\brief Contact generation between a Sphere and a Box.
 *
 * \param s The colliding sphere.
 * \param b The colliding box.
 * \param container Contact container for the generated contacts.
 * \return void
 *
 * In the case of a contact between a sphere and a box, only a single contact point is generated.
 * This contact point is calculated by a projection of the relative distance of the sphere's
 * center of mass from the box's center of mass on the body frame coordinate system of the box.
 * For all three axes, the projection is limited by the lengths of the sides of the box. This
 * projection is transfered to the global world frame, which gives the distance from the box's
 * center of mass to the contact point.
 */
template <typename Container>
inline
bool collide( SphereID s, BoxID b, Container& container )
{
   //WALBERLA_ABORT(s << b);
   const Vec3 l( real_t(0.5)*b->getLengths() );  // Box half side lengths
   const Vec3& spos( s->getPosition() );       // Global position of the sphere
   const Vec3& bpos( b->getPosition() );       // Global position of the box
   const Mat3& R( b->getRotation() );          // Rotation of the box

   bool outside( false );

   // Calculating the distance between the sphere and the box
   const Vec3 d( spos - bpos );

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
      dist -= s->getRadius();

      // Creating a single contact between the sphere and the box
      WALBERLA_LOG_DETAIL(
                "      Contact created between sphere " << s->getID()
             << " and box " << b->getID() << " (dist=" << dist << ")" );

      container.push_back( Contact( s, b, spos, n, dist ) );
      return true;
   }

   const Vec3 q( R * p );  // Transformation from the projection to the global world frame
   const Vec3 n( d - q );  // Normal direction of the contact (pointing from the box to the sphere)

   const real_t dist( n.length() - s->getRadius() );  // Distance between the sphere and the box

   // Creating a single contact between the sphere and the box
   if( dist < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "      Contact created between sphere " << s->getID()
             << " and box " << b->getID() << " (dist=" << dist << ")" );

      container.push_back( Contact( s, b, bpos+q, n.getNormalized(), dist) );
//      WALBERLA_LOG_DEVEL( Contact( s, b, bpos+q, n.getNormalized(), dist) );
      return true;
   }
   return false;
}
//*************************************************************************************************

template <typename Container>
inline
bool collide( BoxID b, SphereID s, Container& container )
{
   return collide(s, b, container);
}

//*************************************************************************************************
/*!\brief Contact generation between two colliding Box primitives.
 *
 * \param b1 The first colliding box.
 * \param b2 The second colliding box.
 * \param contacts Contact container for the generated contacts.
 * \return void
 */
template <typename Container>
inline
bool collide( BoxID b1, BoxID b2, Container& container )
{
   using namespace math;

   bool isContactDetected = false;

   // Force a defined order of collision detection across processes
   if( b2->getSystemID() < b1->getSystemID() )
      std::swap( b1, b2 );

   // TODO: Performance optimization

   WALBERLA_LOG_DETAIL( "      Fine collision detection between box " << b1->getID() << " and box " << b2->getID() );

   real_t maxDepth( -Limits<real_t>::inf() );
   bool invertNormal( false );
   unsigned int contactCase( 0 );
   Vec3 contactNormal;

   const Mat3& R1( b1->getRotation() );
   const Mat3& R2( b2->getRotation() );

   // Calculating the rotation of box 2 relative to the orientation of box 1
   const Mat3 b2_rR( b1->getRotation().getTranspose() * b2->getRotation() );

   // Calculating the absolute values of the relative rotation
   const Mat3 b2_rQ( fabs( b2_rR ) );

   // Calculating the position of box 2 relative to the orientation of box 1
   const Vec3 b2_rPos( b1->getRotation().getTranspose() * ( b2->getPosition() - b1->getPosition() ) );

   // Calculating the half Lengths of both boxes
   const real_t hl1[] = { real_t(0.5) * b1->getLengths()[0],
                        real_t(0.5) * b1->getLengths()[1],
                        real_t(0.5) * b1->getLengths()[2] };
   const real_t hl2[] = { real_t(0.5) * b2->getLengths()[0],
                        real_t(0.5) * b2->getLengths()[1],
                        real_t(0.5) * b2->getLengths()[2] };


   //----- Testing the three axes of box 1 -----

   real_t term1;

   // l = au
   term1 = std::fabs(b2_rPos[0]) - ( hl1[0] + hl2[0]*b2_rQ[0] + hl2[1]*b2_rQ[1] + hl2[2]*b2_rQ[2] );

   if( term1 > contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 1 (l = au) failed!" );
      return isContactDetected;
   }
   else /* if( term1 > maxDepth ) */ {
      maxDepth      = term1;
      contactNormal = Vec3( R1[0], R1[3], R1[6] );
      contactCase   = 1;
      invertNormal  = ( b2_rPos[0] < real_t(0) );

      WALBERLA_LOG_DETAIL(
                "         Contact test 1 succeeded!\n" <<
                "            maxDepth = " << maxDepth );
   }

   // l = av
   term1 = std::fabs(b2_rPos[1]) - ( hl1[1] + hl2[0]*b2_rQ[3] + hl2[1]*b2_rQ[4] + hl2[2]*b2_rQ[5] );

   if( term1 > contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 2 (l = av) failed!" );
      return isContactDetected;
   }
   else if( term1 > maxDepth ) {
      maxDepth      = term1;
      contactNormal = Vec3( R1[1], R1[4], R1[7] );
      contactCase   = 2;
      invertNormal  = ( b2_rPos[1] < real_t(0) );

      WALBERLA_LOG_DETAIL(
                "         Contact test 2 succeeded!\n" <<
                "            maxDepth = " << maxDepth );
   }

   // l = aw
   term1 = std::fabs(b2_rPos[2]) - ( hl1[2] + hl2[0]*b2_rQ[6] + hl2[1]*b2_rQ[7] + hl2[2]*b2_rQ[8] );

   if( term1 > contactThreshold ) {
      WALBERLA_LOG_DETAIL("         Test 3 (l = aw) failed!");
      return isContactDetected;
   }
   else if( term1 > maxDepth ) {
      maxDepth      = term1;
      contactNormal = Vec3( R1[2], R1[5], R1[8] );
      contactCase   = 3;
      invertNormal  = ( b2_rPos[2] < real_t(0) );

      WALBERLA_LOG_DETAIL(
                "         Contact test 3 succeeded!\n" <<
                "            maxDepth = " << maxDepth );
   }


   //----- Testing the three axes of box 2 -----

   real_t term2;

   // l = bu
   term1 = b2_rPos[0]*b2_rR[0] + b2_rPos[1]*b2_rR[3] + b2_rPos[2]*b2_rR[6];
   term2 = std::fabs(term1) - ( hl1[0]*b2_rQ[0] + hl1[1]*b2_rQ[3] + hl1[2]*b2_rQ[6] + hl2[0] );

   if( term2 > contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 4 (l = bu) failed!" );
      return isContactDetected;
   }
   else if( term2 > maxDepth ) {
      maxDepth      = term2;
      contactNormal = Vec3( R2[0], R2[3], R2[6] );
      contactCase   = 4;
      invertNormal  = ( term1 < real_t(0) );

      WALBERLA_LOG_DETAIL(
                "         Contact test 4 succeeded!\n" <<
                "            maxDepth = " << maxDepth );
   }

   // l = bv
   term1 = b2_rPos[0]*b2_rR[1] + b2_rPos[1]*b2_rR[4] + b2_rPos[2]*b2_rR[7];
   term2 = std::fabs(term1) - ( hl1[0]*b2_rQ[1] + hl1[1]*b2_rQ[4] + hl1[2]*b2_rQ[7] + hl2[1] );

   if( term2 > contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 5 (l = bv) failed!" );
      return isContactDetected;
   }
   else if( term2 > maxDepth ) {
      maxDepth      = term2;
      contactNormal = Vec3( R2[1], R2[4], R2[7] );
      contactCase   = 5;
      invertNormal  = ( term1 < real_t(0) );

      WALBERLA_LOG_DETAIL(
                "         Contact test 5 succeeded!\n" <<
                "            maxDepth = " << maxDepth );
   }

   // l = bw
   term1 = b2_rPos[0]*b2_rR[2] + b2_rPos[1]*b2_rR[5] + b2_rPos[2]*b2_rR[8];
   term2 = std::fabs(term1) - ( hl1[0]*b2_rQ[2] + hl1[1]*b2_rQ[5] + hl1[2]*b2_rQ[8] + hl2[2] );

   if( term2 > contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 6 (l = bw) failed!" );
      return isContactDetected;
   }
   else if( term2 > maxDepth ) {
      maxDepth      = term2;
      contactNormal = Vec3( R2[2], R2[5], R2[8] );
      contactCase   = 6;
      invertNormal  = ( term1 < real_t(0) );

      WALBERLA_LOG_DETAIL(
                "         Contact test 6 succeeded!\n" <<
                "            maxDepth = " << maxDepth );
   }


   //----- Testing the all nine combinations of the axes of the two boxes -----

   real_t term3;
   real_t sum;
   real_t length;
   Vec3 normal_c;

   // l = au x bu
   term1 = b2_rPos[2] * b2_rR[3] - b2_rPos[1] * b2_rR[6];
   term2 = hl1[1] * b2_rQ[6] + hl1[2] * b2_rQ[3];
   term3 = hl2[2] * b2_rQ[1] + hl2[1] * b2_rQ[2];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 7 (l = au x bu) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[6]) + sq(b2_rR[3]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( 0, -b2_rR[6]/length, b2_rR[3]/length );
            contactCase  = 7;
            invertNormal = ( term1 < real_t(0) );

            WALBERLA_LOG_DETAIL(
                      "         Contact test 7 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }

   // l = au x bv
   term1 = b2_rPos[2] * b2_rR[4] - b2_rPos[1] * b2_rR[7];
   term2 = hl1[1] * b2_rQ[7] + hl1[2] * b2_rQ[4];
   term3 = hl2[0] * b2_rQ[2] + hl2[2] * b2_rQ[0];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 8 (l = au x bv) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[7]) + sq(b2_rR[4]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( 0, -b2_rR[7]/length, b2_rR[4]/length );
            contactCase  = 8;
            invertNormal = ( term1 < real_t(0) );

            WALBERLA_LOG_DETAIL(
                      "         Contact test 8 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }

   // l = au x bw
   term1 = b2_rPos[2] * b2_rR[5] - b2_rPos[1] * b2_rR[8];
   term2 = hl1[1] * b2_rQ[8] + hl1[2] * b2_rQ[5];
   term3 = hl2[1] * b2_rQ[0] + hl2[0] * b2_rQ[1];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL("         Test 9 (l = au x bw) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[8]) + sq(b2_rR[5]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( 0, -b2_rR[8]/length, b2_rR[5]/length );
            contactCase  = 9;
            invertNormal = ( term1 < real_t(0) );

            WALBERLA_LOG_DETAIL(
                      "         Contact test 9 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }

   // l = av x bu
   term1 = b2_rPos[0] * b2_rR[6] - b2_rPos[2] * b2_rR[0];
   term2 = hl1[2] * b2_rQ[0] + hl1[0] * b2_rQ[6];
   term3 = hl2[2] * b2_rQ[4] + hl2[1] * b2_rQ[5];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 10 (l = av x bu) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[6]) + sq(b2_rR[0]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( b2_rR[6]/length, 0, -b2_rR[0]/length );
            contactCase  = 10;
            invertNormal = ( term1 < real_t(0) ) ;

            WALBERLA_LOG_DETAIL(
                      "         Contact test 10 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }

   // l = av x bv
   term1 = b2_rPos[0] * b2_rR[7] - b2_rPos[2] * b2_rR[1];
   term2 = hl1[2] * b2_rQ[1] + hl1[0] * b2_rQ[7];
   term3 = hl2[0] * b2_rQ[5] + hl2[2] * b2_rQ[3];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 11 (l = av x bv) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[7]) + sq(b2_rR[1]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( b2_rR[7]/length, 0, -b2_rR[1]/length );
            contactCase  = 11;
            invertNormal = ( term1 < real_t(0) );

            WALBERLA_LOG_DETAIL(
                      "         Contact test 11 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }

   // l = av x bw
   term1 = b2_rPos[0] * b2_rR[8] - b2_rPos[2] * b2_rR[2];
   term2 = hl1[2] * b2_rQ[2] + hl1[0] * b2_rQ[8];
   term3 = hl2[1] * b2_rQ[3] + hl2[0] * b2_rQ[4];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 12 (l = av x bw) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[8]) + sq(b2_rR[2]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( b2_rR[8]/length, 0, -b2_rR[2]/length );
            contactCase  = 12;
            invertNormal = ( term1 < real_t(0) );

            WALBERLA_LOG_DETAIL(
                      "         Contact test 12 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }

   // l = aw x bu
   term1 = b2_rPos[1] * b2_rR[0] - b2_rPos[0] * b2_rR[3];
   term2 = hl1[0] * b2_rQ[3] + hl1[1] * b2_rQ[0];
   term3 = hl2[2] * b2_rQ[7] + hl2[1] * b2_rQ[8];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 13 (l = aw x bu) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[3]) + sq(b2_rR[0]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( -b2_rR[3]/length, b2_rR[0]/length, 0 );
            contactCase  = 13;
            invertNormal = ( term1 < real_t(0) );

            WALBERLA_LOG_DETAIL(
                      "         Contact test 13 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }

   // l = aw x bv
   term1 = b2_rPos[1] * b2_rR[1] - b2_rPos[0] * b2_rR[4];
   term2 = hl1[0] * b2_rQ[4] + hl1[1] * b2_rQ[1];
   term3 = hl2[0] * b2_rQ[8] + hl2[2] * b2_rQ[6];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 14 (l = aw x bv) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[4]) + sq(b2_rR[1]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( -b2_rR[4]/length, b2_rR[1]/length, 0 );
            contactCase  = 14;
            invertNormal = ( term1 < real_t(0) );

            WALBERLA_LOG_DETAIL(
                      "         Contact test 14 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }

   // l = aw x bw
   term1 = b2_rPos[1] * b2_rR[2] - b2_rPos[0] * b2_rR[5];
   term2 = hl1[0] * b2_rQ[5] + hl1[1] * b2_rQ[2];
   term3 = hl2[1] * b2_rQ[6] + hl2[0] * b2_rQ[7];
   sum   = std::fabs(term1) - ( term2 + term3 );

   if( std::fabs(term1) > term2 + term3 + contactThreshold ) {
      WALBERLA_LOG_DETAIL( "         Test 15 (l = aw x bw) failed!" );
      return isContactDetected;
   }
   else {
      length = std::sqrt( sq(b2_rR[5]) + sq(b2_rR[2]) );
      if( length > Limits<real_t>::epsilon() ) {
         sum /= length;
         if( sum > maxDepth && std::fabs( sum - maxDepth ) > Limits<real_t>::accuracy() ) {
            maxDepth     = sum;
            normal_c     = Vec3( -b2_rR[5]/length, b2_rR[2]/length, 0 );
            contactCase  = 15;
            invertNormal = ( term1 < real_t(0) );

            WALBERLA_LOG_DETAIL(
                      "         Contact test 15 succeeded!\n" <<
                      "            maxDepth = " << maxDepth );
         }
      }
   }


   if( contactCase == 0 ) {
      return isContactDetected;
   }

   if( contactCase > 6 ) {
      contactNormal =  R1 * normal_c;
   }

   if( invertNormal ) {
      contactNormal = -contactNormal;
   }


   // TEST
   WALBERLA_LOG_DETAIL(
             "         Selected contact case = " << contactCase << "\n"
          << "         Contact normal = " << contactNormal << "\n"
          << "         Normal invertion? " << invertNormal );


   //----- Treating edge/edge collisions -----

   if( contactCase > 6 )
   {
      WALBERLA_LOG_DETAIL(
                "         Treating edge/edge collision between box " << b1->getID()
             << " and " << b2->getID() << "..." );

      Vec3 pB1( b1->getPosition() );
      Vec3 sign = R1.getTranspose() * contactNormal;
      //Vec3 sign = normal_c;
      for( unsigned int i=0; i<3; ++i ) {
         sign[i] = ( sign[i]>0 ) ? ( hl1[i] ) : ( -hl1[i] );
      }
      #ifdef WALBERLA_LOGLEVEL_DETAIL
      const Vec3 tmp1( sign );
      #endif
      pB1 += R1 * sign;

      Vec3 pB2 = b2->getPosition();
      sign = R2.getTranspose() * contactNormal;
      for( size_t i=0; i<3; ++i ) {
         sign[i] = ( sign[i]>0 ) ? ( -hl2[i] ) : ( hl2[i] );
      }
      #ifdef WALBERLA_LOGLEVEL_DETAIL
      const Vec3 tmp2( sign );
      #endif
      pB2 += R2 * sign;

      Vec3 ua, ub;
      for( size_t i=0; i<3; i++ ) {
         ua[i] = R1[(contactCase-7)/3 + i*3];
         ub[i] = R2[(contactCase-7)%3 + i*3];
      }

      real_t s, t;
      geometry::intersectLines( pB1, ua, pB2, ub, s, t );
      pB1 += s * ua;
      pB2 += t * ub;
      Vec3 gpos = real_t(0.5) * ( pB1 + pB2 );
      Vec3 normal = ( ua % ub ).getNormalized();

      WALBERLA_LOG_DETAIL(
                "            box A (b1) = " << b1->getID() << "\n"
             << "            box B (b2) = " << b2->getID() << "\n"
             << "            normal_c   = " << normal_c << "\n"
             << "            tmp1 (for pB1) = " << tmp1 << "\n"
             << "            tmp2 (for pB2) = " << tmp2 << "\n"
             << "            pB1  = " << pB1 << "\n"
             << "            pB2  = " << pB2 << "\n"
             << "            gpos = " << gpos << "\n"
             << "            contactNormal from A to B = " << contactNormal << "\n"
             << "            contactNormal from B to A = " << -contactNormal << "\n"
             << "            ua = " << ua << "\n"
             << "            ub = " << ub << "\n"
             << "            ua x ub = " << ua % ub << "\n"
             << "            normal = " << normal << "\n\n" );

      // TEST
      if( normal*contactNormal < real_t(0) ) {
         WALBERLA_LOG_DETAIL(
                   "         Inverting ub!\n"
                << "         ua = " << ua << "\n"
                << "         ub = " << ub );
         ub = -ub;
      }

      WALBERLA_LOG_DETAIL(
                "      Edge/edge contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << maxDepth << ")" );

      container.push_back( Contact( b2, b1, gpos, contactNormal, maxDepth ) ); // CONTACT GENERATION CHANGED!!!
      return isContactDetected;
   }


   //----- Treating vertex/face collisions -----

   WALBERLA_LOG_DETAIL( "         Treating vertex/face collision..." );

   const real_t* hla( hl1 );
   const real_t* hlb( hl2 );

   if( contactCase > 3 ) {
      std::swap( b1, b2 );
      std::swap( hla, hlb );
      contactNormal = -contactNormal;
   }

   WALBERLA_LOG_DETAIL(
             "            Box A = " << b1->getID() << "\n"
          << "            Box B = " << b2->getID() );

   const Mat3& Ra( b1->getRotation() );
   const Mat3& Rb( b2->getRotation() );


   // Calculating the relative contact normal in the body frame of bx A
   const Vec3 normala( Ra.getTranspose() * contactNormal );

   WALBERLA_LOG_DETAIL( "            Relative contact normal in the body frame of A = " << normala );


   // Calculating the relative contact normal in the body frame of box B
   const Vec3 normalb( Rb.getTranspose() * contactNormal );

   WALBERLA_LOG_DETAIL( "            Relative contact normal in the body frame of B = " << normalb );


   // Estimating the collides face of box B
   unsigned int xb(0), yb(0), zb(0);

   if( std::fabs( normalb[0] ) > std::fabs( normalb[1] ) ) {
      if( std::fabs( normalb[0] ) > std::fabs( normalb[2] ) ) {
         xb = 0; yb = 1; zb = 2;
      }
      else {
         xb = 2; yb = 0; zb = 1;
      }
   }
   else if( std::fabs( normalb[1] ) > std::fabs( normalb[2] ) ) {
      xb = 1; yb = 2; zb = 0;
   }
   else {
      xb = 2; yb = 0; zb = 1;
   }


   // Estimating the colliding face of box A
   unsigned int xa(0), ya(0), za(0);

   if( contactCase < 4 ) {
      xa = contactCase - 1;
   }
   else {
      xa = contactCase - 4;
   }
   if( xa == 0 ) {
      ya = 1; za = 2;
   }
   else if( xa == 1 ) {
      ya = 2; za = 0;
   }
   else {
      ya = 0; za = 1;
   }

   WALBERLA_LOG_DETAIL(
             "            Colliding face of box A:  xa = " << xa << " , ya = " << ya << " , za = " << za
          << "            Colliding face of box B:  xb = " << xb << " , yb = " << yb << " , zb = " << zb );


   // Calculating the four vertices of the colliding face of box B in the frame of box B
   Vec3 vertexb[4];

   vertexb[0][xb] = ( normalb[xb] > real_t(0) )?( -hlb[xb] ):( hlb[xb] );
   vertexb[0][yb] = -hlb[yb];
   vertexb[0][zb] = -hlb[zb];

   vertexb[1][xb] = ( normalb[xb] > real_t(0) )?( -hlb[xb] ):( hlb[xb] );
   vertexb[1][yb] = hlb[yb];
   vertexb[1][zb] = -hlb[zb];

   vertexb[2][xb] = ( normalb[xb] > real_t(0) )?( -hlb[xb] ):( hlb[xb] );
   vertexb[2][yb] = -hlb[yb];
   vertexb[2][zb] = hlb[zb];

   vertexb[3][xb] = ( normalb[xb] > real_t(0) )?( -hlb[xb] ):( hlb[xb] );
   vertexb[3][yb] = hlb[yb];
   vertexb[3][zb] = hlb[zb];


   // Translating the four vertices to the body frame of box A
   const Vec3 ab( b2->getPosition() - b1->getPosition() );

   Vec3 vertexba[] = { Ra.getTranspose() * ( ab + Rb*vertexb[0] ),
                       Ra.getTranspose() * ( ab + Rb*vertexb[1] ),
                       Ra.getTranspose() * ( ab + Rb*vertexb[2] ),
                       Ra.getTranspose() * ( ab + Rb*vertexb[3] ) };


   //----- Calculating the line/line intersections between the two colliding faces -----

   const real_t offseta( ( normala[xa] > real_t(0) )?( hla[xa] ):( -hla[xa] ) );
   real_t s, dist, tmp;

   // Intersection between v0--v1 with hla[ya]
   if( ( vertexba[0][ya] > hla[ya] ) ^ ( vertexba[1][ya] > hla[ya] ) )
   {
      s = ( hla[ya] - vertexba[0][ya] ) / ( vertexba[1][ya] - vertexba[0][ya] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 1\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[0][za] + s*( vertexba[1][za] - vertexba[0][za] ) ) < hla[za] )
      {
         dist = std::fabs( vertexba[0][xa] + s*( vertexba[1][xa] - vertexba[0][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = hla[ya];
            posa[za] = tmp;

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v0--v1 with -hla[ya]
   if( ( vertexba[0][ya] > -hla[ya] ) ^ ( vertexba[1][ya] > -hla[ya] ) )
   {
      s = ( -hla[ya] - vertexba[0][ya] ) / ( vertexba[1][ya] - vertexba[0][ya] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 2\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[0][za] + s*( vertexba[1][za] - vertexba[0][za] ) ) < hla[za] )
      {
         dist = std::fabs( vertexba[0][xa] + s*( vertexba[1][xa] - vertexba[0][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = -hla[ya];
            posa[za] = tmp;

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v0--v1 with hla[za]
   if( ( vertexba[0][za] > hla[za] ) ^ ( vertexba[1][za] > hla[za] ) )
   {
      s = ( hla[za] - vertexba[0][za] ) / ( vertexba[1][za] - vertexba[0][za] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 3\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[0][ya] + s*( vertexba[1][ya] - vertexba[0][ya] ) ) < hla[ya] )
      {
         dist = std::fabs( vertexba[0][xa] + s*( vertexba[1][xa] - vertexba[0][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = tmp;
            posa[za] = hla[za];

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v0--v1 with -hla[za]
   if( ( vertexba[0][za] > -hla[za] ) ^ ( vertexba[1][za] > -hla[za] ) )
   {
      s = ( -hla[za] - vertexba[0][za] ) / ( vertexba[1][za] - vertexba[0][za] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 4\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[0][ya] + s*( vertexba[1][ya] - vertexba[0][ya] ) ) < hla[ya] )
      {
         dist = std::fabs( vertexba[0][xa] + s*( vertexba[1][xa] - vertexba[0][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = tmp;
            posa[za] = -hla[za];

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }


   // Intersection between v0--v2 with hla[ya]
   if( ( vertexba[0][ya] > hla[ya] ) ^ ( vertexba[2][ya] > hla[ya] ) )
   {
      s = ( hla[ya] - vertexba[0][ya] ) / ( vertexba[2][ya] - vertexba[0][ya] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 5\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[0][za] + s*( vertexba[2][za] - vertexba[0][za] ) ) < hla[za] )
      {
         dist = std::fabs( vertexba[0][xa] + s*( vertexba[2][xa] - vertexba[0][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = hla[ya];
            posa[za] = tmp;

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v0--v2 with -hla[ya]
   if( ( vertexba[0][ya] > -hla[ya] ) ^ ( vertexba[2][ya] > -hla[ya] ) )
   {
      s = ( -hla[ya] - vertexba[0][ya] ) / ( vertexba[2][ya] - vertexba[0][ya] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 6\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[0][za] + s*( vertexba[2][za] - vertexba[0][za] ) ) < hla[za] )
      {
         dist = std::fabs( vertexba[0][xa] + s*( vertexba[2][xa] - vertexba[0][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = -hla[ya];
            posa[za] = tmp;

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v0--v2 with hla[za]
   if( ( vertexba[0][za] > hla[za] ) ^ ( vertexba[2][za] > hla[za] ) )
   {
      s = ( hla[za] - vertexba[0][za] ) / ( vertexba[2][za] - vertexba[0][za] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 7\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[0][ya] + s*( vertexba[2][ya] - vertexba[0][ya] ) ) < hla[ya] )
      {
         dist = std::fabs( vertexba[0][xa] + s*( vertexba[2][xa] - vertexba[0][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = tmp;
            posa[za] = hla[za];

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v0--v2 with -hla[za]
   if( ( vertexba[0][za] > -hla[za] ) ^ ( vertexba[2][za] > -hla[za] ) )
   {
      s = ( -hla[za] - vertexba[0][za] ) / ( vertexba[2][za] - vertexba[0][za] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 8\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[0][ya] + s*( vertexba[2][ya] - vertexba[0][ya] ) ) < hla[ya] )
      {
         dist = std::fabs( vertexba[0][xa] + s*( vertexba[2][xa] - vertexba[0][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = tmp;
            posa[za] = -hla[za];

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }


   // Intersection between v3--v1 with hla[ya]
   if( ( vertexba[3][ya] > hla[ya] ) ^ ( vertexba[1][ya] > hla[ya] ) )
   {
      s = ( hla[ya] - vertexba[3][ya] ) / ( vertexba[1][ya] - vertexba[3][ya] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 9\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[3][za] + s*( vertexba[1][za] - vertexba[3][za] ) ) < hla[za] )
      {
         dist = std::fabs( vertexba[3][xa] + s*( vertexba[1][xa] - vertexba[3][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = hla[ya];
            posa[za] = tmp;

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v3--v1 with -hla[ya]
   if( ( vertexba[3][ya] > -hla[ya] ) ^ ( vertexba[1][ya] > -hla[ya] ) )
   {
      s = ( -hla[ya] - vertexba[3][ya] ) / ( vertexba[1][ya] - vertexba[3][ya] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 10\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[3][za] + s*( vertexba[1][za] - vertexba[3][za] ) ) < hla[za] )
      {
         dist = std::fabs( vertexba[3][xa] + s*( vertexba[1][xa] - vertexba[3][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = -hla[ya];
            posa[za] = tmp;

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v3--v1 with hla[za]
   if( ( vertexba[3][za] > hla[za] ) ^ ( vertexba[1][za] > hla[za] ) )
   {
      s = ( hla[za] - vertexba[3][za] ) / ( vertexba[1][za] - vertexba[3][za] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 11\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[3][ya] + s*( vertexba[1][ya] - vertexba[3][ya] ) ) < hla[ya] )
      {
         dist = std::fabs( vertexba[3][xa] + s*( vertexba[1][xa] - vertexba[3][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = tmp;
            posa[za] = hla[za];

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v3--v1 with -hla[za]
   if( ( vertexba[3][za] > -hla[za] ) ^ ( vertexba[1][za] > -hla[za] ) )
   {
      s = ( -hla[za] - vertexba[3][za] ) / ( vertexba[1][za] - vertexba[3][za] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 12\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[3][ya] + s*( vertexba[1][ya] - vertexba[3][ya] ) ) < hla[ya] )
      {
         dist = std::fabs( vertexba[3][xa] + s*( vertexba[1][xa] - vertexba[3][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = tmp;
            posa[za] = -hla[za];

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }


   // Intersection between v3--v2 with hla[ya]
   if( ( vertexba[3][ya] > hla[ya] ) ^ ( vertexba[2][ya] > hla[ya] ) )
   {
      s = ( hla[ya] - vertexba[3][ya] ) / ( vertexba[2][ya] - vertexba[3][ya] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 13\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[3][za] + s*( vertexba[2][za] - vertexba[3][za] ) ) < hla[za] )
      {
         dist = std::fabs( vertexba[3][xa] + s*( vertexba[2][xa] - vertexba[3][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = hla[ya];
            posa[za] = tmp;

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v3--v2 with -hla[ya]
   if( ( vertexba[3][ya] > -hla[ya] ) ^ ( vertexba[2][ya] > -hla[ya] ) )
   {
      s = ( -hla[ya] - vertexba[3][ya] ) / ( vertexba[2][ya] - vertexba[3][ya] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 14\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[3][za] + s*( vertexba[2][za] - vertexba[3][za] ) ) < hla[za] )
      {
         dist = std::fabs( vertexba[3][xa] + s*( vertexba[2][xa] - vertexba[3][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = -hla[ya];
            posa[za] = tmp;

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v3--v2 with hla[za]
   if( ( vertexba[3][za] > hla[za] ) ^ ( vertexba[2][za] > hla[za] ) )
   {
      s = ( hla[za] - vertexba[3][za] ) / ( vertexba[2][za] - vertexba[3][za] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 15\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[3][ya] + s*( vertexba[2][ya] - vertexba[3][ya] ) ) < hla[ya] )
      {
         dist = std::fabs( vertexba[3][xa] + s*( vertexba[2][xa] - vertexba[3][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = tmp;
            posa[za] = hla[za];

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }

   // Intersection between v3--v2 with -hla[za]
   if( ( vertexba[3][za] > -hla[za] ) ^ ( vertexba[2][za] > -hla[za] ) )
   {
      s = ( -hla[za] - vertexba[3][za] ) / ( vertexba[2][za] - vertexba[3][za] );

      WALBERLA_LOG_DETAIL(
                "            Treating case 16\n" <<
                "               s = " << s );

      if( s > real_t(0) && s < real_t(1) &&
          std::fabs( tmp = vertexba[3][ya] + s*( vertexba[2][ya] - vertexba[3][ya] ) ) < hla[ya] )
      {
         dist = std::fabs( vertexba[3][xa] + s*( vertexba[2][xa] - vertexba[3][xa] ) ) - hla[xa];
         if( dist < contactThreshold )
         {
            WALBERLA_LOG_DETAIL(
                      "      Vertex/face contact created between box " << b1->getID()
                   << " and box " << b2->getID() << " (dist=" << dist << ")" );

            Vec3 posa;
            posa[xa] = offseta;
            posa[ya] = tmp;
            posa[za] = -hla[za];

            const Vec3 gpos( b1->pointFromBFtoWF( posa ) );
            container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
            isContactDetected = true;
         }
      }
   }


   //----- Calculating contact points for the vertices of box B -----

   if( std::fabs(vertexba[0][ya]) <= hla[ya] && std::fabs(vertexba[0][za]) <= hla[za] &&
       ( dist = std::fabs(vertexba[0][xa]) - hla[xa] ) < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "            Treating case 17" );

      WALBERLA_LOG_DETAIL(
                "      Vertex/face contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << dist << ")" );

      const Vec3 gpos( b1->pointFromBFtoWF( vertexba[0] ) );
      container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
      isContactDetected = true;
   }

   if( std::fabs(vertexba[1][ya]) <= hla[ya] && std::fabs(vertexba[1][za]) <= hla[za] &&
       ( dist = std::fabs(vertexba[1][xa]) - hla[xa] ) < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "            Treating case 18" );

      WALBERLA_LOG_DETAIL(
                "      Vertex/face contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << dist << ")" );

      const Vec3 gpos( b1->pointFromBFtoWF( vertexba[1] ) );
      container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
      isContactDetected = true;
   }

   if( std::fabs(vertexba[2][ya]) <= hla[ya] && std::fabs(vertexba[2][za]) <= hla[za] &&
       ( dist = std::fabs(vertexba[2][xa]) - hla[xa] ) < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "            Treating case 19" );

      WALBERLA_LOG_DETAIL(
                "      Vertex/face contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << dist << ")" );

      const Vec3 gpos( b1->pointFromBFtoWF( vertexba[2] ) );
      container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
      isContactDetected = true;
   }

   if( std::fabs(vertexba[3][ya]) <= hla[ya] && std::fabs(vertexba[3][za]) <= hla[za] &&
       ( dist = std::fabs(vertexba[3][xa]) - hla[xa] ) < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "            Treating case 20" );

      WALBERLA_LOG_DETAIL(
                "      Vertex/face contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << dist << ")" );

      const Vec3 gpos( b1->pointFromBFtoWF( vertexba[3] ) );
      container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
      isContactDetected = true;
   }


   //----- Calculating contact points for the vertices of box A -----

   // Calculating the four vertices of the colliding face of box B in the frame of box B
   Vec3 vertexa[4];

   vertexa[0][xa] = ( normala[xa] > real_t(0) )?( hla[xa] ):( -hla[xa] );
   vertexa[0][ya] = -hla[ya];
   vertexa[0][za] = -hla[za];

   vertexa[1][xa] = ( normala[xa] > real_t(0) )?( hla[xa] ):( -hla[xa] );
   vertexa[1][ya] = hla[ya];
   vertexa[1][za] = -hla[za];

   vertexa[2][xa] = ( normala[xa] > real_t(0) )?( hla[xa] ):( -hla[xa] );
   vertexa[2][ya] = -hla[ya];
   vertexa[2][za] = hla[za];

   vertexa[3][xa] = ( normala[xa] > real_t(0) )?( hla[xa] ):( -hla[xa] );
   vertexa[3][ya] = hla[ya];
   vertexa[3][za] = hla[za];


   // Translating the four vertices to the body frame of box B
   const Vec3 ba( b1->getPosition() - b2->getPosition() );

   Vec3 vertexab[] = { Rb.getTranspose() * ( ba + Ra*vertexa[0] ),
                       Rb.getTranspose() * ( ba + Ra*vertexa[1] ),
                       Rb.getTranspose() * ( ba + Ra*vertexa[2] ),
                       Rb.getTranspose() * ( ba + Ra*vertexa[3] ) };


   // In order to avoid vertex/vertex-contacts to be generated twice, the evaluation whether
   // a contact point is generated for a vertex of A has the additional requirement, that
   // the vertex of A does not coincide with a vertex of B.

   if( ( ( std::fabs(vertexab[0][yb]) <= hlb[yb] && std::fabs(vertexab[0][zb]) <  hlb[zb] ) ||
         ( std::fabs(vertexab[0][yb]) <  hlb[yb] && std::fabs(vertexab[0][zb]) <= hlb[zb] ) ) &&
       ( dist = std::fabs(vertexab[0][xb]) - hlb[xb] ) < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "            Treating case 21" );

      WALBERLA_LOG_DETAIL(
                "      Vertex/face contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << dist << ")" );

      const Vec3 gpos( b2->pointFromBFtoWF( vertexab[0] ) );
      container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
      isContactDetected = true;
   }

   if( ( ( std::fabs(vertexab[1][yb]) <= hlb[yb] && std::fabs(vertexab[1][zb]) <  hlb[zb] ) ||
         ( std::fabs(vertexab[1][yb]) <  hlb[yb] && std::fabs(vertexab[1][zb]) <= hlb[zb] ) ) &&
       ( dist = std::fabs(vertexab[1][xb]) - hlb[xb] ) < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "            Treating case 22" );
      WALBERLA_LOG_DETAIL(
                "      Vertex/face contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << dist << ")" );

      const Vec3 gpos( b2->pointFromBFtoWF( vertexab[1] ) );
      container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
      isContactDetected = true;
   }

   if( ( ( std::fabs(vertexab[2][yb]) <= hlb[yb] && std::fabs(vertexab[2][zb]) <  hlb[zb] ) ||
         ( std::fabs(vertexab[2][yb]) <  hlb[yb] && std::fabs(vertexab[2][zb]) <= hlb[zb] ) ) &&
       ( dist = std::fabs(vertexab[2][xb]) - hlb[xb] ) < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "            Treating case 23" );

      WALBERLA_LOG_DETAIL(
                "      Vertex/face contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << dist << ")" );

      const Vec3 gpos( b2->pointFromBFtoWF( vertexab[2] ) );
      container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
      isContactDetected = true;
   }

   if( ( ( std::fabs(vertexab[3][yb]) <= hlb[yb] && std::fabs(vertexab[3][zb]) <  hlb[zb] ) ||
         ( std::fabs(vertexab[3][yb]) <  hlb[yb] && std::fabs(vertexab[3][zb]) <= hlb[zb] ) ) &&
       ( dist = std::fabs(vertexab[3][xb]) - hlb[xb] ) < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "            Treating case 24" );

      WALBERLA_LOG_DETAIL(
                "      Vertex/face contact created between box " << b1->getID()
             << " and box " << b2->getID() << " (dist=" << dist << ")" );

      const Vec3 gpos( b2->pointFromBFtoWF( vertexab[3] ) );
      container.push_back( Contact( b2, b1, gpos, contactNormal, dist ) );
      isContactDetected = true;
   }
   return isContactDetected;
}

//*************************************************************************************************
/*!\brief Contact generation between a Box and a Plane.
 *
 * \param b The colliding box.
 * \param p The colliding plane.
 * \param container Contact container for the generated contacts.
 * \return true if a contact is detected, false otherwise
 *
 * For each corner of the box, the distance to the plane is calculated with a projection on
 * the plane's normal. For each corner, which lies on or inside the plane, a contact point is
 * generated.
 */
template <typename Container>
inline
bool collide( BoxID b, PlaneID p, Container& container )
{
   bool collision = false;

   real_t dist, k;
   Vec3 pos;
   const Vec3& bpos( b->getPosition() );
   const Vec3& lengths( b->getLengths() );
   const Mat3& R( b->getRotation() );
   const real_t xrange( static_cast<real_t>( 0.5 )*lengths[0] );
   const real_t yrange( static_cast<real_t>( 0.5 )*lengths[1] );
   const real_t zrange( static_cast<real_t>( 0.5 )*lengths[2] );
   const Vec3 nrel( b->getRotation().getTranspose() * p->getNormal() );
   const real_t minlength( -static_cast<real_t>( 0.99 ) * std::min( xrange, std::min(yrange, zrange) ) );

   // Test of lower-left-front corner
   if( -xrange*nrel[0] - yrange*nrel[1] - zrange*nrel[2] < minlength )
   {
      pos = bpos + R*Vec3( -xrange, -yrange, -zrange );
      k = pos * p->getNormal();
      dist = k - p->getDisplacement();
      if( dist < contactThreshold )
      {
         WALBERLA_LOG_DETAIL( "      Contact created between box " << b->getID() << " and plane " << p->getID() << " (dist=" << dist << ")" );

         container.push_back( Contact(b, p, pos, p->getNormal(), dist) );
         collision = true;
      }
   }

   // Test of the lower-right-front corner
   if( xrange*nrel[0] - yrange*nrel[1] - zrange*nrel[2] < minlength )
   {
      pos = bpos + R*Vec3( xrange, -yrange, -zrange );
      k = pos * p->getNormal();
      dist = k - p->getDisplacement();
      if( dist < contactThreshold )
      {
         WALBERLA_LOG_DETAIL(
                   "      Contact created between box " << b->getID()
                << " and plane " << p->getID() << " (dist=" << dist << ")" );

         container.push_back( Contact(b, p, pos, p->getNormal(), dist) );
         collision = true;
      }
   }

   // Test of the lower-right-back corner
   if( xrange*nrel[0] + yrange*nrel[1] - zrange*nrel[2] < minlength )
   {
      pos = bpos + R*Vec3( xrange, yrange, -zrange );
      k = pos * p->getNormal();
      dist = k - p->getDisplacement();
      if( dist < contactThreshold )
      {
         WALBERLA_LOG_DETAIL(
                   "      Contact created between box " << b->getID()
                << " and plane " << p->getID() << " (dist=" << dist << ")");

         container.push_back( Contact(b, p, pos, p->getNormal(), dist) );
         collision = true;
      }
   }

   // Test of the lower-left-back corner
   if( -xrange*nrel[0] + yrange*nrel[1] - zrange*nrel[2] < minlength )
   {
      pos = bpos + R*Vec3( -xrange, yrange, -zrange );
      k = pos * p->getNormal();
      dist = k - p->getDisplacement();
      if( dist < contactThreshold )
      {
         WALBERLA_LOG_DETAIL(
                   "      Contact created between box " << b->getID()
                << " and plane " << p->getID() << " (dist=" << dist << ")");

         container.push_back( Contact(b, p, pos, p->getNormal(), dist) );
         collision = true;
      }
   }

   // Test of the upper-left-front corner
   if( -xrange*nrel[0] - yrange*nrel[1] + zrange*nrel[2]  < minlength )
   {
      pos = bpos + R*Vec3( -xrange, -yrange, zrange );
      k = pos * p->getNormal();
      dist = k - p->getDisplacement();
      if( dist < contactThreshold )
      {
         WALBERLA_LOG_DETAIL(
                   "      Contact created between box " << b->getID()
                << " and plane " << p->getID() << " (dist=" << dist << ")");

         container.push_back( Contact(b, p, pos, p->getNormal(), dist) );
         collision = true;
      }
   }

   // Test of the upper-right-front corner
   if( xrange*nrel[0] - yrange*nrel[1] + zrange*nrel[2] < minlength )
   {
      pos = bpos + R*Vec3( xrange, -yrange, zrange );
      k = pos * p->getNormal();
      dist = k - p->getDisplacement();
      if( dist < contactThreshold )
      {
         WALBERLA_LOG_DETAIL(
                   "      Contact created between box " << b->getID()
                << " and plane " << p->getID() << " (dist=" << dist << ")");

         container.push_back( Contact(b, p, pos, p->getNormal(), dist) );
         collision = true;
      }
   }

   // Test of the upper-right-back corner
   if( xrange*nrel[0] + yrange*nrel[1] + zrange*nrel[2] < minlength )
   {
      pos = bpos + R*Vec3( xrange, yrange, zrange );
      k = pos * p->getNormal();
      dist = k - p->getDisplacement();
      if( dist < contactThreshold )
      {
         WALBERLA_LOG_DETAIL(
                   "      Contact created between box " << b->getID()
                << " and plane " << p->getID() << " (dist=" << dist << ")");

         container.push_back( Contact(b, p, pos, p->getNormal(), dist) );
         collision = true;
      }
   }

   // Test of the upper-left-back corner
   if( -xrange*nrel[0] + yrange*nrel[1] + zrange*nrel[2] < minlength )
   {
      pos = bpos + R*Vec3( -xrange, yrange, zrange );
      k = pos * p->getNormal();
      dist = k - p->getDisplacement();
      if( dist < contactThreshold )
      {
         WALBERLA_LOG_DETAIL(
                   "      Contact created between box " << b->getID()
                << " and plane " << p->getID() << " (dist=" << dist << ")");

         container.push_back( Contact(b, p, pos, p->getNormal(), dist) );
         collision = true;
      }
   }
   return collision;
}
//*************************************************************************************************

template <typename Container>
inline
bool collide( PlaneID p, BoxID b, Container& container )
{
   return collide(b, p, container);
}


/*!\brief Contact generation between two colliding Capsule primitives.
 *
 * \param c1 The first colliding capsule.
 * \param c2 The second colliding capsule.
 * \param contacts Contact container for the generated contacts.
 * \return void
 *
 * In case of two colliding capsules one or two contact points may be generated. In order to
 * estimate if the two capsules are colliding, the closest points on both centerlines are
 * calculated (see pe::getClosestLineSegmentPoints). If the distance of these two points is
 * smaller than the sum of their radii, the two capsules are colliding. In case they are
 * lying parallel to each other and touch each other along their cylindrical part, two contact
 * points are generated. In all other cases, a single contact point between the calculated
 * closest points is generated.
 */
template <typename Container>
inline
bool collide( CapsuleID c1, CapsuleID c2, Container& container )
{
   bool isContactDetected = false;

   const Mat3& R1( c1->getRotation() );
   const Mat3& R2( c2->getRotation() );

   // Calculating the "up" direction for both capsules that points from
   // the center of mass to the center of the upper cap of each capsule
   const Vec3 c1_up_dir( R1[0], R1[3], R1[6] );
   const Vec3 c2_up_dir( R2[0], R2[3], R2[6] );

   // Calculating the displacement of the center of the upper cap-spheres in world space coordinates
   const Vec3 c1_up( real_t(0.5) * c1->getLength() * c1_up_dir );
   const Vec3 c2_up( real_t(0.5) * c2->getLength() * c2_up_dir );

   // calculate the closest points of the two center lines
   Vec3 cp1, cp2;
   geometry::getClosestLineSegmentPoints( c1->getPosition()+c1_up, c1->getPosition()-c1_up,
                                          c2->getPosition()+c2_up, c2->getPosition()-c2_up, cp1, cp2);

   Vec3 normal( cp1 - cp2 );
   const real_t dist( normal.length() - c1->getRadius() - c2->getRadius() );

   if( dist < contactThreshold )
   {
      normalize(normal);

      // Calculating the relative x-position of the second capsule in the frame of the first capsule
      const real_t c2x( c1_up_dir * ( c2->getPosition() - c1->getPosition() ) );

      // Creating two contact points if the capsules are parallel to each other
      if ( floatIsEqual(std::fabs( c1_up_dir * c2_up_dir ), real_t(1) ) &&
          c2x > -c1->getLength() && c2x < c1->getLength() )
      {
         const real_t k( c1->getRadius() + real_t(0.5) * dist );
         const real_t hl1( real_t(0.5) * c1->getLength() );
         const real_t hl2( real_t(0.5) * c2->getLength() );

         WALBERLA_LOG_DETAIL(
                   "      First contact created between capsule " << c1->getID()
                << " and capsule " << c2->getID() << " (dist=" << dist << ")" );

         // Creating the "upper" contact point in world coordinates
         if( hl1 < c2x + hl2 ) {
            isContactDetected = true;
            container.push_back(  Contact( c1, c2,
                                           c1->getPosition() + c1_up - k*normal,
                                           normal, dist ) );
         }
         else {
            isContactDetected = true;
            container.push_back(  Contact( c1, c2,
                                           c2->getPosition() + hl2 * c1_up_dir + k*normal,
                                           normal, dist ) );
         }

         WALBERLA_LOG_DETAIL(
                   "      Second contact created between capsule " << c1->getID()
                << " and capsule " << c2->getID() << " (dist=" << dist << ")" );

         // Creating the "lower" contact point in world coordinates
         if( -hl1 > c2x - hl2 ) {
            isContactDetected = true;
            container.push_back(  Contact( c1, c2,
                                           c1->getPosition() - c1_up - k*normal,
                                           normal, dist ) );
         }
         else {
            isContactDetected = true;
            container.push_back(  Contact( c1, c2,
                                           c2->getPosition() - hl2 * c1_up_dir + k*normal,
                                           normal, dist ) );
         }
      }

      // Creating a single contact point
      else
      {
         const real_t k( c2->getRadius() + real_t(0.5) * dist );
         const Vec3 gPos( cp2 + normal * k );

         WALBERLA_LOG_DETAIL(
                   "      Single contact created between capsule " << c1->getID()
                << " and capsule " << c2->getID() << " (dist=" << dist << ")" );

         isContactDetected = true;
         container.push_back(  Contact( c1, c2, gPos, normal, dist ) );
      }
   }

   return isContactDetected;
}


//*************************************************************************************************
/*!\brief Contact generation between a Capsule and a Plane.
 *
 * \param c The colliding capsule.
 * \param p The colliding plane.
 * \param container Contact container for the generated contacts.
 * \return true if a contact is detected, false otherwise
 *
 * The collision between a capsule and a plane is handled as two collisions between the two
 * cap spheres and the plane: the contact points are calculated with a projection of the
 * global positions of the spheres' center of mass on the plane's normal. If the length of
 * the projection is smaller than the radius of the cap sphere, a contact point is generated,
 * which results in a maximum of two contact points.
 */
template <typename Container>
inline
bool collide( CapsuleID c, PlaneID p, Container& container )
{
   bool collision = false;

   const Mat3& R( c->getRotation() );

   // Computing the location of the sphere caps of the capsule
   const Vec3 c_up ( real_t(0.5) * c->getLength() * Vec3( R[0], R[3], R[6] ) );
   const Vec3 posUp( c->getPosition() + c_up );
   const Vec3 posDn( c->getPosition() - c_up );

   // Computing the distance between the sphere caps and the plane
   const real_t distUp( posUp * p->getNormal() - p->getDisplacement() - c->getRadius() );
   const real_t distDn( posDn * p->getNormal() - p->getDisplacement() - c->getRadius() );

   // Collision of the upper sphere with the plane
   if( distUp < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "      Contact created between capsule " << c->getID()
             << " and plane " << p->getID() << " (dist=" << distUp << ")");

      container.push_back( Contact(c, p, posUp - c->getRadius()*p->getNormal(), p->getNormal(), distUp) );
      collision = true;
   }

   // Collision of the lower sphere with the plane
   if( distDn < contactThreshold )
   {
      WALBERLA_LOG_DETAIL(
                "      Contact created between capsule " << c->getID()
             << " and plane " << p->getID() << " (dist=" << distUp << ")");

      container.push_back( Contact(c, p, posDn - c->getRadius()*p->getNormal(), p->getNormal(), distDn) );
      collision = true;
   }
   return collision;
}
//*************************************************************************************************

template <typename Container>
inline
bool collide( PlaneID p, CapsuleID c, Container& container )
{
   return collide(c, p, container);
}

template <typename Container>
inline
bool collide( SphereID s, CapsuleID c, Container& container )
{
   const real_t length( real_t(0.5)*c->getLength() );  // Half cylinder length
   const Vec3& spos( s->getPosition() );             // Global position of the sphere
   const Vec3& cpos( c->getPosition() );             // Global position of the capsule
   const Mat3& R( c->getRotation() );                // Rotation of the capsule

   // Calculating the relative x-position of the sphere in the frame of the capsule
   real_t sx( R[0]*(spos[0]-cpos[0]) + R[3]*(spos[1]-cpos[1]) + R[6]*(spos[2]-cpos[2]) );

   // Calculation the center of the sphere representing the capsule
   if( sx > length ) {
      sx = length;
   }
   else if( sx < -length ) {
      sx = -length;
   }

   const Vec3 spos2( sx*R[0]+cpos[0], sx*R[3]+cpos[1], sx*R[6]+cpos[2] );

   // Performing a sphere-sphere collision between the colliding sphere and the
   // capsule representation
   Vec3 normal( spos - spos2 );
   const real_t dist( normal.length() - s->getRadius() - c->getRadius() );

   if( dist < contactThreshold ) {
      normal = normal.getNormalized();
      const real_t k( c->getRadius() + real_t(0.5) * dist );
      const Vec3 gPos( spos2 + normal * k );

      WALBERLA_LOG_DETAIL(
                "      Contact created between sphere " << s->getID()
             << " and capsule " << c->getID() << " (dist=" << dist << ")" );

      container.push_back( Contact(s, c, gPos, normal, dist ) );
      return true;
   }
   return false;
}

template <typename Container>
inline
bool collide( CapsuleID c, SphereID s, Container& container )
{
   return collide(s, c, container);
}


template <typename Container>
inline
bool collide( BoxID b, CapsuleID c, Container& container )
{
   const Mat3& R( c->getRotation() );


   //----- Calculating the first contact point between the capsule and the box -----

   // Computing the displacement of the spherical caps of the capsule in world space coordinates
   const Vec3 c_up( real_t(0.5)*c->getLength()*R[0],
                    real_t(0.5)*c->getLength()*R[3],
                    real_t(0.5)*c->getLength()*R[6] );

   // Computing the centers of the spherical caps in world space coordinates
   const Vec3 up  ( c->getPosition()+c_up );
   const Vec3 down( c->getPosition()-c_up );

   // Computing the closest points on the up-down-axis of the cylinder and the box
   Vec3 cp1, bp1;
   geometry::getClosestLineBoxPoints( up, down, b->getPosition(), b->getRotation(), b->getLengths(), cp1, bp1 );

   Vec3 normal1( bp1 - cp1 );
   real_t dist1( normal1.length() - c->getRadius() );

   // Checking the distance between the capsule and the box
   if( dist1 > contactThreshold ) return false;

   // Calculating the contact normal and the position of the closest contact point
   normalize(normal1);
   const Vec3 gpos1( bp1 - real_t(0.5)*dist1*normal1 );

   WALBERLA_LOG_DETAIL(
             "      Contact created between box " << b->getID() << " and capsule " << c->getID()
          << " (dist=" << dist1 << ", gpos=" << gpos1 << ", normal=" << normal1 << ")" );
   container.push_back( Contact( b, c, gpos1, normal1, dist1 ) );


   //----- Calculating the second contact point between the capsule and the box -----

   // Computing the closest points on the down-up-axis of the cylinder and the box
   Vec3 cp2, bp2;
   geometry::getClosestLineBoxPoints( down, up, b->getPosition(), b->getRotation(), b->getLengths(), cp2, bp2 );

   // Early exit in case the same contact point is found again
   if( cp1 == cp2 ) return true;

   Vec3 normal2( bp2 - cp2 );
   real_t dist2( normal2.length() - c->getRadius() );

   // Calculating the contact normal and the position of the closest contact point
   normalize(normal2);
   const Vec3 gpos2( bp2 - real_t(0.5)*dist2*normal2 );

   WALBERLA_LOG_DETAIL(
             "      Contact created between box " << b->getID() << " and capsule " << c->getID()
          << " (dist=" << dist2 << ", gpos=" << gpos2 << ", normal=" << normal2 << ")" );
   container.push_back( Contact( b, c, gpos2, normal2, dist2 ) );
   return true;
}

template <typename Container>
inline
bool collide( CapsuleID c, BoxID b, Container& container )
{
   return collide(b, c, container);
}

template <typename... BodyTypes, typename BodyB, typename Container>
inline
bool collide( Union<BodyTypes...>* bd1, BodyB* bd2, Container& container )
{
   AnalyticSingleCollideFunctor<BodyB, Container> func(bd2, container);
   bool collision = false;
   for( auto it=bd1->begin(); it!=bd1->end(); ++it )
   {
      collision |= SingleCast<std::tuple<BodyTypes...>, AnalyticSingleCollideFunctor<BodyB, Container>, bool>::execute(it.getBodyID(), func);
   }
   return collision;
}

template <typename BodyA, typename... BodyTypes, typename Container>
inline
bool collide( BodyA* bd1, Union<BodyTypes...>* bd2, Container& container )
{
   return collide (bd2, bd1, container);
}

template <typename... BodyTypesA, typename... BodyTypesB, typename Container>
inline
bool collide( Union<BodyTypesA...>* bd1, Union<BodyTypesB...>* bd2, Container& container )
{
   AnalyticCollideFunctor<Container> func(container);
   bool collision = false;
   for( auto it1=bd1->begin(); it1!=bd1->end(); ++it1 )
   {
      for( auto it2=bd2->begin(); it2!=bd2->end(); ++it2 )
      {
         collision |= DoubleCast<std::tuple<BodyTypesA...>, std::tuple<BodyTypesB...>, AnalyticCollideFunctor<Container>, bool>::execute(it1.getBodyID(), it2.getBodyID(), func);
      }
   }
   return collision;
}

} //namespace analytic
} //namespace fcd
} //namespace pe
} //namespace walberla
