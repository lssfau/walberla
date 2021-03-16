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
//! \file GeometricalFunctions.cpp
//! \ingroup geometry
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Utility functions for geometrical calculations.
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "GeometricalFunctions.h"
#include "core/math/Shims.h"
#include "core/math/Utility.h"


namespace walberla {
namespace geometry {

//=================================================================================================
//
//  GEOMETRY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Find the closest points on a box and a line segment.
 *
 * \param p1 The first end point of the line segment.
 * \param p2 The second end point of the line segment.
 * \param c The center position of the box.
 * \param R The rotation of the box.
 * \param side The box's sides lengths.
 * \param [out] lret The closest point on the line.
 * \param [out] bret The closest point on the box.
 * \return void
 *
 * This function calculates the closest points between a box and a line segment. The variable
 * \a lret will be set to the closest point on the line, the variable \a bret to the closest
 * point on the box. In case \a p1 lies outside the box and the line intersects the box, the
 * intersection point is returned (in both variables). In case \a p1 lies inside the box, both
 * \a lret and \a bret will refer to \a p1.
 *
 * \image html lineBox.png
 * \image latex lineBox.eps "Calculation of the closest points between a line and a box" width=455pt
 */
void getClosestLineBoxPoints( const Vector3<real_t>& p1, const Vector3<real_t>& p2, const Vector3<real_t>& c, const Matrix3<real_t>& R,
                              const Vector3<real_t>& side, Vector3<real_t>& lret, Vector3<real_t>& bret )
{
   WALBERLA_ABORT("This function produces incorrect results! See https://i10git.cs.fau.de/walberla/walberla/-/issues/25");
   using namespace walberla::math;
   //----- Note: All computations will be performed in the box-relative coordinate-system -----


   // Calculating the global and relative direction of the line p1--p2
   const Vector3<real_t> l( p2 - p1 );
   Vector3<real_t> tmp( R * l );

   // Saving the sign of the direction p1--p2
   const real_t sign[] = { math::sign(tmp[0]), math::sign(tmp[1]), math::sign(tmp[2]) };

   // Calculating the absolute values of the direction direction
   const real_t v [] = { sign[0]*tmp[0], sign[1]*tmp[1], sign[2]*tmp[2] };
   const real_t v2[] = { sq( v[0] )   , sq( v[1] )   , sq( v[2] )    };

   // Calculating the starting point of the line p1--p2 in box-relative coordinates
   tmp = p1 - c;
   const real_t s[] = { sign[0]*( R[0]*tmp[0] + R[3]*tmp[1] + R[6]*tmp[2] ),
                      sign[1]*( R[1]*tmp[0] + R[4]*tmp[1] + R[7]*tmp[2] ),
                      sign[2]*( R[2]*tmp[0] + R[5]*tmp[1] + R[8]*tmp[2] ) };

   // Calculating the half lengths of the box
   const real_t h[] = { real_t(0.5)*side[0], real_t(0.5)*side[1], real_t(0.5)*side[2] };


   // Estimating the region of the starting point depending on which side of the
   // box planes the coordinates are (-1,0,1). 'tanchor' stores the next t value
   // causing a transition from one to another region, or the last one if there
   // are no more. Additionally, d|d|^2/dt is computed for t=0. If it is >= 0
   // then p1 is the closest point since the line points away from the box.
   int  region [] = { 0, 0, 0 };
   real_t tanchor[] = { 2, 2, 2 };  // Invalid t values; t cannot be greater than 1
   real_t dd2dt( 0 );

   if( v[0] > Limits<real_t>::epsilon() )
   {
      if( s[0] < -h[0] ) {
         region[0]  = -1;
         tanchor[0] = ( -h[0]-s[0] ) / v[0];
         dd2dt -= v2[0] * tanchor[0];
      }
      else if( s[0] > h[0] ) {
         region[0] = 1;
         tanchor[0] = ( h[0]-s[0] ) / v[0];
         dd2dt -= v2[0] * tanchor[0];
      }
      else {
         tanchor[0] = ( h[0]-s[0] ) / v[0];
      }
   }

   if( v[1] > Limits<real_t>::epsilon() )
   {
      if( s[1] < -h[1] ) {
         region[1]  = -1;
         tanchor[1] = ( -h[1]-s[1] ) / v[1];
         dd2dt -= v2[1] * tanchor[1];
      }
      else if( s[1] > h[1] ) {
         region[1] = 1;
         tanchor[1] = ( h[1]-s[1] ) / v[1];
         dd2dt -= v2[1] * tanchor[1];
      }
      else {
         tanchor[1] = ( h[1]-s[1] ) / v[1];
      }
   }

   if( v[2] > Limits<real_t>::epsilon() )
   {
      if( s[2] < -h[2] ) {
         region[2]  = -1;
         tanchor[2] = ( -h[2]-s[2] ) / v[2];
         dd2dt -= v2[2] * tanchor[2];
      }
      else if( s[2] > h[2] ) {
         region[2] = 1;
         tanchor[2] = ( h[2]-s[2] ) / v[2];
         dd2dt -= v2[2] * tanchor[2];
      }
      else {
         tanchor[2] = ( h[2]-s[2] ) / v[2];
      }
   }


   // Calculating the value t for the closest point on the line
   real_t t( 0 );

   if( dd2dt < -Limits<real_t>::accuracy() )
   {
      do {
         // Finding the t value for the next clip plane/line intersection
         real_t next_t( 1 );

         if( ( tanchor[0] > t ) && ( tanchor[0] < real_t(1) ) && ( tanchor[0] < next_t ) ) {
            next_t = tanchor[0];
         }
         if( ( tanchor[1] > t ) && ( tanchor[1] < real_t(1) ) && ( tanchor[1] < next_t ) ) {
            next_t = tanchor[1];
         }
         if( ( tanchor[2] > t ) && ( tanchor[2] < real_t(1) ) && ( tanchor[2] < next_t ) ) {
            next_t = tanchor[2];
         }

         // Computing d|d|^2/dt for the next t
         real_t next_dd2dt( 0 );

         if( region[0] != 0 ) {
            next_dd2dt += v2[0] * ( next_t - tanchor[0] );
         }
         if( region[1] != 0 ) {
            next_dd2dt += v2[1] * ( next_t - tanchor[1] );
         }
         if( region[2] != 0 ) {
            next_dd2dt += v2[2] * ( next_t - tanchor[2] );
         }

         // if the sign of d|d|^2/dt has changed, solution = the crossover point
         if( next_dd2dt >= 0 ) {
            t -= dd2dt / ( ( next_dd2dt - dd2dt ) / ( next_t - t ) );
            break;
         }

         // Advancing to the next anchor point / region
         if( floatIsEqual(tanchor[0], next_t) ) {
            tanchor[0] = ( h[0] - s[0] ) / v[0];
            ++region[0];
         }
         if( floatIsEqual(tanchor[1], next_t) ) {
            tanchor[1] = ( h[1] - s[1] ) / v[1];
            ++region[1];
         }
         if( floatIsEqual(tanchor[2], next_t) ) {
            tanchor[2] = ( h[2] - s[2] ) / v[2];
            ++region[2];
         }

         t = next_t;
         dd2dt = next_dd2dt;
      }
      while( t < real_t(1) );
   }

   WALBERLA_ASSERT_GREATER_EQUAL( t, real_t(0), "Invalid line point" );
   WALBERLA_ASSERT_LESS_EQUAL( t, real_t(1), "Invalid line point" );

   // Computing the closest point on the line
   lret = p1 + t * l;

   // Computing the closest point on the box
   tmp[0] = sign[0] * ( s[0] + t * v[0] );
   if     ( tmp[0] < -h[0] ) tmp[0] = -h[0];
   else if( tmp[0] >  h[0] ) tmp[0] =  h[0];

   tmp[1] = sign[1] * ( s[1] + t * v[1] );
   if     ( tmp[1] < -h[1] ) tmp[1] = -h[1];
   else if( tmp[1] >  h[1] ) tmp[1] =  h[1];

   tmp[2] = sign[2] * ( s[2] + t * v[2] );
   if     ( tmp[2] < -h[2] ) tmp[2] = -h[2];
   else if( tmp[2] >  h[2] ) tmp[2] =  h[2];

   bret = ( R * tmp ) + c;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Find the closest points on two line segments.
 *
 * \param a1 The first end point of the first line segment.
 * \param a2 The second end point of the first line segment.
 * \param b1 The first end point of the second line segment.
 * \param b2 The second end point of the second line segment.
 * \param [out] cp1 The closest point on line one.
 * \param [out] cp2 The closest point on line two.
 * \return void
 *
 * Given the two line segments A and B, this function returns the points that are closest to
 * each other. In the case of parallel lines, where multiple solutions are possible, a solution
 * involving the endpoint of at least one line will be returned. This will also work correctly
 * for zero length lines, e.g. if \a a1 = \a a2 and/or \a b1 = \a b2.
 */
void getClosestLineSegmentPoints( const Vector3<real_t>& a1, const Vector3<real_t>& a2, const Vector3<real_t>& b1, const Vector3<real_t>& b2,
                                  Vector3<real_t>& cp1, Vector3<real_t>& cp2 )
{
   using namespace walberla::math;
   //----- Checking the vertex-vertex features -----

   const Vector3<real_t> a1a2( a2 - a1 );
   const Vector3<real_t> b1b2( b2 - b1 );
   const Vector3<real_t> a1b1( b1 - a1 );
   const real_t da1 ( a1a2 * a1b1 );
   const real_t db1 ( b1b2 * a1b1 );
   if( da1 <= +Limits<real_t>::fpuAccuracy() && db1 >= -Limits<real_t>::fpuAccuracy() ) {
      cp1 = a1;
      cp2 = b1;
      return;
   }

   const Vector3<real_t> a1b2( b2 - a1 );
   const real_t da2 ( a1a2 * a1b2 );
   const real_t db2 ( b1b2 * a1b2 );
   if( da2 <= +Limits<real_t>::fpuAccuracy() && db2 <= +Limits<real_t>::fpuAccuracy() ) {
      cp1 = a1;
      cp2 = b2;
      return;
   }

   const Vector3<real_t> a2b1( b1 - a2 );
   const real_t da3 ( a1a2 * a2b1 );
   const real_t db3 ( b1b2 * a2b1 );
   if( da3 >= -Limits<real_t>::fpuAccuracy() && db3 >= -Limits<real_t>::fpuAccuracy() ) {
      cp1 = a2;
      cp2 = b1;
      return;
   }

   const Vector3<real_t> a2b2( b2 - a2 );
   const real_t da4 ( a1a2 * a2b2 );
   const real_t db4 ( b1b2 * a2b2 );
   if( da4 >= -Limits<real_t>::fpuAccuracy() && db4 <= +Limits<real_t>::fpuAccuracy() ) {
      cp1 = a2;
      cp2 = b2;
      return;
   }


   //----- Checking the edge-vertex features -----
   // If one or both of the line segments have zero length, we will never get here. Therefore
   // we don't have to worry about possible divisions by zero in the following calculations.

   Vector3<real_t> n;
   Vector3<real_t> k;

   const real_t la( a1a2 * a1a2 );
   if( da1 >= -Limits<real_t>::fpuAccuracy() && da3 <= +Limits<real_t>::fpuAccuracy() ) {
      k = (da1 / la) * a1a2;
      n = a1b1 - k;
      if( b1b2 * n >= -Limits<real_t>::fpuAccuracy() ) {
         cp1 = a1 + k;
         cp2 = b1;
         return;
      }
   }

   if( da2 >= -Limits<real_t>::fpuAccuracy() && da4 <= +Limits<real_t>::fpuAccuracy() ) {
      k = (da2 / la) * a1a2;
      n = a1b2 - k;
      if( b1b2 * n <= +Limits<real_t>::fpuAccuracy() ) {
         cp1 = a1 + k;
         cp2 = b2;
         return;
      }
   }

   const real_t lb( b1b2 * b1b2 );
   if( db1 <= +Limits<real_t>::fpuAccuracy() && db2 >= -Limits<real_t>::fpuAccuracy() ) {
      k = (-db1 / lb) * b1b2;
      n = -a1b1 - k;
      if( a1a2 * n >= -Limits<real_t>::fpuAccuracy() ) {
         cp1 = a1;
         cp2= b1 + k;
         return;
      }
   }

   if( db3 <= +Limits<real_t>::fpuAccuracy() && db4 >= -Limits<real_t>::fpuAccuracy() ) {
      k = (-db3 / lb) * b1b2;
      n = -a2b1 - k;
      if( a1a2 * n <= +Limits<real_t>::fpuAccuracy() ) {
         cp1 = a2;
         cp2 = b1 + k;
         return;
      }
   }


   //----- Checking the edge-edge features -----

   const real_t scale( a1a2 * b1b2 );
   real_t div( la * lb - math::sq(scale) );

   WALBERLA_ASSERT_GREATER( div, real_t(0), std::setprecision(16) << "Invalid division\n" << a1 << "\n" << a2 << "\n" << b1 << "\n" << b2 );
   div = real_t(1) / div;

   const real_t s( ( lb    * da1 - scale * db1 ) * div );
   const real_t t( ( scale * da1 - la    * db1 ) * div );

   cp1 = a1 + s * a1a2;
   cp2 = b1 + t * b1b2;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Find the intersection point or the point of the closest approach between two straight
 * \brief lines \f$ l_1(s) = o_1+sd_1 \f$ and \f$ l_2(t) = o_2+td_2 \f$.
 * \ingroup geometry
 *
 * \param o1 A point on line \f$ l_1 \f$.
 * \param d1 The direction of line \f$ l_1 \f$.
 * \param o2 A point on line \f$ l_2 \f$.
 * \param d2 The direction of line \f$ l_2 \f$.
 * \param [out] s The resolved parameter for line \f$ l_1 \f$.
 * \param [out] t The resolved parameter for line \f$ l_2 \f$.
 * \return void

   \f[  o_1+sd_1 = o_2+td_2  \f]
   \f[  \Longleftrightarrow  \f]
   \f[  s = \frac{\det(o_2-o_1,d_2,d_1 \times d_2)}{\left \Vert d_1 \times d_2 \right \| ^2 }  \f]
   \f[  t = \frac{\det(o_2-o_1,d_1,d_1 \times d_2)}{\left \Vert d_1 \times d_2 \right \| ^2 }  \f]
 */
void intersectLines( const Vector3<real_t>& o1, const Vector3<real_t>& d1, const Vector3<real_t>& o2, const Vector3<real_t>& d2,
                     real_t& s, real_t& t )
{
   using namespace walberla::math;

   const real_t sqrlen( ( d1 % d2 ).sqrLength() );

   if( floatIsEqual(sqrlen, 0) )
   {
      s = real_t(0);
      t = real_t(0);
   }
   else
   {
      const real_t isqrlen( real_t(1) / sqrlen );
      const Vector3<real_t> p( o2 - o1 );
      const real_t dot(  d1 * d2 );
      const real_t a  (  d1 * p  );
      const real_t b  ( -d2 * p  );

      s = ( a * ( d2 * d2 ) + dot*b ) * isqrlen;
      t = ( b * ( d1 * d1 ) + dot*a ) * isqrlen;
   }
}
//*************************************************************************************************

} // namespace math
} // namespace walberla
