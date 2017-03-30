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
//! \file GeometricalFunctions.h
//! \ingroup geometry
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Utility functions for geometrical calculations.
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "core/DataTypes.h"

#include "core/math/Limits.h"
#include "core/math/Matrix3.h"
#include "core/math/Vector3.h"

namespace walberla {
namespace geometry {

//=================================================================================================
//
//  GEOMETRY FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\name Geometry functions */
//@{
void getClosestLineBoxPoints( const Vector3<real_t>& p1, const Vector3<real_t>& p2, const Vector3<real_t>& c, const Matrix3<real_t>& R,
                              const Vector3<real_t>& side, Vector3<real_t>& lret, Vector3<real_t>& bret);

void getClosestLineSegmentPoints( const Vector3<real_t>& a1, const Vector3<real_t>& a2, const Vector3<real_t>& b1, const Vector3<real_t>& b2,
                                  Vector3<real_t>& cp1, Vector3<real_t>& cp2 );

void intersectLines( const Vector3<real_t>& o1, const Vector3<real_t>& d1, const Vector3<real_t>& o2, const Vector3<real_t>& d2,
                     real_t& s, real_t& t );

inline bool originInTetrahedron( const Vector3<real_t>& A, const Vector3<real_t>& B, const Vector3<real_t>& C, const Vector3<real_t>& D );
inline bool pointInTetrahedron ( const Vector3<real_t>& A, const Vector3<real_t>& B, const Vector3<real_t>& C, const Vector3<real_t>& D,
                                 const Vector3<real_t>& point );

inline bool pointInFrontOfPlane( const Vector3<real_t>& normal, const Vector3<real_t>& pointOnPlane, const Vector3<real_t>& point );
//@}
//*************************************************************************************************



//=================================================================================================
//
//  GEOMETRY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Estimates whether or not the origin lies within the given tetrahedron.
 * \ingroup geometry
 *
 * \param A Vertex that sees the triangle BCD in counterclockwise order.
 * \param B Vertex that sees the triangle ADC in counterclockwise order.
 * \param C Vertex that sees the triangle ABD in counterclockwise order.
 * \param D Vertex that sees the triangle ACB in counterclockwise order.
 * \return \a true if the origin lies within the tetrahedron, otherwise \a false.
 *
 * \note Points on the surface of the tetrahedron are considered inside.
 * 
 * \todo Review documentation
 */
inline bool originInTetrahedron( const Vector3<real_t>& A, const Vector3<real_t>& B, const Vector3<real_t>& C, const Vector3<real_t>& D ) {
   using namespace walberla::math;

   Vector3<real_t> aoT = A;

   //using fpuAccuracy instead of real(0.0) to avoid numeric problems
   if((aoT * (B % C)) < -Limits<real_t>::fpuAccuracy()) {
      //if volume of ABC and Origin <0.0 than the origin is on the wrong side of ABC
      //http://mathworld.wolfram.com/Tetrahedron.html volume formula
      return false;
   }
   if((aoT * (C % D)) < -Limits<real_t>::fpuAccuracy()) {
      return false;
   }
   if((aoT * (D % B)) < -Limits<real_t>::fpuAccuracy()) {
      return false;
   }
   if((B * (D % C)) < -Limits<real_t>::fpuAccuracy()) {
      return false;
   }
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates whether or not a given point lies within the given tetrahedron.
 * \ingroup geometry
 *
 * \param A Vertex that sees the triangle BCD in counterclockwise order.
 * \param B Vertex that sees the triangle ADC in counterclockwise order.
 * \param C Vertex that sees the triangle ABD in counterclockwise order.
 * \param D Vertex that sees the triangle ACB in counterclockwise order.
 * \param point The point whose position is check for being inside the tetrahedron.
 * \return \a true if the origin lies within the tetrahedron, otherwise \a false.
 *
 * \note Points on the surface of the tetrahedron are considered inside.
 * \todo Review documentation
 */
inline bool pointInTetrahedron( const Vector3<real_t>& A, const Vector3<real_t>& B, const Vector3<real_t>& C, const Vector3<real_t>& D, const Vector3<real_t>& point ) {
   return originInTetrahedron( A-point, B-point, C-point, D-point );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates whether a given point is in front of a plane.
 * \ingroup geometry
 *
 * \param normal The normal of the plane, does not have to be normalised.
 * \param pointOnPlane Any point on the Plane.
 * \param point The point whose position is check for being in front of the plane.
 * \return \a true if the origin lies in front of the plane, otherwise \a false.
 *
 * \note Points on the surface of the plane are considered not in front of the plane.
 * \todo Review documentation
 */
inline bool pointInFrontOfPlane( const Vector3<real_t>& normal, const Vector3<real_t>& pointOnPlane, const Vector3<real_t>& point ) {
   return (normal * (point - pointOnPlane)) > 0;
}
//*************************************************************************************************


} // namespace geometry
} // namespace walberla
