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
//! \file PeIntersectionRatio.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Specialization of the intersectionRatio function for certain pe bodies
//
//======================================================================================================================

#include "pe_coupling/geometry/PeIntersectionRatio.h"

#include "core/math/Limits.h"
#include "core/math/Matrix3.h"

#include <algorithm>

namespace walberla {
namespace lbm {

real_t intersectionRatioSpherePe( const pe::Sphere & sphere,
                                  const Vector3<real_t> & fluidPoint,
                                  const Vector3<real_t> & direction )
{
   WALBERLA_ASSERT( !walberla::geometry::contains( sphere, fluidPoint ), "fluidPoint: " << fluidPoint );
   WALBERLA_ASSERT(  walberla::geometry::contains( sphere, fluidPoint + direction ), "fluidPoint + direction: " << fluidPoint + direction );

   // get the physical sphere center
   const Vector3< real_t >& sphereCenter( sphere.getPosition() );

   real_t dirLength = direction.length();
   Vector3< real_t > l = sphereCenter - fluidPoint;
   real_t s = l * direction / dirLength;
   real_t lsqr = l.sqrLength();
   real_t msqr = lsqr - s * s;
   real_t rsqr = real_c( sphere.getRadius() * sphere.getRadius() );
   real_t q = std::sqrt( rsqr - msqr );
   real_t delta = ( lsqr > rsqr ) ? s - q : s + q;
   delta /= dirLength;

   WALBERLA_ASSERT_GREATER_EQUAL( delta, real_t( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( delta - real_t(1), math::Limits<real_t>::accuracy(), delta );

   return std::min(std::max(delta, real_c(0)), real_c(1));

}

real_t intersectionRatioPlanePe( const pe::Plane & plane,
                                 const Vector3<real_t> & fluidPoint,
                                 const Vector3<real_t> & direction )
{
   WALBERLA_ASSERT( !walberla::geometry::contains<pe::RigidBody>( plane, fluidPoint ), "fluidPoint: " << fluidPoint );
   WALBERLA_ASSERT(  walberla::geometry::contains<pe::RigidBody>( plane, fluidPoint + direction ), "fluidPoint + direction: " << fluidPoint + direction );


   const Vector3<real_t>& planeCenter( plane.getPosition() );
   const Vector3<real_t>& planeNormal( plane.getNormal() );

   real_t denom = planeNormal * direction;

   auto diff = planeCenter - fluidPoint;

   WALBERLA_ASSERT_FLOAT_UNEQUAL(denom, real_t(0));

   real_t delta = diff * planeNormal / denom;

   WALBERLA_ASSERT_GREATER_EQUAL( delta, real_t( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( delta - real_t(1), math::Limits<real_t>::accuracy(), delta );

   return std::min(std::max(delta, real_c(0)), real_c(1));

}

real_t intersectionRatioEllipsoidPe( const pe::Ellipsoid & ellipsoid,
                                     const Vector3<real_t> & fluidPoint,
                                     const Vector3<real_t> & direction )
{
   WALBERLA_ASSERT( !walberla::geometry::contains<pe::RigidBody>( ellipsoid, fluidPoint ), "fluidPoint: " << fluidPoint );
   WALBERLA_ASSERT(  walberla::geometry::contains<pe::RigidBody>( ellipsoid, fluidPoint + direction ), "fluidPoint + direction: " << fluidPoint + direction );

   Vector3<real_t> transformedP = ellipsoid.pointFromWFtoBF(fluidPoint);
   Vector3<real_t> transformedDir = ellipsoid.vectorFromWFtoBF(direction);

   Vector3<real_t> semiAxes = ellipsoid.getSemiAxes();

   Matrix3<real_t> M = Matrix3<real_t>::makeDiagonalMatrix(real_t(1)/semiAxes[0], real_t(1)/semiAxes[1], real_t(1)/semiAxes[2]);

   Vector3<real_t> P_M = M*transformedP;
   Vector3<real_t> d_M = M*transformedDir;

   const real_t a = d_M*d_M;
   const real_t b = real_t(2)*P_M*d_M;
   const real_t c = P_M*P_M - real_t(1);

   const real_t discriminant = b*b - real_t(4)*a*c;

   WALBERLA_ASSERT_GREATER_EQUAL(discriminant, real_t(0), "No intersection possible!");
   WALBERLA_ASSERT_FLOAT_UNEQUAL(a, real_t(0));

   const real_t root = std::sqrt(discriminant);
   real_t delta = (-b - root) / (real_t(2) * a);

   WALBERLA_ASSERT_GREATER_EQUAL( delta, real_t( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( delta - real_t(1), math::Limits<real_t>::accuracy(), delta );

   return std::min(std::max(delta, real_c(0)), real_c(1));

}


} // namespace lbm
} // namespace walberla
