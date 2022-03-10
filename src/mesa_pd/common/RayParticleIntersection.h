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
//! \file
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/common/Contains.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>
#include <mesa_pd/kernel/SingleCast.h>

namespace walberla {
namespace mesa_pd {

/*
 * ray - particle intersection ratio functionality
 * can either be formulated in world frame coordinates (then the rotation of the geometry is not taken into account)
 * or in body frame coordinates (BF) which requires the point to be first transformed
 */

real_t raySphereIntersectionRatio( const Vec3& rayOrigin, const Vec3& rayDirection,
                                   const Vec3& spherePosition, const real_t sphereRadius )
{
   WALBERLA_ASSERT( !isPointInsideSphere(rayOrigin, spherePosition, sphereRadius ), "rayOrigin: " << rayOrigin );
   WALBERLA_ASSERT(  isPointInsideSphere(rayOrigin + rayDirection, spherePosition, sphereRadius ), "rayOrigin + rayDirection: " << rayOrigin + rayDirection );

   // get the physical sphere center

   real_t dirLength = rayDirection.length();
   Vector3< real_t > l = spherePosition - rayOrigin;
   real_t s = l * rayDirection / dirLength;
   real_t lsqr = l.sqrLength();
   real_t msqr = lsqr - s * s;
   real_t rsqr = real_c( sphereRadius * sphereRadius );
   real_t q = std::sqrt( rsqr - msqr );
   real_t delta = ( lsqr > rsqr ) ? s - q : s + q;
   delta /= dirLength;

   WALBERLA_ASSERT_GREATER_EQUAL( delta, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( delta, 1_r );

   return delta;
}

real_t rayHalfSpaceIntersectionRatio( const Vec3& rayOrigin, const Vec3& rayDirection,
                                      const Vec3& halfSpacePosition, const Vec3& halfSpaceNormal)
{
   WALBERLA_ASSERT( !isPointInsideHalfSpace( rayOrigin, halfSpacePosition, halfSpaceNormal ), "rayOrigin: " << rayOrigin );
   WALBERLA_ASSERT(  isPointInsideHalfSpace( rayOrigin + rayDirection, halfSpacePosition, halfSpaceNormal ), "rayOrigin + rayDirection: " << rayOrigin + rayDirection );

   real_t denom = halfSpaceNormal * rayDirection;

   auto diff = halfSpacePosition - rayOrigin;

   WALBERLA_ASSERT_FLOAT_UNEQUAL(denom, 0_r);

   real_t delta = diff * halfSpaceNormal / denom;

   WALBERLA_ASSERT_GREATER_EQUAL( delta, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( delta, 1_r );

   return delta;
}

real_t rayEllipsoidIntersectionRatioBF( const Vec3& rayOriginBF, const Vec3& rayDirectionBF,
                                        const Vec3& ellipsoidSemiAxes)
{
   WALBERLA_ASSERT( !isPointInsideEllipsoidBF( rayOriginBF, ellipsoidSemiAxes ), "rayOriginBF: " << rayOriginBF );
   WALBERLA_ASSERT(  isPointInsideEllipsoidBF( rayOriginBF + rayDirectionBF, ellipsoidSemiAxes ), "rayOriginBF + rayDirectionBF: " << rayOriginBF + rayDirectionBF );

   Matrix3<real_t> M = Matrix3<real_t>::makeDiagonalMatrix(1_r/ellipsoidSemiAxes[0], 1_r/ellipsoidSemiAxes[1], 1_r/ellipsoidSemiAxes[2]);

   Vec3 P_M = M*rayOriginBF;
   Vec3 d_M = M*rayDirectionBF;

   const real_t a = d_M*d_M;
   const real_t b = 2_r*P_M*d_M;
   const real_t c = P_M*P_M - 1_r;

   const real_t discriminant = b*b - 4_r*a*c;

   WALBERLA_ASSERT_GREATER_EQUAL(discriminant, 0_r, "No intersection possible!");
   WALBERLA_ASSERT_FLOAT_UNEQUAL(a, 0_r);

   const real_t root = std::sqrt(discriminant);
   real_t delta = (-b - root) / (2_r * a);

   WALBERLA_ASSERT_GREATER_EQUAL( delta, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( delta - 1_r, math::Limits<real_t>::accuracy(), delta );

   return std::min(std::max(delta, real_c(0)), real_c(1));

}

template <typename ParticleAccessor_T>
real_t intersectionRatioBisection( const size_t particleIdx, const ParticleAccessor_T& ac,
                                   const Vec3& rayOrigin, const Vec3& rayDirection,
                                   real_t epsilon)
{
   mesa_pd::kernel::SingleCast singleCast;
   ContainsPointFunctor containsPointFctr;

   WALBERLA_ASSERT( !singleCast(particleIdx, ac, containsPointFctr, ac, rayOrigin), "rayOrigin: " << rayOrigin );
   WALBERLA_ASSERT(  singleCast(particleIdx, ac, containsPointFctr, ac, rayOrigin + rayDirection), "rayOrigin + rayDirection: " << rayOrigin + rayDirection );

   const real_t sqEpsilon         = epsilon * epsilon;
   const real_t sqDirectionLength = rayDirection.sqrLength();

   real_t q( 0.5 );
   real_t qDelta( 0.5 );

   while( qDelta * qDelta * sqDirectionLength >= sqEpsilon )
   {
      qDelta *= 0.5_r;
      Vec3 p = rayOrigin + q * rayDirection;
      if( singleCast(particleIdx, ac, containsPointFctr, ac, p) )
      {
         q -= qDelta;
      }
      else
      {
         q += qDelta;
      }
   }

   WALBERLA_ASSERT_GREATER_EQUAL( q, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( q, 1_r );

   return q;
}

struct RayParticleIntersectionRatioFunctor
{
   template<typename ParticleAccessor_T, typename Shape_T>
   real_t operator()(const size_t particleIdx, const Shape_T& /*shape*/, const ParticleAccessor_T& ac, const Vec3& rayOrigin, const Vec3& rayDirection, real_t epsilon )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return intersectionRatioBisection(particleIdx, ac, rayOrigin, rayDirection, epsilon );
   }

   template<typename ParticleAccessor_T>
   real_t operator()(const size_t particleIdx, const mesa_pd::data::Sphere& sphere, const ParticleAccessor_T& ac, const Vec3& rayOrigin, const Vec3& rayDirection, real_t /*epsilon*/ )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return raySphereIntersectionRatio(rayOrigin, rayDirection, ac.getPosition(particleIdx), sphere.getRadius() );
   }

   template<typename ParticleAccessor_T>
   real_t operator()(const size_t particleIdx, const mesa_pd::data::HalfSpace& halfSpace, const ParticleAccessor_T& ac, const Vec3& rayOrigin, const Vec3& rayDirection, real_t /*epsilon*/ )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return rayHalfSpaceIntersectionRatio(rayOrigin, rayDirection, ac.getPosition(particleIdx), halfSpace.getNormal() );
   }

   template<typename ParticleAccessor_T>
   real_t operator()(const size_t particleIdx, const mesa_pd::data::Ellipsoid& ellipsoid, const ParticleAccessor_T& ac, const Vec3& rayOrigin, const Vec3& rayDirection, real_t /*epsilon*/ )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return rayEllipsoidIntersectionRatioBF(mesa_pd::transformPositionFromWFtoBF(particleIdx, ac, rayOrigin), mesa_pd::transformVectorFromWFtoBF(particleIdx, ac, rayDirection), ellipsoid.getSemiAxes() );
   }

};

} //namespace mesa_pd
} //namespace walberla
