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
//! \file RayParticleIntersection.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/common/Contains.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>
#include <mesa_pd/kernel/SingleCast.h>

namespace walberla {
namespace mesa_pd {

/*
 * ray - particle intersection ratio functionality
 */

real_t raySphereIntersectionRatio( const Vector3<real_t> & rayOrigin, const Vector3<real_t> & rayDirection,
                                   const Vector3<real_t> & spherePosition, const real_t sphereRadius )
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

   WALBERLA_ASSERT_GREATER_EQUAL( delta, real_t( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( delta, real_t( 1 ) );

   return delta;
}

real_t rayHalfSpaceIntersectionRatio( const Vector3<real_t> & rayOrigin, const Vector3<real_t> & rayDirection,
                                      const Vector3<real_t> & halfSpacePosition, const Vector3<real_t> & halfSpaceNormal, const math::Matrix3<real_t>& halfSpaceRotationMatrix)
{
   WALBERLA_ASSERT( !isPointInsideHalfSpace( rayOrigin, halfSpacePosition, halfSpaceNormal, halfSpaceRotationMatrix ), "rayOrigin: " << rayOrigin );
   WALBERLA_ASSERT(  isPointInsideHalfSpace( rayOrigin + rayDirection, halfSpacePosition, halfSpaceNormal, halfSpaceRotationMatrix ), "rayOrigin + rayDirection: " << rayOrigin + rayDirection );

   const Vector3<real_t> planeNormal( halfSpaceRotationMatrix * halfSpaceNormal );

   real_t denom = planeNormal * rayDirection;

   auto diff = halfSpacePosition - rayOrigin;

   WALBERLA_ASSERT_FLOAT_UNEQUAL(denom, real_t(0));

   real_t delta = diff * planeNormal / denom;

   WALBERLA_ASSERT_GREATER_EQUAL( delta, real_t( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( delta, real_t( 1 ) );

   return delta;
}

//TODO add ellipsoids from pe_coupling/geometry/PeIntersectionRatio.cpp

template <typename ParticleAccessor_T>
real_t intersectionRatioBisection( const size_t particleIdx, const ParticleAccessor_T& ac,
                                   const Vector3<real_t>& rayOrigin, const Vector3<real_t>& rayDirection,
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
      qDelta *= real_t( 0.5 );
      Vector3<real_t> p = rayOrigin + q * rayDirection;
      if( singleCast(particleIdx, ac, containsPointFctr, ac, p) )
      {
         q -= qDelta;
      }
      else
      {
         q += qDelta;
      }
   }

   WALBERLA_ASSERT_GREATER_EQUAL( q, real_t( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( q, real_t( 1 ) );

   return q;
}

struct RayParticleIntersectionRatioFunctor
{
   template<typename ParticleAccessor_T, typename Shape_T>
   real_t operator()(const size_t particleIdx, const Shape_T& /*shape*/, const ParticleAccessor_T& ac, const Vector3<real_t>& rayOrigin, const Vector3<real_t>& rayDirection, real_t epsilon )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return intersectionRatioBisection(particleIdx, ac, rayOrigin, rayDirection, epsilon );
   }

   template<typename ParticleAccessor_T>
   real_t operator()(const size_t particleIdx, const mesa_pd::data::Sphere& sphere, const ParticleAccessor_T& ac, const Vector3<real_t>& rayOrigin, const Vector3<real_t>& rayDirection, real_t /*epsilon*/ )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return raySphereIntersectionRatio(rayOrigin, rayDirection, ac.getPosition(particleIdx), sphere.getRadius() );
   }

   template<typename ParticleAccessor_T>
   real_t operator()(const size_t particleIdx, const mesa_pd::data::HalfSpace& halfSpace, const ParticleAccessor_T& ac, const Vector3<real_t>& rayOrigin, const Vector3<real_t>& rayDirection, real_t /*epsilon*/ )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return rayHalfSpaceIntersectionRatio(rayOrigin, rayDirection, ac.getPosition(particleIdx), halfSpace.getNormal(), ac.getRotation(particleIdx).getMatrix() );
   }

};

} //namespace mesa_pd
} //namespace walberla
