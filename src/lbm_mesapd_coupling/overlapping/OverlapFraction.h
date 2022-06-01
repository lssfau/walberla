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
//! \file OverlapFraction.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \brief Functor that provides overlap fraction computations for different MESA-PD shapes (used for SingleCast)
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "geometry/bodies/BodyOverlapFunctions.h"

#include "mesa_pd/common/Contains.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/shape/Sphere.h"

#include "shapes/BoxWithOverlap.h"
#include "shapes/CylindricalBoundaryWithOverlap.h"
#include "shapes/EllipsoidWithOverlap.h"
#include "shapes/HalfSpaceWithOverlap.h"
#include "shapes/SphereWithOverlap.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{

struct OverlapFractionFunctor
{
   template< typename ParticleAccessor_T, typename Shape_T >
   real_t operator()(const size_t /*particleIdx*/, const Shape_T& /*shape*/,
                     const shared_ptr< ParticleAccessor_T >& /*ac*/, const Vector3< real_t >& /*point*/,
                     const Vector3< real_t >& /*dxVec*/, uint_t /*superSamplingDepth*/)
   {
      WALBERLA_ABORT("OverlapFraction not implemented!");
   }

   template< typename ParticleAccessor_T >
   real_t operator()(const size_t particleIdx, const mesa_pd::data::Sphere& sphere,
                     const shared_ptr< ParticleAccessor_T >& ac, const Vector3< real_t >& point,
                     const Vector3< real_t >& dxVec, uint_t superSamplingDepth)
   {
      WALBERLA_STATIC_ASSERT((std::is_base_of< mesa_pd::data::IAccessor, ParticleAccessor_T >::value));

      SphereWithOverlap sphereWithOverlap(particleIdx, ac, sphere);
      return geometry::overlapFraction< geometry::AbstractBody >(sphereWithOverlap, point, dxVec, superSamplingDepth);
   }

   template< typename ParticleAccessor_T >
   real_t operator()(const size_t particleIdx, const mesa_pd::data::HalfSpace& halfSphere,
                     const shared_ptr< ParticleAccessor_T >& ac, const Vector3< real_t >& point, real_t dxVec,
                     uint_t superSamplingDepth)
   {
      WALBERLA_STATIC_ASSERT((std::is_base_of< mesa_pd::data::IAccessor, ParticleAccessor_T >::value));

      HalfSpaceWithOverlap halfSpaceWithOverlap(particleIdx, ac, halfSphere);
      return geometry::overlapFraction< geometry::AbstractBody >(halfSpaceWithOverlap, point, dxVec,
                                                                 superSamplingDepth);
   }

   template< typename ParticleAccessor_T >
   real_t operator()(const size_t particleIdx, const mesa_pd::data::CylindricalBoundary& cylindricalBoundary,
                     const shared_ptr< ParticleAccessor_T >& ac, const Vector3< real_t >& point, real_t dxVec,
                     uint_t superSamplingDepth)
   {
      WALBERLA_STATIC_ASSERT((std::is_base_of< mesa_pd::data::IAccessor, ParticleAccessor_T >::value));

      CylindricalBoundaryWithOverlap cylindricalBoundaryWithOverlap(particleIdx, ac, cylindricalBoundary);
      return geometry::overlapFraction< geometry::AbstractBody >(cylindricalBoundaryWithOverlap, point, dxVec,
                                                                 superSamplingDepth);
   }

   template< typename ParticleAccessor_T >
   real_t operator()(const size_t particleIdx, const mesa_pd::data::Box& box,
                     const shared_ptr< ParticleAccessor_T >& ac, const Vector3< real_t >& point, real_t dxVec,
                     uint_t superSamplingDepth)
   {
      WALBERLA_STATIC_ASSERT((std::is_base_of< mesa_pd::data::IAccessor, ParticleAccessor_T >::value));

      BoxWithOverlap boxWithOverlap(particleIdx, ac, box);
      return geometry::overlapFraction< geometry::AbstractBody >(boxWithOverlap, point, dxVec, superSamplingDepth);
   }

   template< typename ParticleAccessor_T >
   real_t operator()(const size_t particleIdx, const mesa_pd::data::Ellipsoid& ellipsoid,
                     const shared_ptr< ParticleAccessor_T >& ac, const Vector3< real_t >& point, real_t dxVec,
                     uint_t superSamplingDepth)
   {
      WALBERLA_STATIC_ASSERT((std::is_base_of< mesa_pd::data::IAccessor, ParticleAccessor_T >::value));

      EllipsoidWithOverlap ellipsoidWithOverlap(particleIdx, ac, ellipsoid);
      return geometry::overlapFraction< geometry::AbstractBody >(ellipsoidWithOverlap, point, dxVec,
                                                                 superSamplingDepth);
   }
};

} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
