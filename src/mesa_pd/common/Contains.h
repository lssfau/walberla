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
//! \file Contains.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>

namespace walberla {
namespace mesa_pd {

/*
 * contains functionality
 */

bool isPointInsideSphere(const Vector3<real_t>& point,
                         const Vector3<real_t>& spherePosition, const real_t sphereRadius )
{
   return !((point - spherePosition).sqrLength() > sphereRadius * sphereRadius);
}

bool isPointInsideHalfSpace(const Vector3<real_t>& point,
                            const Vector3<real_t>& halfSpacePosition, const Vector3<real_t>& halfSpaceNormal )
{
   return !((point - halfSpacePosition) * halfSpaceNormal > real_t(0));
}

//TODO add ellipsoids

struct ContainsPointFunctor
{
   template<typename ParticleAccessor_T, typename Shape_T>
   bool operator()(const size_t /*particleIdx*/, const Shape_T& /*shape*/, const ParticleAccessor_T& /*ac*/, const Vector3<real_t>& /*point*/ )
   {
      WALBERLA_ABORT("ContainsPoint not implemented!");
   }

   template<typename ParticleAccessor_T>
   bool operator()(const size_t particleIdx, const mesa_pd::data::Sphere& sphere, const ParticleAccessor_T& ac, const Vector3<real_t>& point )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return isPointInsideSphere(point, ac.getPosition(particleIdx), sphere.getRadius() );
   }

   template<typename ParticleAccessor_T>
   bool operator()(const size_t particleIdx, const mesa_pd::data::HalfSpace& halfSpace, const ParticleAccessor_T& ac, const Vector3<real_t>& point )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      return isPointInsideHalfSpace(point, ac.getPosition(particleIdx), halfSpace.getNormal() );
   }

};


} //namespace mesa_pd
} //namespace walberla
