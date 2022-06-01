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
//! \file SphereWithOverlap.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \brief Wrapper class that provides a "contains" function for MESA-PD spheres
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "geometry/bodies/BodyOverlapFunctions.h"

#include "mesa_pd/common/Contains.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/shape/Sphere.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{

template< typename ParticleAccessor_T >
class SphereWithOverlap : public geometry::AbstractBody
{
 public:
   SphereWithOverlap(const size_t idx, const shared_ptr< ParticleAccessor_T >& ac, const mesa_pd::data::Sphere& sphere)
      : idx_(idx), ac_(ac), sphere_(sphere)
   {}

   bool contains(const Vector3< real_t >& point) const
   {
      return mesa_pd::isPointInsideSphere(point, ac_->getPosition(idx_), sphere_.getRadius());
   }

 private:
   size_t idx_;
   shared_ptr< ParticleAccessor_T > ac_;
   mesa_pd::data::Sphere sphere_;
};

} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
