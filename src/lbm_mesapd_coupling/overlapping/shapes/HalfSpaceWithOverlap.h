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
//! \file HalfSpaceWithOverlap.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \brief Wrapper class that provides a "contains" function for MESA-PD half spaces
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "geometry/bodies/BodyOverlapFunctions.h"

#include "mesa_pd/common/Contains.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/shape/HalfSpace.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{

template< typename ParticleAccessor_T >
class HalfSpaceWithOverlap : public geometry::AbstractBody
{
 public:
   HalfSpaceWithOverlap(const size_t idx, const shared_ptr< ParticleAccessor_T >& ac,
                        const mesa_pd::data::HalfSpace& halfSpace)
      : idx_(idx), ac_(ac), halfSpace_(halfSpace)
   {}

   bool contains(const Vector3< real_t >& point) const
   {
      return mesa_pd::isPointInsideHalfSpace(point, ac_->getPosition(idx_), halfSpace_.getNormal());
   }

 private:
   size_t idx_;
   shared_ptr< ParticleAccessor_T > ac_;
   mesa_pd::data::HalfSpace halfSpace_;
};

} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
