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
//! \file EllipsoidWithOverlap.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \brief Wrapper class that provides a "contains" function for MESA-PD ellipsoids
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "geometry/bodies/BodyOverlapFunctions.h"

#include "mesa_pd/common/Contains.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/shape/Ellipsoid.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{

template< typename ParticleAccessor_T >
class EllipsoidWithOverlap : public geometry::AbstractBody
{
 public:
   EllipsoidWithOverlap(const size_t idx, const shared_ptr< ParticleAccessor_T >& ac,
                        const mesa_pd::data::Ellipsoid& ellipsoid)
      : idx_(idx), ac_(ac), ellipsoid_(ellipsoid)
   {}

   bool contains(const Vector3< real_t >& point) const
   {
      return mesa_pd::isPointInsideEllipsoidBF(mesa_pd::transformPositionFromWFtoBF(idx_, ac_, point),
                                               ellipsoid_.getSemiAxes());
   }

 private:
   size_t idx_;
   shared_ptr< ParticleAccessor_T > ac_;
   mesa_pd::data::Ellipsoid ellipsoid_;
};

} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
