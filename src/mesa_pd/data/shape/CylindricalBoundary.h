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
//! \file CylindricalBoundary.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/shape/BaseShape.h>

#include <core/math/Constants.h>

namespace walberla {
namespace mesa_pd {
namespace data {

class CylindricalBoundary : public BaseShape
{
public:
   explicit CylindricalBoundary(const real_t& radius, const Vec3& axis)
      : BaseShape(CylindricalBoundary::SHAPE_TYPE)
      , radius_(radius)
      , axis_(axis)
   {
      WALBERLA_CHECK_FLOAT_EQUAL(axis.sqrLength(), real_t(1.0), "Axis has to be normalized!");
   }

   const real_t& getRadius() const { return radius_; }
   const Vec3&   getAxis() const   { return axis_; }

   real_t getVolume() const override { return std::numeric_limits<real_t>::infinity(); };
   void   updateMassAndInertia(const real_t density) override;

   constexpr static int SHAPE_TYPE = 2; ///< Unique shape type identifier for cylindrical boundaries.\ingroup mesa_pd_shape

private:
      real_t radius_; ///< radius of the cylinder
      Vec3   axis_;   ///< axis of the cylinder
};

inline
void CylindricalBoundary::updateMassAndInertia(const real_t /*density*/)
{
   invMass_      = real_t(0.0);
   invInertiaBF_ = Mat3(real_t(0.0));
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
