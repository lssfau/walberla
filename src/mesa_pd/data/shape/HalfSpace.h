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
//! \file HalfSpace.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/shape/BaseShape.h>

namespace walberla {
namespace mesa_pd {
namespace data {

/**
 * Half space shape class.
 *
 * \attention
 * The HalfSpace does not obay any rotation. The rotation is given purely be the normal of the half space!
 */
class HalfSpace : public BaseShape
{
public:
   explicit HalfSpace(const Vec3& normal) : BaseShape(HalfSpace::SHAPE_TYPE), normal_(normal) { updateMassAndInertia(real_t(1.0)); }

   void updateMassAndInertia(const real_t density) override;

   real_t getVolume() const override { return std::numeric_limits<real_t>::infinity(); }

   const Vec3& getNormal() const { return normal_; }

   static const int SHAPE_TYPE = 0; ///< Unique shape type identifier for planes.\ingroup mesa_pd_shape
private:
   /**
    * Normal of the plane in reference to the global world frame.
    *
    * The normal of the plane is always pointing towards the halfspace outside the plane.
   **/
   Vec3   normal_;
};

inline
void HalfSpace::updateMassAndInertia(const real_t /*density*/)
{
   getInvMass()      = real_t(0);
   getInvInertiaBF() = Mat3(real_t(0));
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
