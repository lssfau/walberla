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
//! \file Sphere.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/shape/BaseShape.h>

#include <core/math/Constants.h>

namespace walberla {
namespace mesa_pd {
namespace data {

class Sphere : public BaseShape
{
public:
   explicit Sphere(const real_t& radius) : BaseShape(Sphere::SHAPE_TYPE), radius_(radius) {}

   const real_t& getRadius() const { return radius_; }

   void updateMassAndInertia(const real_t density) override;

   real_t getVolume() const override { return (real_t(4) / real_t(3)) * math::pi * getRadius() * getRadius() * getRadius(); }

   Vec3 support( const Vec3& d ) const override;

   constexpr static int SHAPE_TYPE = 1; ///< Unique shape type identifier for spheres.\ingroup mesa_pd_shape

private:
      real_t radius_; ///< radius of the sphere
};

inline
void Sphere::updateMassAndInertia(const real_t density)
{
   const real_t m = (real_c(4.0)/real_c(3.0) * math::pi) * getRadius() * getRadius() * getRadius() * density;
   const Mat3   I = Mat3::makeDiagonalMatrix( real_c(0.4) * m * getRadius() * getRadius() );

   mass_         = m;
   invMass_      = real_t(1.0) / m;

   inertiaBF_    = I;
   invInertiaBF_ = I.getInverse();
}

inline
Vec3 Sphere::support( const Vec3& d ) const
{
   return radius_ * d;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
