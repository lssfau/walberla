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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/shape/BaseShape.h>

#include <core/math/Constants.h>

namespace walberla {
namespace mesa_pd {
namespace data {

class Ellipsoid : public BaseShape
{
public:
   explicit Ellipsoid(const Vec3& semiAxes = Vec3(real_t(1)))
      : BaseShape(Ellipsoid::SHAPE_TYPE)
      , semiAxes_(semiAxes)
   {}

   const Vec3&   getSemiAxes() const   { return semiAxes_; }

   real_t getVolume() const override { return real_c(4.0/3.0) * math::pi * semiAxes_[0] * semiAxes_[1] * semiAxes_[2]; }
   void   updateMassAndInertia(const real_t density) override;

   Vec3 support( const Vec3& d_loc ) const override;

   void pack(walberla::mpi::SendBuffer& buf) override;
   void unpack(walberla::mpi::RecvBuffer& buf) override;

   constexpr static int SHAPE_TYPE = 4; ///< Unique shape type identifier for boxes.\ingroup mesa_pd_shape

private:
   Vec3   semiAxes_;   ///< edge length of the box
};

inline
void Ellipsoid::updateMassAndInertia(const real_t density)
{
   const real_t m = getVolume() * density;
   const Mat3   I = Mat3::makeDiagonalMatrix(
                       real_c(0.2) * m * (semiAxes_[1] * semiAxes_[1] + semiAxes_[2] * semiAxes_[2]),
         real_c(0.2) * m * (semiAxes_[2] * semiAxes_[2] + semiAxes_[0] * semiAxes_[0]),
         real_c(0.2) * m * (semiAxes_[0] * semiAxes_[0] + semiAxes_[1] * semiAxes_[1]));

   mass_         = m;
   invMass_      = real_t(1.0) / m;
   inertiaBF_    = I;
   invInertiaBF_ = I.getInverse();
}

inline
Vec3 Ellipsoid::support( const Vec3& d_loc ) const
{
   Vec3 norm_vec(d_loc[0] * semiAxes_[0], d_loc[1] * semiAxes_[1], d_loc[2] * semiAxes_[2]);
   real_t norm_length = norm_vec.length();
   Vec3 local_support = (real_t(1.0) / norm_length) * Vec3(semiAxes_[0] * semiAxes_[0] * d_loc[0],
         semiAxes_[1] * semiAxes_[1] * d_loc[1],
         semiAxes_[2] * semiAxes_[2] * d_loc[2]);
   return local_support;
}

inline
void Ellipsoid::pack(walberla::mpi::SendBuffer& buf)
{
   BaseShape::pack(buf);
   buf << semiAxes_;
}
inline
void Ellipsoid::unpack(walberla::mpi::RecvBuffer& buf)
{
   BaseShape::unpack(buf);
   buf >> semiAxes_;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
