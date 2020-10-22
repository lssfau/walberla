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
 * The HalfSpace does not obey any rotation. The rotation is given purely by the normal of the half space!
 */
class HalfSpace : public BaseShape
{
public:
   explicit HalfSpace(const Vec3& normal = Vec3(real_t(1), real_t(0), real_t(0)))
      : BaseShape(HalfSpace::SHAPE_TYPE), normal_(normal.getNormalized())
   { updateMassAndInertia(real_t(1.0)); }

   void updateMassAndInertia(const real_t density) override;

   real_t getVolume() const override { return std::numeric_limits<real_t>::infinity(); }

   const Vec3& getNormal() const { return normal_; }

   void pack(walberla::mpi::SendBuffer& buf) override;
   void unpack(walberla::mpi::RecvBuffer& buf) override;

   constexpr static int SHAPE_TYPE = 0; ///< Unique shape type identifier for half spaces.\ingroup mesa_pd_shape
private:
   /**
    * Normal of the half space in reference to the global world frame.
    *
    * The normal of the half space is always pointing from the occupied half space towards the open half space.
   **/
   Vec3   normal_;
};

inline
void HalfSpace::updateMassAndInertia(const real_t /*density*/)
{
   mass_         = std::numeric_limits<real_t>::infinity();
   invMass_      = real_t(0);

   inertiaBF_    = Mat3(std::numeric_limits<real_t>::infinity());
   invInertiaBF_ = Mat3(real_t(0));
}

inline
void HalfSpace::pack(walberla::mpi::SendBuffer& buf)
{
   BaseShape::pack(buf);
   buf << normal_;
}
inline
void HalfSpace::unpack(walberla::mpi::RecvBuffer& buf)
{
   BaseShape::unpack(buf);
   buf >> normal_;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
