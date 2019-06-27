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
//! \file Support.h
//! \author Tobias Scharpff
//! \author Tobias Leemann
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/shape/BaseShape.h>

namespace walberla {
namespace mesa_pd {
namespace collision_detection {

class Support
{
public:
   /// \attention Lifetime of shape has to surpass lifetime of Support!!!
   /// No copy is stored.
   Support(const Vec3& pos, const Rot3& rot, const data::BaseShape& shape)
      : pos_(pos)
      , rot_(rot)
      , shape_(shape)
   {}

   /**
    * \brief Estimates the point which is farthest in direction \a d.
    *
    * \param d The normalized search direction in world-frame coordinates.
    * \return The support point in world-frame coordinates in direction a\ d.
    */
   Vec3 support( Vec3 d ) const;
public:
   Vec3 pos_;
   Rot3 rot_;
   const data::BaseShape& shape_;
};

inline
Vec3 Support::support( Vec3 d ) const
{
   auto len = d.sqrLength();
   if (math::equal(len, real_t(0)))
      return Vec3(0,0,0);

   d = rot_.getMatrix().getTranspose() * d.getNormalized(); //vectorFromWFtoBF
   d = shape_.support(d);
   return pos_ + rot_.getMatrix() * d; //vectorFromBFtoWF
}

} // namespace collision_detection
} // namespace mesa_pd
} // namespace walberla
