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
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

#include <core/math/Constants.h>
#include <core/logging/Logging.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {
namespace VBondModel {

/**
 * Repulsive wall interaction kernel.
 *
 * Implementation follows:
 * Wittmaack, Volkov, Zhigilei, "Phase transformation as the mechanism of mechanical deformation of vertically aligned CNT arrays" - Carbon, 2019
 * https://doi.org/10.1016/j.carbon.2018.11.066
 *
 * The force is divided into three areas.
 * ========================== z=0, position of the wall
 * close to wall
 * -------------------------- z=z1+r
 * spline interpolation
 * -------------------------- z=z2+r
 * far away from wall
 */
class WallContact
{
public:
   explicit WallContact(real_t zPos) : zPos_(zPos) {}

   template<typename Accessor>
   void operator()(const size_t p_idx1,
                   Accessor &ac);

   static constexpr real_t r = 6.78_r; ///< A
   static constexpr real_t eps = 0.254e-3_r; ///< eV/amu
   static constexpr real_t m = 2648.8_r; ///< amu
   static constexpr real_t s = 3.6_r; ///< A
   static constexpr real_t s12 = ((s * s) * (s * s) * (s * s)) * ((s * s) * (s * s) * (s * s));

   static constexpr real_t z1 = 10_r; ///< A
   static constexpr real_t z2 = 12_r; ///< A

   auto getLastForce() const {return F_;}

   void setZPos(const real_t& zPos) {zPos_ = zPos;}
   auto getZPos() const {return zPos_;}

private:
   real_t zPos_ = 0_r;
   real_t F_ = 0_r; ///< resulting force from the last interaction
};

template<typename Accessor>
inline void WallContact::operator()(const size_t p_idx,
                                    Accessor &ac)
{
   auto dz = ac.getPosition(p_idx)[2] - zPos_;
   F_ = std::copysign(1_r, dz);
   dz = std::abs(dz);

   if (dz < r + z1)
   {
      //close to wall
      const auto tmp = dz - r;
      const auto pow = ((tmp * tmp) * (tmp * tmp) * (tmp * tmp)) * ((tmp * tmp) * (tmp * tmp) * (tmp * tmp)) * tmp;
      F_ *= m * eps * s12 * 12_r / pow;
   } else if (dz < r + z2)
   {
      //cubic spline interpolation
      auto S = [](real_t x) { return (3_r * (12_r + r) * (14_r + r) - 6_r * (13_r + r) * x + 3_r * x * x) * 5e-14_r; };
      F_ *= m * eps * s12 * S(dz);
   } else
   {
      //far away from wall
      F_ = 0_r;
   }

   addForceAtomic( p_idx, ac, Vec3(0_r, 0_r, F_) );
}

} //namespace VBondModel
} //namespace kernel
} //namespace mesa_pd
} //namespace walberla