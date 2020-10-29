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
//! \file QCriterion.h
//! \ingroup lbm
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"


// Back-end for calculating macroscopic values
// You should never use these functions directly, always refer to the member functions
// of PdfField or the free functions that can be found in MacroscopicValueCalculation.h

namespace walberla {
namespace lbm {

struct Vorticity {
   template< typename VelocityField_T, typename Filter_T >
   static inline Vector3<real_t> get(const VelocityField_T & velocityField, const Filter_T & filter,
           const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
           real_t dx = real_t(1), real_t dy = real_t(1), real_t dz = real_t(1)) {
      const auto one = cell_idx_t(1);

      if(filter(x,y,z) && filter(x+one,y,z) && filter(x-one,y,z) && filter(x,y+one,z)
         && filter(x,y-one,z) && filter(x,y,z+one) && filter(x,y,z-one)) {
         const Vector3<real_t> xa = velocityField.get(x+one,y,z);
         const Vector3<real_t> xb = velocityField.get(x-one,y,z);
         const Vector3<real_t> ya = velocityField.get(x,y+one,z);
         const Vector3<real_t> yb = velocityField.get(x,y-one,z);
         const Vector3<real_t> za = velocityField.get(x,y,z+one);
         const Vector3<real_t> zb = velocityField.get(x,y,z-one);

         return calculate(xa, xb, ya, yb, za, zb, dx, dy, dz);
      }

      return Vector3<real_t>(real_t(0));
   }

   static inline Vector3<real_t> calculate(const Vector3<real_t> xa, const Vector3<real_t> xb,
           const Vector3<real_t> ya, const Vector3<real_t> yb,
           const Vector3<real_t> za, const Vector3<real_t> zb,
           const real_t dx, const real_t dy, const real_t dz) {

      const auto halfInvDx = real_t(0.5) / dx;
      const auto halfInvDy = real_t(0.5) / dy;
      const auto halfInvDz = real_t(0.5) / dz;

      const real_t duxdy = (ya[0] - yb[0]) * halfInvDy;
      const real_t duxdz = (za[0] - zb[0]) * halfInvDz;

      const real_t duydx = (xa[1] - xb[1]) * halfInvDx;
      const real_t duydz = (za[1] - zb[1]) * halfInvDz;

      const real_t duzdx = (xa[2] - xb[2]) * halfInvDx;
      const real_t duzdy = (ya[2] - yb[2]) * halfInvDy;

      const Vector3< real_t > curl( duzdy - duydz, duxdz - duzdx, duydx - duxdy );

      return curl;
   }
};


} // namespace lbm
} // namespace walberla
