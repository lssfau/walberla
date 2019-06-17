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
//! \file LinearizedCompareFunctor.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "LinearizedCompareFunctor.h"

#include <core/DataTypes.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/math/Vector3.h>
#include <core/mpi/MPIManager.h>

#include <stack>
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace sorting {

LinearizedCompareFunctor::LinearizedCompareFunctor(const math::AABB& domain, const Vector3<uint_t> cells)
   : domain_(domain)
   , cells_(cells)
{
   inverse_dx[0] = real_c(cells_[0]) / domain_.xSize();
   inverse_dx[1] = real_c(cells_[1]) / domain_.ySize();
   inverse_dx[2] = real_c(cells_[2]) / domain_.zSize();
}

bool LinearizedCompareFunctor::operator()(const data::Particle p1, const data::Particle p2) const
{
   const auto hash1 = discretize(p1.getPosition() - domain_.minCorner());
   WALBERLA_ASSERT_LESS( hash1, cells_[0]*cells_[1]*cells_[2]);
   const auto hash2 = discretize(p2.getPosition() - domain_.minCorner());
   WALBERLA_ASSERT_LESS( hash2, cells_[0]*cells_[1]*cells_[2]);

   return hash1 < hash2;
}

uint_t LinearizedCompareFunctor::discretize(const Vec3& pos) const
{
   int x = int_c(pos[0] * inverse_dx[0]);
   int y = int_c(pos[1] * inverse_dx[1]);
   int z = int_c(pos[2] * inverse_dx[2]);
   if (x<0) x=0;
   if (y<0) y=0;
   if (z<0) z=0;
   if (x>=int_c(cells_[0])) x=int_c(cells_[0])-1;
   if (y>=int_c(cells_[1])) y=int_c(cells_[1])-1;
   if (z>=int_c(cells_[2])) z=int_c(cells_[2])-1;

   return uint_c(z) * cells_[1] * cells_[0] + uint_c(y) * cells_[0] + uint_c(x);
}

} //namespace sorting
} //namespace mesa_pd
} //namespace walberla
