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
//! \file Overlap.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "Overlap.h"

#include "core/debug/Debug.h"
#include <core/math/Constants.h>

namespace walberla {
namespace pe {

real_t getSphereSphereOverlap(const real_t d, const real_t r1, const real_t r2)
{
   WALBERLA_ASSERT_GREATER(d, r1);
   WALBERLA_ASSERT_GREATER(d, r2);

   //http://math.stackexchange.com/questions/297751/overlapping-spheres
   if (d > r1+r2)
      return 0; else
      return math::pi / (real_c(12.0)  * d) * (r1 + r2 - d) * (r1 + r2 - d) * (d*d + 2*d*(r1+r2) - 3*(r1-r2)*(r1-r2));
}

}  // namespace pe
}  // namespace walberla
