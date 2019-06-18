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
//! \file   Contact.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/DataTypes.h>

namespace walberla {
namespace mesa_pd {

struct Contact
{
   Contact(const size_t idx1,
           const size_t idx2,
           const Vec3   contactPoint,
           const Vec3   contactNormal,
           const real_t penetrationDepth)
      : idx1_(idx1)
      , idx2_(idx2)
      , contactPoint_(contactPoint)
      , contactNormal_(contactNormal)
      , penetrationDepth_(penetrationDepth) {}

   size_t idx1_;
   size_t idx2_;
   Vec3   contactPoint_;
   Vec3   contactNormal_;
   real_t penetrationDepth_;
};

} // namespace mesa_pd
} // namespace walberla
