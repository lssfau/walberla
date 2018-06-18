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
//! \file PeriodicIntersect.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "PeriodicIntersect.h"
#include "MapPointToPeriodicDomain.h"

namespace walberla {
namespace domain_decomposition {

bool periodicIntersect( const std::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB & box1, const math::AABB & box2 )
{
   auto min1 = box1.minCorner();
   auto min2 = box2.minCorner();
   auto max1 = box1.maxCorner();
   auto max2 = box2.maxCorner();

   mapPointToPeriodicDomain(periodic, domain, min1);
   mapPointToPeriodicDomain(periodic, domain, min2);
   mapPointToPeriodicDomain(periodic, domain, max1);
   mapPointToPeriodicDomain(periodic, domain, max2);

   bool flag = false;

   if (max1[0] > max2[0])
   {
      if ( (max1[0]-max2[0]) < box1.xSize() ) flag = true;
   } else
   {
      if ( (max2[0]-max1[0]) < box2.xSize() ) flag = true;
   }
   if (min1[0] > min2[0])
   {
      if ( (min1[0]-min2[0]) < box2.xSize() ) flag = true;
   } else
   {
      if ( (min2[0]-min1[0]) < box1.xSize() ) flag = true;
   }

   if (!flag) return false;
   flag = false;

   if (max1[1] > max2[1])
   {
      if ( (max1[1]-max2[1]) < box1.ySize() ) flag = true;
   } else
   {
      if ( (max2[1]-max1[1]) < box2.ySize() ) flag = true;
   }
   if (min1[1] > min2[1])
   {
      if ( (min1[1]-min2[1]) < box2.ySize() ) flag = true;
   } else
   {
      if ( (min2[1]-min1[1]) < box1.ySize() ) flag = true;
   }

   if (!flag) return false;
   flag = false;

   if (max1[2] > max2[2])
   {
      if ( (max1[2]-max2[2]) < box1.zSize() ) flag = true;
   } else
   {
      if ( (max2[2]-max1[2]) < box2.zSize() ) flag = true;
   }
   if (min1[2] > min2[2])
   {
      if ( (min1[2]-min2[2]) < box2.zSize() ) flag = true;
   } else
   {
      if ( (min2[2]-min1[2]) < box1.zSize() ) flag = true;
   }

   return flag;
}

bool periodicIntersect( const std::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB& box1, const math::AABB & box2, const real_t dx )
{
   return periodicIntersect(periodic, domain, box1.getExtended( dx ), box2);
}


} // namespace domain_decomposition
} // namespace walberla
