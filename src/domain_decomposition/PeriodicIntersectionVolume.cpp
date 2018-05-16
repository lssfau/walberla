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
//! \file PeriodicIntersectionVolume.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/logging/Logging.h"
#include "MapPointToPeriodicDomain.h"

namespace walberla {
namespace domain_decomposition {

real_t periodicIntersectionVolume( const std::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB & box1, const math::AABB & box2 )
{
   const auto diagonal = (domain.maxCorner() - domain.minCorner());
   const auto halfDiagonal = real_t(0.5) * diagonal;
   const auto center1 = box1.center();
   auto center2 = box2.center();

   for (size_t dim = 0; dim < 3; ++dim)
   {
      if (periodic[dim])
      {
         while ((center2[dim]-center1[dim])>halfDiagonal[dim]) center2[dim] -= diagonal[dim];
         while ((center2[dim]-center1[dim])<-halfDiagonal[dim]) center2[dim] += diagonal[dim];
      }
   }

   return box1.intersectionVolume(box2.getTranslated(center2 - box2.center()));
}

real_t periodicIntersectionVolume( const std::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB& box1, const math::AABB & box2, const real_t dx )
{
   return periodicIntersectionVolume(periodic, domain, box1.getExtended( dx ), box2);
}


} // namespace domain_decomposition
} // namespace walberla
