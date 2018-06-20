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
//! \file Ray.h
//! \author Lukas Werner
//
//======================================================================================================================

#include <pe/raytracing/Ray.h>

namespace walberla {
namespace pe {
namespace raytracing {

/*!\brief Global output operator for rays.
 *
 * \param os Reference to the output stream.
 * \param ray Reference to a constant ray object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const Ray& ray) {
   return os << "<o: " << ray.getOrigin()
   << ", d: " << ray.getDirection()
   << ", c: (" << ray.getImageX() << "/" << ray.getImageY() << ")>";
}

} //namespace raytracing
} //namespace pe
} //namespace walberla
