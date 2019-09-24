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
//! \file   Color.cpp
//! \author Lukas Werner
//
//======================================================================================================================

#include "Color.h"

namespace walberla {
namespace pe {
namespace raytracing {

/*!\brief Create a Color object from HSV values.
 * \param hue Hue value in degrees from 0-360
 * \param saturation Saturation value from 0-1
 * \param value Value from 0-1
 */
Color Color::colorFromHSV(real_t hue, real_t saturation, real_t value) {
   // based on Max K. Agoston: Computer Graphics and Geometric Modeling - Implementation and Algorithms
   real_t r;
   real_t g;
   real_t b;
   
   if (realIsEqual(hue, real_t(360))) {
      hue = real_t(0);
   } else {
      hue /= real_t(60);
   }
   real_t fract = hue - std::floor(hue);
   
   real_t P = value*(real_t(1) - saturation);
   real_t Q = value*(real_t(1) - saturation*fract);
   real_t T = value*(real_t(1) - saturation*(real_t(1) - fract));
   
   if (real_t(0) <= hue && hue < real_t(1)) {
      r = value;
      g = T;
      b = P;
   } else if (real_t(1) <= hue && hue < real_t(2)) {
      r = Q;
      g = value,
      b = P;
   } else if (real_t(2) <= hue && hue < real_t(3)) {
      r = P;
      g = value;
      b = T;
   } else if (real_t(3) <= hue && hue < real_t(4)) {
      r = P;
      g = Q;
      b = value;
   } else if (real_t(4) <= hue && hue < real_t(5)) {
      r = T;
      g = P;
      b = value;
   } else if (real_t(5) <= hue && hue < real_t(6)) {
      r = value;
      g = P;
      b = Q;
   } else {
      r = g = b = real_t(0);
   }
   
   return Color(r, g, b);
}

}
}
}
