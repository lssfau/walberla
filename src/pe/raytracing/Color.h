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
//! \file   Color.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/Types.h>
#include "core/math/Vector3.h"

namespace walberla {
namespace pe {
namespace raytracing {

class Color: public Vector3<real_t> {
public:
   /*!\name Constructors */
   //@{
   /*!\brief Instantiation constructor for the Color class. Defaults to white.
    */
   Color () : Color(1, 1, 1) {

   }
   
   /*!\brief Instantiation constructor for the Color class.
   * \param r Red component
   * \param g Green component
   * \param b Blue component
   * Instantiation constructor for the Color class with RGB components. Each value should be between 0 and 1 (soft limits)
   */
   Color (real_t r, real_t g, real_t b) : Vector3<real_t>(r, g, b) {
      
   }
   
   /*!\brief Instantiation constructor for the Color class.
    * \param vector RGB vector
    * Instantiation constructor for the Color class with RGB components. Each value should be between 0 and 1 (soft limits)
    */
   Color (const Vec3& vector) : Color(vector[0], vector[1], vector[2]) {
      
   }
   //@}
   
   /*!\brief Multiply this color with another component wise.
    * \return Color with components of this and other multiplied.
    */
   inline Color mulComponentWise(const Color& other) const {
      return Color((*this)[0]*other[0],
                   (*this)[1]*other[1],
                   (*this)[2]*other[2]);
   }
   
   /*!\brief Clamps this colors component values between 0 and 1.
    */
   inline void clamp() {
      (*this)[0] = std::min(std::max((*this)[0], real_t(0)), real_t(1));
      (*this)[1] = std::min(std::max((*this)[1], real_t(0)), real_t(1));
      (*this)[2] = std::min(std::max((*this)[2], real_t(0)), real_t(1));
   }
   
   static Color colorFromHSV(real_t hue, real_t saturation, real_t value);
};

} //namespace raytracing
} //namespace pe
} //namespace walberla
