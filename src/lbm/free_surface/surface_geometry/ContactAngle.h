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
//! \file ContactAngle.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Class to avoid re-computing sine and cosine of the contact angle.
//
//======================================================================================================================

#pragma once

#include "core/math/Constants.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Class to avoid re-computing sine and cosine of the contact angle.
 **********************************************************************************************************************/
class ContactAngle
{
 public:
   ContactAngle(real_t angleInDegrees)
      : angleDegrees_(angleInDegrees), angleRadians_(math::pi / real_c(180) * angleInDegrees),
        sinAngle_(std::sin(angleRadians_)), cosAngle_(std::cos(angleRadians_))
   {}

   inline real_t getInDegrees() const { return angleDegrees_; }
   inline real_t getInRadians() const { return angleRadians_; }
   inline real_t getSin() const { return sinAngle_; }
   inline real_t getCos() const { return cosAngle_; }

 private:
   real_t angleDegrees_;
   real_t angleRadians_;
   real_t sinAngle_;
   real_t cosAngle_;
}; // class ContactAngle
} // namespace free_surface
} // namespace walberla