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
//! \file RefinementScaling.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

namespace walberla::lbm_generated
{

inline real_t relaxationRateScaling( real_t relaxationRate, uint_t refinementLevel )
{
   const real_t levelScaleFactor = real_c(uint_c(1) << refinementLevel);
   const real_t one                = real_c(1.0);
   const real_t half               = real_c(0.5);

   return real_c(relaxationRate / (levelScaleFactor * (-relaxationRate * half + one) + relaxationRate * half));
}

} // namespace walberla::lbm_generated