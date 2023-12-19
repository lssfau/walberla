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
//! \file DataTypes.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "DataTypes.h"


namespace walberla {

namespace real_comparison
{
   #ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
//   const    bfloat16 Epsilon<    bfloat16 >::value = static_cast<    bfloat16 >(1e-2); // machine eps is 2^-7
   const     float16 Epsilon<     float16 >::value = static_cast<     float16 >(1e-3); // machine eps is 2^-10
   // Note, depending on the kind of float16 <bfloat, float16> another Epsilon must be used.
   #endif
   const       float Epsilon<       float >::value = static_cast<       float >(1e-4);
   const      double Epsilon<      double >::value = static_cast<      double >(1e-8);
   const long double Epsilon< long double >::value = static_cast< long double >(1e-10);
}

} // namespace walberla
