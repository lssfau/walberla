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
//! \file StreamOutput.h
//! \ingroup simd
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "SIMD.h"
#include <type_traits>


namespace walberla {
namespace simd {


// enable output operator only for vector4 types
template<typename VEC4>
typename std::enable_if<is_vector4_type<VEC4>::value, std::ostream&>::type
operator<< ( std::ostream & os, VEC4 v ) {
   os << "[" << getComponent(v,3) <<", "
             << getComponent(v,2) <<", "
             << getComponent(v,1) <<", "
             << getComponent(v,0) <<"]";

   return os;
}



} // namespace simd
} // namespace walberla


