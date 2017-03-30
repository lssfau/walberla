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
//! \file IntegerFactorization.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include <vector>


namespace walberla {
namespace math {

std::vector< uint_t > getFactors( const uint_t number, const uint_t numberOfFactors );
std::vector< uint_t > getFactors( const uint_t number, const uint_t numberOfFactors, const std::vector< real_t >& weights );

Vector3<uint_t> getFactors3D( const uint_t number );
Vector3<uint_t> getFactors3D( const uint_t number, const Vector3< real_t >& weights );

} // namespace math
} // namespace walberla
