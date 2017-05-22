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
//! \file Primes.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Declares functions dealing with prime numbers.
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include <set>
#include <vector>


namespace walberla {
namespace math {

bool                isPrime        ( const uint_t n );
std::vector<uint_t> getPrimes      ( const uint_t n );
std::vector<uint_t> getPrimeFactors( const uint_t n );
std::set<uint_t>    getDevisors    ( const uint_t n );

} // namespace math
} // namespace walberla
