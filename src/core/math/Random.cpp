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
//! \file Random.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Random.h"


namespace walberla {
namespace math {



namespace internal {

static std::mt19937 generator; // static std::mt19937_64 generator;

std::mt19937 & getGenerator() // std::mt19937_64
{
   return generator;
}

}



void seedRandomGenerator( const std::mt19937::result_type & seed )
{
#ifdef _OPENMP
   #pragma omp critical (Random_random)
#endif
   internal::getGenerator().seed( seed );
}



} // namespace math
} // namespace walberla
