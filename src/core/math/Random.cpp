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

static boost::mt19937 generator; // static boost::random::mt19937_64 generator;

boost::mt19937 & getGenerator() // boost::random::mt19937_64
{
   return generator;
}

}



void seedRandomGenerator( const boost::mt19937::result_type & seed )
{
#ifdef _OPENMP
   #pragma omp critical (random)
#endif
   internal::getGenerator().seed( seed );
}



} // namespace math
} // namespace walberla
