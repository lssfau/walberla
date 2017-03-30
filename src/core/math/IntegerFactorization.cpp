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
//! \file IntegerFactorization.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "IntegerFactorization.h"
#include "Primes.h"
#include "core/debug/Debug.h"


namespace walberla {
namespace math {



//**********************************************************************************************************************
/*!
*   Computes and returns an integer factorization of any arbitrary integer
*
*   Goal: returned factors should roughly be of equal size
*   Typical use case: uniform process distribution
*
*   \param  number           integer to factorize
*   \param  numberOfFactors  number of factors to create (= size of returned vector)
*/
//**********************************************************************************************************************

std::vector< uint_t > getFactors( const uint_t number, const uint_t numberOfFactors )
{
   WALBERLA_ASSERT( numberOfFactors > 0 );

   std::vector< uint_t > factors( numberOfFactors, 1 );

   std::vector< uint_t > primes = getPrimeFactors( number );

   for( auto prime = primes.rbegin(); prime != primes.rend(); ++prime )
   {
      auto smallestCurrentFactor = std::min_element( factors.begin(), factors.end() );
      WALBERLA_ASSERT( smallestCurrentFactor != factors.end(), "vector\"factors\" is empty" );
      *smallestCurrentFactor *= *prime;
   }

   return factors;
}



std::vector< uint_t > getFactors( const uint_t number, const uint_t numberOfFactors, const std::vector< real_t >& weight )
{
   WALBERLA_ASSERT( numberOfFactors > 0 );
   WALBERLA_ASSERT( weight.size() == numberOfFactors );

   std::vector< uint_t > factors( numberOfFactors, 1 );

   std::vector< uint_t > primes = getPrimeFactors( number );

   for( auto prime = primes.rbegin(); prime != primes.rend(); ++prime )
   {
      uint_t smallestWeightedFactor = 0;
      for( uint_t factor = 1; factor < numberOfFactors; ++factor )
      {
         if( ( real_c(factors[factor]) / weight[factor] ) < ( real_c(factors[smallestWeightedFactor]) / weight[smallestWeightedFactor] ) )
            smallestWeightedFactor = factor;
      }
      factors[smallestWeightedFactor] *= *prime;
   }

   return factors;
}

Vector3<uint_t> getFactors3D( const uint_t number )
{
   std::vector<uint_t> result = getFactors( number, 3u );

   return Vector3<uint_t>( result[0], result[1], result[2] );
}

Vector3<uint_t> getFactors3D( const uint_t number, const Vector3< real_t >& weights )
{
   std::vector<real_t> weights2( 3u );
   weights2[0] = weights[0];
   weights2[1] = weights[1];
   weights2[2] = weights[2];

   std::vector<uint_t> result = getFactors( number, 3u, weights2 );
   return Vector3<uint_t>( result[0], result[1], result[2] );
}



} // namespace math
} // namespace walberla
