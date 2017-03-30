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
//! \file PrimesTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Implements test for prime number functions.
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Primes.h"
#include "core/mpi/MPIManager.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <algorithm>
#include <numeric>


using namespace walberla;

template <typename InputIterator>
bool is_sorted( InputIterator begin, InputIterator end )
{
   return std::adjacent_find( begin, end, std::greater<uint_t>() ) == end;
}

void runTests( const uint_t n )
{
   auto primes = math::getPrimes(n);
   WALBERLA_CHECK( ::is_sorted( primes.begin(), primes.end() ) );
   for(uint_t i = 0; i <= n; ++i)
      WALBERLA_CHECK_EQUAL( walberla::math::isPrime(i), std::binary_search( primes.begin(), primes.end(), i ) );

   if( n != 0 )
   {
      auto primeFactors = math::getPrimeFactors(n);
      WALBERLA_CHECK( ::is_sorted( primeFactors.begin(), primeFactors.end() ) );
      WALBERLA_CHECK_EQUAL( n, std::accumulate( primeFactors.begin(), primeFactors.end(), uint_c(1), std::multiplies<uint_t>() ) );
      for( auto it = primeFactors.begin(); it != primeFactors.end(); ++it )
         WALBERLA_CHECK( math::isPrime(*it) );
   }
}

int main(int argc, char * argv[])
{
   debug::enterTestMode();

   MPIManager::instance()->initializeMPI(&argc, &argv);

   // Verify corner cases
   WALBERLA_CHECK( !math::isPrime( uint_c(0) ) );
   WALBERLA_CHECK( !math::isPrime( uint_c(1) ) );
   WALBERLA_CHECK(  math::isPrime( uint_c(2) ) );

   for(uint_t n = 0; n < 100; ++n)
   {
      runTests( n );
   }

   boost::mt11213b rng;
   boost::uniform_int<uint_t> dist( 100, 10000 );
   for(int i = 0; i < 100; ++i)
   {
      const uint_t n = dist(rng);
      runTests( n );
   }
}


