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

#include <random>
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
   auto devisors = math::getDevisors(n);

   WALBERLA_CHECK( ::is_sorted( primes.begin(), primes.end() ) );
   for(uint_t i = 1; i <= n; ++i)
   {
      WALBERLA_CHECK_EQUAL( walberla::math::isPrime(i), std::binary_search( primes.begin(), primes.end(), i ) );
      WALBERLA_CHECK_EQUAL( n % i == 0, devisors.find( i ) != devisors.end() );
   }

   for( auto it = primes.begin(); it != primes.end(); ++it )
      WALBERLA_CHECK( math::isPrime( *it ) );

   for( auto it = devisors.begin(); it != devisors.end(); ++it )
      WALBERLA_CHECK( n % *it == 0 );

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

   typedef std::mersenne_twister_engine< walberla::uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 0xa37d3c92 > mt11213b;
   mt11213b rng;
   std::uniform_int_distribution<uint_t> dist( 100, 10000 );
   for(int i = 0; i < 100; ++i)
   {
      const uint_t n = dist(rng);
      runTests( n );
   }
}


