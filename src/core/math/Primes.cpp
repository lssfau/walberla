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
//! \file Primes.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Implements functions dealing with prime numbers.
//
//======================================================================================================================

#include "Primes.h"
#include "core/debug/Debug.h"

#include <algorithm>
#include <cmath>


namespace walberla {
namespace math {

/*******************************************************************************************************************//**
 * \brief   Query if n is prime.
 *
 * \param   n  The number to be checked.
 *
 * \return  true if n is prime, false if not.
 **********************************************************************************************************************/
bool isPrime( const uint_t n )
{
   switch(n)
   {
   case 0:
   case 1:
      return false;
   case 2:
      return true;
   default:
      if( n % 2 == 0 )
         return false;
      uint_t sqrtN = uint_c( std::sqrt( real_c(n) ) );
      for( uint_t i = 3; i <= sqrtN; i += 2 )
         if( n % i == 0 )
            return false;
      break;
   }

   return true;
}

/*******************************************************************************************************************//**
 * \brief   Determine prime numbers up to an upper bound.
 *
 * Uses the 'Sieve of Eratosthenes' algorithm.
 * See http://en.wikipedia.org/w/index.php?title=Sieve_of_Eratosthenes&oldid=527676398.
 *
 * \param   n  The maximum prime number.
 *
 * \return  All prime numbers <= n in ascending order.
 **********************************************************************************************************************/
std::vector<uint_t> getPrimes( const uint_t n )
{
   if( n < 2 )
      return std::vector<uint_t>();

   std::vector<bool> markers( n+1, false );
   std::vector<uint_t> primes;

   size_t p = 2;
   while( p <= n / 2 )
   {
      primes.push_back(p);
      for( size_t i = p + p; i <= n; i += p )
         markers[i] = true;

      for( p = p + 1; p <= n / 2 && markers[p]; ++p )
         ; // empty body
   }

  for( uint_t i = n / 2 + 1; i <= n; ++i )
     if( !markers[i] )
        primes.push_back(i);

   return primes;
}

/*******************************************************************************************************************//**
 * \brief   Gets all prime factors of a number.
 *
 * Uses trial division algorithm.
 * See http://en.wikipedia.org/w/index.php?title=Trial_division&oldid=518625973.
 *
 * \param   n  The number to be factorized.
 *
 * \pre     n > 0
 *
 * \return  The prime factors in ascending order.
 **********************************************************************************************************************/
std::vector<uint_t> getPrimeFactors( const uint_t n )
{
   WALBERLA_ASSERT( n != 0 );

   auto primes = getPrimes(n);
   std::vector<uint_t> primeFactors;

   uint_t n_rest = n;
   for(auto primeIt = primes.begin(); primeIt != primes.end(); ++primeIt)
   {
      if( *primeIt * *primeIt > n )
         break;
      while( n_rest % *primeIt == 0)
      {
         n_rest /= *primeIt;
         primeFactors.push_back(*primeIt);
      }
   }

   if( n_rest != 1 )
      primeFactors.push_back(n_rest);

   return primeFactors;
}


/*******************************************************************************************************************//**
 * \brief   Gets all devisors of a number n.
 *
 * Computes all unique products of n's prime factors.
 *
 * \param   n  The number the devisors are computed for.
 *
 * \pre     n > 0
 *
 * \return  A set of the devisors including 1 and n
 **********************************************************************************************************************/
std::set<uint_t> getDevisors( const uint_t n )
{
   if( n == uint_t(0) )
      return std::set<uint_t>();

   std::vector<uint_t> factors = getPrimeFactors( n );

   std::set<uint_t> devisors;
   std::vector<uint_t> tmpDevisors;

   devisors.insert( uint_t(1) );
   tmpDevisors.reserve( ( size_t(1) << factors.size() ) - size_t(1) );

   for( auto fIt = factors.begin(); fIt != factors.end(); ++fIt )
   {
      tmpDevisors.clear();
      for(auto pIt = devisors.begin(); pIt != devisors.end(); ++pIt)
      {
         tmpDevisors.push_back( *pIt * *fIt );
      }
      devisors.insert( tmpDevisors.begin(), tmpDevisors.end() );
   }

   return devisors;
}

} // namespace math
} // namespace walberla
