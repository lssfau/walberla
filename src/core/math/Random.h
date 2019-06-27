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
//! \file Random.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include <random>
#include <limits>


namespace walberla {
namespace math {



namespace internal {
std::mt19937 & getGenerator(); // std::mt19937_64
}



void seedRandomGenerator( const std::mt19937::result_type & seed ); // std::mt19937_64



//**********************************************************************************************************************
/*!
*   \brief Returns a random integral number of type INT in the range [min,max] (max included!)
*/
//**********************************************************************************************************************
template< typename INT >
INT intRandom( const INT min, const INT max, std::mt19937 & generator )
{
   static_assert_int_t< INT >();
   static_assert(sizeof(INT) > sizeof(char), "cannot use char");
   WALBERLA_ASSERT_LESS_EQUAL( min, max );

   std::uniform_int_distribution< INT > distribution( min, max );

   INT value;
#ifdef _OPENMP
   #pragma omp critical (Random_random)
#endif
   { value = distribution( generator ); }

   return value;
}

template<>
inline char intRandom<char>( const char min, const char max, std::mt19937 & generator )
{
   static_assert_int_t< char >();
   WALBERLA_ASSERT_LESS_EQUAL( min, max );

   std::uniform_int_distribution< int16_t > distribution( min, max );

   char value;
#ifdef _OPENMP
   #pragma omp critical (Random_random)
#endif
   { value = static_cast<char>( distribution( generator ) ); }

   return value;
}

template<>
inline unsigned char intRandom<unsigned char>( const unsigned char min, const unsigned char max, std::mt19937 & generator )
{
   static_assert_int_t< unsigned char >();
   WALBERLA_ASSERT_LESS_EQUAL( min, max );

   std::uniform_int_distribution< int16_t > distribution( min, max );

   unsigned char value;
#ifdef _OPENMP
   #pragma omp critical (Random_random)
#endif
   { value = static_cast<unsigned char>( distribution( generator ) ); }

   return value;
}

template<>
inline signed char intRandom<signed char>( const signed char min, const signed char max, std::mt19937 & generator )
{
   static_assert_int_t< signed char >();
   WALBERLA_ASSERT_LESS_EQUAL( min, max );

   std::uniform_int_distribution< int16_t > distribution( min, max );

   signed char value;
#ifdef _OPENMP
   #pragma omp critical (Random_random)
#endif
   { value = static_cast<signed char>( distribution( generator ) ); }

   return value;
}


template< typename INT >
INT intRandom( const INT min )
{
   return intRandom( min, std::numeric_limits<INT>::max(), internal::getGenerator() );
}

template< typename INT >
INT intRandom( const INT min, const INT max )
{
   return intRandom( min, max, internal::getGenerator() );
}

template< typename INT >
INT intRandom()
{
   return intRandom( std::numeric_limits<INT>::min(), std::numeric_limits<INT>::max() );
}



template< typename INT >
class IntRandom
{
public:
   IntRandom( const std::mt19937::result_type & seed = std::mt19937::result_type() ) { generator_.seed( seed ); }
   INT operator()( const INT min = std::numeric_limits<INT>::min(), const INT max = std::numeric_limits<INT>::max() )
   {
      return intRandom( min, max, generator_ );
   }
private:
   std::mt19937 generator_;
};



//**********************************************************************************************************************
/*!
*   \brief Returns a random floating point number of type REAL in the range [min,max) (max excluded!)
*/
//**********************************************************************************************************************
template< typename REAL = real_t >
REAL realRandom( const REAL min = REAL(0), const REAL max = REAL(1), std::mt19937 & generator = internal::getGenerator() )
{
   static_assert( std::numeric_limits<REAL>::is_specialized && !std::numeric_limits<REAL>::is_integer, "Floating point type required/expected!" );
   WALBERLA_ASSERT_LESS( min, max );

   std::uniform_real_distribution< REAL > distribution( min, max );

   REAL value;
#ifdef _OPENMP
   #pragma omp critical (Random_random)
#endif
   { value = distribution( generator ); }

   return value;
}



template< typename REAL >
class RealRandom
{
public:
   RealRandom( const std::mt19937::result_type & seed = std::mt19937::result_type() ) { generator_.seed( seed ); }
   REAL operator()( const REAL min = REAL(0), const REAL max = REAL(1) )
   {
      return realRandom( min, max, generator_ );
   }
private:
   std::mt19937 generator_;
};



inline bool boolRandom() { ///< Randomly returns 'true' or 'false'
   return realRandom<real_t>() >= real_t(0.5);
}



class BoolRandom
{
public:
   BoolRandom( const std::mt19937::result_type & seed = std::mt19937::result_type() ) { generator_.seed( seed ); }
   bool operator()() { return boolRandom(); }
private:
   std::mt19937 generator_;
};



} // namespace math
} // namespace walberla
