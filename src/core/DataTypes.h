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
//! \file DataTypes.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>

namespace walberla {


#define WALBERLA_STATIC_ASSERT(x) static_assert(x, "Assertion failed")


template <typename> struct never_true : std::false_type {};

template< typename T > bool isIdentical( const T a, const T b );


// shared ptr

using std::shared_ptr;
using std::weak_ptr;
using std::make_shared;
using std::dynamic_pointer_cast;


// numeric cast (performs range checks in debug mode)

template< typename S, typename T >
inline S numeric_cast( T t ) {
#ifndef NDEBUG
   if( std::is_integral<S>::value && std::is_integral<T>::value && !std::is_same<S,T>::value )
        // integer to different integer: check that forward and back conversion does not change value
   {
      if( !isIdentical( static_cast<T>( static_cast<S>(t) ), t ) )
      {
         throw std::range_error("out of range");
      }
   }
   else if( !std::is_integral<S>::value && !std::is_integral<T>::value && sizeof(S) < sizeof(T) )
       // float to shorter float: check that value within limits of shorter type
   {
      using H = typename std::conditional< !std::is_integral<S>::value && !std::is_integral<T>::value && (sizeof(S) < sizeof(T)), T, long double >::type; // always true, but makes Intel's overflow check happy
      H h = static_cast<H>(t);
      if( h < static_cast<H>(std::numeric_limits<S>::lowest()) || h > static_cast<H>(std::numeric_limits<S>::max()) ) {
         throw std::range_error("out of range");
      }
   }
   else if( std::is_integral<S>::value && !std::is_integral<T>::value )
       // float to integer: check that value within limits of integer
   {
      using H = typename std::conditional< std::is_integral<S>::value && !std::is_integral<T>::value, T, long double >::type; // always true, but makes Intel's overflow check happy
      H h = static_cast<H>(t);
      if( h < static_cast<H>(std::numeric_limits<S>::lowest()) || h > static_cast<H>(std::numeric_limits<S>::max()) ) {
         throw std::range_error("out of range");
      }
   }
#endif
   return static_cast< S >(t);
}


// fixed size signed integral types
typedef std::int8_t   int8_t;    ///<  8 bit signed integer
typedef std::int16_t  int16_t;   ///< 16 bit signed integer
typedef std::int32_t  int32_t;   ///< 32 bit signed integer
typedef std::int64_t  int64_t;   ///< 64 bit signed integer

template< typename T > inline int8_t   int8_c( T t ) { return numeric_cast< int8_t  >(t); } ///< cast to type int8_t  using "int8_c(x)"
template< typename T > inline int16_t int16_c( T t ) { return numeric_cast< int16_t >(t); } ///< cast to type int16_t using "int16_c(x)"
template< typename T > inline int32_t int32_c( T t ) { return numeric_cast< int32_t >(t); } ///< cast to type int32_t using "int32_c(x)"
template< typename T > inline int64_t int64_c( T t ) { return numeric_cast< int64_t >(t); } ///< cast to type int64_t using "int64_c(x)"



// fixed size unsigned integral types

typedef std::uint8_t  uint8_t;    ///<  8 bit unsigned integer
typedef std::uint16_t uint16_t;   ///< 16 bit unsigned integer
typedef std::uint32_t uint32_t;   ///< 32 bit unsigned integer
typedef std::uint64_t uint64_t;   ///< 64 bit unsigned integer
typedef uint8_t byte_t;
typedef uint64_t id_t;            //sid datatype for pe

template< typename T > inline uint8_t   uint8_c( T t ) { return numeric_cast< uint8_t  >(t); } ///< cast to type uint8_t  using "uint8_c(x)"
template< typename T > inline uint16_t uint16_c( T t ) { return numeric_cast< uint16_t >(t); } ///< cast to type uint16_t using "uint16_c(x)"
template< typename T > inline uint32_t uint32_c( T t ) { return numeric_cast< uint32_t >(t); } ///< cast to type uint32_t using "uint32_c(x)"
template< typename T > inline uint64_t uint64_c( T t ) { return numeric_cast< uint64_t >(t); } ///< cast to type uint64_t using "uint64_c(x)"



// signed integral type

template< typename T > inline int int_c( T t ) { return numeric_cast< int >(t); } ///< cast to type int using "int_c(x)"

template< typename INT >
inline void static_assert_int_t() {
   static_assert( std::numeric_limits<INT>::is_specialized && std::numeric_limits<INT>::is_integer, "Integer type required/expected!" );
}



// unsigned integral type

typedef size_t uint_t;

static_assert( std::numeric_limits<uint_t>::is_specialized &&
               std::numeric_limits<uint_t>::is_integer &&
              !std::numeric_limits<uint_t>::is_signed, "Type \"uint_t\" must be an unsigned integer!" );

template< typename T > inline uint_t uint_c( T t ) { return numeric_cast< uint_t >(t); } ///< cast to type uint_t using "uint_c(x)"

template< typename UINT >
inline void static_assert_uint_t() {
   static_assert( std::numeric_limits<UINT>::is_specialized &&
                  std::numeric_limits<UINT>::is_integer     &&
                 !std::numeric_limits<UINT>::is_signed, "Unsigned integer type required/expected!" );
}

// data structure specific data types

typedef int cell_idx_t;
//typedef int64_t cell_idx_t;

WALBERLA_STATIC_ASSERT( std::numeric_limits<cell_idx_t>::is_specialized &&
                     std::numeric_limits<cell_idx_t>::is_integer &&
                     std::numeric_limits<cell_idx_t>::is_signed );

template< typename T > inline cell_idx_t cell_idx_c( T t ) { return numeric_cast< cell_idx_t >(t); } ///< cast to type cell_idx_t using "cell_idx_c(x)"



// floating point type

#ifdef WALBERLA_DOUBLE_ACCURACY
typedef double real_t;
#else
typedef float  real_t;
#endif

inline real_t operator"" _r( long double t ) { return static_cast< real_t >(t); }
inline real_t operator"" _r( unsigned long long int t ) { return static_cast< real_t >(t); }
template< typename T > inline real_t real_c  ( T t ) { return numeric_cast< real_t >(t); } ///< cast to type real_t using "real_c(x)"
template< typename T > inline double double_c( T t ) { return numeric_cast< double >(t); } ///< cast to type double
template< typename T > inline float  float_c ( T t ) { return numeric_cast< float > (t); } ///< cast to type float

/// If you want to compare two reals using operator == and you really know what you are doing, you can use the following function:


template <typename T>
inline bool isIdentical( const T a, const T b )
{
#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
   return a == b;
#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif
}

inline bool realIsIdentical( const real_t a, const real_t b )
{
   return isIdentical( a, b );
}

// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
// conclusion: Comparing Floating Point Numbers - Know what you’re doing. There is no silver bullet. You have to choose wisely.

/// \cond internal
namespace real_comparison
{
   template< class T > struct Epsilon;
   template<> struct Epsilon<       float > { static const       float value; };
   template<> struct Epsilon<      double > { static const      double value; };
   template<> struct Epsilon< long double > { static const long double value; };
}
/// \endcond

inline bool realIsEqual( const real_t a, const real_t b, const real_t eps = real_comparison::Epsilon<real_t>::value ) {
   return std::fabs( a - b ) < eps;
}


inline bool floatIsEqual( long double lhs, long double rhs, const long double epsilon = real_comparison::Epsilon<long double>::value )
{
   return std::fabs( lhs - rhs ) < epsilon;
}

inline bool floatIsEqual( double lhs, double rhs, const double epsilon = real_comparison::Epsilon<double>::value )
{
   return std::fabs( lhs - rhs ) < epsilon;
}

inline bool floatIsEqual( float lhs, float rhs, const float epsilon = real_comparison::Epsilon<float>::value )
{
   return std::fabs( lhs - rhs ) < epsilon;
}

} // namespace walberla

#define WALBERLA_UNUSED(x)  (void)(x)

