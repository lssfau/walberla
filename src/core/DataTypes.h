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
#include <string>
#include <type_traits>

#ifdef __GLIBCXX__ 
#define HAVE_CXXABI_H
#include <cxxabi.h>
#else
#ifdef __has_include
#if __has_include(<cxxabi.h>)
#define HAVE_CXXABI_H
#include <cxxabi.h>
#endif
#endif
#endif


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



template<typename S>
inline S string_to_num( std::string & t );
template <> inline float              string_to_num( std::string & t ) { return std::stof(t); }
template <> inline double             string_to_num( std::string & t ) { return std::stod(t); }
template <> inline long double        string_to_num( std::string & t ) { return std::stold(t); }
template <> inline int                string_to_num( std::string & t ) { return std::stoi(t); }
template <> inline long               string_to_num( std::string & t ) { return std::stol(t); }
template <> inline long long          string_to_num( std::string & t ) { return std::stoll(t); }
template <> inline unsigned long      string_to_num( std::string & t ) { return std::stoul(t); }
template <> inline unsigned long long string_to_num( std::string & t ) { return std::stoull(t); }



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

template< uint_t size > struct uintFromBitWidth;
template<> struct uintFromBitWidth<  8 > { typedef uint8_t  type; };
template<> struct uintFromBitWidth< 16 > { typedef uint16_t type; };
template<> struct uintFromBitWidth< 32 > { typedef uint32_t type; };
template<> struct uintFromBitWidth< 64 > { typedef uint64_t type; };

constexpr uint_t leastUnsignedIntegerBitWidth( uint_t width )
{
   if ( width <=  8 ) return  8;
   if ( width <= 16 ) return 16;
   if ( width <= 32 ) return 32;
   if ( width <= 64 ) return 64;
   return width;
}

/// \brief Provides the smallest unsigned integer type that has at least minSize bits.
///
/// Example:
///
///   leastUnsignedInteger< 5 >::type a; // a is an 8-bit unsigned integer
///   leastUnsignedInteger< 9 >::type b; // b is a 16-bit unsigned integer
///
template< uint_t minSize >
struct leastUnsignedInteger
{
   typedef typename uintFromBitWidth< leastUnsignedIntegerBitWidth( minSize ) >::type type;
};

/// \cond internal
static const uint_t UINT_BITS  = static_cast< uint_t >( std::numeric_limits< uint_t >::digits );
static const uint_t UINT_BYTES = static_cast< uint_t >( std::numeric_limits< uint_t >::digits ) >> 3;

static_assert( !(UINT_BITS & (UINT_BITS - 1)), "Type \"uint_t\" must consist of 2^x Bits!" ); // power of two

template< int N >
struct int_ld
{
   static_assert( N >= 1 && !(N & (N - 1)), "Calculating log_2(N) -> \"N\" must be a power of two!" );
   static const uint_t exp = 1 + int_ld< (N >> 1) >::exp;
};

template< int N > const uint_t int_ld<N>::exp;

template<>
struct int_ld<1>
{
   static const uint_t exp = 0;
};

static const uint_t UINT_BITS_LD = int_ld< std::numeric_limits< uint_t >::digits >::exp;
/// \endcond



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



// data type to string conversion

template< typename T > inline const char* typeToString();

#define TypeToString(X) template<> inline const char* typeToString< X >() { \
   static char string[] = #X; \
   return string; \
}

TypeToString(bool)
TypeToString(char)
TypeToString(short)
TypeToString(int)
TypeToString(long)
TypeToString(long long)
TypeToString(unsigned char)
TypeToString(unsigned short)
TypeToString(unsigned int)
TypeToString(unsigned long)
TypeToString(unsigned long long)
TypeToString(float)
TypeToString(double)

#undef TypeToString

template< typename T > inline const char* typeToString( T ) { return typeToString<T>(); }

// type info demangling

inline std::string demangle( const std::string & name )
{
#ifdef HAVE_CXXABI_H
   int status = 0;
   std::size_t size = 0;
   const char * demangled = abi::__cxa_demangle( name.c_str(), NULL, &size, &status );
   if( demangled == nullptr )
   {
      return name;
   }
   std::string demangled_str(demangled);
   std::free( const_cast<char*>(demangled) );
   return demangled_str;
#else
   return name;
#endif
}



} // namespace walberla

#define WALBERLA_UNUSED(x)  (void)(x)

