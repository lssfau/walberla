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
//! \file Uint.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include <type_traits>


namespace walberla {
namespace math {



//**********************************************************************************************************************
/*!
*   \brief Returns true only if 'value' is an even number
*/
//**********************************************************************************************************************
template< typename UINT >
inline bool uintIsEven( const UINT value ) {
   static_assert_uint_t< UINT >();
   return (value & static_cast< UINT >(1)) == static_cast< UINT >(0);
}



//**********************************************************************************************************************
/*!
*   \brief Returns true only if 'value' (must be an unsigned integer) is a power of 2
*          (http://graphics.stanford.edu/~seander/bithacks.html)
*/
//**********************************************************************************************************************
template< typename UINT >
inline bool uintIsPowerOfTwo( const UINT value ) {
   static_assert_uint_t< UINT >();
   return value && !(value & (value - static_cast< UINT >(1)));
}



//**********************************************************************************************************************
/*!
*   \brief Returns the result of 2^exp
*/
//**********************************************************************************************************************
template< typename UINT >
inline UINT uintPow2( const UINT exp ) {

   static_assert_uint_t< UINT >();
   WALBERLA_ASSERT_LESS( exp, static_cast< UINT >( std::numeric_limits< UINT >::digits ) );

   return static_cast< UINT >( static_cast< UINT >(1) << exp );
}



//**********************************************************************************************************************
/*!
*   \brief Returns the result of 4^exp
*/
//**********************************************************************************************************************
template< typename UINT >
inline UINT uintPow4( UINT exp ) {

   exp *= static_cast< UINT >(2);

   static_assert_uint_t< UINT >();
   WALBERLA_ASSERT_LESS( exp, static_cast< UINT >( std::numeric_limits< UINT >::digits ) );

   return static_cast< UINT >( static_cast< UINT >(1) << exp );
}



//**********************************************************************************************************************
/*!
*   \brief Returns the result of 8^exp
*/
//**********************************************************************************************************************
template< typename UINT >
inline UINT uintPow8( UINT exp ) {

   exp *= static_cast< UINT >(3);

   static_assert_uint_t< UINT >();
   WALBERLA_ASSERT_LESS( exp, static_cast< UINT >( std::numeric_limits< UINT >::digits ) );

   return static_cast< UINT >( static_cast< UINT >(1) << exp );
}



//**********************************************************************************************************************
/*!
*   \brief Returns the result of 16^exp
*/
//**********************************************************************************************************************
template< typename UINT >
inline UINT uintPow16( UINT exp ) {

   exp *= static_cast< UINT >(4);

   static_assert_uint_t< UINT >();
   WALBERLA_ASSERT_LESS( exp, static_cast< UINT >( std::numeric_limits< UINT >::digits ) );

   return static_cast< UINT >( static_cast< UINT >(1) << exp );
}



// http://graphics.stanford.edu/~seander/bithacks.html

static const uint8_t msbLookupTable[256] =
{
#define msbLT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
      0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
      msbLT(5), msbLT(6), msbLT(6), msbLT(7), msbLT(7), msbLT(7), msbLT(7),
      msbLT(8), msbLT(8), msbLT(8), msbLT(8), msbLT(8), msbLT(8), msbLT(8), msbLT(8)
#undef msbLT
};

//**********************************************************************************************************************
/*!
*   \brief Calculation of the position of the most significant bit of a variable of type UINT
*
*   Examples: "[...] 1001" -> 4, "[...] 0011" -> 2, "[...] 0110" -> 3, "[...] 0001" -> 1, "[...] 0000" -> 0, etc.
*/
//**********************************************************************************************************************
template< typename UINT > uint_t uintMSBPosition( UINT value )
{
   static_assert( std::is_unsigned< UINT >::value, "uintMSBPosition can only be used with unsigned integer types!" );

   switch( std::numeric_limits<UINT>::digits )
   {
      case  8: return uintMSBPosition(  uint8_c( value ) );
      case 16: return uintMSBPosition( uint16_c( value ) );
      case 32: return uintMSBPosition( uint32_c( value ) );
      case 64: return uintMSBPosition( uint64_c( value ) );
   }

   WALBERLA_ABORT( "Unsupported unsigned int data type for function \"math::uintMSBPosition\"!\n"
                   "Data type uses " << std::numeric_limits<UINT>::digits << " digits, waLBerla only supports 8, 16, 32, and 64." );
}

template<> uint_t uintMSBPosition< uint64_t >( uint64_t value ); // -> Uint.cpp

template<> inline uint_t uintMSBPosition< uint32_t >( uint32_t value ) {

   uint32_t i, j;

   j = value >> 16;
   if( j != 0 ) {
      i = value >> 24;
      return ( i != 0 ) ? (24 + msbLookupTable[i]) : (16 + msbLookupTable[j]);
   }
   i = value >> 8;
   return ( i != 0 ) ? (8 + msbLookupTable[i]) : msbLookupTable[value];
}

template<> inline uint_t uintMSBPosition< uint16_t >( uint16_t value ) {

// uint16_t i = value >> 8; // strange (?) g++ error: conversion to "uint16_t {aka short unsigned int}"
                            //                        from "int" may alter its value [-Werror=conversion]
   uint32_t i = uint32_c(value) >> 8;

   return ( i != 0 ) ? (8 + msbLookupTable[i]) : msbLookupTable[value];
}

template<> inline uint_t uintMSBPosition< uint8_t >( uint8_t value ) {

   return msbLookupTable[value];
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
static constexpr uint_t UINT_BITS  = static_cast< uint_t >( std::numeric_limits< uint_t >::digits );
static constexpr uint_t UINT_BYTES = static_cast< uint_t >( std::numeric_limits< uint_t >::digits ) >> 3;

static_assert( !(UINT_BITS & (UINT_BITS - 1)), "Type \"uint_t\" must consist of 2^x Bits!" ); // power of two

template< int N >
struct int_ld
{
   static_assert( N >= 1 && !(N & (N - 1)), "Calculating log_2(N) -> \"N\" must be a power of two!" );
   static constexpr uint_t exp = 1 + int_ld< (N >> 1) >::exp;
   static_assert( exp > 0 );
};

template<>
struct int_ld<1>
{
   static constexpr uint_t exp = 0;
};

static constexpr uint_t UINT_BITS_LD = int_ld< std::numeric_limits< uint_t >::digits >::exp;
/// \endcond

} // namespace math
} // namespace walberla
