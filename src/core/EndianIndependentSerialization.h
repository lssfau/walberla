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
//! \file EndianIndependentSerialization.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "DataTypes.h"
#include "core/debug/Debug.h"

#include <cmath>
#include <limits>
#include <vector>



namespace walberla {



//**********************************************************************************************************************
/*!
*   \brief Stores the 'bytes' least significant bytes of 'value' in 'array' starting at 'array[offset]'
*          (complementary operation/function: byteArrayToUint)
*
*   \code
*     uint_t number = uint_c(199994); // '[...] 0000 0011 0000 1101 0011 1010'
*     std::vector< uint8_t > byteStream( 10, uint8_c(0) );
*     uintToByteArray( number, byteStream, 1, 2 );
*     for( int i = 0; i != 10; ++i )
*        std::cout << byteStream[i] << " ";
*     // prints: "0 58 13 0 0 0 0 0 0 0" (58 <-> 0011 1010 -- 13 <-> 0000 1101)
*   \endcode
*/
//**********************************************************************************************************************
inline void uintToByteArray( uint_t value, std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes )
{
   WALBERLA_ASSERT_LESS_EQUAL( offset + bytes, array.size() );
   WALBERLA_ASSERT_LESS_EQUAL( bytes, UINT_BYTES );

   for( uint_t i = 0; i != bytes; ++i ) {

      array[ offset + i ] = uint8_c( value & uint_c(255) );
      value >>= 8;
   }
}

//**********************************************************************************************************************
/*!
*   \brief Converts 'bytes' bytes stored in 'array' starting at 'array[offset]' into a value of type uint_t
*          (complementary operation/function: uintToByteArray)
*
*   \code
*     std::vector< uint8_t > byteStream( 10, uint8_c(0) );
*     byteStream[0] = uint8_c(114); // '0111 0010'
*     byteStream[1] = uint8_c( 85); // '0101 0101'
*     byteStream[2] = uint8_c(213); // '1101 0101'
*     uint_t value = byteArrayToUint( byteStream, 1, 2 );
*     std::cout << value << std::endl; // prints: "54613" ('[...] 0000 1101 0101 0101 0101')
*   \endcode
*/
//**********************************************************************************************************************
inline uint_t byteArrayToUint( const std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes )
{
   WALBERLA_ASSERT_LESS_EQUAL( offset + bytes, array.size() );
   WALBERLA_ASSERT_LESS_EQUAL( bytes, UINT_BYTES );

   uint_t value = 0;
   for( uint_t i = 0; i != bytes; ++i )
      value |= uint_c( array[ offset + i ] ) << ( i * 8 );

   return value;
}



template< typename REAL_T >
uint_t realToByteArray( const REAL_T value, std::vector< uint8_t >& array, const uint_t offset )
{
   static_assert( sizeof( uint64_t ) >= sizeof( REAL_T ), "type uint64_t must have at least as many bits as type REAL_T" );

   const uint_t size = uint_c( 1 + 2 + sizeof( REAL_T ) );
   WALBERLA_ASSERT_LESS_EQUAL( offset + size, array.size() );

   int exp;
   const REAL_T x = std::frexp( value, &exp );
   const REAL_T sign = ( x < REAL_T(0) ) ? REAL_T(-1) : REAL_T(1);

   uint_t signByte = ( ( exp >= 0 ) ? uint_c(0) : uint_c(1) ) + ( ( sign < REAL_T(0) ) ? uint_c(2) : uint_c(0) );
   array[ offset ] = uint8_c( signByte );

   uint32_t uexp = ( exp >= 0 ) ? uint32_c( exp ) : uint32_c( -1 * exp );
   for( uint_t i = 0; i != 2; ++i )
   {
      array[ offset + 1 + i ] = uint8_c( uexp & uint32_c(255) );
      uexp >>= 8;
   }

   uint64_t mant = uint64_c( real_c( uint64_c(1) << uint64_c( 1 + std::numeric_limits<REAL_T>::digits ) ) * sign * x );
   for( uint_t i = 0; i != sizeof( REAL_T ); ++i )
   {
      array[ offset + 3 + i ] = uint8_c( mant & uint64_c(255) );
      mant >>= 8;
   }

   return size;
}

template< typename REAL_T >
REAL_T byteArrayToReal( const std::vector< uint8_t >& array, const uint_t offset )
{
   static_assert( sizeof( uint64_t ) >= sizeof( REAL_T ), "type uint64_t must have at least as many bits as type REAL_T" );

   WALBERLA_ASSERT_LESS_EQUAL( offset + uint_c( 1 + 2 + sizeof( REAL_T ) ), array.size() );

   uint32_t uexp = 0;
   for( uint_t i = 0; i != 2; ++i )
      uexp |= uint32_c( array[ offset + 1 + i ] ) << ( i * 8 );

   const int exp = ( ( uint_c( array[ offset ] ) & uint_c(1) ) == uint_c(0) ) ? int_c( uexp ) : ( -1 * int_c( uexp ) );

   uint64_t mant = 0;
   for( uint_t i = 0; i != sizeof( REAL_T ); ++i )
      mant |= uint64_c( array[ offset + 3 + i ] ) << ( i * 8 );

   const REAL_T sign = ( ( uint_c( array[ offset ] ) & uint_c(2) ) == uint_c(0) ) ? REAL_T(1) : REAL_T(-1);
   return std::ldexp( sign * numeric_cast<REAL_T>( mant ) / numeric_cast<REAL_T>( uint64_c(1) << uint64_c( 1 + std::numeric_limits<REAL_T>::digits ) ), exp );
}


inline void boolVectorToByteArray( const std::vector< bool >& boolVec, std::vector< uint8_t >& array, const uint_t offset )
{
   static const uint_t bit[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

   WALBERLA_ASSERT_EQUAL( boolVec.size() % 8, uint_t(0) );
   const uint_t bytes =  uint_c(boolVec.size() / 8);

   WALBERLA_ASSERT_LESS_EQUAL( offset + bytes, array.size() );

   for( uint_t i = 0; i != bytes; ++i ) {
      WALBERLA_ASSERT_LESS_EQUAL(offset + i, array.size());
      array[ offset + i ] = uint8_c( boolVec[i*8+0]*bit[0] + boolVec[i*8+1]*bit[1] + boolVec[i*8+2]*bit[2] + boolVec[i*8+3]*bit[3] +
                                     boolVec[i*8+4]*bit[4] + boolVec[i*8+5]*bit[5] + boolVec[i*8+6]*bit[6] + boolVec[i*8+7]*bit[7] );
   }
}

inline std::vector< bool > byteArrayToBoolVector( const std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes )
{
   static const uint8_t bit[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

   WALBERLA_ASSERT_LESS_EQUAL( offset + bytes, array.size() );

   std::vector< bool > boolVec( 8 * bytes );

   for( uint_t i = 0; i != bytes; ++i )
      for( uint_t j = 0; j != 8; ++j )
         boolVec[i * 8 + j] = (array[offset + i] & bit[j]) != uint8_t(0);

   return boolVec;
}



} // namespace walberla
