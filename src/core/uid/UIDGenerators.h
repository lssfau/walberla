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
//! \file UIDGenerators.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Feichtinger
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"

#include <cstdint>

#include <limits>


namespace walberla {
namespace uid {



//**********************************************************************************************************************
/*
*   Concept for UID generators
*   --------------------------
*
*   Every UID generator class must use two template parameters ( -> "template< typename T, typename UINT >" )
*
*   Every UID generator must define 2 types:
*
*      - 'uint_type' (must be equal to the second template parameter)
*      - 'generator_type' (must be unique for each UID generator class)
*
*   Every UID generator must implement 4 static member functions:
*
*      - uint_type generateUID()
*      - uint_type toIndex( const uint_type uid )
*      - uint_type toBitMask( const uint_type uid )
*      - const char* getType()
*
*   These 4 functions must implement the following behavior:
*
*      - generateUID(): every time this function is called it must return a new, unique identifier of type uint_type
*      - toIndex(): this function must return '0' for the first UID created by generateUID(), '1' for the second, '2'
*                   for the third, '3' for the fourth, etc.
*      - toBitMask(): this function must return '[...] 0001' for the first UID created by generateUID(), '[...] 0010'
*                     for the second, '[...] 0100' for the third, etc.
*      - getType(): for each UID generator, this function should return a characteristic/unique string
*/
//**********************************************************************************************************************

// currently available UID generators:

template< typename T, typename UINT > class IndexGenerator;
template< typename T, typename UINT > class StandardGenerator;
template< typename T, typename UINT > class BitGenerator;
template< typename T, typename UINT > class SingletonGenerator;





struct index_generated_tag {};

//**********************************************************************************************************************
/*!
*   \brief Generates UIDs (unique IDs) per class T starting from 0. Increment is 1.
*/
//**********************************************************************************************************************

template< typename T, typename UINT >
class IndexGenerator {
public:

   typedef UINT uint_type;
   typedef index_generated_tag generator_type;

   static uint_type generateUID() {
      static uint_type uid(0);
      return uid++;
   }

   static uint_type firstUID() { return 0; }

   static uint_type nextUID( const uint_type uid ) { return uid+1; }

   static uint_type toIndex( const uint_type uid ) { return uid; }

   static uint_type toBitMask( const uint_type uid ) { return static_cast< uint_type >(1) << uid; }

   static const char* getType() { static const char* const type = "index generator"; return type; }

   WALBERLA_STATIC_ASSERT( std::numeric_limits<UINT>::is_specialized && std::numeric_limits<UINT>::is_integer &&
                       !std::numeric_limits<UINT>::is_signed );
};



struct standard_generated_tag {};

//**********************************************************************************************************************
/*!
*   \brief Generates UIDs (unique IDs) per class T starting from 1. Increment is 1.
*/
//**********************************************************************************************************************

template< typename T, typename UINT >
class StandardGenerator {
public:

   typedef UINT uint_type;
   typedef standard_generated_tag generator_type;

   static uint_type generateUID() {
      static uint_type uid(0);
      return ++uid;
   }

   static uint_type firstUID() { return 1; }

   static uint_type nextUID( const uint_type uid ) { return uid+1; }

   static uint_type toIndex( const uint_type uid ) { return uid - 1; }

   static uint_type toBitMask( const uint_type uid ) { return static_cast< uint_type >(1) << ( uid - 1 ); }

   static const char* getType() { static const char* const type = "standard generator"; return type; }

   WALBERLA_STATIC_ASSERT( std::numeric_limits<UINT>::is_specialized && std::numeric_limits<UINT>::is_integer &&
                       !std::numeric_limits<UINT>::is_signed );
};



///////////////////////////////////////
// LOG BASE 2 OF AN UNSIGNED INTEGER //
///////////////////////////////////////

// http://graphics.stanford.edu/~seander/bithacks.html

static const uint8_t logBase2LookupTable[256] =
{
#define logBase2LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
      0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
      logBase2LT(4), logBase2LT(5), logBase2LT(5), logBase2LT(6), logBase2LT(6), logBase2LT(6), logBase2LT(6),
      logBase2LT(7), logBase2LT(7), logBase2LT(7), logBase2LT(7), logBase2LT(7), logBase2LT(7), logBase2LT(7), logBase2LT(7)
#undef logBase2LT
};

template< typename UINT > size_t logBase2( UINT value ); // 1000 -> 3, 0010 -> 1, 0001 -> 0, etc.

template<> inline size_t logBase2< uint64_t >( uint64_t value ) {

   uint64_t i, j;

   i = value >> 32;
   if( i != 0 ) {
      j = value >> 48;
      if( j != 0 ) {
         i = value >> 56;
         return ( i != 0 ) ? (56 + logBase2LookupTable[i]) : (48 + logBase2LookupTable[j]);
      }
      j = value >> 40;
      return ( j != 0 ) ? (40 + logBase2LookupTable[j]) : (32 + logBase2LookupTable[i]);
   }
   j = value >> 16;
   if( j != 0 ) {
      i = value >> 24;
      return ( i != 0 ) ? (24 + logBase2LookupTable[i]) : (16 + logBase2LookupTable[j]);
   }
   i = value >> 8;
   return ( i != 0 ) ? (8 + logBase2LookupTable[i]) : logBase2LookupTable[value];
}

template<> inline size_t logBase2< uint32_t >( uint32_t value ) {

   uint32_t i, j;

   j = value >> 16;
   if( j != 0 ) {
      i = value >> 24;
      return ( i != 0 ) ? (24 + logBase2LookupTable[i]) : (16 + logBase2LookupTable[j]);
   }
   i = value >> 8;
   return ( i != 0 ) ? (8 + logBase2LookupTable[i]) : logBase2LookupTable[value];
}

struct bit_generated_tag {};

//**********************************************************************************************************************
/*!
*   \brief Generates UIDs (unique IDs) per class T starting from 1. This generator uses a bitwise increment (0001 ->
*          0010 -> 0100 -> 1000).
*
*   In debug mode, generating more UIDs than there are bits available for an unsigned int of type "uint_type" will
*   trigger an assertion and fail.
*/
//**********************************************************************************************************************

template< typename T, typename UINT >
class BitGenerator {
public:

   typedef UINT uint_type;
   typedef bit_generated_tag generator_type;

   static uint_type generateUID() {
      static uint_type uid(1);
      //WALBERLA_ASSERT_LESS( uid, ( static_cast< uint_type >(1) << ( static_cast< uint_type >( std::numeric_limits< uint_type >::digits - 1 ) ) ) );
      WALBERLA_ASSERT_UNEQUAL( uid, 0 );
      uint_type ret = uid;
      uid = static_cast< uint_type >(uid << 1);
      return ret;
   }

   static uint_type firstUID() { return 1; }

   static uint_type nextUID( const uint_type uid ) {
      WALBERLA_ASSERT( uid && !(uid & (uid - 1)) ); // assert( uid is power of two )
      WALBERLA_ASSERT_LESS( uid, ( static_cast< uint_type >(1) << ( static_cast< uint_type >( std::numeric_limits< uint_type >::digits - 1 ) ) ) );
      return uid << 1;
   }

   static uint_type toIndex( const uint_type uid ) { return static_cast< uint_type >( logBase2( uid ) ); }

   static uint_type toBitMask( const uint_type uid ) { return uid; }

   static const char* getType() { static const char* const type = "bit generator"; return type; }

   WALBERLA_STATIC_ASSERT( std::numeric_limits<UINT>::is_specialized && std::numeric_limits<UINT>::is_integer &&
                       !std::numeric_limits<UINT>::is_signed );
};



struct singleton_generated_tag {};

//**********************************************************************************************************************
/*!
*   \brief Generates UIDs (unique IDs) per class T starting from 1. This generator only allows one UID to be created.
*
*   In debug mode, generating more than one UID for class T will trigger an assertion and fail.
*/
//**********************************************************************************************************************

template< typename T, typename UINT >
class SingletonGenerator {
public:

   typedef UINT uint_type;
   typedef singleton_generated_tag generator_type;

   static uint_type generateUID() {
      static uint_type uid(0);
      WALBERLA_ASSERT_EQUAL( uid, 0 );
      return ++uid;
   }

   static uint_type firstUID() { return 1; }

   static uint_type nextUID( const uint_type uid ) { WALBERLA_ASSERT( false ); return 1; }

   static uint_type toIndex( const uint_type uid ) { WALBERLA_ASSERT_EQUAL( uid, 1 ); return 0; }

   static uint_type toBitMask( const uint_type uid ) { WALBERLA_ASSERT_EQUAL( uid, 1 ); return 1; }

   static const char* getType() { static const char* const type = "singleton generator"; return type; }

   WALBERLA_STATIC_ASSERT( std::numeric_limits<UINT>::is_specialized && std::numeric_limits<UINT>::is_integer &&
                       !std::numeric_limits<UINT>::is_signed );
};



} // namespace uid
} // namespace walberla
