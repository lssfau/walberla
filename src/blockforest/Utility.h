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
//! \file Utility.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Types.h"

#include "core/EndianIndependentSerialization.h"
#include "core/math/KahanSummation.h"
#include "core/math/Uint.h"

#include <limits>
#include <ostream>
#include <string>
#include <vector>


namespace walberla {
namespace blockforest {


using math::uintPow2;
using math::uintMSBPosition;


//**********************************************************************************************************************
/*!
*   \brief Returns a string that stores the bitwise representation of 'value' (must be an unsigned integer)
*
*   \code{.unparsed}
*      8bit display: 0101_1101
*     16bit display: 1110_0101.1100_0001
*     32bit display: 1000_0011.0110_1101.0000_0001.1010_0110
*     64bit analogously ...
*   \endcode
*/
//**********************************************************************************************************************
template< typename UINT > std::string uintToBitString( const UINT value ) {

   static_assert_uint_t< UINT >();

   std::string bitsString;

   UINT exp = static_cast< UINT >( std::numeric_limits< UINT >::digits ) - 1;

   static const uint_t bytes = uint_c( std::numeric_limits< UINT >::digits ) >> 3;

   for( uint_t i = 0; i != bytes; ++i ) {

      for( uint_t j = 0; j != 4; ++j, --exp )
         bitsString.push_back( (value & uintPow2(exp)) ? '1' : '0' );

      bitsString.push_back('_');

      for( uint_t j = 0; j != 4; ++j, --exp )
         bitsString.push_back( (value & uintPow2(exp)) ? '1' : '0' );

      if( i != bytes-1 ) bitsString.append(".");
   }

   return bitsString;
}



template< typename T >
workload_t workloadSum( const T& array )
{
   math::KahanAccumulator< workload_t > acc;

   for( uint_t i = 0; i < array.size(); ++i )
      acc += array[i]->getWorkload();

   return acc.get();
}

template< typename T >
memory_t memorySum( const T& array )
{
   math::KahanAccumulator< memory_t > acc;

   for( uint_t i = 0; i < array.size(); ++i )
      acc += array[i]->getMemory();

   return acc.get();
}



inline memory_t bytesToMiB( memory_t bytes ) {
   return bytes / static_cast< memory_t >( 1024.0 * 1024.0 );
}



//**********************************************************************************************************************
/*!
*   \brief Returns the string representation of 'number', every three digits the character 'separator' is inserted
*          (172408725 -> "172 408 725")
*/
//**********************************************************************************************************************
std::string naturalNumberToGroupedThousandsString( const uint_t number, const char separator = ' ' );

inline std::string naturalNumberToGroupedThousandsString( const real_t number, const char separator = ' ' ) {

   return naturalNumberToGroupedThousandsString( static_cast< uint_t >( 0.5 + number ), separator );
}



inline void fillStream( std::ostream& ostream, const char fill, uint_t length ) {
   while( length-- > 0 ) ostream << fill;
}



} // namespace blockforest
} // namespace walberla
