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
//! \file Uint.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Uint.h"


namespace walberla {
namespace math {

template<> uint_t uintMSBPosition< uint64_t >( uint64_t value ) { // for the documentation see the header file

   uint64_t i;
   uint64_t j;

   i = value >> 32;
   if( i != 0 ) {
      j = value >> 48;
      if( j != 0 ) {
         i = value >> 56;
         return ( i != 0 ) ? (56 + msbLookupTable[i]) : (48 + msbLookupTable[j]);
      }
      j = value >> 40;
      return ( j != 0 ) ? (40 + msbLookupTable[j]) : (32 + msbLookupTable[i]);
   }
   j = value >> 16;
   if( j != 0 ) {
      i = value >> 24;
      return ( i != 0 ) ? (24 + msbLookupTable[i]) : (16 + msbLookupTable[j]);
   }
   i = value >> 8;
   return ( i != 0 ) ? (8 + msbLookupTable[i]) : msbLookupTable[value];
}

} // namespace math
} // namespace walberla
