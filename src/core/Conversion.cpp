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
//! \file Conversion.cpp
//! \ingroup core
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "Conversion.h"

namespace walberla {

void convert(const std::vector<bool> & from, std::vector<uint8_t> & to)
{
   to.clear();
   to.reserve( from.size() / 8 + 1);

   auto it = from.begin();
   while( it != from.end() )
   {
      uint8_t container = 0;

      for( size_t bitPos = 0; bitPos < 8 && it != from.end(); ++bitPos, ++it )
         container = numeric_cast<uint8_t>( container | numeric_cast<uint8_t>( ( static_cast<bool>( *it ) ? ( 1U << bitPos ) : 0U ) ) );

      to.push_back(container);
   }
}

void convert(const std::vector<uint8_t> & from, std::vector<bool> & to)
{
   to.clear();
   to.resize(from.size() * sizeof(uint8_t) * 8);

   auto unpackIt = from.begin();
   auto it = to.begin();
   while( it != to.end() )
   {
      uint8_t container;
      container = *unpackIt;
      ++unpackIt;

      for( size_t bitPos = 0; bitPos < 8 && it != to.end(); ++bitPos, ++it )
         *it = ( container & ( static_cast<uint8_t>(1) << bitPos ) ) > 0U;
   }
}

}
