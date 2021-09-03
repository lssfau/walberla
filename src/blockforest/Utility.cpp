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
//! \file Utility.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Utility.h"

#include <iomanip>
#include <sstream>
#include <stack>


namespace walberla {
namespace blockforest {



std::string naturalNumberToGroupedThousandsString( uint_t number, const char separator ) { // for the documentation see the header file

   std::ostringstream oss;

   std::stack< uint_t > numbers;

   while( true ) {

      numbers.push( number % 1000 );

      if( number <= 999 ) break;

      number /= 1000;
   }

   oss << numbers.top();
   numbers.pop();

   while( !numbers.empty() ) {

      oss << separator << std::setfill( '0' ) << std::setw( 3 ) << numbers.top();
      oss << std::setw( 0 );
      numbers.pop();
   }

   return oss.str();
}



} // namespace blockforest
} // namespace walberla
