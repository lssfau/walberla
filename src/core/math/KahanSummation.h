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
//! \file KahanSummation.h
//! \ingroup math
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include <iterator>


namespace walberla {
namespace math {



template< typename T >
class KahanAccumulator
{
public:

   KahanAccumulator() : sum( T(0) ), c( T(0) ) {}

   void reset() { sum = T(0); c = T(0); }
   
   void operator+=( const T & number )
   {
      const T y = number - c;
      const T t = sum + y;

      c   = ( t - sum ) - y;
      sum = t;   
   }
   
   T get() const { return sum; }
   
private:

   T sum;
   T c;
};



template< typename Iterator >
typename std::iterator_traits< Iterator >::value_type
kahanSummation( const Iterator & begin, const Iterator & end )
{
   KahanAccumulator< typename std::iterator_traits< Iterator >::value_type > acc;

   for( auto element = begin; element != end; ++element )
      acc += *element;

   return acc.get();
}



} // namespace math
} // namespace walberla
