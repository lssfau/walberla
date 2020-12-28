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
//! \file VariantTest.cpp
//! \ingroup core
//! \author Step
//! \brief Test if there are no typos in the collective header files "all.h"
//
//======================================================================================================================

#include "core/Variant.h"
#include "core/debug/all.h"

#include <string>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std::literals; // C++ 14
using var_t = walberla::variant<int, long, double, std::string>;

int main( int /*argc*/, char** /*argv*/ )
{
   walberla::debug::enterTestMode();



   // Example from: https://en.cppreference.com/w/cpp/utility/variant (license (CC-BY-SA) and GPL)
   walberla::variant<int, float> v, w;
   v = 12; // v contains int
   int i = walberla::get<int>( v );
   WALBERLA_CHECK( 0 == 12 - i );
   w = walberla::get<int>( v );
   w = v;
   WALBERLA_CHECK( 0 == 12 - i );

   //  walberla::get<double>(v); // error: no double in [int, float]
   //  walberla::get<3>(v);      // error: valid index values are 0 and 1

   try {
      float f = walberla::get<float>( w ); // w contains int, not float: will throw
      std::cout << f << std::endl;
   } catch ( const walberla::bad_variant_access& ) {}

   walberla::variant<std::string> x( "abc" ); // converting constructors work when unambiguous
   x = "def"; // converting assignment also works when unambiguous

   std::cout << "hallo" << std::endl;
   walberla::variant<std::string, bool> y( true );
   std::cout << "eoo" << std::endl;
   WALBERLA_CHECK( walberla::holds_alternative<bool>( y ) ); // succeeds
   y = "xyz"s;
   WALBERLA_CHECK( walberla::holds_alternative<std::string>( y ) ); //succeeds

   std::cout << "bye" << std::endl;
   std::vector<var_t> vec = {10, 15l, 1.5, "hello"};

   for ( auto& z : vec ) {
      // 1. void visitor, only called for side-effects (here, for I/O)
      walberla::visit( []( auto && arg ) {std::cout << arg;}, z );
   }

   return 0;
}
