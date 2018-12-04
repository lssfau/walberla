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
//! \file SetSelectableObjectTest.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/selectable/SetSelectableObject.h"

#include <string>


using namespace walberla;



inline static void add( selectable::SetSelectableObject< std::string, size_t >& container, const std::string& function,
                        const Set<size_t>& required, const Set<size_t>& incompatible = Set<size_t>::emptySet() ) {

   container.add( function, required, incompatible, function );
}



int main( int /*argc*/, char** /*argv*/ ) {

   debug::enterTestMode();

   using A = Set<size_t>;

   selectable::SetSelectableObject< std::string, size_t > container;

   add( container, "function_1", A(1)+A(2)+A(3) );
   add( container, "function_2", A(5)+A(6)+A(7) );
   add( container, "function_3", A(1)+A(2) );
   add( container, "function_4", A(1)+A(2)+A(3)+A(4) );
   add( container, "function_5", A(1)+A(2)+A(3), A(4) );
   add( container, "function_6", A::emptySet(), A(1)+A(2) );
   add( container, "function_7", A::emptySet(), A(4) );

   std::string function;

   std::vector< std::string > functions;
   std::vector< std::string > expected;

   WALBERLA_CHECK_EQUAL( container.get( function, A(1)+A(2)+A(3) ), static_cast< size_t >(2) );

   functions.clear();
   container.get( functions, A(1)+A(2)+A(3) );

   expected.clear();
   expected.emplace_back("function_1");
   expected.emplace_back("function_5");

   WALBERLA_CHECK_EQUAL( functions, expected );

   WALBERLA_CHECK_EQUAL( container.get( function, A(1)+A(2)+A(3)+A(4) ), static_cast< size_t >(1) );

   expected.clear();
   expected.emplace_back("function_4");
   WALBERLA_CHECK_EQUAL( function, expected[0] );

   WALBERLA_CHECK_EQUAL( container.get( function, A(1)+A(2) ), static_cast< size_t >(1) );

   expected.clear();
   expected.emplace_back("function_3");
   WALBERLA_CHECK_EQUAL( function, expected[0] );

   WALBERLA_CHECK_EQUAL( container.get( function, A(1)+A(2)+A(3)+A(5)+A(6)+A(7) ), static_cast< size_t >(3) );

   functions.clear();
   container.get( functions, A(1)+A(2)+A(3)+A(5)+A(6)+A(7) );

   expected.clear();
   expected.emplace_back("function_1");
   expected.emplace_back("function_2");
   expected.emplace_back("function_5");

   WALBERLA_CHECK_EQUAL( functions, expected );

   WALBERLA_CHECK_EQUAL( container.get( function, A(1)+A(2)+A(3)+A(4)+A(5)+A(6)+A(7) ), static_cast< size_t >(1) );

   expected.clear();
   expected.emplace_back("function_4");
   WALBERLA_CHECK_EQUAL( function, expected[0] );

   WALBERLA_CHECK_EQUAL( container.get( function, A(3)+A(5)+A(6) ), static_cast< size_t >(2) );

   functions.clear();
   container.get( functions, A(3)+A(5)+A(6) );

   expected.clear();
   expected.emplace_back("function_6");
   expected.emplace_back("function_7");

   WALBERLA_CHECK_EQUAL( functions, expected );

   WALBERLA_CHECK_EQUAL( container.get( function, A(1)+A(5)+A(6) ), static_cast< size_t >(1) );

   expected.clear();
   expected.emplace_back("function_7");
   WALBERLA_CHECK_EQUAL( function, expected[0] );

   WALBERLA_CHECK_EQUAL( container.get( function, A(4)+A(5)+A(6) ), static_cast< size_t >(1) );

   expected.clear();
   expected.emplace_back("function_6");
   WALBERLA_CHECK_EQUAL( function, expected[0] );

   WALBERLA_CHECK_EQUAL( container.get( function, A(7)+A(5)+A(6) ), static_cast< size_t >(1) );

   expected.clear();
   expected.emplace_back("function_2");
   WALBERLA_CHECK_EQUAL( function, expected[0] );

   functions.clear();
   for( auto it = container.begin(); it != container.end(); ++it )
      functions.push_back( *it );

   expected.clear();
   expected.emplace_back("function_1");
   expected.emplace_back("function_2");
   expected.emplace_back("function_3");
   expected.emplace_back("function_4");
   expected.emplace_back("function_5");
   expected.emplace_back("function_6");
   expected.emplace_back("function_7");

   WALBERLA_CHECK_EQUAL( functions, expected );

   functions.clear();
   for( auto it = container.begin(); it != container.end(); ++it )
      functions.push_back( it.identifier() );

   WALBERLA_CHECK_EQUAL( functions, expected );

   return EXIT_SUCCESS;
}
