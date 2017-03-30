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
//! \file SetTest.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "core/Set.h"
#include "core/debug/TestSubsystem.h"

#include <cstdlib>
#include <iostream>


using namespace walberla;



int main( int /*argc*/, char** /*argv*/ ) {

   debug::enterTestMode();

   Set<int> A = Set<int>(1) + Set<int>(2) + Set<int>(3);
   A.insert(4);
   Set<int> B = Set<int>(5) + Set<int>(6) + Set<int>(7);
   Set<int> C = Set<int>(4) + Set<int>(5) + Set<int>(6);
   Set<int> D = Set<int>(1) + Set<int>(3) + Set<int>(5);
   Set<int> E = Set<int>(4) + Set<int>(3) + Set<int>(2);

   WALBERLA_CHECK( !A.contains(5) && A.contains(1) && A.contains(3) && A.contains(4) );

   WALBERLA_CHECK( !A.contains(B) && !A.contains(D) && A.contains(E) && A.contains(A) && D.contains(D) );

   WALBERLA_CHECK( !A.intersects(B) && A.intersects(C) && A.intersects(D) && A.intersects(A) );

   WALBERLA_CHECK( !A.isEmpty() );

   WALBERLA_CHECK( !A.equalSize(B) && B.equalSize(C) );

   WALBERLA_CHECK_GREATER( A, B );       // A >  B
   WALBERLA_CHECK_GREATER_EQUAL( A, B ); // A >= B
   WALBERLA_CHECK_LESS( B, A );          // B <  A
   WALBERLA_CHECK_LESS_EQUAL( B, A );    // B <= A

   WALBERLA_CHECK_LESS_EQUAL( B, C );    // B <= C
   WALBERLA_CHECK_LESS_EQUAL( C, B );    // C <= B
   WALBERLA_CHECK_GREATER_EQUAL( B, C ); // B >= C
   WALBERLA_CHECK_GREATER_EQUAL( C, B ); // C >= B

   E.clear();
   WALBERLA_CHECK( E.isEmpty() );

   WALBERLA_CHECK( Set<int>::emptySet().isEmpty() );

   WALBERLA_CHECK( A.contains( Set<int>::emptySet() ) && Set<int>::emptySet().contains( Set<int>::emptySet() ) );

   WALBERLA_CHECK( !A.equalSize( Set<int>::emptySet() ) );

   WALBERLA_CHECK( A > Set<int>::emptySet() && A >= Set<int>::emptySet() );

   E = Set<int>(4) + Set<int>(3) + Set<int>(2) + Set<int>(1);

   WALBERLA_CHECK_EQUAL( A, E );   // A == E
   WALBERLA_CHECK_UNEQUAL( C, D ); // C != D

   // A and A
   Set<int>   result = A + A;
   Set<int> expected = Set<int>(1) + Set<int>(2) + Set<int>(3) + Set<int>(4);
   WALBERLA_CHECK_EQUAL( result, expected );

   result = A;
   result.insert( A.begin(), A.end() );
   WALBERLA_CHECK_EQUAL( result, expected );

   result   = A & A;
   expected = Set<int>(1) + Set<int>(2) + Set<int>(3) + Set<int>(4);
   WALBERLA_CHECK_EQUAL( result, expected );

   result   = A - A;
   expected = Set<int>::emptySet();
   WALBERLA_CHECK_EQUAL( result, expected );

   // A and B
   result    = A + B;
   expected = Set<int>(1) + Set<int>(2) + Set<int>(3) + Set<int>(4) + Set<int>(5) + Set<int>(6) + Set<int>(7);
   WALBERLA_CHECK_EQUAL( result, expected );

   result = A;
   result.insert( B.begin(), B.end() );
   WALBERLA_CHECK_EQUAL( result, expected );

   result   = A & B;
   expected = Set<int>::emptySet();
   WALBERLA_CHECK_EQUAL( result, expected );

   result   = A - B;
   expected = A;
   WALBERLA_CHECK_EQUAL( result, expected );

   // A and C
   result   = A + C;
   expected = Set<int>(1) + Set<int>(2) + Set<int>(3) + Set<int>(4) + Set<int>(5) + Set<int>(6);
   WALBERLA_CHECK_EQUAL( result, expected );

   result = A;
   result.insert( C.begin(), C.end() );
   WALBERLA_CHECK_EQUAL( result, expected );

   result   = A & C;
   expected = Set<int>(4);
   WALBERLA_CHECK_EQUAL( result, expected );

   result   = A - C;
   expected = Set<int>(1) + Set<int>(2) + Set<int>(3);
   WALBERLA_CHECK_EQUAL( result, expected );

   // A and D
   result   = A + D;
   expected = Set<int>(1) + Set<int>(2) + Set<int>(3) + Set<int>(4) + Set<int>(5);
   WALBERLA_CHECK_EQUAL( result, expected );

   result = A;
   result.insert( D.begin(), D.end() );
   WALBERLA_CHECK_EQUAL( result, expected );

   result   = A & D;
   expected = Set<int>(1) + Set<int>(3);
   WALBERLA_CHECK_EQUAL( result, expected );

   result   = A - D;
   expected = Set<int>(2) + Set<int>(4);
   WALBERLA_CHECK_EQUAL( result, expected );

   return EXIT_SUCCESS;
}
