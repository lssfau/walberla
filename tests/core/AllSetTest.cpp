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
//! \file AllSetTest.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "core/AllSet.h"
#include "core/debug/TestSubsystem.h"

#include <cstdlib>
#include <iostream>


using namespace walberla;



int normalAndNormal() {

   AllSet<int> A = AllSet<int>(1) + AllSet<int>(2) + AllSet<int>(3);
   A.insert(4);
   AllSet<int> B = AllSet<int>(5) + AllSet<int>(6) + AllSet<int>(7);
   AllSet<int> C = AllSet<int>(4) + AllSet<int>(5) + AllSet<int>(6);
   AllSet<int> D = AllSet<int>(1) + AllSet<int>(3) + AllSet<int>(5);
   AllSet<int> E = AllSet<int>(4) + AllSet<int>(3) + AllSet<int>(2);

   WALBERLA_CHECK( !A.contains(5) && A.contains(1) && A.contains(3) && A.contains(4) );

   WALBERLA_CHECK( !A.contains(B) && !A.contains(D) && A.contains(E) && A.contains(A) && D.contains(D) );

   WALBERLA_CHECK( !A.intersects(B) && A.intersects(C) && A.intersects(D) && A.intersects(A) );

   WALBERLA_CHECK( !A.isEmpty() && A.isCountable() && !A.isAll() && !A.isUniverse() );

   WALBERLA_CHECK( !A.equalSize(B) && B.equalSize(C) );

   WALBERLA_CHECK( A > B && A >= B && B < A && B <= A );

   WALBERLA_CHECK( B <= C && C <= B && B >= C && C >= B );

   E.clear();
   WALBERLA_CHECK( E.isEmpty() && E.isCountable() && !E.isAll() && !E.isUniverse() );

   WALBERLA_CHECK( AllSet<int>::emptySet().isEmpty() && AllSet<int>::emptySet().isCountable() &&
         !AllSet<int>::emptySet().isAll() && !AllSet<int>::emptySet().isUniverse() );

   WALBERLA_CHECK( !AllSet<int>::all().isEmpty() && !AllSet<int>::all().isCountable() &&
           AllSet<int>::all().isAll() && AllSet<int>::all().isUniverse() );

   WALBERLA_CHECK( A.contains( AllSet<int>::emptySet() ) && AllSet<int>::emptySet().contains( AllSet<int>::emptySet() ) );
   WALBERLA_CHECK( AllSet<int>::all().contains( AllSet<int>::emptySet() ) );

   WALBERLA_CHECK( !A.equalSize( AllSet<int>::all() ) && !A.equalSize( AllSet<int>::emptySet() ) );

   WALBERLA_CHECK( AllSet<int>::all() > A && AllSet<int>::all() >= A && A < AllSet<int>::all() && A <= AllSet<int>::all() );

   WALBERLA_CHECK( A > AllSet<int>::emptySet() && A >= AllSet<int>::emptySet() );

   WALBERLA_CHECK( AllSet<int>::all() > AllSet<int>::emptySet() && AllSet<int>::all() >= AllSet<int>::emptySet() &&
         AllSet<int>::emptySet() < AllSet<int>::all() && AllSet<int>::emptySet() <= AllSet<int>::all() );

   E = AllSet<int>(4) + AllSet<int>(3) + AllSet<int>(2) + AllSet<int>(1);

   WALBERLA_CHECK( A == E && C != D );

   // A and A
   AllSet<int>   result = A + A;
   AllSet<int> expected = AllSet<int>(1) + AllSet<int>(2) + AllSet<int>(3) + AllSet<int>(4);
   WALBERLA_CHECK( result == expected );

   result   = A & A;
   expected = AllSet<int>(1) + AllSet<int>(2) + AllSet<int>(3) + AllSet<int>(4);
   WALBERLA_CHECK( result == expected );

   result   = A - A;
   expected = AllSet<int>::emptySet();
   WALBERLA_CHECK( result == expected );

   // A and B
   result    = A + B;
   expected = AllSet<int>(1) + AllSet<int>(2) + AllSet<int>(3) + AllSet<int>(4) + AllSet<int>(5) + AllSet<int>(6) + AllSet<int>(7);
   WALBERLA_CHECK( result == expected );

   result   = A & B;
   expected = AllSet<int>::emptySet();
   WALBERLA_CHECK( result == expected );

   result   = A - B;
   expected = A;
   WALBERLA_CHECK( result == expected );

   // A and C
   result   = A + C;
   expected = AllSet<int>(1) + AllSet<int>(2) + AllSet<int>(3) + AllSet<int>(4) + AllSet<int>(5) + AllSet<int>(6);
   WALBERLA_CHECK( result == expected );

   result   = A & C;
   expected = AllSet<int>(4);
   WALBERLA_CHECK( result == expected );

   result   = A - C;
   expected = AllSet<int>(1) + AllSet<int>(2) + AllSet<int>(3);
   WALBERLA_CHECK( result == expected );

   // A and D
   result   = A + D;
   expected = AllSet<int>(1) + AllSet<int>(2) + AllSet<int>(3) + AllSet<int>(4) + AllSet<int>(5);
   WALBERLA_CHECK( result == expected );

   result   = A & D;
   expected = AllSet<int>(1) + AllSet<int>(3);
   WALBERLA_CHECK( result == expected );

   result   = A - D;
   expected = AllSet<int>(2) + AllSet<int>(4);
   WALBERLA_CHECK( result == expected );

   return EXIT_SUCCESS;
}



int normalAndAll() {
/*
   AllSet<int> A = AllSet<int>(1) + AllSet<int>(2) + AllSet<int>(3) + AllSet<int>(4);

   // A and all()
   AllSet<int> result   = A + AllSet<int>::all();
   AllSet<int> expected = AllSet<int>::all();
   WALBERLA_CHECK( result == expected );

   result   = A & AllSet<int>::all();
   expected = A;
   WALBERLA_CHECK( result == expected );

   result   = A - AllSet<int>::all();
   expected = AllSet<int>::emptySet();
   WALBERLA_CHECK( result == expected );

   // A and B (all)
   AllSet<int> B = AllSet<int>::all() - AllSet<int>(5);

   result   = A + B;
   expected = AllSet<int>::all() - AllSet<int>(5);
   WALBERLA_CHECK( result == expected );

   result   = A & B;
   expected = A;
   WALBERLA_CHECK( result == expected );

   result   = A - B;
   expected = AllSet<int>::emptySet();
   WALBERLA_CHECK( result == expected );

   WALBERLA_CHECK( B.contains(A) && B.intersects(A) && B > A && B >= A && !B.equalSize(A) && B != A )

   B = AllSet<int>::all() - AllSet<int>(2);

   result   = A + B;
   expected = AllSet<int>::all();
   WALBERLA_CHECK( result == expected );

   result   = A & B;
   expected = AllSet<int>(1) + AllSet<int>(3) + AllSet<int>(4);
   WALBERLA_CHECK( result == expected );

   result   = A - B;
   expected = AllSet<int>(2);
   WALBERLA_CHECK( result == expected );

   WALBERLA_CHECK( !B.contains(A) && B.intersects(A) && B > A && B >= A && !B.equalSize(A) && B != A )

   B = AllSet<int>::all() - AllSet<int>(2) - AllSet<int>(3);

   result   = A + B;
   expected = AllSet<int>::all();
   WALBERLA_CHECK( result == expected );

   result   = A & B;
   expected = AllSet<int>(1) + AllSet<int>(4);
   WALBERLA_CHECK( result == expected );

   result   = A - B;
   expected = AllSet<int>(2) + AllSet<int>(3);
   WALBERLA_CHECK( result == expected );

   WALBERLA_CHECK( !B.contains(A) && B.intersects(A) && B > A && B >= A && !B.equalSize(A) && B != A )

   B = AllSet<int>::all() - A;

   result   = A + B;
   expected = AllSet<int>::all();
   WALBERLA_CHECK( result == expected );

   result   = A & B;
   expected = AllSet<int>::emptySet();
   WALBERLA_CHECK( result == expected );

   result   = A - B;
   expected = A;
   WALBERLA_CHECK( result == expected );

   WALBERLA_CHECK( !B.contains(A) && !B.intersects(A) && B > A && B >= A && !B.equalSize(A) && B != A )
*/
   return EXIT_SUCCESS;
}



int allAndNormal() {

   // TODO

   return EXIT_SUCCESS;
}



int allAndAll() {

   // TODO

   return EXIT_SUCCESS;
}



int main( int /*argc*/, char** /*argv*/ ) {

   debug::enterTestMode();

   if( normalAndNormal() == EXIT_FAILURE )
      return EXIT_FAILURE;

   if( normalAndAll() == EXIT_FAILURE )
      return EXIT_FAILURE;

   if( allAndNormal() == EXIT_FAILURE )
      return EXIT_FAILURE;

   if( allAndAll() == EXIT_FAILURE )
      return EXIT_FAILURE;

   return EXIT_SUCCESS;
}
