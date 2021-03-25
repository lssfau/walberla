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
//! \file MultiArrayIOTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/MultiArrayIO.h"
#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

#include <sstream>
#include <iostream>

using namespace walberla;


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   using namespace std;

   string test1 = "[[ 0.2   0.24  0.2   0.17]\n[ 0.24  0.2   0.2   0.2 ]\n\t[ 0.2   0.2   0.2   0.2 ]\n[ 0.17  0.2   0.2   0.2 ]]";
   stringstream ss1 ( test1 );
   //string betweenBrackets;
   //WALBERLA_CHECK( readContentBetweenBrackets( ss1, betweenBrackets  ) );
   //cout << betweenBrackets;

   boost::multi_array<real_t,2> arr1;
   bool res = !( ss1 >> arr1 ).fail();
   WALBERLA_CHECK( res );
   cout << arr1 << endl;


   string test2 = " [1 2 3,4,5,6,7 ]";
   stringstream ss2 ( test2 );

   boost::multi_array<real_t,1> arr2;
   bool res2 = !(ss2 >> arr2).fail() ;
   WALBERLA_CHECK( res2 );
   cout << arr2 << endl;


   return 0;
}