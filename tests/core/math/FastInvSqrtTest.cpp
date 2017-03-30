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
//! \file FastInvSqrtTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/math/FastInvSqrt.h"
#include "core/logging/Logging.h"

#include <cmath>

using namespace walberla;


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   WALBERLA_LOG_RESULT( "Testing double" );
   for ( double testNr = 1e-12; testNr < 1e12; testNr = testNr * 2.1 )
   {
      double error = math::fastInvSqrt<3>( testNr ) - ( 1.0 / std::sqrt( testNr ) );
      WALBERLA_CHECK( floatIsEqual(error, 0.0, 1e-4 ) );
   }

   for ( double testNr = 1e-3; testNr < 1e12; testNr = testNr * 2.1 )
   {
      double error = math::fastInvSqrt<3>( testNr ) - ( 1.0 / std::sqrt( testNr ) );
      WALBERLA_CHECK( floatIsEqual(error, 0.0, 1e-9 ) );
   }

   for ( double testNr = 1e-3; testNr < 1e12; testNr = testNr * 2.1 )
   {
      double error = math::fastInvSqrt<2>( testNr ) - ( 1.0 / std::sqrt( testNr ) );
      WALBERLA_CHECK( floatIsEqual(error, 0.0, 1e-4 ) );
   }

   WALBERLA_LOG_RESULT( "Testing float" );
   for ( float testNr = 1e-4f; testNr < 1e12f; testNr = testNr * 2.1f )
   {
      double error = math::fastInvSqrt<2>( testNr ) - ( 1.0 / std::sqrt( testNr ) );
      WALBERLA_CHECK( floatIsEqual(error, 0.0, 1e-2 ) );
   }


   return 0;
}
