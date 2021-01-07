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
//! \file 
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/math/Angles.h"
#include "core/debug/TestSubsystem.h"

namespace walberla{

int main( int /*argc*/, char** /*argv*/ )
{
   using namespace walberla::math;

   debug::enterTestMode();

   WALBERLA_CHECK_FLOAT_EQUAL( radToDeg(half_pi), 90_r );
   WALBERLA_CHECK_FLOAT_EQUAL( half_pi, degToRad(90_r) );

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char** argv )
{
   return walberla::main(argc, argv);
}
