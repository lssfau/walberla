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
//! \file DataTypesTest.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Uint.h"

#include <cstdlib>
#include <iostream>


using namespace walberla;

uint_t uintBitsLd( const uint_t bits )
{
   WALBERLA_ASSERT( bits >= 1 );
   return ( bits == 1 ) ? 0 : ( 1 + uintBitsLd( bits >> 1 ) );
}

int main( int /*argc*/, char** /*argv*/ )
{
   debug::enterTestMode();

   WALBERLA_CHECK_EQUAL( math::UINT_BITS_LD, uintBitsLd( math::UINT_BITS ) );

   WALBERLA_CHECK_EQUAL( math::int_ld< 2048 >::exp, 11 );

   WALBERLA_CHECK_EQUAL( math::int_ld< 64 >::exp, 6 );
   WALBERLA_CHECK_EQUAL( math::int_ld< 32 >::exp, 5 );
   WALBERLA_CHECK_EQUAL( math::int_ld< 16 >::exp, 4 );
   WALBERLA_CHECK_EQUAL( math::int_ld<  8 >::exp, 3 );
   WALBERLA_CHECK_EQUAL( math::int_ld<  4 >::exp, 2 );
   WALBERLA_CHECK_EQUAL( math::int_ld<  2 >::exp, 1 );
   WALBERLA_CHECK_EQUAL( math::int_ld<  1 >::exp, 0 );

   WALBERLA_CHECK_IDENTICAL(1.23456_r, real_c(1.23456));
   WALBERLA_CHECK_IDENTICAL(1_r, real_c(1));
   WALBERLA_CHECK_IDENTICAL(-1_r, real_c(-1));

   return EXIT_SUCCESS;
}
