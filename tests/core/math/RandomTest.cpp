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
//! \file RandomTest.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"


using namespace walberla;

template< typename INT >
void runInt()
{
   math::seedRandomGenerator( 23 );

   std::vector< INT > v;
   for( int i = 0; i != 100; ++i )
      v.push_back( math::intRandom< INT >() );

   math::seedRandomGenerator( 23 );

   for( int i = 0; i != 100; ++i )
      WALBERLA_CHECK_EQUAL( v[ uint_c(i) ], math::intRandom< INT >() );
}

template< typename REAL >
void runReal()
{
   math::seedRandomGenerator( 23 );

   std::vector< REAL > v;
   for( int i = 0; i != 100; ++i )
      v.push_back( math::realRandom< REAL >() );

   math::seedRandomGenerator( 23 );

   for( int i = 0; i != 100; ++i )
      WALBERLA_CHECK_IDENTICAL( v[uint_c(i)], math::realRandom< REAL >() );
}

int main( int /*argc*/, char** /*argv*/ )
{
   debug::enterTestMode();

   runInt< walberla::uint8_t >();
   runInt< walberla::uint16_t >();
   runInt< walberla::uint32_t >();
   runInt< walberla::uint64_t >();

   runInt< walberla::uint_t >();

   runInt< walberla::int8_t >();
   runInt< walberla::int16_t >();
   runInt< walberla::int32_t >();
   runInt< walberla::int64_t >();

   runInt< int >();

   runReal< float >();
   runReal< double >();
   runReal< walberla::real_t >();

   return EXIT_SUCCESS;
}
