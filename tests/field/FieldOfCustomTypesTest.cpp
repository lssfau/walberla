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
//! \file FieldOfCustomTypesTest.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "field/GhostLayerField.h"
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include <iostream>
#include <set>


using namespace walberla;


using std::cout;
using std::endl;




struct MyClass
{
   MyClass() {
      ++constructorCalls;
   }

   ~MyClass() {
      ++destructorCalls;
   }

   static uint_t constructorCalls;
   static uint_t destructorCalls;
};
uint_t MyClass::constructorCalls = 0;
uint_t MyClass::destructorCalls  = 0;



void testCorrectNumberOfConstructorCalls()
{
   MyClass::constructorCalls = 0;
   MyClass::destructorCalls  = 0;

   const uint_t xs = 7;
   const uint_t ys = 5;
   const uint_t zs = 3;
   const uint_t fs = 2;
   const uint_t gl = 2;
   const uint_t expectedCalls = ( xs + 2*gl) * ( ys + 2*gl ) * (zs + 2*gl ) * fs;

   {
      GhostLayerField<MyClass,fs>  field( xs, ys, zs, gl, field::fzyx,
                                          shared_ptr<field::AllocateAligned<MyClass,32> >() );

      WALBERLA_CHECK_EQUAL( expectedCalls, MyClass::constructorCalls );
      auto clonedField = field.cloneUninitialized();
      WALBERLA_CHECK_EQUAL( 2* expectedCalls, MyClass::constructorCalls );

      delete clonedField;
   }
   WALBERLA_CHECK_EQUAL( 2* expectedCalls, MyClass::destructorCalls );

   MyClass::constructorCalls = 0;
   MyClass::destructorCalls  = 0;
}



int main( int argc, char**argv )
{
   walberla::Environment walberlaEnv( argc, argv );
   debug::enterTestMode();

   testCorrectNumberOfConstructorCalls();

   return 0;
}
