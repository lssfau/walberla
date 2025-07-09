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
//! \file FunctionTraits.h
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/FunctionTraits.h"

#include <type_traits>

using namespace walberla;

// FunctionTraits are used in a similar way in gpu/Kernel.h. As explained below, special attention is required.
template< typename F>
struct SomeClass
{
   template< typename T>
   requires( FunctionTraits<F>::arity > 1 )
   static bool checkParameter1() {

      // The keyword "template" before "argument<1>" is crucial when these functions are inside a class. If the
      // keyword is dropped, compilers interpret "<1" as an arithmetic expression and therefore throw errors.
      using argType = typename FunctionTraits<F>::template argument<1>::type;
      return std::is_same_v< T, argType >;
   }

   template< typename T>
   requires( FunctionTraits<F>::arity <= 1 )
   static bool checkParameter1() {
      return false;
   }

};


int main( int /*argc*/, char** /*argv*/ )
{
   debug::enterTestMode();

   // obtain a function pointer's type
   using FuncType = std::remove_pointer_t<bool (*)(int, float, double)>;

   // check return type
   constexpr auto retTypeEqual = std::is_same_v< FunctionTraits<FuncType>::return_type, bool>;
   WALBERLA_CHECK( retTypeEqual );

   // check arity
   constexpr auto argNumber = FunctionTraits<FuncType>::arity;
   WALBERLA_CHECK_EQUAL( argNumber, 3 );

   // check argument types
   constexpr auto arg0TypeEqual = std::is_same_v< FunctionTraits<FuncType>::argument<0>::type, int >;
   WALBERLA_CHECK ( arg0TypeEqual );

   constexpr auto arg1TypeEqual = std::is_same_v< FunctionTraits<FuncType>::argument<1>::type, float >;
   WALBERLA_CHECK ( arg1TypeEqual );

   constexpr auto arg2TypeEqual = std::is_same_v< FunctionTraits<FuncType>::argument<2>::type, double >;
   WALBERLA_CHECK ( arg2TypeEqual );

   // check usage inside class
   WALBERLA_CHECK( SomeClass<FuncType>::checkParameter1<float>() );
   WALBERLA_CHECK( ! (SomeClass<FuncType>::checkParameter1<double>()) );

   return EXIT_SUCCESS;
}
