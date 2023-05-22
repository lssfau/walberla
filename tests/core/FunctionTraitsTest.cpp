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
   static bool checkParameter1( typename std::enable_if< (FunctionTraits<F>::arity > 1 ), T >::type * = nullptr) {

      // The keyword "template" before "argument<1>" is crucial when these functions are inside a class. If the
      // keyword is dropped, compilers interpret "<1" as an arithmetic expression and therefore throw errors.
      typedef typename FunctionTraits<F>::template argument<1>::type argType;
      return std::is_same< T, argType >::value;
   }

   template< typename T>
   static bool checkParameter1( typename std::enable_if< (FunctionTraits<F>::arity <= 1 ), T >::type * = 0) {
      return false;
   }

};


int main( int /*argc*/, char** /*argv*/ )
{
   debug::enterTestMode();

   // obtain a function pointer's type
   typedef typename std::remove_pointer< bool (*)(int, float, double) >::type FuncType;

   // check return type
   constexpr auto retTypeEqual = std::is_same< FunctionTraits<FuncType>::return_type, bool>::value;
   WALBERLA_CHECK( retTypeEqual );

   // check arity
   constexpr auto argNumber = FunctionTraits<FuncType>::arity;
   WALBERLA_CHECK_EQUAL( argNumber, 3 );

   // check argument types
   constexpr auto arg0TypeEqual = std::is_same< FunctionTraits<FuncType>::argument<0>::type, int >::value;
   WALBERLA_CHECK ( arg0TypeEqual );

   constexpr auto arg1TypeEqual = std::is_same< FunctionTraits<FuncType>::argument<1>::type, float >::value;
   WALBERLA_CHECK ( arg1TypeEqual );

   constexpr auto arg2TypeEqual = std::is_same< FunctionTraits<FuncType>::argument<2>::type, double >::value;
   WALBERLA_CHECK ( arg2TypeEqual );

   // check usage inside class
   WALBERLA_CHECK( SomeClass<FuncType>::checkParameter1<float>() );
   WALBERLA_CHECK( ! (SomeClass<FuncType>::checkParameter1<double>()) );

   return EXIT_SUCCESS;
}