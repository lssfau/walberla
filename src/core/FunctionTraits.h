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
//! \ingroup core
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#pragma once

#include <cstddef>
#include <tuple>

namespace walberla {

//**********************************************************************************************************************
/*! Gives a function's return type as well as the number and type of arguments accepted by a function.
*
* This variadic template substitutes <boost/function_traits>.
*
* Usage:
* FunctionTraits<F>::return_type       Type returned by function type F.
* FunctionTraits<F>::arity             Number of arguments accepted by function type F.
* FunctionTraits<F>::argument<N>       Type of the Nth argument of function type F with 0 <= N < arity of F.
*
*/
//**********************************************************************************************************************

template< typename F >
struct FunctionTraits;

template< typename R, typename ...Args >
struct FunctionTraits< R( Args... ) >
{
   using return_type = R;

   static constexpr std::size_t arity = sizeof...(Args);

   template< std::size_t N >
   struct argument
   {
      static_assert(N < arity, "Error: Parameter index is not valid!");
      using type = typename std::tuple_element< N, std::tuple<Args...> >::type;
   };
};

} // namespace walberla