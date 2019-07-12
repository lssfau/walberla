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
//! \file Debug.h
//! \ingroup core
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Check whether a type implements a certain operator.
//
//======================================================================================================================

#pragma once

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <utility>

namespace walberla {
namespace debug {

//**********************************************************************************************************************
/*!
* The implementation is directly taken from https://stackoverflow.com/a/39348287
*/
//**********************************************************************************************************************

template<class X, class Y, class Op>
struct op_valid_impl
{
   template<class U, class L, class R>
   static auto test(int) -> decltype(std::declval<U>()(std::declval<L>(), std::declval<R>()), void(), std::true_type());

   template<class U, class L, class R>
   static auto test(...) -> std::false_type;

   using type = decltype(test<Op, X, Y>(0));

};

template<class X, class Y, class Op> using op_valid = typename op_valid_impl<X, Y, Op>::type;

namespace notstd {

struct left_shift {

   template <class L, class R>
   constexpr auto operator()(L&& l, R&& r) const
   noexcept(noexcept(std::forward<L>(l) << std::forward<R>(r)))
   -> decltype(std::forward<L>(l) << std::forward<R>(r))
   {
      return std::forward<L>(l) << std::forward<R>(r);
   }
};

struct right_shift {

   template <class L, class R>
   constexpr auto operator()(L&& l, R&& r) const
   noexcept(noexcept(std::forward<L>(l) >> std::forward<R>(r)))
   -> decltype(std::forward<L>(l) >> std::forward<R>(r))
   {
      return std::forward<L>(l) >> std::forward<R>(r);
   }
};

}

template<class X, class Y> using has_equality = op_valid<X, Y, std::equal_to<>>;
template<class X, class Y> using has_inequality = op_valid<X, Y, std::not_equal_to<>>;
template<class X, class Y> using has_less_than = op_valid<X, Y, std::less<>>;
template<class X, class Y> using has_less_equal = op_valid<X, Y, std::less_equal<>>;
template<class X, class Y> using has_greater_than = op_valid<X, Y, std::greater<>>;
template<class X, class Y> using has_greater_equal = op_valid<X, Y, std::greater_equal<>>;
template<class X, class Y> using has_bit_xor = op_valid<X, Y, std::bit_xor<>>;
template<class X, class Y> using has_bit_or = op_valid<X, Y, std::bit_or<>>;
template<class X, class Y> using has_left_shift = op_valid<X, Y, notstd::left_shift>;
template<class X, class Y> using has_right_shift = op_valid<X, Y, notstd::right_shift>;


} // namespace debug
} // namespace walberla
