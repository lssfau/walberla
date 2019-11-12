//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) OPTIONAL later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  OPTIONAL WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Variant.h
//! \ingroup core
//! \author Stephan Seitz <stephan.seitz@fau.de>
//
//======================================================================================================================

#pragma once


#if defined(WALBERLA_USE_STD_VARIANT)
#include <variant>
#else
#include <boost/variant.hpp>
#endif



namespace walberla
{

#if defined(WALBERLA_USE_STD_VARIANT)
using std::variant;
using std::visit;
using std::get;
using std::holds_alternative;
using std::bad_variant_access;
#else
using boost::variant;
using boost::get;
template <class T, class... Types>
constexpr bool holds_alternative( const boost::variant<Types...>& v ) noexcept
{
   return v.type() == typeid( T );
}
using bad_variant_access = boost::bad_get;
template<typename Visitor, typename... Variant>
decltype( auto ) visit( Visitor&& visitor, Variant&& ... variant )
{
   return boost::apply_visitor( visitor, variant... );
}
#endif

}
