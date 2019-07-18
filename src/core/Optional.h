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
//! \file Optional.h
//! \ingroup core
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once


#if defined(WALBERLA_USE_STD_OPTIONAL)
#include <optional>
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_OPTIONAL)
#undef _LIBCPP_WARN_ON_DEPRECATED_EXPERIMENTAL_HEADER
#include <experimental/optional>
#else
#include <boost/optional.hpp>
#endif



namespace walberla {

#if defined(WALBERLA_USE_STD_OPTIONAL)
using std::optional;
using std::nullopt;
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_OPTIONAL)
using std::experimental::optional;
using std::experimental::nullopt;
#else
using boost::optional;
const boost::none_t nullopt = boost::none;
#endif

}
