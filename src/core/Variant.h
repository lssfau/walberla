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

// suppress conversion warning that is emitted by icpc compiler with `-Werror`
#include "waLBerlaDefinitions.h"
#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning(push)
#pragma warning disable 68
#endif

#include <variant>

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning(pop)
#endif


namespace walberla
{

using std::variant;
using std::visit;
using std::get;
using std::holds_alternative;
using std::bad_variant_access;

}
