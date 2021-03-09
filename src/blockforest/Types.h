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
//! \file types.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include <limits>


namespace walberla {
namespace blockforest {



// data structure specific types

using workload_t = real_t;
using memory_t = real_t;

WALBERLA_STATIC_ASSERT( sizeof( workload_t ) == 4 || sizeof( workload_t ) == 8 );
WALBERLA_STATIC_ASSERT( sizeof( memory_t   ) == 4 || sizeof( memory_t   ) == 8 );

template< typename T > inline workload_t workload_c( T t ) { return numeric_cast< workload_t >(t); } ///< cast to type workload_t using "workload_c(x)"
template< typename T > inline memory_t   memory_c  ( T t ) { return numeric_cast< memory_t >(t); }   ///< cast to type memory_t   using "memory_c(x)"



} // namespace blockforest

using blockforest::workload_t;
using blockforest::memory_t;

} // namespace walberla


