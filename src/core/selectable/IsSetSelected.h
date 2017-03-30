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
//! \file IsSetSelected.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Set.h"

namespace walberla {
namespace selectable {

template <typename T>
bool isSetSelected( const Set<T> & state, const Set<T> & requiredSUIDs, const Set<T> & incompatibleSUIDSs )
{
   return state.contains( requiredSUIDs ) && !state.intersects( incompatibleSUIDSs );
}

} // namespace selectable
} // namespace walberla
