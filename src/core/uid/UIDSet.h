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
//! \file UIDSet.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Feichtinger
//
//======================================================================================================================

#pragma once

#include "UID.h"
#include "core/Set.h"


namespace walberla {
namespace uid {



template< typename T >
inline Set< UID<T> > operator+( const UID<T>& a, const UID<T>& b ) {

   return setUnion( Set< UID<T> >(a), Set< UID<T> >(b) );
}



template< typename T >
inline Set< UID<T> > operator&( const UID<T>& a, const UID<T>& b ) {

   return setIntersection( Set< UID<T> >(a), Set< UID<T> >(b) );
}



} // namespace uid
} // namespace walberla
