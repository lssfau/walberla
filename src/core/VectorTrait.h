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
//! \file VectorTrait.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include <boost/type_traits/is_arithmetic.hpp>

#include <boost/utility/enable_if.hpp>


namespace walberla {


//**********************************************************************************************************************
/*! Provides information on how to serialize (=extract components) from a compound data type
*
* This is the general implementation for arithmetic data-types, for a specialization example see Vector3.h
*/
//**********************************************************************************************************************

template< typename T >
struct VectorTrait
{
   typedef T OutputType;

   static const uint_t F_SIZE = 1u;
   static T get   ( T   value, uint_t /*f*/ )        { return value; }
   static void set( T & value, uint_t /*f*/, T val ) { value = val;  }
   static_assert( boost::is_arithmetic<T>::value, "Specialize OutputTrait for your type!" );
};


} // namespace walberla
