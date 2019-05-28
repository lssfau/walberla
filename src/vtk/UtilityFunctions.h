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
//! \file UtilityFunctions.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "VTKTrait.h"

#include "core/DataTypes.h"

#include <ostream>
#include <string>


namespace walberla {
namespace vtk {



// "toStream" should be used in order to prevent 8bit signed/unsigned data types to be interpreted as a character rather than a number
template< typename T > inline void toStream( std::ostream& os, const T       value ) { os << value; }
template<>             inline void toStream( std::ostream& os, const  int8_t value ) { os <<  int_c( value ); }
template<>             inline void toStream( std::ostream& os, const uint8_t value ) { os << uint_c( value ); }



template< typename T >
inline std::string typeToString()
{
   return VTKTrait<T>::type_string;
}

} // namespace vtk
} // namespace walberla
