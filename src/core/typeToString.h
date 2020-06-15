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
//! \file typeToString.h
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

namespace walberla {

// data type to string conversion

template< typename T > inline const char* typeToString();

#define TypeToString(X) template<> inline const char* typeToString< X >() { \
   static char string[] = #X; \
   return string; \
}

TypeToString(bool)
TypeToString(char)
TypeToString(short)
TypeToString(int)
TypeToString(long)
TypeToString(long long)
TypeToString(unsigned char)
TypeToString(unsigned short)
TypeToString(unsigned int)
TypeToString(unsigned long)
TypeToString(unsigned long long)
TypeToString(float)
TypeToString(double)

#undef TypeToString

template< typename T > inline const char* typeToString( T ) { return typeToString<T>(); }

} // namespace walberla
