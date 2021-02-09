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
//! \file demangle.h
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include <string>

#ifdef __GLIBCXX__ 
#define HAVE_CXXABI_H
#include <cxxabi.h>
#else
#ifdef __has_include
#if __has_include(<cxxabi.h>)
#define HAVE_CXXABI_H
#include <cxxabi.h>
#endif
#endif
#endif

namespace walberla {
namespace debug {

inline std::string demangle( const std::string & name )
{
#ifdef HAVE_CXXABI_H
   int status = 0;
   std::size_t size = 0;
   const char * demangled = abi::__cxa_demangle( name.c_str(), nullptr, &size, &status );
   if( demangled == nullptr )
   {
      return name;
   }
   std::string demangled_str(demangled);
   std::free( const_cast<char*>(demangled) );
   return demangled_str;
#else
   return name;
#endif
}

} // namespace debug
} // namespace walberla
