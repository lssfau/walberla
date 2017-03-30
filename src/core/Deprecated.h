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
//! \file Deprecated.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once



/*******************************************************************************************************************//**
 * \brief   Macro to mark a function as deprecated
 *
 * To mark a function as deprecated wrap its signature with the WALBERLA_DEPRECATED macro at the definition.
 * E.g. to mark int f(const & foo) const; as deprecated use WALBERLA_DEPRECATED(int f(const & foo) const);
 * Note that this won't work for functions with a comma in their return type definition! E.g. to mark
 * std::pair<int,int> f(); as deprecated add a typedef std::pair<int, int> PairInt; first an then deprecate
 * via WALBERLA_DEPRECATE(PairInt f());
 * Deprecation of a function will produce a compiler warning if the function is used despite it's deprecation.
 *
 * \param   func  The function to be marked deprecated.
 **********************************************************************************************************************/

#ifdef __GNUC__
#  define WALBERLA_DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#  define WALBERLA_DEPRECATED(func) __declspec(deprecated) func
#elif defined(__IBMCPP__)
#  define WALBERLA_DEPRECATED(func) func // xlc++ (12.1) has no obvious possibility for deprecation
#elif defined(_SX)
#  define WALBERLA_DEPRECATED(func) __declspec(deprecated) func
#else
#  pragma message("WARNING: You need to implement WALBERLA_DEPRECATED for this compiler!")
#  define WALBERLA_DEPRECATED(func) func
#endif
