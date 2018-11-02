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
//! \file Macros.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

// restrict keyword

#ifdef __GNUC__
#define WALBERLA_RESTRICT __restrict__
#elif _MSC_VER
#define WALBERLA_RESTRICT __restrict
#else
#define WALBERLA_RESTRICT
#endif

// forced inline for different compilers
#if defined(__INTEL_COMPILER)
# define  WALBERLA_FORCE_INLINE(func) __forceinline func
#elif defined(_MSC_VER)
#  define WALBERLA_FORCE_INLINE(func) __forceinline func
#elif defined(__GNUC__)
#  define WALBERLA_FORCE_INLINE(func) inline func __attribute__ ((always_inline))
#else
#  pragma message("WARNING: You need to implement WALBERLA_FORCE_INLINE for this compiler!")
#  define WALBERLA_FORCE_INLINE(func) inline func
#endif

// pragma in macros (-> https://stackoverflow.com/a/3030312)

#ifdef WALBERLA_CXX_COMPILER_IS_MSVC
#define WALBERLA_PRAGMA(x) __pragma(x)
#else
#define WALBERLA_PRAGMA(x) _Pragma(#x)
#endif

#ifndef WALBERLA_CXX_COMPILER_IS_IBM
#define WALBERLA_OVERRIDE
#else
#define WALBERLA_OVERRIDE override
#endif

// macro overloading (-> https://stackoverflow.com/a/24028231)

#define WALBERLA_GLUE(x, y) x y

#define WALBERLA_RETURN_ARG_COUNT(_1_, _2_, _3_, _4_, _5_, _6_, _7_, _8_, _9_, _10_, _11_, _12_, _13_, _14_, _15_, _16_, _17_, _18_, _19_, _20_, count, ...) count
#define WALBERLA_EXPAND_ARGS(args) WALBERLA_RETURN_ARG_COUNT args
#define WALBERLA_COUNT_ARGS_MAX20(...) WALBERLA_EXPAND_ARGS((__VA_ARGS__, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0))

#define WALBERLA_OVERLOAD_MACRO2(name, count) name##count
#define WALBERLA_OVERLOAD_MACRO1(name, count) WALBERLA_OVERLOAD_MACRO2(name, count)
#define WALBERLA_OVERLOAD_MACRO(name, count) WALBERLA_OVERLOAD_MACRO1(name, count)

#define WALBERLA_MACRO_OVERLOAD(name, ...) WALBERLA_GLUE(WALBERLA_OVERLOAD_MACRO(name, WALBERLA_COUNT_ARGS_MAX20(__VA_ARGS__)), (__VA_ARGS__))
