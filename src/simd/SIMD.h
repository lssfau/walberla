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
//! \file SIMD.h
//! \ingroup simd
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h" // for WALBERLA_SIMD_FORCE_SCALAR


//===================================================================================================================
//
//  Macros Testing which instruction sets are available
//  WALBERLA_SIMD_*_AVAILABLE
//
//===================================================================================================================


#ifdef __AVX2__
#define WALBERLA_SIMD_AVX2_AVAILABLE 1
#endif


#ifdef __AVX__
#define WALBERLA_SIMD_AVX_AVAILABLE 1
#endif


#ifdef  __SSE4_2__
#define WALBERLA_SIMD_SSE4_AVAILABLE 1
#endif

#ifdef  __SSE2__
#define WALBERLA_SIMD_SSE2_AVAILABLE 1
#endif




//===================================================================================================================
//
//  Selection of best available instruction set
//  WALBERLA_USE_*
//
//===================================================================================================================

#ifndef WALBERLA_SIMD_FORCE_SCALAR

// AVX2 ( Intel Haswell )
#if defined( WALBERLA_SIMD_AVX2_AVAILABLE )
#include "AVX2.h"
#define WALBERLA_USE_AVX2 1
#define WALBERLA_USE_SIMD 1
#endif

// AVX ( Intel Sandybridge )
#if defined( WALBERLA_SIMD_AVX_AVAILABLE ) && !defined( WALBERLA_USE_SIMD )
#include "AVX.h"
#define WALBERLA_USE_AVX  1
#define WALBERLA_USE_SIMD 1
#endif

// SSE4 ( Intel )
#if defined( WALBERLA_SIMD_SSE4_AVAILABLE ) && !defined( WALBERLA_USE_SIMD )
#include "SSE4.h"
#define WALBERLA_USE_SSE4  1
#define WALBERLA_USE_SIMD 1
#endif

// SSE2 ( Intel )
#if defined( WALBERLA_SIMD_SSE2_AVAILABLE ) && !defined( WALBERLA_USE_SIMD )
#include "SSE2.h"
#define WALBERLA_USE_SSE2  1
#define WALBERLA_USE_SIMD 1
#endif

#endif //WALBERLA_SIMD_FORCE_SCALAR



// Fallback if no instruction set extension is available ( or WALBERLA_SIMD_FORCE_SCALAR is active )
#if !defined( WALBERLA_USE_SIMD )
#include "Scalar.h"
#define WALBERLA_USE_SCALAR_SIMD_EMULATION 1
#define WALBERLA_USE_SIMD 1
#endif





//===================================================================================================================
//
//  Including the instruction set header and import its symbols into simd namespace
//
//===================================================================================================================

namespace walberla {
namespace simd {

template<typename T> struct is_vector4_type {  static const bool value = false; };

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wignored-attributes"
#endif


#ifdef WALBERLA_USE_AVX2
  using namespace avx2;
  template<> struct is_vector4_type<avx2::double4_t> {  static const bool value = true; };
#endif


#ifdef WALBERLA_USE_AVX
   using namespace avx;
   template<> struct is_vector4_type<avx::double4_t> {  static const bool value = true; };
#endif


#ifdef WALBERLA_USE_SSE4
   using namespace sse4;
   template<> struct is_vector4_type<sse4::double4_t> {  static const bool value = true; };
#endif

#ifdef WALBERLA_USE_SSE2
   using namespace sse2;
   template<> struct is_vector4_type<sse2::double4_t> {  static const bool value = true; };
#endif


#ifdef WALBERLA_USE_SCALAR_SIMD_EMULATION
   using namespace scalar;
   template<> struct is_vector4_type<scalar::double4_t> {  static const bool value = true; };
#endif


#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#   pragma GCC diagnostic pop
#endif


inline bool available() {
#ifdef WALBERLA_USE_SIMD
   return true;
#else
   return false;
#endif
}




} // namespace simd
} // namespace walberla

