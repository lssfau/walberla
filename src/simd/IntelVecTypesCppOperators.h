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
//! \file IntelVecTypesCppOperators.h
//! \ingroup avx
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief C++ Operators for the Intel intrinsics vector types
//
//======================================================================================================================

#pragma once

#include "immintrin.h"

#include "waLBerlaDefinitions.h"

// Intel Compiler has no overloaded operators for vector type
#ifdef WALBERLA_CXX_COMPILER_IS_INTEL
inline __m256d operator+( __m256d a, __m256d b ) { return _mm256_add_pd ( a, b); }
inline __m256d operator-( __m256d a, __m256d b ) { return _mm256_sub_pd ( a, b); }
inline __m256d operator*( __m256d a, __m256d b ) { return _mm256_mul_pd ( a, b); }
inline __m256d operator/( __m256d a, __m256d b ) { return _mm256_div_pd ( a, b); }
#endif


