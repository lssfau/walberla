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
//! \file Constants.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for mathematical constants
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <core/DataTypes.h>

#include <cmath>

// Disable false warnings in GCC 5
#if ( defined __GNUC__ ) && ( __GNUC__ == 5 ) && ( __GNUC_MINOR__ == 1 )
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace walberla {
namespace math {

# undef M_E
# undef M_LOG2E
# undef M_LOG10E
# undef M_LN2
# undef M_LN10
# undef M_PI
# undef M_PI_2
# undef M_PI_4
# undef M_1_PI
# undef M_2_PI
# undef M_2_SQRTPI
# undef M_SQRT2
# undef M_SQRT1_2

constexpr real_t M_E          = real_t( 2.7182818284590452354  );  /* e */
constexpr real_t M_LOG2E      = real_t( 1.4426950408889634074  );  /* log_2 e */
constexpr real_t M_LOG10E     = real_t( 0.43429448190325182765 );  /* log_10 e */
constexpr real_t M_LN2        = real_t( 0.69314718055994530942 );  /* log_e 2 */
constexpr real_t M_LN10       = real_t( 2.30258509299404568402 );  /* log_e 10 */
constexpr real_t M_PI         = real_t( 3.14159265358979323846 );  /* pi */
constexpr real_t M_PI_2       = real_t( 1.57079632679489661923 );  /* pi/2 */
constexpr real_t M_PI_4       = real_t( 0.78539816339744830962 );  /* pi/4 */
constexpr real_t M_1_PI       = real_t( 0.31830988618379067154 );  /* 1/pi */
constexpr real_t M_2_PI       = real_t( 0.63661977236758134308 );  /* 2/pi */
constexpr real_t M_2_SQRTPI   = real_t( 1.12837916709551257390 );  /* 2/sqrt(pi) */
constexpr real_t M_SQRT2      = real_t( 1.41421356237309504880 );  /* sqrt(2) */
constexpr real_t M_SQRT1_2    = real_t( 0.70710678118654752440 );  /* 1/sqrt(2) */

} // namespace math
}
