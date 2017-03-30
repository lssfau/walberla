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
//! \ingroup core
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
#include <boost/math/constants/constants.hpp>

// Disable false warnings in GCC 5
#if ( defined __GNUC__ ) && ( __GNUC__ == 5 ) && ( __GNUC_MINOR__ == 1 )
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace walberla {
namespace math {

//=================================================================================================
//
//  MATHEMATICAL CONSTANT E
//
//=================================================================================================

#ifdef M_E
#  undef M_E
#endif

//*************************************************************************************************
/*!\brief Definition of the mathematical constant \f$ e \f$.
 * \ingroup math
 */
const real_t M_E = boost::math::constants::e<real_t>();
//*************************************************************************************************




#if 0
//=================================================================================================
//
//  MATHEMATICAL CONSTANT LOG2E
//
//=================================================================================================

#ifdef M_LOG2E
#  undef M_LOG2E
#endif

//*************************************************************************************************
/*!\brief Definition of the mathematical constant \f$ \log_2 e \f$.
 * \ingroup math
 */
const real_t M_LOG2E = 1.4426950408889634073599246810018921;
//*************************************************************************************************
#endif




//=================================================================================================
//
//  MATHEMATICAL CONSTANT LOG10E
//
//=================================================================================================

#if 0
#ifdef M_LOG10E
#  undef M_LOG10E
#endif

//*************************************************************************************************
/*!\brief Definition of the mathematical constant \f$ \log_{10} e \f$.
 * \ingroup math
 */
const real_t M_LOG10E = boost::math::constants::log10_e<real_t>();
//*************************************************************************************************
#endif




//=================================================================================================
//
//  MATHEMATICAL CONSTANT LN2
//
//=================================================================================================

#ifdef M_LN2
#  undef M_LN2
#endif

//*************************************************************************************************
/*!\brief Definition of the mathematical constant \f$ \ln 2 \f$.
 * \ingroup math
 */
const real_t M_LN2 = boost::math::constants::ln_two<real_t>();
//*************************************************************************************************




//=================================================================================================
//
//  MATHEMATICAL CONSTANT LN10
//
//=================================================================================================

#if 0
#ifdef M_LN10
#  undef M_LN10
#endif

//*************************************************************************************************
/*!\brief Definition of the mathematical constant \f$ \ln 10 \f$.
 * \ingroup math
 */
const real_t M_LN10 = boost::math::constants::one_div_log10_e<real_t>();
//*************************************************************************************************
#endif




//=================================================================================================
//
//  MATHEMATICAL CONSTANT PI
//
//=================================================================================================

#ifdef M_PI
#  undef M_PI
#endif

//*************************************************************************************************
/*!\brief Definition of the mathematical constant \f$ \pi \f$.
 * \ingroup math
 */
const real_t M_PI = boost::math::constants::pi<real_t>();
//*************************************************************************************************




//=================================================================================================
//
//  MATHEMATICAL CONSTANT SQRT2
//
//=================================================================================================

#ifdef M_SQRT2
#  undef M_SQRT2
#endif

//*************************************************************************************************
/*!\brief Definition of the mathematical constant \f$ \sqrt{2} \f$.
 * \ingroup math
 */
const real_t M_SQRT2 = boost::math::constants::root_two<real_t>();
//*************************************************************************************************




//=================================================================================================
//
//  MATHEMATICAL CONSTANT SQRT3
//
//=================================================================================================

#if 0
#ifdef M_SQRT3
#  undef M_SQRT3
#endif

//*************************************************************************************************
/*!\brief Definition of the mathematical constant \f$ \sqrt{3} \f$.
 * \ingroup math
 */
const real_t M_SQRT3 = boost::math::constants::root_three<real_t>();
//*************************************************************************************************
#endif

} // namespace math
}
