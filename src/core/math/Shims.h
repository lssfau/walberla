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
//! \file Sample.h
//! \ingroup core
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for shims
//
//======================================================================================================================

#pragma once

#include <type_traits>

namespace walberla {
namespace math {

//=================================================================================================
//
//  CLEAR SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Clearing the given value/object to the default state.
 * \ingroup math_shims
 *
 * \param clearable The value/object to be cleared.
 * \return void
 *
 * The clear shim represents an abstract interface for clearing a value/object of any given
 * data type to its default state. Values of built-in data type are reset to zero.
 */
template< typename Type >
inline void clear( Type& clearable )
{
   clearable = Type(0);
}
//*************************************************************************************************

//=================================================================================================
//
//  INVERT SHIMS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Inverting the given single precision value.
 * \ingroup math_shims
 *
 * \param a The single precision value to be inverted.
 * \return The inverse of the given value.
 *
 * The invert shim represents an abstract interface for inverting a value/object of any given
 * data type. For single precision floating point values this results in \f$ \frac{1}{a} \f$.
 */
inline float inv( float a )
{
   return ( 1.0F / a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given double precision value.
 * \ingroup math_shims
 *
 * \param a The double precision value to be inverted.
 * \return The inverse of the given value.
 *
 * The invert shim represents an abstract interface for inverting a value/object of any given
 * data type. For double precision floating point values this results in \f$ \frac{1}{a} \f$.
 */
inline double inv( double a )
{
   return ( 1.0 / a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given long double value.
 * \ingroup math_shims
 *
 * \param a The long double value to be inverted.
 * \return The inverse of the given value.
 *
 * The invert shim represents an abstract interface for inverting a value/object of any given
 * data type. For long double floating point values this results in \f$ \frac{1}{a} \f$.
 */
inline long double inv( long double a )
{
   return ( 1.0L / a );
}
//*************************************************************************************************

//=================================================================================================
//
//  ISDEFAULT SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the given value/object is in default state.
 * \ingroup math_shims
 *
 * \param v The value/object to be tested for its default state.
 * \return \a true in case the given value/object is in its default state, \a false otherwise.
 *
 * The isDefault shim represents an abstract interface for testing a value/object whether
 * it is in its default state or not. In case the value/object is in its default state, the
 * function returns \a true, otherwise it returns \a false. For built-in data types, the
 * function returns \a true in case the current value is zero.

   \code
   const int i = 0;          // isDefault( i ) returns true
   double    d = 2.0;        // isDefault( d ) returns false
   Vec3      v1;             // isDefault( v1 ) returns true
   Vec3      v2( 0, 0, 0 );  // isDefault( v2 ) returns true since (0,0,0) is the default state
   Vec3      v3( 1, 2, 3 );  // isDefault( v3 ) returns false
   \endcode
 */
template< typename Type >
inline bool isDefault( const Type& v )
{
   return v == Type(0);
}
//*************************************************************************************************

//=================================================================================================
//
//  RESET SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Resetting the given value/object to the default value.
 * \ingroup math_shims
 *
 * \param resettable The value/object to be resetted.
 * \return void
 *
 * The reset shim represents an abstract interface for the resetting of a value/object of
 * any given data type to its default value. Values of built-in data type are reset to zero.
 */
template< typename Type >
inline void reset( Type& resettable )
{
   resettable = Type(0);
}
//*************************************************************************************************

//=================================================================================================
//
//  SQUARE SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Squaring the given value/object.
 * \ingroup math_shims
 *
 * \param a The value/object to be squared.
 * \return The result of the square operation.
 *
 * The square shim represents an abstract interface for squaring a value/object of any given
 * data type. For values of built-in data type this results in a plain multiplication.
 */
template< typename T >
inline const T sq( const T& a )
{
   return a*a;
}
//*************************************************************************************************

} // namespace math
}
