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
//! \file Utility.h
//! \ingroup core
//! \author Klaus Iglberger
//! \brief Header file for mathematical functions and constants
//
//======================================================================================================================

#pragma once

#include "MathTrait.h"
#include "core/DataTypes.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>


namespace walberla {
namespace math {

//======================================================================================================================
//
//  MATHEMATICAL UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\name Mathematical utility functions */
//@{
template< typename T >
inline const T sign( T a );

template< typename T >
inline const typename std::enable_if<  std::is_unsigned<T>::value, T >::type abs( T a );
template< typename T >
inline const typename std::enable_if< ! std::is_unsigned<T>::value, T >::type abs( T a );


template< typename T1, typename T2 >
inline const typename MathTrait<T1,T2>::High min( const T1& a, const T2& b );

template< typename T1, typename T2, typename T3 >
inline const typename MathTrait< typename MathTrait<T1,T2>::High, T3 >::High min( const T1& a, const T2& b, const T3& c );

template< typename T1, typename T2 >
inline const typename MathTrait<T1,T2>::High max( const T1& a, const T2& b );

template< typename T1, typename T2, typename T3 >
inline const typename MathTrait< typename MathTrait<T1,T2>::High, T3 >::High max( const T1& a, const T2& b, const T3& c );

template< typename T >
inline const T sqr( const T& a );

template< typename T1, typename T2 >
inline bool equal( T1 a, T2 b );

inline real_t round( real_t a );
//@}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Sign function.
// \ingroup math
//
// \param a The signed value.
// \return 1 if the value is greater than or equal to zero, -1 if the value is smaller than zero.
//
// The sign function only works for signed built-in data types. The attempt to use unsigned data
// types or user-defined class types will result in a compile time error.
//
// http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
 */
template< typename T >
inline const T sign( T a )
{
   WALBERLA_STATIC_ASSERT( std::numeric_limits<T>::is_signed );
   return ( a < T(0) )?( T(-1) ):( T(1) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Absolute value function.
// \ingroup math
//
// \param a The value.
// \return The value if it is greater than or equal to zero, -1 times the value if the value is smaller than zero.
 */
template< typename T >
inline const typename std::enable_if<  std::is_unsigned<T>::value, T >::type abs( T a )
{
   return a;
}

template< typename T >
inline const typename std::enable_if< !std::is_unsigned<T>::value, T >::type abs( T a )
{
   return std::abs( a );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Minimum function for two arguments.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return The minimum of the two values.
//
// This function returns the minimum of the two given data values. The return type of the
// function is determined by the data types of the given arguments (for further detail see
// the MathTrait class description).
 */
template< typename T1, typename T2 >
inline const typename MathTrait<T1,T2>::High min( const T1& a, const T2& b )
{
   // The return type of the function is only a copy of the one of the arguments for two reasons:
   //  - in case the data types T1 and T2 are equal, a reference return type could result in a
   //    bug if combined with literals
   //  - in case the two data types are unequal, the result of the comparison could be converted
   //    to the more significant data type, which results in a local temporary value
   // The copy operation might cause a performance decrease for class types, which is probably
   // avoided if the function is inlined.
   return ( a < b )?( a ):( b );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Minimum function for three arguments.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \param c Third value
// \return The minimum of the three values.
//
// This function returns the minimum of the three given data values. The return type of the
// function is determined by the data types of the given arguments (for further detail see
// the MathTrait class description).
 */
template< typename T1, typename T2, typename T3 >
inline const typename MathTrait< typename MathTrait<T1,T2>::High, T3 >::High min( const T1& a, const T2& b, const T3& c )
{
   // The return type of the function is only a copy of the one of the arguments for two reasons:
   //  - in case the data types T1, T2, and T3 are equal, a reference return type could result in
   //    a bug if combined with literals
   //  - in case the three data types are unequal, the result of the comparison could be converted
   //    to the more significant data type, which results in a local temporary value
   // The copy operation might cause a performance decrease for class types, which is probably
   // avoided if the function is inlined.
   return ( a < b )?( ( a < c )?( a ):( c ) ):( ( b < c )?( b ):( c ) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Maximum function for two arguments.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return The maximum of the two values.
//
// This function returns the maximum of the two given data values. The return type of the
// function is determined by the data types of the given arguments (for further detail see
// the MathTrait class description).
 */
template< typename T1, typename T2 >
inline const typename MathTrait<T1,T2>::High max( const T1& a, const T2& b )
{
   // The return type of the function is only a copy of the one of the arguments for two reasons:
   //  - in case the data types T1 and T2 are equal, a reference return type could result in a
   //    bug if combined with literals
   //  - in case the two data types are unequal, the result of the comparison could be converted
   //    to the more significant data type, which results in a local temporary value
   // The copy operation might cause a performance decrease for class types, which is probably
   // avoided if the function is inlined.
   return ( a > b )?( a ):( b );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Maximum function for three arguments.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \param c Third value.
// \return The maximum of the three values.
//
// This function returns the maximum of the three given data values. The return type of the
// function is determined by the data types of the given arguments (for further detail see
// the MathTrait class description).
 */
template< typename T1, typename T2, typename T3 >
inline const typename MathTrait< typename MathTrait<T1,T2>::High, T3 >::High max( const T1& a, const T2& b, const T3& c )
{
   // The return type of the function is only a copy of the one of the arguments for two reasons:
   //  - in case the data types T1, T2, and T3 are equal, a reference return type could result in
   //    a bug if combined with literals
   //  - in case the three data types are unequal, the result of the comparison could be converted
   //    to the more significant data type, which results in a local temporary value
   // The copy operation might cause a performance decrease for class types, which is probably
   // avoided if the function is inlined.
   return ( a > b )?( ( a > c )?( a ):( c ) ):( ( b > c )?( b ):( c ) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Square function.
// \ingroup math
//
// \param a Value to be squared.
// \return The square of the input value.
 */
template< typename T >
inline const T sqr( const T& a )
{
   return a*a;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*! \cond internal */
/*!\brief Equality check for integral data types.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of two integral values.
 */
template< typename T >
inline bool equal_backend( T a, T b )
{
   return a == b;
}
/*! \endcond */
//**********************************************************************************************************************


//**********************************************************************************************************************
/*! \cond internal */
/*!\brief Equality check for single precision floating point values.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of two single precision floating point numbers. Due to the
// limited machine accuracy, a direct comparison of two floating point numbers should be avoided.
// This functions offers the possibility to compare two floating-point values with a certain
// accuracy margin.
 */
template<>
inline bool equal_backend<float>( float a, float b )
{
   return std::fabs( a - b ) < 1E-8;
}
/*! \endcond */
//**********************************************************************************************************************


//**********************************************************************************************************************
/*! \cond internal */
/*!\brief Equality check for double precision floating point values.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of two double precision floating point numbers. Due to the
// limited machine accuracy, a direct comparison of two floating point numbers should be avoided.
// This functions offers the possibility to compare two floating-point values with a certain
// accuracy margin.
 */
template<>
inline bool equal_backend<double>( double a, double b )
{
   return std::fabs( a - b ) < 1E-8;
}
/*! \endcond */
//**********************************************************************************************************************


//**********************************************************************************************************************
/*! \cond internal */
/*!\brief Equality check for long double precision floating point values.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the two values are equal, \a false if not.
//
// Equal function for the comparison of two long double precision floating point numbers. Due
// to the limited machine accuracy, a direct comparison of two floating point numbers should be
// avoided. This functions offers the possibility to compare two floating-point values with a
// certain accuracy margin.
 */
template<>
inline bool equal_backend<long double>( long double a, long double b )
{
   return std::fabs( a - b ) < 1E-10;
}
/*! \endcond */
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Generic equality check.
// \ingroup math
//
// \param a First value.
// \param b Second value.
// \return \a true if the two values are equal, \a false if not.
//
// Generic equal function for the comparison of two numeric values. Depending on the types of
// the two arguments, a special comparison for floating point values is selected that takes
// the limited machine accuracy into account.
 */
template< typename T1, typename T2 >
inline bool equal( T1 a, T2 b )
{
   typedef typename MathTrait<T1,T2>::High High;
   return equal_backend<High>( a, b );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Rounding function.
// \ingroup math
//
// \param  a Floating point value.
// \return The floating point value rounded towards the next whole number.

 */
inline real_t round( real_t a )
{
   return std::floor( a + real_t(0.5) );
}
//**********************************************************************************************************************



} // namespace math
} // namespace walberla
