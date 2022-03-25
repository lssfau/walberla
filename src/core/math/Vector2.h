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
//! \file Vector2.h
//! \ingroup core
//! \author Klaus Iglberger
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Header file for the implementation of a 2D vector
//
//======================================================================================================================

#pragma once

#include "FPClassify.h"
#include "MathTrait.h"
#include "SqrtTrait.h"
#include "Utility.h"

#include "core/DataTypes.h"
#include "core/VectorTrait.h"
#include "core/debug/Debug.h"
#include "core/mpi/Datatype.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <type_traits>

namespace walberla {
namespace math {


//**********************************************************************************************************************
// Definitions
//**********************************************************************************************************************

//! High-order return value.
/*! Abbreviation for the evaluation of the higher-order data type in a numerical operation. */
#define HIGH typename MathTrait<Type,Other>::High



//======================================================================================================================
//
//  CLASS DEFINITION
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Efficient, generic implementation of a 2-dimensional vector.
// \ingroup math
//
// The Vector2 class is the representation of a 2D vector with a total of 2 statically allocated
// elements of arbitrary type. The naming convention of the elements is as follows:

                             \f[\left(\begin{array}{*{2}{c}}
                             x & y \\
                             \end{array}\right)\f]

// These elements can be accessed directly with the subscript operator. The numbering of the
// vector elements is

                             \f[\left(\begin{array}{*{2}{c}}
                             0 & 1 \\
                             \end{array}\right)\f]
*/
template< typename Type >
class Vector2
{
   static_assert( std::is_arithmetic<Type>::value, "Vector2 only accepts arithmetic data types" );

private:
   //**Friend declarations*************************************************************************
   /*! \cond internal */
   template< typename Other > friend class Vector2;
   /*! \endcond */
   //*******************************************************************************************************************

public:
   //**Type definitions****************************************************************************
   using Length = typename SqrtTrait<Type>::Type;  //!< Vector length return type.
                                                   /*!< Return type of the Vector2<Type>::length
                                                        function. */
   using value_type = Type;
   //*******************************************************************************************************************

   //**Constructors*****************************************************************************************************
                              explicit inline Vector2() = default;
                              explicit inline Vector2( Type init );
   template< typename Other > explicit inline Vector2( Other init );
                              explicit inline Vector2( Type x, Type y );
                              explicit inline Vector2( const Type* init );
                                       inline Vector2( const Vector2& v ) = default;

   template< typename Other >
   inline Vector2( const Vector2<Other>& v );
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   // No explicitly declared destructor.
   //*******************************************************************************************************************

   //**Operators********************************************************************************************************
   /*!\name Operators */
   //@{
   inline Vector2&                              operator= ( const Vector2& v ) = default;
   template< typename Other > inline Vector2&   operator= ( const Vector2<Other>& v );
   template< typename Other > inline bool       operator==( Other rhs )                 const;
   template< typename Other > inline bool       operator==( const Vector2<Other>& rhs ) const;
   template< typename Other > inline bool       operator!=( Other rhs )                 const;
   template< typename Other > inline bool       operator!=( const Vector2<Other>& rhs ) const;
   inline Type&                                 operator[]( uint_t index );
   inline const Type&                           operator[]( uint_t index )                const;
   //@}
   //*******************************************************************************************************************

   //**Arithmetic operators************************************************************************
   /*!\name Arithmetic operators
   // \brief The return type of the arithmetic operators depends on the involved data types of
   // \brief the vectors. HIGH denotes the more significant data type of the arithmetic operation
   // \brief (for further detail see the MathTrait class description).
   */
   //@{
                              inline Vector2       operator-()                             const;
   template< typename Other > inline Vector2&      operator+=( const Vector2<Other>& rhs );
   template< typename Other > inline Vector2&      operator-=( const Vector2<Other>& rhs );
   template< typename Other > inline Vector2&      operator*=( Other rhs );
   template< typename Other > inline Vector2&      operator/=( Other rhs );
   template< typename Other > inline Vector2<HIGH> operator+ ( const Vector2<Other>& rhs ) const;
   template< typename Other > inline Vector2<HIGH> operator- ( const Vector2<Other>& rhs ) const;
   template< typename Other > inline Vector2<HIGH> operator* ( Other rhs )                 const;
   template< typename Other > inline HIGH          operator* ( const Vector2<Other>& rhs ) const;
   template< typename Other > inline Vector2<HIGH> operator/ ( Other rhs )                 const;
   //@}
   //*******************************************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline uint_t          indexOfMax( )                  const;
   inline uint_t          indexOfMin( )                  const;
   inline void            set( Type x, Type y );
   inline Length          length()                       const;
   inline Type            sqrLength()                    const;
   inline Vector2<Length> getNormalized()                const;
   inline Type*           data()                         {return v_;}
   inline Type const *    data()                         const {return v_;}
   //@}
   //*******************************************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   /**
    * The two statically allocated vector elements.
    *
    * Access to the vector values is gained via the subscript operator.
    * The order of the elements is
    * \f[\left(\begin{array}{*{2}{c}}
    * 0 & 1 \\
    * \end{array}\right)\f]
   **/
   Type v_[2] = {Type(), Type()};
   //@}
   //*******************************************************************************************************************
};
static_assert( std::is_trivially_copyable<Vector2<real_t>>::value, "Vector2<real_t> has to be trivially copyable!");
//**********************************************************************************************************************

template<typename T>
Vector2<T> & normalize( Vector2<T> & v );


//======================================================================================================================
//
//  CONSTRUCTORS
//
//======================================================================================================================


//**********************************************************************************************************************
/*!\fn Vector2<Type>::Vector2( Type init )
// \brief Constructor for a homogenous initialization of all elements.
//
// \param init Initial value for all vector elements.
*/
template< typename Type >
inline Vector2<Type>::Vector2( Type init )
{
   v_[0] = v_[1] = init;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Type>::Vector2( Type init )
// \brief Constructor for a homogenous initialization of all elements.
//
// \param init Initial value for all vector elements.
*/
template< typename Type >
template< typename Other >
inline Vector2<Type>::Vector2( Other init )
{
   static_assert( std::is_arithmetic<Other>::value, "Vector2 only accepts arithmetic data types in Vector2( Other init )");

   v_[0] = v_[1] = init;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Type>::Vector2( Type x, Type y )
// \brief Constructor for a direct initialization of all vector elements.
//
// \param x The initial value for the x-component.
// \param y The initial value for the y-component.
*/
template< typename Type >
inline Vector2<Type>::Vector2( Type x, Type y )
{
   v_[0] = x;
   v_[1] = y;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Type>::Vector2( const Type* init )
// \brief Constructor for an array initializer.
//
// \param init Pointer to the initialization array.
//
// The array is assumed to have at least three valid elements.
*/
template< typename Type >
inline Vector2<Type>::Vector2( const Type* init )
{
   v_[0] = init[0];
   v_[1] = init[1];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Type>::Vector2( const Vector2<Other>& v )
// \brief Conversion constructor from different Vector2 instances.
//
// \param v Vector to be copied.
*/
template< typename Type >
template< typename Other >
inline Vector2<Type>::Vector2( const Vector2<Other>& v )
{
   v_[0] = v.v_[0];
   v_[1] = v.v_[1];
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  OPERATORS
//
//======================================================================================================================


//**********************************************************************************************************************
/*!\fn Vector2<Type>& Vector2<Type>::operator=( const Vector2<Other>& v )
// \brief Assignment operator for different Vector2 instances.
//
// \param v Vector to be copied.
// \return Reference to the assigned vector.
*/
template< typename Type >
template< typename Other >
inline Vector2<Type>& Vector2<Type>::operator=( const Vector2<Other>& v )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = v.v_[0];
   v_[1] = v.v_[1];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Vector2<Type>::operator==( Other rhs ) const
// \brief Equality operator for the comparison of a vector and a scalar value.
//
// \param rhs The right-hand-side scalar value for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
template< typename Other >
inline bool Vector2<Type>::operator==( Other rhs ) const
{
   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   return equal( v_[0], rhs ) && equal( v_[1], rhs );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Vector2<Type>::operator==( const Vector2<Other>& rhs ) const
// \brief Equality operator for the comparison of two vectors.
//
// \param rhs The right-hand-side vector for the comparison.
// \return bool
*/
template< typename Type >
template< typename Other >
inline bool Vector2<Type>::operator==( const Vector2<Other>& rhs ) const
{
   // In order to compare the two vectors, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   return equal( v_[0], rhs.v_[0] ) && equal( v_[1], rhs.v_[1] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Vector2<Type>::operator!=( Other rhs ) const
// \brief Inequality operator for the comparison of a vector and a scalar value.
//
// \param rhs The right-hand-side scalar value for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
template< typename Other >
inline bool Vector2<Type>::operator!=( Other rhs ) const
{
   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   return !(*this == rhs);
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Vector2<Type>::operator!=( const Vector2<Other>& rhs ) const
// \brief Inequality operator for the comparison of two vectors.
//
// \param rhs The right-hand-side vector for the comparison.
// \return bool
*/
template< typename Type >
template< typename Other >
inline bool Vector2<Type>::operator!=( const Vector2<Other>& rhs ) const
{
   // In order to compare the two vectors, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   return !(*this == rhs);
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type& Vector2<Type>::operator[]( uint_t index )
// \brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..2]\f$.
// \return Reference to the accessed value.
*/
template< typename Type >
inline Type& Vector2<Type>::operator[]( uint_t index )
{
   WALBERLA_ASSERT_LESS( index, 2 , "Invalid vector access index" );
   return v_[index];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Type& Vector2<Type>::operator[]( uint_t index ) const
// \brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..2]\f$.
// \return Reference-to-const to the accessed value.
*/
template< typename Type >
inline const Type& Vector2<Type>::operator[]( uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, 2, "Invalid vector access index" );
   return v_[index];
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  ARITHMETIC OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Vector2<Type> Vector2<Type>::operator-() const
// \brief Unary minus operator for the inversion of a vector (\f$ \vec{a} = -\vec{b} \f$).
//
// \return The inverse of the vector.
*/
template< typename Type >
inline Vector2<Type> Vector2<Type>::operator-() const
{
   return Vector2( -v_[0], -v_[1] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Type>& Vector2<Type>::operator+=( const Vector2<Other>& rhs )
// \brief Addition assignment operator for the addition of two vectors (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand-side vector to be added to the vector.
// \return Reference to the vector.
*/
template< typename Type >
template< typename Other >
inline Vector2<Type>& Vector2<Type>::operator+=( const Vector2<Other>& rhs )
{
   v_[0] += rhs.v_[0];
   v_[1] += rhs.v_[1];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Type>& Vector2<Type>::operator-=( const Vector2<Other>& rhs )
// \brief Subtraction assignment operator for the subtraction of two vectors
// \brief (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand-side vector to be subtracted from the vector.
// \return Reference to the vector.
*/
template< typename Type >
template< typename Other >
inline Vector2<Type>& Vector2<Type>::operator-=( const Vector2<Other>& rhs )
{
   v_[0] -= rhs.v_[0];
   v_[1] -= rhs.v_[1];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Type>& Vector2<Type>::operator*=( Other rhs )
// \brief Multiplication assignment operator for the multiplication between a vector and
// \brief a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand-side scalar value for the multiplication.
// \return Reference to the vector.
*/
template< typename Type >
template< typename Other >
inline Vector2<Type>& Vector2<Type>::operator*=( Other rhs )
{
   v_[0] *= rhs;
   v_[1] *= rhs;
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Type>& Vector2<Type>::operator/=( Other rhs )
// \brief Division assignment operator for the division of a vector by a scalar value
// \brief (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand-side scalar value for the division.
// \return Reference to the vector.
//
// \b Note: No check for 0 is applied.
*/
template< typename Type >
template< typename Other >
inline Vector2<Type>& Vector2<Type>::operator/=( Other rhs )
{
   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( std::numeric_limits<HIGH>::is_integer ) {
      v_[0] /= rhs;
      v_[1] /= rhs;
      return *this;
   }
   else {
      const HIGH tmp( 1/static_cast<HIGH>( rhs ) );
      v_[0] = static_cast<Type>( static_cast<HIGH>( v_[0] ) * tmp );
      v_[1] = static_cast<Type>( static_cast<HIGH>( v_[1] ) * tmp );
      return *this;
   }
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<HIGH> Vector2<Type>::operator+( const Vector2<Other>& rhs ) const
// \brief Addition operator for the addition of two vectors (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
//
// \param rhs The right-hand-side vector to be added to the vector.
// \return The sum of the two vectors.
//
// The operator returns a vector of the higher-order data type of the two involved vector
// data types (in fact the two template arguments \a Type and \a Other ).
*/
template< typename Type >
template< typename Other >
inline Vector2<HIGH> Vector2<Type>::operator+( const Vector2<Other>& rhs ) const
{
   return Vector2<HIGH>( v_[0]+rhs.v_[0], v_[1]+rhs.v_[1] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<HIGH> Vector2<Type>::operator-( const Vector2<Other>& rhs ) const
// \brief Subtraction operator for the subtraction of two vectors (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
//
// \param rhs The right-hand-side vector to be subtracted from the vector.
// \return The difference of the two vectors.
//
// The operator returns a vector of the higher-order data type of the two involved vector
// data types (in fact the two template arguments \a Type and \a Other ).
*/
template< typename Type >
template< typename Other >
inline Vector2<HIGH> Vector2<Type>::operator-( const Vector2<Other>& rhs ) const
{
   return Vector2<HIGH>( v_[0]-rhs.v_[0], v_[1]-rhs.v_[1] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<HIGH> Vector2<Type>::operator*( Other rhs ) const
// \brief Multiplication operator for the multiplication of a vector and a scalar value
// \brief (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param rhs The right-hand-side scalar value for the multiplication.
// \return The scaled result vector.
//
// The operator returns a vector of the higher-order data type of the two involved data types
// (in fact the two template arguments \a Type and \a Other ).
*/
template< typename Type >
template< typename Other >
inline Vector2<HIGH> Vector2<Type>::operator*( Other rhs ) const
{
   return Vector2<HIGH>( v_[0]*rhs, v_[1]*rhs );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn HIGH Vector2<Type>::operator*( const Vector2<Other>& rhs ) const
// \brief Multiplication operator for the scalar product (inner product) of two vectors
// \brief (\f$ s=\vec{a}*\vec{b} \f$).
//
// \param rhs The right-hand-side vector for the inner product.
// \return The scalar product.
//
// The operator returns a scalar value of the higher-order data type of the two involved data
// types (in fact the two template arguments \a Type and \a Other ).
*/
template< typename Type >
template< typename Other >
inline HIGH Vector2<Type>::operator*( const Vector2<Other>& rhs ) const
{
   return ( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[1]  );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<HIGH> Vector2<Type>::operator/( Other rhs ) const
// \brief Division operator for the division of a vector by a scalar value
// \brief (\f$ \vec{a}=\vec{b}/s \f$).
//
// \param rhs The right-hand-side scalar value for the division.
// \return The scaled result vector.
//
// The operator returns a vector of the higher-order data type of the two involved data types
// (in fact the two template arguments \a Type and \a Other ).\n
// \b Note: No check for 0 is applied.
*/
template< typename Type >
template< typename Other >
inline Vector2<HIGH> Vector2<Type>::operator/( Other rhs ) const
{
   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( std::numeric_limits<HIGH>::is_integer ) {
      return Vector2<HIGH>( v_[0]/rhs, v_[1]/rhs );
   }
   else {
      const HIGH tmp( 1/static_cast<HIGH>( rhs ) );
      return Vector2<HIGH>( v_[0]*tmp, v_[1]*tmp );
   }
}
//**********************************************************************************************************************


//======================================================================================================================
//
//  UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Vector2<Type>::indexOfMax( )
// \brief Returns index of absolute maximum value
*/
template< typename Type >
inline uint_t Vector2<Type>::indexOfMax( ) const {
   return math::abs( v_[0] ) > math::abs( v_[1] ) ? 0u : 1u;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Vector2<Type>::indexOfMin( )
// \brief Returns index of absolute minimum value
*/
template< typename Type >
inline uint_t Vector2<Type>::indexOfMin( ) const {
   return math::abs( v_[0] ) < math::abs( v_[1] ) ? 0u : 1u;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Vector2<Type>::set( Type x, Type y )
// \brief Set function for a direct assignment of all vector elements.
//
// \param x The initial value for the x-component.
// \param y The initial value for the y-component.
*/
template< typename Type >
inline void Vector2<Type>::set( Type x, Type y )
{
   v_[0] = x;
   v_[1] = y;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Length Vector2<Type>::length() const
// \brief Calculation of the vector length \f$|\vec{a}|\f$.
//
// \return The length of the vector.
//
// The return type of the length function depends on the actual type of the vector instance:
//
// <table border="0" cellspacing="0" cellpadding="1">
//    <tr>
//       <td width="250px"> \b Type </td>
//       <td width="100px"> \b Length </td>
//    </tr>
//    <tr>
//       <td>float</td>
//       <td>float</td>
//    </tr>
//    <tr>
//       <td>integral data types and double</td>
//       <td>double</td>
//    </tr>
//    <tr>
//       <td>long double</td>
//       <td>long double</td>
//    </tr>
// </table>
*/
template< typename Type >
inline typename SqrtTrait<Type>::Type Vector2<Type>::length() const
{
   return std::sqrt( static_cast<typename SqrtTrait<Type>::Type>( v_[0]*v_[0] + v_[1]*v_[1] ) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type Vector2<Type>::sqrLength() const
// \brief Calculation of the vector square length \f$|\vec{a}|^2\f$.
//
// \return The square length of the vector.
*/
template< typename Type >
inline Type Vector2<Type>::sqrLength() const
{
   return ( v_[0]*v_[0] + v_[1]*v_[1] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<Length>& Vector2<Type>::getNormalized() const
// \brief Calculation of the normalized vector (\f$|\vec{a}|=1\f$).
//
// \pre \f$|\vec{a}|\neq0\f$
//
// \return The normalized vector.
//
// The function returns the normalized vector.
*/
template< typename Type >
Vector2<typename Vector2<Type>::Length> Vector2<Type>::getNormalized() const
{
   const Length len( length() );

   WALBERLA_ASSERT_FLOAT_UNEQUAL( len, Length(0) );

   Vector2<Length> result ( static_cast<Length>( v_[0] ) / len,
                            static_cast<Length>( v_[1] ) / len );

   WALBERLA_ASSERT_FLOAT_EQUAL( result.sqrLength(), 1.0 );

   return result;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  GLOBAL OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\name Vector2 operators */
//@{

// The following overloads of the comparison operators are necessary to disambiguate
// comparisons between a scalar value and a vector.
template< typename Type > inline bool operator==( unsigned char  scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( char           scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( signed char    scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( wchar_t        scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( unsigned short scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( short          scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( unsigned int   scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( int            scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( unsigned long  scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( long           scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( float          scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( double         scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator==( long double    scalar, const Vector2<Type>& vec );

template< typename Type > inline bool operator!=( unsigned char  scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( char           scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( signed char    scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( wchar_t        scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( unsigned short scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( short          scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( unsigned int   scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( int            scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( unsigned long  scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( long           scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( float          scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( double         scalar, const Vector2<Type>& vec );
template< typename Type > inline bool operator!=( long double    scalar, const Vector2<Type>& vec );

template< typename Type >
std::ostream& operator<<( std::ostream& os, const Vector2<Type>& v );

template< typename Type >
std::istream& operator>>( std::istream& is, Vector2<Type>& v );

template< typename Type >
inline bool isnan( const Vector2<Type>& v );

template< typename Type >
inline const Vector2<Type> abs( const Vector2<Type>& v );

template< typename Type >
inline const Vector2<Type> fabs( const Vector2<Type>& v );
//@}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( unsigned char scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of an unsigned char scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( unsigned char scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( char scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of a char scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( char scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( signed char scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of a signed char scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( signed char scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( wchar_t scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of a wchar_t scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( wchar_t scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( unsigned short scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of an unsigned short scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( unsigned short scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( short scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of a short scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( short scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( uint_t scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of an uint_t scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( unsigned int scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( int scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of an int scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( int scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( unsigned long scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of an unsigned long scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( unsigned long scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( long scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of a long scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( long scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( float scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of a float scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( float scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( double scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of a double scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( double scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( long double scalar, const Vector2<Type>& vec )
// \brief Equality operator for the comparison of a long double scalar and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If all values of the vector are equal to the scalar value, the equality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator==( long double scalar, const Vector2<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( unsigned char scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of an unsigned char scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( unsigned char scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( char scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of a char scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( char scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( signed char scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of a signed char scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( signed char scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( wchar_t scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of a wchar_t scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( wchar_t scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( unsigned short scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of an unsigned short scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( unsigned short scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( short scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of a short scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( short scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( uint_t scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of an uint_t scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( unsigned int scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( int scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of an int scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( int scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( unsigned long scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of an unsigned long scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( unsigned long scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( long scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of a long scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( long scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( float scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of a float scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( float scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( double scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of a double scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( double scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( long double scalar, const Vector2<Type>& vec )
// \brief Inequality operator for the comparison of a long double scalar value and a vector.
//
// \param scalar The left-hand-side scalar value for the comparison.
// \param vec The right-hand-side vector for the comparison.
// \return bool
//
// If one value of the vector is inequal to the scalar value, the inequality test returns true,
// otherwise false.
*/
template< typename Type >
inline bool operator!=( long double scalar, const Vector2<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector2<HIGH> operator*( Other scalar, const Vector2<Type>& vec )
// \brief Multiplication operator for the multiplication of a scalar value and a vector
// \brief (\f$ \vec{a}=s*\vec{b} \f$).
//
// \param scalar The left-hand-side scalar value for the multiplication.
// \param vec The right-hand-side vector for the multiplication.
// \return The scaled result vector.
*/
template< typename Type, typename Other >
inline typename std::enable_if< std::is_fundamental<Other>::value, Vector2<HIGH> >::type
   operator*( Other scalar, const Vector2<Type>& vec )
{
   return vec * scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::ostream& operator<<( std::ostream& os, const Vector2<Type>& v )
// \brief Global output operator for 2-dimensional vectors.
//
// \param os Reference to the output stream.
// \param v Reference to a constant vector object.
// \return Reference to the output stream.
*/
template< typename Type >
std::ostream& operator<<( std::ostream& os, const Vector2<Type>& v )
{
   return os << "<" << v[0] << "," << v[1] << ">";
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::istream& operator>>( std::istream& is, Vector2<Type>& v )
// \brief Global input operator for 2-dimensional vectors.
//
// \param is Reference to the input stream.
// \param v Reference to a vector object.
// \return The input stream.
*/
template< typename Type >
std::istream& operator>>( std::istream& is, Vector2<Type>& v )
{
   if( !is ) return is;

   char bracket1, bracket2, comma1;
   Type x(0), y(0);
   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );

   // Setting the 'skip whitespaces' flag
   is >> std::skipws;

   // Extracting the vector
   if( !(is >> bracket1 >> x >> comma1 >> y >> bracket2) ||
       bracket1 != '<' || comma1 != ',' || bracket2 != '>' ) {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      is.flags( oldFlags );
      return is;
   }

   // Transferring the input to the vector values
   v[0] = x; v[1] = y;

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::istream& operator>>( std::istream& is, Vector2<bool>& v )
// \brief Specialization for input operator for 2-dimensional vectors of bool.
//
// The operator can parse e.g. <1, 0, 1> and <true, false, true> from the istream
//
// \param is Reference to the input stream.
// \param v Reference to a bool vector object.
// \return The input stream.
*/
template<>
inline std::istream& operator>>( std::istream& is, Vector2<bool>& v )
{
   if( !is ) return is;

   char bracket1, bracket2, comma1;
   bool x(false), y(false);
   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );


   // try parsing e.g. <true, false, true>
   // Setting the 'skip whitespaces' flag
   is >> std::skipws >> std::boolalpha;

   // Extracting the vector
   if( !(is >> bracket1 >> x >> comma1 >> y >> bracket2) ||
       bracket1 != '<' || comma1 != ','  || bracket2 != '>' ) {
      is.clear();
      is.seekg( pos );
      is.flags( oldFlags );

      // try parsing e.g. <1, 0, 1>
      // Setting the 'skip whitespaces' flag
      is >> std::skipws >> std::noboolalpha;

      // Extracting the vector
      if( !(is >> bracket1 >> x >> comma1 >> y >> bracket2) ||
         bracket1 != '<' || comma1 != ','  || bracket2 != '>' ) {
            is.clear();
            is.seekg( pos );
            is.setstate( std::istream::failbit );
            is.flags( oldFlags );

            return is;
      }
   }

   // Transferring the input to the vector values
   v[0] = x; v[1] = y;

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool isnan( const Vector2<Type>& v )
// \brief Checks the given vector for not-a-number elements.
//
// \param v The vector to be checked for not-a-number elements.
// \return \a true if at least one element of the vector is not-a-number, \a false otherwise.
*/
template< typename Type >
inline bool isnan( const Vector2<Type>& v )
{
   return walberla::math::isnan( v[0] ) || walberla::math::isnan( v[1] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector2<Type> abs( const Vector2<Type>& v )
// \brief Returns a vector containing the absolute values of each single element of \a v.
//
// \param v The integral input vector.
// \return The absolute value of each single element of \a v.
//
// The \a abs function calculates the absolute value of each element of the input vector \a v.
// This function can only be applied to vectors of integral data type. For floating point vectors,
// the pe::fabs( const Vector2& ) function can be used.
*/
template< typename Type >
inline const Vector2<Type> abs( const Vector2<Type>& v )
{
   static_assert( std::numeric_limits<Type>::is_integer, "v has to be of integral type!" );
   return Vector2<Type>( std::abs(v[0]), std::abs(v[1]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector2<Type> fabs( const Vector2<Type>& v )
// \brief Returns a vector containing the absolute values of each single element of \a v.
//
// \param v The floating point input vector.
// \return The absolute value of each single element of \a v.
//
// The \a fabs function calculates the absolute value of each element of the input vector \a v.
// This function can only be applied to floating point vectors. For vectors of integral data
// type, the pe::abs( const Vector2& ) function can be used.
*/
template< typename Type >
inline const Vector2<Type> fabs( const Vector2<Type>& v )
{
   static_assert( !std::numeric_limits<Type>::is_integer, "v has to be of floating point type!" );
   return Vector2<Type>( std::fabs(v[0]), std::fabs(v[1]) );
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/**
// \brief Normalization of the vector (\f$|\vec{a}|=1\f$).
//
// \pre \f$|\vec{a}|\neq0\f$
//
// \return Reference to the normalized vector.
//
// Normalization of the vector to a length of 1. If you want to normalize a vector of integral T
// use member function getNormalized() instead
*/
template<typename T>
Vector2<T> & normalize( Vector2<T> & v )
{
   static_assert( std::is_floating_point<T>::value,
      "You can only normalize floating point vectors in-place!");
   static_assert( (std::is_same<T, typename Vector2<T>::Length>::value),
      "The type of your Vector2's length does not match its element type!" );

   const T len( v.length() );

   WALBERLA_ASSERT_FLOAT_UNEQUAL( len, typename Vector2<T>::Length(0) );

   v[0] /= len;
   v[1] /= len;

   WALBERLA_ASSERT_FLOAT_EQUAL( v.sqrLength(), T(1.0) );

   return v;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/**
// \brief Functor providing a lexicographical ordering for Vector2.
//
// \tparam T Datatype of the compared Vector2's elements.
*/
template< typename T >
struct Vector2LexicographicalyLess
{
   /**
   // \brief Provides a lexicographical less-than-operator for Vector2.
   //
   // \param lhs left hand side of less-than-operator.
   // \param rhs right hand side of less-than-operator.
   // \returns true if lhs < rhs (lexicographically), else returns false.
   */
   bool operator()( const Vector2<T> & lhs, const Vector2<T> & rhs ) const
   {
      return std::lexicographical_compare( &lhs[0], &lhs[0] + 2, &rhs[0], &rhs[0] + 2 );
   }
};
//**********************************************************************************************************************


//**********************************************************************************************************************
/**
// \brief Function providing a hash value for Vector2.
//
// \tparam  T Datatype of the Vector2's elements (only integers are supported).
// \param   v The vector the hash is computed for.
// \returns   A hash for the entire Vector2.
*/
template< typename T, typename Enable = std::enable_if_t<std::is_integral<T>::value> >
std::size_t hash_value( const Vector2<T> & v )
{
   std::size_t seed;

   if( sizeof(std::size_t) >= 8 )
   {
      seed = (static_cast<std::size_t>(v[0]) << 42) +
             (static_cast<std::size_t>(v[1]) << 21);
   }
   else
   {
      seed = (static_cast<std::size_t>(v[0]) << 21) +
             (static_cast<std::size_t>(v[1]) << 10);
   }

   return seed;
}
//**********************************************************************************************************************


} // namespace math

using math::Vector2;

} // namespace walberla

#undef HIGH


//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

   template< typename T,    // Element type of SendBuffer
             typename G,    // Growth policy of SendBuffer
             typename VT >  // Element type of vector
   mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const Vector2<VT> & vec )
   {
      buf.addDebugMarker( "v2" );
      static_assert ( std::is_trivially_copyable< Vector2<VT> >::value,
                      "type has to be trivially copyable for the memcpy to work correctly" );
      auto pos = buf.forward(sizeof(Vector2<VT>));
      std::memcpy(pos, &vec, sizeof(Vector2<VT>));
      return buf;
   }

   template< typename T,    // Element type  of RecvBuffer
             typename VT >  // Element type of vector
   mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, Vector2<VT> & vec )
   {
      buf.readDebugMarker( "v2" );
      static_assert ( std::is_trivially_copyable< Vector2<VT> >::value,
                      "type has to be trivially copyable for the memcpy to work correctly" );
      auto pos = buf.skip(sizeof(Vector2<VT>));
      //suppress https://gcc.gnu.org/onlinedocs/gcc/C_002b_002b-Dialect-Options.html#index-Wclass-memaccess
      std::memcpy(static_cast<void*>(&vec), pos, sizeof(Vector2<VT>));
      return buf;
   }

   template<typename VT>
   struct BufferSizeTrait< walberla::math::Vector2<VT> > {
      static const bool constantSize = true;
      static const uint_t size = 2 * BufferSizeTrait<VT>::size + mpi::BUFFER_DEBUG_OVERHEAD;
   };


}
}

//======================================================================================================================
//
//  MPI Datatype
//
//======================================================================================================================

namespace walberla {

   template< typename T>
   struct MPITrait< Vector2<T> >
   {
      static inline MPI_Datatype type()
      {
         // cannot use mpi::Datatype here because its destructor calls MPI_Type_free and static variables are destroyed after the MPI_Finalize
         static MPI_Datatype datatype;
         static bool initialized = false;

         if( ! initialized ) {
            MPI_Type_contiguous(2, MPITrait<T>::type(), &datatype );
            MPI_Type_commit( &datatype );
            initialized = true;
         }
         return datatype;
      }
   };
} // namespace walberla



//======================================================================================================================
//
//  Vector Trait Specialization
//
//======================================================================================================================

namespace walberla {

   // Specialization of VectorTrait for Vector2's
   template<typename T>
   struct VectorTrait< Vector2<T> >
   {
      using OutputType = T;

      static const uint_t F_SIZE =  2u;
      static   T  get( const Vector2<T> & v, uint_t f )       { return v[f]; }
      static void set(       Vector2<T> & v, uint_t f, T val) { v[f] = val;  }
   };

} // namespace walberla

namespace std
{
    template<typename T>
    struct hash< walberla::Vector2<T> >
    {
        std::size_t operator()( walberla::Vector2<T> const & v ) const noexcept
        {
            return walberla::math::hash_value( v );
        }
    };
} // namespace std
