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
//! \file Vector3.h
//! \ingroup core
//! \author Klaus Iglberger
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Header file for the implementation of a 3D vector
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
#include "core/debug/CheckFunctions.h"

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
/*!\brief Efficient, generic implementation of a 3-dimensional vector.
// \ingroup math
//
// The Vector3 class is the representation of a 3D vector with a total of 3 statically allocated
// elements of arbitrary type. The naming convention of the elements is as follows:

                             \f[\left(\begin{array}{*{3}{c}}
                             x & y & z \\
                             \end{array}\right)\f]

// These elements can be accessed directly with the subscript operator. The numbering of the
// vector elements is

                             \f[\left(\begin{array}{*{3}{c}}
                             0 & 1 & 2 \\
                             \end{array}\right)\f]
*/
template< typename Type >
class Vector3
{
   static_assert( std::is_arithmetic<Type>::value, "Vector3 only accepts arithmetic data types" );

private:
   //**Friend declarations*************************************************************************
   /*! \cond internal */
   template< typename Other > friend class Vector3;
   /*! \endcond */
   //*******************************************************************************************************************

public:
   //**Type definitions****************************************************************************
   using Length = typename SqrtTrait<Type>::Type;  //!< Vector length return type.
                                                   /*!< Return type of the Vector3<Type>::length
                                                        function. */
   using value_type = Type;
   //*******************************************************************************************************************

   //**Constructors*****************************************************************************************************
                              explicit inline constexpr Vector3() = default;
                              explicit inline constexpr Vector3( Type init );
   template< typename Other > explicit inline constexpr Vector3( Other init );
                              explicit inline constexpr Vector3( Type x, Type y, Type z );
                              explicit inline constexpr Vector3( const Type* init );
                                       inline constexpr Vector3( const Vector3& v ) = default;

   template< typename Other >
   inline constexpr Vector3( const Vector3<Other>& v );
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   // No explicitly declared destructor.
   //*******************************************************************************************************************

   //**Operators********************************************************************************************************
   /*!\name Operators */
   //@{
   inline Vector3&                              operator= ( const Vector3& v ) = default;
   template< typename Other > inline Vector3&   operator= ( const Vector3<Other>& v );
   template< typename Other > inline bool       operator==( Other rhs )                 const;
   template< typename Other > inline bool       operator==( const Vector3<Other>& rhs ) const;
   template< typename Other > inline bool       operator!=( Other rhs )                 const;
   template< typename Other > inline bool       operator!=( const Vector3<Other>& rhs ) const;
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
                              inline Vector3       operator-()                             const;
   template< typename Other > inline Vector3&      operator%=( const Vector3<Other>& rhs );       //cross product
   template< typename Other > inline Vector3&      operator+=( const Vector3<Other>& rhs );
   template< typename Other > inline Vector3&      operator-=( const Vector3<Other>& rhs );
   template< typename Other > inline Vector3&      operator*=( Other rhs );
   template< typename Other > inline Vector3&      operator/=( Other rhs );
   template< typename Other > inline Vector3<HIGH> operator% ( const Vector3<Other>& rhs ) const; //cross product
   template< typename Other > inline Vector3<HIGH> operator+ ( const Vector3<Other>& rhs ) const;
   template< typename Other > inline Vector3<HIGH> operator- ( const Vector3<Other>& rhs ) const;
   template< typename Other > inline Vector3<HIGH> operator* ( Other rhs )                 const;
   template< typename Other > inline HIGH          operator* ( const Vector3<Other>& rhs ) const;
   template< typename Other > inline Vector3<HIGH> operator/ ( Other rhs )                 const;
   //@}
   //*******************************************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline uint_t          indexOfMax( )                  const;
   inline uint_t          indexOfMin( )                  const;
   inline Type            max( )                         const;
   inline Type            min( )                         const;
   inline void            set( Type x, Type y, Type z );
   inline Length          length()                       const;
   inline Type            sqrLength()                    const;
   inline Vector3<Length> getNormalized()                const;
   inline Vector3<Length> getNormalizedOrZero()          const;
   inline void            reset();
   inline Type*           data()                         {return v_;}
   inline Type const *    data()                         const {return v_;}
   //@}
   //*******************************************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   /**
    * The three statically allocated vector elements.
    *
    * Access to the vector values is gained via the subscript operator.
    * The order of the elements is
    * \f[\left(\begin{array}{*{3}{c}}
    * 0 & 1 & 2 \\
    * \end{array}\right)\f]
   **/
   Type v_[3] = {Type(), Type(), Type()};
   //@}
   //*******************************************************************************************************************
};
static_assert( std::is_trivially_copyable<Vector3<real_t>>::value, "Vector3<real_t> has to be trivially copyable!");
//**********************************************************************************************************************

template<typename T>
Vector3<T> & normalize( Vector3<T> & v );


//======================================================================================================================
//
//  CONSTRUCTORS
//
//======================================================================================================================


//**********************************************************************************************************************
/*!\fn Vector3<Type>::Vector3( Type init )
// \brief Constructor for a homogenous initialization of all elements.
//
// \param init Initial value for all vector elements.
*/
template< typename Type >
inline constexpr Vector3<Type>::Vector3( Type init )
{
   v_[0] = v_[1] = v_[2] = init;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type>::Vector3( Type init )
// \brief Constructor for a homogenous initialization of all elements.
//
// \param init Initial value for all vector elements.
*/
template< typename Type >
template< typename Other >
inline constexpr Vector3<Type>::Vector3( Other init )
{
   static_assert( std::is_arithmetic<Other>::value, "Vector3 only accepts arithmetic data types in Vector3( Other init )");

   v_[0] = v_[1] = v_[2] = numeric_cast<Type>(init);
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type>::Vector3( Type x, Type y, Type z )
// \brief Constructor for a direct initialization of all vector elements.
//
// \param x The initial value for the x-component.
// \param y The initial value for the y-component.
// \param z The initial value for the z-component.
*/
template< typename Type >
inline constexpr Vector3<Type>::Vector3( Type x, Type y, Type z )
{
   v_[0] = x;
   v_[1] = y;
   v_[2] = z;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type>::Vector3( const Type* init )
// \brief Constructor for an array initializer.
//
// \param init Pointer to the initialization array.
//
// The array is assumed to have at least three valid elements.
*/
template< typename Type >
inline constexpr Vector3<Type>::Vector3( const Type* init )
{
   v_[0] = init[0];
   v_[1] = init[1];
   v_[2] = init[2];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type>::Vector3( const Vector3<Other>& v )
// \brief Conversion constructor from different Vector3 instances.
//
// \param v Vector to be copied.
*/
template< typename Type >
template< typename Other >
inline constexpr Vector3<Type>::Vector3( const Vector3<Other>& v )
{
   v_[0] = numeric_cast<Type>( v.v_[0] );
   v_[1] = numeric_cast<Type>( v.v_[1] );
   v_[2] = numeric_cast<Type>( v.v_[2] );
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  OPERATORS
//
//======================================================================================================================


//**********************************************************************************************************************
/*!\fn Vector3<Type>& Vector3<Type>::operator=( const Vector3<Other>& v )
// \brief Assignment operator for different Vector3 instances.
//
// \param v Vector to be copied.
// \return Reference to the assigned vector.
*/
template< typename Type >
template< typename Other >
inline Vector3<Type>& Vector3<Type>::operator=( const Vector3<Other>& v )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = v.v_[0];
   v_[1] = v.v_[1];
   v_[2] = v.v_[2];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Vector3<Type>::operator==( Other rhs ) const
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
inline bool Vector3<Type>::operator==( Other rhs ) const
{
   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   return equal( v_[0], rhs ) && equal( v_[1], rhs ) && equal( v_[2], rhs );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Vector3<Type>::operator==( const Vector3<Other>& rhs ) const
// \brief Equality operator for the comparison of two vectors.
//
// \param rhs The right-hand-side vector for the comparison.
// \return bool
*/
template< typename Type >
template< typename Other >
inline bool Vector3<Type>::operator==( const Vector3<Other>& rhs ) const
{
   // In order to compare the two vectors, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   return equal( v_[0], rhs.v_[0] ) && equal( v_[1], rhs.v_[1] ) && equal( v_[2], rhs.v_[2] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Vector3<Type>::operator!=( Other rhs ) const
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
inline bool Vector3<Type>::operator!=( Other rhs ) const
{
   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   return !(*this == rhs);
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Vector3<Type>::operator!=( const Vector3<Other>& rhs ) const
// \brief Inequality operator for the comparison of two vectors.
//
// \param rhs The right-hand-side vector for the comparison.
// \return bool
*/
template< typename Type >
template< typename Other >
inline bool Vector3<Type>::operator!=( const Vector3<Other>& rhs ) const
{
   // In order to compare the two vectors, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   return !(*this == rhs);
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type& Vector3<Type>::operator[]( uint_t index )
// \brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..2]\f$.
// \return Reference to the accessed value.
*/
template< typename Type >
inline Type& Vector3<Type>::operator[]( uint_t index )
{
   WALBERLA_ASSERT_LESS( index, 3 , "Invalid vector access index" );
   return v_[index];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Type& Vector3<Type>::operator[]( uint_t index ) const
// \brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..2]\f$.
// \return Reference-to-const to the accessed value.
*/
template< typename Type >
inline const Type& Vector3<Type>::operator[]( uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, 3, "Invalid vector access index" );
   return v_[index];
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  ARITHMETIC OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Vector3<Type> Vector3<Type>::operator-() const
// \brief Unary minus operator for the inversion of a vector (\f$ \vec{a} = -\vec{b} \f$).
//
// \return The inverse of the vector.
*/
template< typename Type >
inline Vector3<Type> Vector3<Type>::operator-() const
{
   return Vector3( -v_[0], -v_[1], -v_[2] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type> & Vector3<Type>::operator%=( const Vector3<Other>& rhs )
// \brief Cross product (outer product) of two vectors (\f$ \vec{a}=\vec{b}\times\vec{c} \f$).
//
// \param rhs The right-hand-side vector for the cross product.
// \return The cross product.
*/
template< typename Type >
template< typename Other >
inline Vector3<Type>& Vector3<Type>::operator%=( const Vector3<Other>& rhs )
{
   Type tmp0 = v_[1] * rhs.v_[2] - v_[2] * rhs.v_[1];
   Type tmp1 = v_[2] * rhs.v_[0] - v_[0] * rhs.v_[2];
   v_[2]     = v_[0] * rhs.v_[1] - v_[1] * rhs.v_[0];
   v_[1]     = tmp1;
   v_[0]     = tmp0;
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type>& Vector3<Type>::operator+=( const Vector3<Other>& rhs )
// \brief Addition assignment operator for the addition of two vectors (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand-side vector to be added to the vector.
// \return Reference to the vector.
*/
template< typename Type >
template< typename Other >
inline Vector3<Type>& Vector3<Type>::operator+=( const Vector3<Other>& rhs )
{
   v_[0] += numeric_cast<Type>(rhs.v_[0]);
   v_[1] += numeric_cast<Type>(rhs.v_[1]);
   v_[2] += numeric_cast<Type>(rhs.v_[2]);
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type>& Vector3<Type>::operator-=( const Vector3<Other>& rhs )
// \brief Subtraction assignment operator for the subtraction of two vectors
// \brief (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand-side vector to be subtracted from the vector.
// \return Reference to the vector.
*/
template< typename Type >
template< typename Other >
inline Vector3<Type>& Vector3<Type>::operator-=( const Vector3<Other>& rhs )
{
   v_[0] -= numeric_cast<Type>(rhs.v_[0]);
   v_[1] -= numeric_cast<Type>(rhs.v_[1]);
   v_[2] -= numeric_cast<Type>(rhs.v_[2]);
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type>& Vector3<Type>::operator*=( Other rhs )
// \brief Multiplication assignment operator for the multiplication between a vector and
// \brief a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand-side scalar value for the multiplication.
// \return Reference to the vector.
*/
template< typename Type >
template< typename Other >
inline Vector3<Type>& Vector3<Type>::operator*=( Other rhs )
{
   v_[0] *= numeric_cast<Type>(rhs);
   v_[1] *= numeric_cast<Type>(rhs);
   v_[2] *= numeric_cast<Type>(rhs);
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Type>& Vector3<Type>::operator/=( Other rhs )
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
inline Vector3<Type>& Vector3<Type>::operator/=( Other rhs )
{
   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( std::numeric_limits<HIGH>::is_integer ) {
      v_[0] /= numeric_cast<Type>(rhs);
      v_[1] /= numeric_cast<Type>(rhs);
      v_[2] /= numeric_cast<Type>(rhs);
      return *this;
   }
   else {
      const HIGH tmp( 1/static_cast<HIGH>( rhs ) );
      v_[0] = static_cast<Type>( static_cast<HIGH>( v_[0] ) * tmp );
      v_[1] = static_cast<Type>( static_cast<HIGH>( v_[1] ) * tmp );
      v_[2] = static_cast<Type>( static_cast<HIGH>( v_[2] ) * tmp );
      return *this;
   }
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<HIGH> Vector3<Type>::operator%( const Vector3<Other>& rhs ) const
// \brief Cross product (outer product) of two vectors (\f$ \vec{a}=\vec{b}\times\vec{c} \f$).
//
// \param rhs The right-hand-side vector for the cross product.
// \return The cross product.
*/
template< typename Type >
template< typename Other >
inline Vector3<HIGH> Vector3<Type>::operator%( const Vector3<Other>& rhs ) const
{
   return Vector3<HIGH>( v_[1] * rhs.v_[2] - v_[2] * rhs.v_[1],
                         v_[2] * rhs.v_[0] - v_[0] * rhs.v_[2],
                         v_[0] * rhs.v_[1] - v_[1] * rhs.v_[0] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Cross product (outer product) of two vectors (\f$ \vec{a}=\vec{b}\times\vec{c} \f$).
//
// \param lhs The left-hand-side vector for the cross product.
// \param rhs The right-hand-side vector for the cross product.
// \return The cross product.
*/
template< typename Type >
inline Vector3<Type> cross( const Vector3<Type>& lhs, const Vector3<Type>& rhs )
{
   return lhs % rhs;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<HIGH> Vector3<Type>::operator+( const Vector3<Other>& rhs ) const
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
inline Vector3<HIGH> Vector3<Type>::operator+( const Vector3<Other>& rhs ) const
{
   return Vector3<HIGH>( v_[0]+numeric_cast<Type>(rhs.v_[0]), v_[1]+numeric_cast<Type>(rhs.v_[1]), v_[2]+numeric_cast<Type>(rhs.v_[2]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<HIGH> Vector3<Type>::operator-( const Vector3<Other>& rhs ) const
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
inline Vector3<HIGH> Vector3<Type>::operator-( const Vector3<Other>& rhs ) const
{
   return Vector3<HIGH>( v_[0]-numeric_cast<Type>(rhs.v_[0]), v_[1]-numeric_cast<Type>(rhs.v_[1]), v_[2]-numeric_cast<Type>(rhs.v_[2]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<HIGH> Vector3<Type>::operator*( Other rhs ) const
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
inline Vector3<HIGH> Vector3<Type>::operator*( Other rhs ) const
{
   return Vector3<HIGH>( v_[0]*numeric_cast<Type>(rhs), v_[1]*numeric_cast<Type>(rhs), v_[2]*numeric_cast<Type>(rhs) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn HIGH Vector3<Type>::operator*( const Vector3<Other>& rhs ) const
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
inline HIGH Vector3<Type>::operator*( const Vector3<Other>& rhs ) const
{
   return ( v_[0]*numeric_cast<Type>(rhs.v_[0]) + v_[1]*numeric_cast<Type>(rhs.v_[1]) + v_[2]*numeric_cast<Type>(rhs.v_[2]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<HIGH> Vector3<Type>::operator/( Other rhs ) const
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
inline Vector3<HIGH> Vector3<Type>::operator/( Other rhs ) const
{
   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( std::numeric_limits<HIGH>::is_integer ) {
      return Vector3<HIGH>( v_[0]/rhs, v_[1]/rhs, v_[2]/rhs );
   }
   else {
      const HIGH tmp( 1/static_cast<HIGH>( rhs ) );
      return Vector3<HIGH>( v_[0]*tmp, v_[1]*tmp, v_[2]*tmp );
   }
}
//**********************************************************************************************************************

template< typename Type, typename Other >
inline Vector3<HIGH> operator/( Other lhs, const Vector3<Type>& rhs )
{
   return Vector3<HIGH>( lhs/rhs[0], lhs/rhs[1], lhs/rhs[2] );
}
//**********************************************************************************************************************


//======================================================================================================================
//
//  UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Vector3<Type>::indexOfMax( )
// \brief Returns index of absolute maximum value
*/
template< typename Type >
inline uint_t Vector3<Type>::indexOfMax( ) const {
   if(math::abs(v_[1]) > math::abs(v_[2]))
      return (math::abs(v_[0]) > math::abs(v_[1])) ? 0u : 1u;
   else
      return (math::abs(v_[0]) > math::abs(v_[2])) ? 0u : 2u;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Vector3<Type>::indexOfMin( )
// \brief Returns index of absolute minimum value
*/
template< typename Type >
inline uint_t Vector3<Type>::indexOfMin( ) const {
   if(math::abs(v_[2]) < math::abs(v_[1]))
      return (math::abs(v_[2]) < math::abs(v_[0])) ? 2u : 0u;
   else
      return (math::abs(v_[1]) < math::abs(v_[0])) ? 1u : 0u;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Vector3<Type>::max( )
// \brief Returns maximum value
*/
template< typename Type >
inline Type Vector3<Type>::max( ) const {
   return std::max(v_[0], std::max(v_[1], v_[2]));
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Vector3<Type>::min( )
// \brief Returns minimum value
*/
template< typename Type >
inline Type Vector3<Type>::min( ) const {
   return std::min(v_[0], std::min(v_[1], v_[2]));
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Vector3<Type>::set( Type x, Type y, Type z )
// \brief Set function for a direct assignment of all vector elements.
//
// \param x The initial value for the x-component.
// \param y The initial value for the y-component.
// \param z The initial value for the z-component.
*/
template< typename Type >
inline void Vector3<Type>::set( Type x, Type y, Type z )
{
   v_[0] = x;
   v_[1] = y;
   v_[2] = z;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Length Vector3<Type>::length() const
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
inline typename SqrtTrait<Type>::Type Vector3<Type>::length() const
{
   return std::sqrt( static_cast<typename SqrtTrait<Type>::Type>( v_[0]*v_[0] + v_[1]*v_[1] + v_[2]*v_[2] ) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type Vector3<Type>::sqrLength() const
// \brief Calculation of the vector square length \f$|\vec{a}|^2\f$.
//
// \return The square length of the vector.
*/
template< typename Type >
inline Type Vector3<Type>::sqrLength() const
{
   return ( v_[0]*v_[0] + v_[1]*v_[1] + v_[2]*v_[2] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<Length>& Vector3<Type>::getNormalized() const
// \brief Calculation of the normalized vector (\f$|\vec{a}|=1\f$).
//
// \pre \f$|\vec{a}|\neq0\f$
//
// \return The normalized vector.
//
// The function returns the normalized vector.
*/
template< typename Type >
Vector3<typename Vector3<Type>::Length> Vector3<Type>::getNormalized() const
{
   const Length len( length() );

   WALBERLA_ASSERT_FLOAT_UNEQUAL( len, Length(0) );

   const Length ilen( Length(1) / len );

   Vector3<Length> result ( static_cast<Length>( v_[0] ) * ilen,
                            static_cast<Length>( v_[1] ) * ilen,
                            static_cast<Length>( v_[2] ) * ilen );

   WALBERLA_ASSERT_FLOAT_EQUAL( result.sqrLength(), 1.0 );

   return result;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Vector3<Length>& Vector3<Type>::getNormalized() const
// \brief Calculation of the normalized vector (\f$|\vec{a}|=1\f$) without precondition.
//
// \return The normalized vector or the original vector if vector is too small.
*/
template< typename Type >
Vector3<typename Vector3<Type>::Length> Vector3<Type>::getNormalizedOrZero() const
{
   const Length len( length() );

   if (floatIsEqual( len, Length(0) )) return *this;

   const Length ilen( Length(1) / len );

   Vector3<Length> result ( static_cast<Length>( v_[0] ) * ilen,
                            static_cast<Length>( v_[1] ) * ilen,
                            static_cast<Length>( v_[2] ) * ilen );

   WALBERLA_ASSERT_FLOAT_EQUAL( result.sqrLength(), 1.0, "initial vector: " << result );

   return result;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn void Vector3<Type>::reset()
// \brief Sets all components of the vector to 0.
*/
template< typename Type >
void Vector3<Type>::reset()
{
   v_[0] = 0;
   v_[1] = 0;
   v_[2] = 0;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  GLOBAL OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\name Vector3 operators */
//@{

// The following overloads of the comparison operators are necessary to disambiguate
// comparisons between a scalar value and a vector.
template< typename Type > inline bool operator==( unsigned char  scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( char           scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( signed char    scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( wchar_t        scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( unsigned short scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( short          scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( unsigned int   scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( int            scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( unsigned long  scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( long           scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( float          scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( double         scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator==( long double    scalar, const Vector3<Type>& vec );

template< typename Type > inline bool operator!=( unsigned char  scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( char           scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( signed char    scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( wchar_t        scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( unsigned short scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( short          scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( unsigned int   scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( int            scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( unsigned long  scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( long           scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( float          scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( double         scalar, const Vector3<Type>& vec );
template< typename Type > inline bool operator!=( long double    scalar, const Vector3<Type>& vec );

template< typename Type >
std::ostream& operator<<( std::ostream& os, const Vector3<Type>& v );

template< typename Type >
std::istream& operator>>( std::istream& is, Vector3<Type>& v );

template< typename Type >
inline bool isnan( const Vector3<Type>& v );

template< typename Type >
inline bool isinf( const Vector3<Type>& v );

template< typename Type >
inline bool finite( const Vector3<Type>& v );

template< typename Type >
inline const Vector3<Type> abs( const Vector3<Type>& v );

template< typename Type >
inline const Vector3<Type> fabs( const Vector3<Type>& v );
//@}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( unsigned char scalar, const Vector3<Type>& vec )
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
inline bool operator==( unsigned char scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( char scalar, const Vector3<Type>& vec )
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
inline bool operator==( char scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( signed char scalar, const Vector3<Type>& vec )
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
inline bool operator==( signed char scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( wchar_t scalar, const Vector3<Type>& vec )
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
inline bool operator==( wchar_t scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( unsigned short scalar, const Vector3<Type>& vec )
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
inline bool operator==( unsigned short scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( short scalar, const Vector3<Type>& vec )
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
inline bool operator==( short scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( uint_t scalar, const Vector3<Type>& vec )
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
inline bool operator==( unsigned int scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( int scalar, const Vector3<Type>& vec )
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
inline bool operator==( int scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( unsigned long scalar, const Vector3<Type>& vec )
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
inline bool operator==( unsigned long scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( long scalar, const Vector3<Type>& vec )
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
inline bool operator==( long scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( float scalar, const Vector3<Type>& vec )
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
inline bool operator==( float scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( double scalar, const Vector3<Type>& vec )
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
inline bool operator==( double scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator==( long double scalar, const Vector3<Type>& vec )
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
inline bool operator==( long double scalar, const Vector3<Type>& vec )
{
   return vec == scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( unsigned char scalar, const Vector3<Type>& vec )
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
inline bool operator!=( unsigned char scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( char scalar, const Vector3<Type>& vec )
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
inline bool operator!=( char scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( signed char scalar, const Vector3<Type>& vec )
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
inline bool operator!=( signed char scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( wchar_t scalar, const Vector3<Type>& vec )
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
inline bool operator!=( wchar_t scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( unsigned short scalar, const Vector3<Type>& vec )
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
inline bool operator!=( unsigned short scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( short scalar, const Vector3<Type>& vec )
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
inline bool operator!=( short scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( uint_t scalar, const Vector3<Type>& vec )
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
inline bool operator!=( unsigned int scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( int scalar, const Vector3<Type>& vec )
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
inline bool operator!=( int scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( unsigned long scalar, const Vector3<Type>& vec )
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
inline bool operator!=( unsigned long scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( long scalar, const Vector3<Type>& vec )
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
inline bool operator!=( long scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( float scalar, const Vector3<Type>& vec )
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
inline bool operator!=( float scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( double scalar, const Vector3<Type>& vec )
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
inline bool operator!=( double scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool operator!=( long double scalar, const Vector3<Type>& vec )
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
inline bool operator!=( long double scalar, const Vector3<Type>& vec )
{
   return vec != scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Vector3<HIGH> operator*( Other scalar, const Vector3<Type>& vec )
// \brief Multiplication operator for the multiplication of a scalar value and a vector
// \brief (\f$ \vec{a}=s*\vec{b} \f$).
//
// \param scalar The left-hand-side scalar value for the multiplication.
// \param vec The right-hand-side vector for the multiplication.
// \return The scaled result vector.
*/
template< typename Type, typename Other >
inline typename std::enable_if< std::is_fundamental<Other>::value, Vector3<HIGH> >::type
   operator*( Other scalar, const Vector3<Type>& vec )
{
   return vec * scalar;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::ostream& operator<<( std::ostream& os, const Vector3<Type>& v )
// \brief Global output operator for 3-dimensional vectors.
//
// \param os Reference to the output stream.
// \param v Reference to a constant vector object.
// \return Reference to the output stream.
*/
template< typename Type >
std::ostream& operator<<( std::ostream& os, const Vector3<Type>& v )
{
   return os << "<" << v[0] << "," << v[1] << "," << v[2] << ">";
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::istream& operator>>( std::istream& is, Vector3<Type>& v )
// \brief Global input operator for 3-dimensional vectors.
//
// \param is Reference to the input stream.
// \param v Reference to a vector object.
// \return The input stream.
*/
template< typename Type >
std::istream& operator>>( std::istream& is, Vector3<Type>& v )
{
   if( !is ) return is;

   char bracket1, bracket2, comma1, comma2;
   Type x(0), y(0), z(0);
   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );

   // Setting the 'skip whitespaces' flag
   is >> std::skipws;

   // Extracting the vector
   if( !(is >> bracket1 >> x >> comma1 >> y >> comma2 >> z >> bracket2) ||
       bracket1 != '<' || comma1 != ',' || comma2 != ',' || bracket2 != '>' ) {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      is.flags( oldFlags );
      return is;
   }

   // Transferring the input to the vector values
   v[0] = x; v[1] = y; v[2] = z;

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::istream& operator>>( std::istream& is, Vector3<bool>& v )
// \brief Specialization for input operator for 3-dimensional vectors of bool.
//
// The operator can parse e.g. <1, 0, 1> and <true, false, true> from the istream
//
// \param is Reference to the input stream.
// \param v Reference to a bool vector object.
// \return The input stream.
*/
template<>
inline std::istream& operator>>( std::istream& is, Vector3<bool>& v )
{
   if( !is ) return is;

   char bracket1, bracket2, comma1, comma2;
   bool x(false), y(false), z(false);
   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );


   // try parsing e.g. <true, false, true>
   // Setting the 'skip whitespaces' flag
   is >> std::skipws >> std::boolalpha;

   // Extracting the vector
   if( !(is >> bracket1 >> x >> comma1 >> y >> comma2 >> z >> bracket2) ||
       bracket1 != '<' || comma1 != ',' || comma2 != ',' || bracket2 != '>' ) {
      is.clear();
      is.seekg( pos );
      is.flags( oldFlags );

      // try parsing e.g. <1, 0, 1>
      // Setting the 'skip whitespaces' flag
      is >> std::skipws >> std::noboolalpha;

      // Extracting the vector
      if( !(is >> bracket1 >> x >> comma1 >> y >> comma2 >> z >> bracket2) ||
         bracket1 != '<' || comma1 != ',' || comma2 != ',' || bracket2 != '>' ) {
            is.clear();
            is.seekg( pos );
            is.setstate( std::istream::failbit );
            is.flags( oldFlags );

            return is;
      }
   }

   // Transferring the input to the vector values
   v[0] = x; v[1] = y; v[2] = z;

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool isnan( const Vector3<Type>& v )
// \brief Checks the given vector for not-a-number elements.
//
// \param v The vector to be checked for not-a-number elements.
// \return \a true if at least one element of the vector is not-a-number, \a false otherwise.
*/
template< typename Type >
inline bool isnan( const Vector3<Type>& v )
{
   return walberla::math::isnan( v[0] ) || walberla::math::isnan( v[1] ) || walberla::math::isnan( v[2] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool isinf( const Vector3<Type>& v )
// \brief Checks the given vector for infinite elements.
//
// \param v The vector to be checked for infinite elements.
// \return \a true if at least one element of the vector is infinite, \a false otherwise.
*/
template< typename Type >
inline bool isinf( const Vector3<Type>& v )
{
   return walberla::math::isinf( v[0] ) || walberla::math::isinf( v[1] ) || walberla::math::isinf( v[2] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool finite( const Vector3<Type>& v )
// \brief Checks if the given vector has only finite elements.
//
// \param v The vector to be checked for only having finite elements.
// \return \a true if all elements of the vector are finite, \a false otherwise.
*/
template< typename Type >
inline bool finite( const Vector3<Type>& v )
{
   return walberla::math::finite( v[0] ) && walberla::math::finite( v[1] ) && walberla::math::finite( v[2] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector3<Type> abs( const Vector3<Type>& v )
// \brief Returns a vector containing the absolute values of each single element of \a v.
//
// \param v The integral input vector.
// \return The absolute value of each single element of \a v.
//
// The \a abs function calculates the absolute value of each element of the input vector \a v.
// This function can only be applied to vectors of integral data type. For floating point vectors,
// the pe::fabs( const Vector3& ) function can be used.
*/
template< typename Type >
inline const Vector3<Type> abs( const Vector3<Type>& v )
{
   static_assert( std::numeric_limits<Type>::is_integer, "v has to be of integral type!" );
   return Vector3<Type>( std::abs(v[0]), std::abs(v[1]), std::abs(v[2]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector3<Type> fabs( const Vector3<Type>& v )
// \brief Returns a vector containing the absolute values of each single element of \a v.
//
// \param v The floating point input vector.
// \return The absolute value of each single element of \a v.
//
// The \a fabs function calculates the absolute value of each element of the input vector \a v.
// This function can only be applied to floating point vectors. For vectors of integral data
// type, the pe::abs( const Vector3& ) function can be used.
*/
template< typename Type >
inline const Vector3<Type> fabs( const Vector3<Type>& v )
{
   static_assert( !std::numeric_limits<Type>::is_integer, "v has to be of floating point type!" );
   return Vector3<Type>( std::fabs(v[0]), std::fabs(v[1]), std::fabs(v[2]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void normals(Vector3<Type>& v, Vector3<Type>& defNor, Vector3<Type>& comNor)
// \brief Computes two normal values to a vector.
//
// The first normal is defined by switching the absolute highest and absolute lowest value of the
// vector and negating the higher value. The second normal is computed by the cross product of
// the input vector and the defined normal.
//
// \param v The vector for the normal computation.
// \param defNor The defined normal from the vector.
// \param comNor The computed normal from the vector.
*/
template< typename Type >
inline void normals(const Vector3<Type>& v, Vector3<Type>& defNor, Vector3<Type>& comNor)
{
   uint_t iMax = v.indexOfMax( );
   defNor.set( Type(0), Type(0), Type(0) );

   if ( equal( v[iMax], Type(0) ) )
   {
      comNor.set( Type(0), Type(0), Type(0) );
   }
   else
   {
      uint_t iMin = v.indexOfMin( );

      Type f = Type(1) / std::sqrt(v[iMax]*v[iMax] + v[iMin]*v[iMin]);
      defNor[iMin] = -v[iMax] * f;
      defNor[iMax] =  v[iMin] * f;

      comNor = v;
      comNor %= defNor;
      normalize(comNor);

      WALBERLA_ASSERT_FLOAT_EQUAL( v * defNor, 0.0 );
      WALBERLA_ASSERT_FLOAT_EQUAL( v * comNor, 0.0 );
      WALBERLA_ASSERT_FLOAT_EQUAL( defNor.sqrLength(), typename Vector3<Type>::Length(1.0) );
      WALBERLA_ASSERT_FLOAT_EQUAL( comNor.sqrLength(), typename Vector3<Type>::Length(1.0) );
   }
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
Vector3<T> & normalize( Vector3<T> & v )
{
   static_assert( std::is_floating_point<T>::value,
      "You can only normalize floating point vectors in-place!");
   static_assert( (std::is_same<T, typename Vector3<T>::Length>::value),
      "The type of your Vector3's length does not match its element type!" );

   const T len( v.length() );

   const T ilen( T(1) / len );

   v[0] *= ilen;
   v[1] *= ilen;
   v[2] *= ilen;

   WALBERLA_ASSERT_FLOAT_EQUAL( v.sqrLength(), T(1.0), "Problem after normalization of vector " << v <<
                                                       " (length prior to normalization was: " << len << ")" );

   return v;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/**
// \brief Length of the vector.
//
// \return Length of the vector.
*/
template<typename T>
real_t length( const Vector3<T> & v )
{
   return v.length();
}
//**********************************************************************************************************************
//**********************************************************************************************************************
/**
// \brief Length of the vector squared.
//
// \return Length of the vector squared.
*/
template<typename T>
real_t sqrLength( const Vector3<T> & v )
{
   return v.sqrLength();
}
//**********************************************************************************************************************
//**********************************************************************************************************************
/**
// \brief Dot product of two vectors.
*/
template<typename T>
real_t dot( const Vector3<T> & v1, const Vector3<T> & v2 )
{
   return v1*v2;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/**
// \brief Functor providing a lexicographical ordering for Vector3.
//
// \tparam T Datatype of the compared Vector3's elements.
*/
template< typename T >
struct Vector3LexicographicalyLess
{
   /**
   // \brief Provides a lexicographical less-than-operator for Vector3.
   //
   // \param lhs left hand side of less-than-operator.
   // \param rhs right hand side of less-than-operator.
   // \returns true if lhs < rhs (lexicographically), else returns false.
   */
   bool operator()( const Vector3<T> & lhs, const Vector3<T> & rhs ) const
   {
      return std::lexicographical_compare( &lhs[0], &lhs[0] + 3, &rhs[0], &rhs[0] + 3 );
   }
};
//**********************************************************************************************************************


//**********************************************************************************************************************
/**
// \brief Function providing a hash value for Vector3.
//
// \tparam  T Datatype of the Vector3's elements (only integers are supported).
// \param   v The vector the hash is computed for.
// \returns   A hash for the entire Vector3.
*/
template< typename T, typename Enable = std::enable_if_t<std::is_integral_v<T>> >
std::size_t hash_value( const Vector3<T> & v )
{
   std::size_t seed;

   if constexpr( sizeof(std::size_t) >= 8 )
   {
      seed = (static_cast<std::size_t>(v[0]) << 42) +
             (static_cast<std::size_t>(v[1]) << 21) +
             (static_cast<std::size_t>(v[2]) << 0);
   }
   else
   {
      seed = (static_cast<std::size_t>(v[0]) << 21) +
             (static_cast<std::size_t>(v[1]) << 10) +
             (static_cast<std::size_t>(v[2]) << 0);
   }

   return seed;
}
} // namespace math

using math::Vector3;

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
   mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const Vector3<VT> & vec )
   {
      buf.addDebugMarker( "v3" );
      static_assert ( std::is_trivially_copyable< Vector3<VT> >::value,
                      "type has to be trivially copyable for the memcpy to work correctly" );
      auto pos = buf.forward(sizeof(Vector3<VT>));
      std::memcpy(pos, &vec, sizeof(Vector3<VT>));
      return buf;
   }

   template< typename T,    // Element type  of RecvBuffer
             typename VT >  // Element type of vector
   mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, Vector3<VT> & vec )
   {
      buf.readDebugMarker( "v3" );
      static_assert ( std::is_trivially_copyable< Vector3<VT> >::value,
                      "type has to be trivially copyable for the memcpy to work correctly" );
      auto pos = buf.skip(sizeof(Vector3<VT>));
      //suppress https://gcc.gnu.org/onlinedocs/gcc/C_002b_002b-Dialect-Options.html#index-Wclass-memaccess
      std::memcpy(static_cast<void*>(&vec), pos, sizeof(Vector3<VT>));
      return buf;
   }

   template<typename VT>
   struct BufferSizeTrait< walberla::math::Vector3<VT> > {
      static const bool constantSize = true;
      static const uint_t size = 3 * BufferSizeTrait<VT>::size + mpi::BUFFER_DEBUG_OVERHEAD;
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
   struct MPITrait< Vector3<T> >
   {
      static inline MPI_Datatype type()
      {
         // cannot use mpi::Datatype here because its destructor calls MPI_Type_free and static variables are destroyed after the MPI_Finalize
         static MPI_Datatype datatype;
         static bool initialized = false;

         if( ! initialized ) {
            MPI_Type_contiguous(3, MPITrait<T>::type(), &datatype );
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

// Specialization of VectorTrait for Vector3's
template<typename T>
struct VectorTrait< Vector3<T> >
{
   using OutputType = T;

   static const uint_t F_SIZE =  3u;
   static T    get( const Vector3<T> & v, uint_t f )       { return v[f]; }
   static void set(       Vector3<T> & v, uint_t f, T val) { v[f] = val;  }
};

} // namespace walberla


//======================================================================================================================
//
//  comparison backend for Vector3<real_t>
//
//======================================================================================================================

namespace walberla {
namespace debug {
namespace check_functions_detail {

template< >
inline bool check_float_equal( const math::Vector3<real_t> & lhs, const math::Vector3<real_t> & rhs )
{
   return floatIsEqual( lhs[0], rhs[0] ) && floatIsEqual( lhs[1], rhs[1] ) && floatIsEqual( lhs[2], rhs[2] );
}

template< >
inline bool check_float_equal_eps( const math::Vector3<real_t> & lhs, const math::Vector3<real_t> & rhs, const real_t epsilon )
{
   return floatIsEqual( lhs[0], rhs[0], epsilon ) && floatIsEqual( lhs[1], rhs[1], epsilon ) && floatIsEqual( lhs[2], rhs[2], epsilon );
}

}
}
}

namespace std
{
    template<typename T>
    struct hash< walberla::Vector3<T> >
    {
        std::size_t operator()( walberla::Vector3<T> const & v ) const noexcept
        {
            return walberla::math::hash_value( v );
        }
    };
} // namespace std
