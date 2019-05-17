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
//! \file Matrix2.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Implementation of a 3x3 matrix
//
//======================================================================================================================

#pragma once

#include "FPClassify.h"
#include "MathTrait.h"
#include "Vector2.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>


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
/*!\brief Efficient, generic implementation of a 2x2 matrix.
// \ingroup math
//
// The Matrix2 class is the representation of a 2x2 matrix with a total of 9 statically allocated
// elements of arbitrary type. The naming convention of the elements is as follows:

                          \f[\left(\begin{array}{*{2}{c}}
                          xx & xy  \\
                          yx & yy  \\
                          \end{array}\right)\f]\n

// These elements can be accessed directly with the 1D subscript operator or with the 2D function
// operator. The numbering of the matrix elements is

                          \f[\left(\begin{array}{*{3}{c}}
                          0 & 1  \\
                          2 & 3  \\
                          \end{array}\right)\f]
*/
template< typename Type >
class Matrix2
{
private:
   //**Friend declarations*************************************************************************
   /*! \cond internal */
   template< typename Other > friend class Matrix2;
   /*! \endcond */
   //*******************************************************************************************************************

public:
   //**Constructors*****************************************************************************************************
   explicit inline Matrix2();
   explicit inline Matrix2( Type init );
   explicit inline Matrix2( Type xx, Type xy, Type yx, Type yy );
   explicit inline Matrix2( const Type* init );


   inline Matrix2( const Matrix2& m );

   template< typename Other >
   inline Matrix2( const Matrix2<Other>& m );
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   // No explicitly declared destructor.
   //*******************************************************************************************************************

   //**Operators********************************************************************************************************
   /*!\name Operators */
   //@{
                              inline Matrix2&    operator= ( Type set );
                              inline Matrix2&    operator= ( const Matrix2& set );
   template< typename Other > inline Matrix2&    operator= ( const Matrix2<Other>& set );
   template< typename Other > inline bool        operator==( const Matrix2<Other>& rhs )   const;
   template< typename Other > inline bool        operator!=( const Matrix2<Other>& rhs )   const;
                              inline Type&       operator[]( uint_t index );
                              inline const Type& operator[]( uint_t index )                const;
                              inline Type&       operator()( uint_t i, uint_t j );
                              inline const Type& operator()( uint_t i, uint_t j )          const;
   //@}
   //*******************************************************************************************************************

   //**Arithmetic operators************************************************************************
   /*!\name Arithmetic operators
   // \brief The return type of the arithmetic operators depends on the involved data types of
   // \brief the matrices. HIGH denotes the more significant data type of the arithmetic operation
   // \brief (for further detail see the MathTrait class description).
   */
   //@{
   template< typename Other > inline Matrix2&            operator+=( const Matrix2<Other>& rhs );
   template< typename Other > inline Matrix2&            operator-=( const Matrix2<Other>& rhs );
   template< typename Other > inline Matrix2&            operator*=( Other rhs );
   template< typename Other > inline Matrix2&            operator*=( const Matrix2<Other>& rhs );
   template< typename Other > inline const Matrix2<HIGH> operator+ ( const Matrix2<Other>& rhs ) const;
   template< typename Other > inline const Matrix2<HIGH> operator- ( const Matrix2<Other>& rhs ) const;
   template< typename Other > inline const Matrix2<HIGH> operator* ( Other rhs )                 const;
   template< typename Other > inline const Vector2<HIGH> operator* ( const Vector2<Other>& rhs ) const;
   template< typename Other > inline const Matrix2<HIGH> operator* ( const Matrix2<Other>& rhs ) const;
   //@}
   //*******************************************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions
   // \brief The return type of the utility functions depends on the involved data types of the
   // \brief matrices. HIGH denotes the more significant data type of the utility operations
   // \brief (for further detail see the MathTrait class description).
   */
   //@{
   inline Type                getDeterminant()                           const;
   inline Matrix2&            transpose();
   inline const Matrix2       getTranspose()                             const;
   inline Matrix2&            invert();
   inline const Matrix2       getInverse()                               const;
   inline bool                isSingular()                               const;
   inline bool                isSymmetric()                              const;
   inline Type*               data()                                     {return v_;}
   //@}
   //*******************************************************************************************************************


private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Type v_[4];  //!< The nine statically allocated matrix elements.
                /*!< Access to the matrix elements is gained via the subscript or function call
                     operator. The order of the elements is
                     \f[\left(\begin{array}{*{2}{c}}
                     0 & 1  \\
                     2 & 3  \\
                     \end{array}\right)\f] */
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************




//======================================================================================================================
//
//  CONSTRUCTORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Matrix2<Type>::Matrix2()
// \brief The default constructor for Matrix2.
//
// The diagonal matrix elements are initialized with 1, all other elements are initialized
// with 0.
*/
template< typename Type >
inline Matrix2<Type>::Matrix2()
{
   v_[0] = v_[3] = Type(1);
   v_[1] = v_[2] = Type(0);
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>::Matrix2( Type init )
// \brief Constructor for a homogeneous initialization of all elements.
//
// \param init Initial value for all matrix elements.
*/
template< typename Type >
inline Matrix2<Type>::Matrix2( Type init )
{
   v_[0] = v_[1] = v_[2] = v_[3] = init;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>::Matrix2( Type xx, Type xy,  Type yx, Type yy )
// \brief Constructor for a direct initialization of all matrix elements
//
// \param xx The initial value for the xx-component.
// \param xy The initial value for the xy-component.
// \param yx The initial value for the yx-component.
// \param yy The initial value for the yy-component.
*/
template< typename Type >
inline Matrix2<Type>::Matrix2( Type xx, Type xy,
                               Type yx, Type yy  )
{
   v_[0] = xx; v_[1] = xy;
   v_[2] = yx; v_[3] = yy;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>::Matrix2( const Type* init )
// \brief Constructor for an array initializer
//
// \param init Pointer to the initialization array.
//
// The array is assumed to have at least nine valid elements.
*/
template< typename Type >
inline Matrix2<Type>::Matrix2( const Type* init )
{
   v_[0] = init[0];
   v_[1] = init[1];
   v_[2] = init[2];
   v_[3] = init[3];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>::Matrix2( const Matrix2& m )
// \brief The copy constructor for Matrix2.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined in order to enable/facilitate NRV optimization.
*/
template< typename Type >
inline Matrix2<Type>::Matrix2( const Matrix2& m )
{
   v_[0] = m.v_[0];
   v_[1] = m.v_[1];
   v_[2] = m.v_[2];
   v_[3] = m.v_[3];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>::Matrix2( const Matrix2<Other>& m )
// \brief Conversion constructor from different Matrix2 instances.
//
// \param m Matrix to be copied.
*/
template< typename Type >
template< typename Other >
inline Matrix2<Type>::Matrix2( const Matrix2<Other>& m )
{
   v_[0] = m.v_[0];
   v_[1] = m.v_[1];
   v_[2] = m.v_[2];
   v_[3] = m.v_[3];
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::operator=( Type set )
// \brief Homogenous assignment to all matrix elements.
//
// \param set Scalar value to be assigned to all matrix elements.
// \return Reference to the assigned matrix.
*/
template< typename Type >
inline Matrix2<Type>& Matrix2<Type>::operator=( Type set )
{
   v_[0] = set;
   v_[1] = set;
   v_[2] = set;
   v_[3] = set;
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::operator=( const Matrix2& set )
// \brief Copy assignment operator for Matrix2.
//
// \param set Matrix to be copied.
// \return Reference to the assigned matrix.
//
// Explicit definition of a copy assignment operator for performance reasons.
*/
template< typename Type >
inline Matrix2<Type>& Matrix2<Type>::operator=( const Matrix2& set )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = set.v_[0];
   v_[1] = set.v_[1];
   v_[2] = set.v_[2];
   v_[3] = set.v_[3];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::operator=( const Matrix2<Other>& set )
// \brief Assignment operator for different Matrix2 instances.
//
// \param set Matrix to be copied.
// \return Reference to the assigned matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix2<Type>& Matrix2<Type>::operator=( const Matrix2<Other>& set )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = set.v_[0];
   v_[1] = set.v_[1];
   v_[2] = set.v_[2];
   v_[3] = set.v_[3];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Matrix2<Type>::operator==( const Matrix2<Other>& rhs ) const
// \brief Equality operator for the comparison of two matrices.
//
// \param rhs The right-hand-side matrix for the comparison.
// \return bool
*/
template< typename Type >
template< typename Other >
inline bool Matrix2<Type>::operator==( const Matrix2<Other>& rhs ) const
{
   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type.
   if( !equal( v_[0], rhs.v_[0] ) ||
       !equal( v_[1], rhs.v_[1] ) ||
       !equal( v_[2], rhs.v_[2] ) ||
       !equal( v_[3], rhs.v_[3] ) )
      return false;
   else return true;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Matrix2<Type>::operator!=( const Matrix2<Other>& rhs ) const
// \brief Inequality operator for the comparison of two matrices.
//
// \param rhs The right-hand-side matrix for the comparison.
// \return bool
*/
template< typename Type >
template< typename Other >
inline bool Matrix2<Type>::operator!=( const Matrix2<Other>& rhs ) const
{
   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type.
   if( !equal( v_[0], rhs.v_[0] ) ||
       !equal( v_[1], rhs.v_[1] ) ||
       !equal( v_[2], rhs.v_[2] ) ||
       !equal( v_[3], rhs.v_[3] ) )
      return true;
   else return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type& Matrix2<Type>::operator[]( uint_t index )
// \brief 1D-access to the matrix elements.
//
// \param index Access index. The index has to be in the range \f$[0..8]\f$.
// \return Reference to the accessed value.
*/
template< typename Type >
inline Type& Matrix2<Type>::operator[]( uint_t index )
{
   WALBERLA_ASSERT_LESS( index, 4 ,"Invalid matrix access index" );
   return v_[index];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Type& Matrix2<Type>::operator[]( uint_t index ) const
// \brief 1D-access to the matrix elements.
//
// \param index Access index. The index has to be in the range \f$[0..8]\f$.
// \return Reference-to-const to the accessed value.
*/
template< typename Type >
inline const Type& Matrix2<Type>::operator[]( uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, 4 ,"Invalid matrix access index" );
   return v_[index];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type& Matrix2<Type>::operator()( uint_t i, uint_t j )
// \brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..2].
// \param j Access index for the column. The index has to be in the range [0..2].
// \return Reference to the accessed value.
*/
template< typename Type >
inline Type& Matrix2<Type>::operator()( uint_t i, uint_t j )
{
   WALBERLA_ASSERT( i<2 && j<2,"Invalid matrix access index" );
   return v_[ i*2 +j ];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Type& Matrix2<Type>::operator()( uint_t i, uint_t j ) const
// \brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..2].
// \param j Access index for the column. The index has to be in the range [0..2].
// \return Reference-to-const to the accessed value.
*/
template< typename Type >
inline const Type& Matrix2<Type>::operator()( uint_t i, uint_t j ) const
{
   WALBERLA_ASSERT( i<2 && j<2 ,"Invalid matrix access index" );
   return v_[ i*2 +j ];
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  ARITHMETIC OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::operator+=( const Matrix2<Other>& rhs )
// \brief Addition assignment operator for the addition of two matrices (\f$ A+=B \f$).
//
// \param rhs The right-hand-side matrix to be added to the matrix.
// \return Reference to the matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix2<Type>& Matrix2<Type>::operator+=( const Matrix2<Other>& rhs )
{
   v_[0] += rhs.v_[0];
   v_[1] += rhs.v_[1];
   v_[2] += rhs.v_[2];
   v_[3] += rhs.v_[3];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::operator-=( const Matrix2<Other>& rhs )
// \brief Subtraction assignment operator for the subtraction of two matrices (\f$ A-=B \f$).
//
// \param rhs The right-hand-side matrix to be subtracted from the matrix.
// \return Reference to the matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix2<Type>& Matrix2<Type>::operator-=( const Matrix2<Other>& rhs )
{
   v_[0] -= rhs.v_[0];
   v_[1] -= rhs.v_[1];
   v_[2] -= rhs.v_[2];
   v_[3] -= rhs.v_[3];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::operator*=( Other rhs )
// \brief Multiplication assignment operator for the multiplication between a matrix and
// \brief a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand-side scalar value for the multiplication.
// \return Reference to the matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix2<Type>& Matrix2<Type>::operator*=( Other rhs )
{
   v_[0] *= rhs;
   v_[1] *= rhs;
   v_[2] *= rhs;
   v_[3] *= rhs;
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::operator*=( const Matrix2<Other>& rhs )
// \brief Multiplication assignment operator for the multiplication between two matrices
// \brief (\f$ A*=B \f$).
//
// \param rhs The right-hand-side matrix for the multiplication.
// \return Reference to the matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix2<Type>& Matrix2<Type>::operator*=( const Matrix2<Other>& rhs )
{
   // Creating a temporary due to data dependencies
   Matrix2 tmp( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[2],
                v_[0]*rhs.v_[1] + v_[1]*rhs.v_[3],
                v_[2]*rhs.v_[0] + v_[3]*rhs.v_[2],
                v_[2]*rhs.v_[1] + v_[3]*rhs.v_[3] );

   return this->operator=( tmp );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<HIGH> Matrix2<Type>::operator+( const Matrix2<Other>& rhs ) const
// \brief Addition operator for the addition of two matrices (\f$ A=B+C \f$).
//
// \param rhs The right-hand-side matrix to be added to the matrix.
// \return The sum of the two matrices.
*/
template< typename Type >
template< typename Other >
inline const Matrix2<HIGH> Matrix2<Type>::operator+( const Matrix2<Other>& rhs ) const
{
   return Matrix2<HIGH>( v_[0] + rhs.v_[0],
                         v_[1] + rhs.v_[1],
                         v_[2] + rhs.v_[2],
                         v_[3] + rhs.v_[3] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<HIGH> Matrix2<Type>::operator-( const Matrix2<Other>& rhs ) const
// \brief Subtraction operator for the subtraction of two matrices (\f$ A=B-C \f$).
//
// \param rhs The right-hand-side matrix to be subtracted from the matrix.
// \return The difference of the two matrices.
*/
template< typename Type >
template< typename Other >
inline const Matrix2<HIGH> Matrix2<Type>::operator-( const Matrix2<Other>& rhs ) const
{
   return Matrix2<HIGH>( v_[0] - rhs.v_[0],
                         v_[1] - rhs.v_[1],
                         v_[2] - rhs.v_[2],
                         v_[3] - rhs.v_[3] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<HIGH> Matrix2<Type>::operator*( Other rhs ) const
// \brief Multiplication operator for the multiplication of a matrix and a scalar value
// \brief (\f$ A=B*s \f$).
//
// \param rhs The right-hand-side scalar value for the multiplication.
// \return The scaled result matrix.
*/
template< typename Type >
template< typename Other >
inline const Matrix2<HIGH> Matrix2<Type>::operator*( Other rhs ) const
{
   return Matrix2<HIGH>( v_[0]*rhs, v_[1]*rhs,
                         v_[2]*rhs, v_[3]*rhs  );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector3<HIGH> Matrix2<Type>::operator*( const Vector3<Other>& rhs ) const
// \brief Multiplication operator for the multiplication of a matrix and a vector
// \brief (\f$ \vec{a}=B*\vec{c} \f$).
//
// \param rhs The right-hand-side vector for the multiplication.
// \return The resulting vector.
*/
template< typename Type >
template< typename Other >
inline const Vector2<HIGH> Matrix2<Type>::operator*( const Vector2<Other>& rhs ) const
{
   return Vector2<HIGH>( v_[0]*rhs[0] + v_[1]*rhs[1],
                         v_[2]*rhs[0] + v_[3]*rhs[1] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<HIGH> Matrix2<Type>::operator*( const Matrix2<Other>& rhs ) const
// \brief Multiplication operator for the multiplication of two matrices (\f$ A=B*C \f$).
//
// \param rhs The right-hand-side matrix for the multiplication.
// \return The resulting matrix.
*/
template< typename Type >
template< typename Other >
inline const Matrix2<HIGH> Matrix2<Type>::operator*( const Matrix2<Other>& rhs ) const
{
   return Matrix2<HIGH>( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[2],
                         v_[0]*rhs.v_[1] + v_[1]*rhs.v_[3],
                         v_[2]*rhs.v_[0] + v_[3]*rhs.v_[2],
                         v_[2]*rhs.v_[1] + v_[3]*rhs.v_[3] );

}
//**********************************************************************************************************************




//======================================================================================================================
//
//  UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Type Matrix2<Type>::getDeterminant() const
// \brief Calculation of the determinant of the matrix.
//
// \return The determinant of the matrix.
*/
template< typename Type >
inline Type Matrix2<Type>::getDeterminant() const
{
   return   v_[0] * v_[3]
          - v_[2] * v_[1] ;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::transpose()
// \brief Transposing the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type >
inline Matrix2<Type>& Matrix2<Type>::transpose()
{
   std::swap( v_[1], v_[2] );
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<Type> Matrix2<Type>::getTranspose() const
// \brief Calculation of the transpose of the matrix.
//
// \return The transpose of the matrix.
*/
template< typename Type >
inline const Matrix2<Type> Matrix2<Type>::getTranspose() const
{
   return Matrix2( v_[0], v_[2], v_[1], v_[3] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix2<Type>& Matrix2<Type>::invert()
// \brief Inverting the matrix.
//
// \return Reference to the inverted matrix.
//
// The calculation is performed with the matrix inversion by Cramer. This function is only
// defined for matrices of floating point type. The attempt to use this function with matrices
// of integral data types will result in a compile time error.
*/
template< typename Type >
inline Matrix2<Type>& Matrix2<Type>::invert()
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   Type det = getDeterminant();

   WALBERLA_ASSERT_NOT_IDENTICAL( det, Type(0), "Matrix is singular and cannot be inverted" );

   det = Type(1) / det;

   // Creating a temporary due to data dependencies
   return *this = Matrix2( det * v_[3],  - det * v_[1], - det * v_[2], det * v_[0] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<Type> Matrix2<Type>::getInverse() const
// \brief Calculation of the inverse of the matrix.
//
// \return The inverse of the matrix.
//
// The calculation is performed with the matrix inversion by Cramer. This function is only
// defined for matrices of floating point type. The attempt to use this function with matrices
// of integral data types will result in a compile time error.
*/
template< typename Type >
inline const Matrix2<Type> Matrix2<Type>::getInverse() const
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   Type det = getDeterminant();

   WALBERLA_ASSERT_NOT_IDENTICAL( det, Type(0), "Matrix is singular and cannot be inverted" );

   det = Type(1) / det;

   // Creating a temporary due to data dependencies
   return Matrix2( det * v_[3],  - det * v_[1],
                 - det * v_[2], det * v_[0]     );
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\fn bool Matrix2<Type>::isSingular() const
// \brief Singularity check for the matrix (det=0).
//
// \return \a true if the matrix is singular, \a false if not.
*/
template< typename Type >
inline bool Matrix2<Type>::isSingular() const
{
   if( equal( getDeterminant(), Type(0) ) )
      return true;
   else return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Matrix2<Type>::isSymmetric() const
// \brief Checks if the matrix is symmetric.
//
// \return \a true if the matrix is symmetric, \a false if not.
*/
template< typename Type >
inline bool Matrix2<Type>::isSymmetric() const
{
   return ( v_[1] == v_[2] );
}
//**********************************************************************************************************************



//======================================================================================================================
//
//  GLOBAL OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\name Matrix2 operators */
//@{
//template< typename Type, typename Other >
//inline const Matrix2<HIGH> operator*( Other scalar, const Matrix2<Type>& matrix );

template< typename Type >
std::ostream& operator<<( std::ostream& os, const Matrix2<Type>& m );

template< typename Type >
inline bool isnan( const Matrix2<Type>& m );

template< typename Type >
inline const Matrix2<Type> abs( const Matrix2<Type>& m );

template< typename Type >
inline const Matrix2<Type> fabs( const Matrix2<Type>& m );
//@}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<Type> operator*( Other scalar, const Matrix2<Type>& matrix )
// \brief Multiplication operator for the multiplication of a scalar value and a matrix
// \brief (\f$ A=s*B \f$).
//
// \param scalar The left-hand-side scalar value for the multiplication.
// \param matrix The right-hand-side matrix for the multiplication.
// \return The scaled result matrix.
*/
//template< typename Type, typename Other >
//inline const Matrix2<HIGH> operator*( Other scalar, const Matrix2<Type>& matrix )
//{
//   static_assert( ! std::is_scalar<Other>::value, "Only scalar types allowed" );
//   return matrix*scalar;
//}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::ostream& operator<<( std::ostream& os, const Matrix2<Type>& m )
// \brief Global output operator for 2x2 matrices.
//
// \param os Reference to the output stream.
// \param m Reference to a constant matrix object.
// \return Reference to the output stream.
*/
template< typename Type >
std::ostream& operator<<( std::ostream& os, const Matrix2<Type>& m )
{
   return os << " ( " << m[0] << " , " << m[1] << " )\n"
             << " ( " << m[2] << " , " << m[3] << " )\n";
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool math::isnan( const Matrix2<Type>& m )
// \brief Checks the given matrix for not-a-number elements.
//
// \param m The matrix to be checked for not-a-number elements.
// \return \a true if at least one element of the matrix is not-a-number, \a false otherwise.
*/
template< typename Type >
inline bool isnan( const Matrix2<Type>& m )
{
   if( math::isnan( m[0] ) || math::isnan( m[1] )||
       math::isnan( m[2] ) || math::isnan( m[3] )  )
      return true;
   else
      return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<Type> abs( const Matrix2<Type>& m )
// \brief Returns a matrix containing the absolute values of each single element of \a m.
//
// \param m The integral input matrix.
// \return The absolute value of each single element of \a m.
//
// The \a abs function calculates the absolute value of each element of the input matrix
// \a m. This function can only be applied to matrices of integral data type. For floating
// point matrices, the pe::fabs( const Matrix2& ) function can be used.
*/
template< typename Type >
inline const Matrix2<Type> abs( const Matrix2<Type>& m )
{
   WALBERLA_STATIC_ASSERT( std::numeric_limits<Type>::is_integer );
   return Matrix2<Type>( std::abs(m[0]), std::abs(m[1]),
                         std::abs(m[2]), std::abs(m[3]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<Type> fabs( const Matrix2<Type>& m )
// \brief Returns a matrix containing the absolute values of each single element of \a m.
//
// \param m The floating point input matrix.
// \return The absolute value of each single element of \a m.
//
// The \a fabs function calculates the absolute value of each element of the input matrix
// \a m. This function can only be applied to  floating point matrices. For matrices of
// integral data type, the pe::abs( const Matrix2& ) function can be used.
*/
template< typename Type >
inline const Matrix2<Type> fabs( const Matrix2<Type>& m )
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );
   return Matrix2<Type>( std::fabs(m[0]), std::fabs(m[1]),
                         std::fabs(m[2]), std::fabs(m[3]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix2<Type> fabs( const Matrix2<Type>& m )
// \brief Returns a matrix containing the absolute values of each single element of \a m.
//
// \param m The floating point input matrix.
// \return The absolute value of each single element of \a m.
//
// The \a fabs function calculates the absolute value of each element of the input matrix
// \a m. This function can only be applied to  floating point matrices. For matrices of
// integral data type, the pe::abs( const Matrix2& ) function can be used.
*/
template< typename T0, typename T1>
Matrix2< typename MathTrait<T0,T1>::High > tensorProduct( Vector2<T0> v0, Vector2<T1> v1 ){
   return Matrix2< typename MathTrait<T0,T1>::High > (
      v0[0] * v1[0], v0[0] * v1[1],
      v0[1] * v1[0], v0[1] * v1[1] );
}
//**********************************************************************************************************************



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename MT >  // Element type of matrix
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const Matrix2<MT> & m )
{
   buf.addDebugMarker( "m2" );
   static_assert ( std::is_trivially_copyable< Matrix2<MT> >::value,
                   "type has to be trivially copyable for the memcpy to work correctly" );
   auto pos = buf.forward(sizeof(Matrix2<MT>));
   std::memcpy(pos, &m, sizeof(Matrix2<MT>));
   return buf;
}

template< typename T,    // Element type  of RecvBuffer
          typename MT >  // Element type of matrix
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, Matrix2<MT> & m )
{
   buf.readDebugMarker( "m2" );
   static_assert ( std::is_trivially_copyable< Matrix2<MT> >::value,
                   "type has to be trivially copyable for the memcpy to work correctly" );
   auto pos = buf.skip(sizeof(Matrix2<MT>));
   //suppress https://gcc.gnu.org/onlinedocs/gcc/C_002b_002b-Dialect-Options.html#index-Wclass-memaccess
   std::memcpy(static_cast<void*>(&m), pos, sizeof(Matrix2<MT>));
   return buf;
}



} // namespace math

using math::Matrix2;

} // namespace walberla

//======================================================================================================================
//
//  Vector Trait Specialization
//
//======================================================================================================================

namespace walberla {

// Specialization of VectorTrait for Matrix2s
template<typename T>
struct VectorTrait< Matrix2<T> >
{
   typedef T OutputType;

   static const uint_t F_SIZE =  4u;
   static T    get( const Matrix2<T> & v, uint_t f )       { return v[f]; }
   static void set(       Matrix2<T> & v, uint_t f, T val) { v[f] = val;  }
};

} // namespace walberla

#undef HIGH
