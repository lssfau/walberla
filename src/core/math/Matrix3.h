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
//! \file Matrix3.h
//! \ingroup core
//! \author Klaus Iglberger
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Implementation of a 3x3 matrix
//
//======================================================================================================================

#pragma once

#include "FPClassify.h"
#include "MathTrait.h"
#include "Vector3.h"

#include "core/debug/Debug.h"
#include "core/mpi/Datatype.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <type_traits>

#include <algorithm>
#include <array>
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
/*!\brief Efficient, generic implementation of a 3x3 matrix.
// \ingroup math
//
// The Matrix3 class is the representation of a 3x3 matrix with a total of 9 statically allocated
// elements of arbitrary type. The naming convention of the elements is as follows:

                          \f[\left(\begin{array}{*{3}{c}}
                          xx & xy & xz \\
                          yx & yy & yz \\
                          zx & zy & zz \\
                          \end{array}\right)\f]\n

// These elements can be accessed directly with the 1D subscript operator or with the 2D function
// operator. The numbering of the matrix elements is

                          \f[\left(\begin{array}{*{3}{c}}
                          0 & 1 & 2 \\
                          3 & 4 & 5 \\
                          6 & 7 & 8 \\
                          \end{array}\right)\f]
*/
template< typename Type >
class Matrix3
{
private:
   //**Friend declarations*************************************************************************
   /*! \cond internal */
   template< typename Other > friend class Matrix3;
   /*! \endcond */
   //*******************************************************************************************************************

public:
   //**Constructors*****************************************************************************************************
   explicit inline constexpr Matrix3() = default;
   explicit inline constexpr Matrix3( Type init );
   explicit inline constexpr Matrix3( const Vector3<Type>& a, const Vector3<Type>& b, const Vector3<Type>& c );
   explicit inline constexpr Matrix3( Type xx, Type xy, Type xz, Type yx, Type yy, Type yz, Type zx, Type zy, Type zz );
   explicit inline constexpr Matrix3( const Type* init );

   template< typename Axis, typename Angle >
   explicit Matrix3( Vector3<Axis> axis, Angle angle );

   inline Matrix3( const Matrix3& m ) = default;

   template< typename Other >
   inline Matrix3( const Matrix3<Other>& m );

   inline static Matrix3 makeDiagonalMatrix( const Type xx, const Type yy, const Type zz );
   inline static Matrix3 makeDiagonalMatrix( const Type d );
   inline static Matrix3 makeIdentityMatrix();
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   // No explicitly declared destructor.
   //*******************************************************************************************************************

   //**Operators********************************************************************************************************
   /*!\name Operators */
   //@{
                              inline Matrix3&    operator= ( Type set );
                              inline Matrix3&    operator= ( const Matrix3& set ) = default;
   template< typename Other > inline Matrix3&    operator= ( const Matrix3<Other>& set );
   template< typename Other > inline bool        operator==( const Matrix3<Other>& rhs )   const;
   template< typename Other > inline bool        operator!=( const Matrix3<Other>& rhs )   const;
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
   template< typename Other > inline Matrix3&            operator+=( const Matrix3<Other>& rhs );
   template< typename Other > inline Matrix3&            operator-=( const Matrix3<Other>& rhs );
   template< typename Other > inline Matrix3&            operator*=( const Matrix3<Other>& rhs );
   template< typename Other > inline const Matrix3<HIGH> operator+ ( const Matrix3<Other>& rhs ) const;
   template< typename Other > inline const Matrix3<HIGH> operator- ( const Matrix3<Other>& rhs ) const;
   template< typename Other > inline const Vector3<HIGH> operator* ( const Vector3<Other>& rhs ) const;
   template< typename Other > inline const Matrix3<HIGH> operator* ( const Matrix3<Other>& rhs ) const;

   template< typename Other > inline typename std::enable_if< std::is_arithmetic<Other>::value, Matrix3&            >::type operator*=( Other rhs );
   template< typename Other > inline typename std::enable_if< std::is_arithmetic<Other>::value, const Matrix3<HIGH> >::type operator* ( Other rhs ) const;
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
                              inline Matrix3&            transpose();
                              inline const Matrix3       getTranspose()                             const;
                              inline Matrix3&            invert();
                              inline const Matrix3       getInverse()                               const;
   template< typename Other > inline const Vector3<HIGH> multTranspose( const Vector3<Other>& rhs ) const;
   template< typename Other > inline const Matrix3<HIGH> rotate( const Matrix3<Other>& m )          const;
   template< typename Other > inline const Matrix3<HIGH> diagRotate( const Matrix3<Other>& m )      const;
                              inline bool                isSingular()                               const;
                              inline bool                isSymmetric()                              const;
                              inline bool                isZero()                                   const;
                              inline const Matrix3       getCholesky()                              const;
   template< typename Other > inline const Vector3<HIGH> solve( const Vector3<Other> &rhs )         const;
                              inline Type                trace()                                    const;
                              inline Type*               data()                                     {return v_.data();}
                              inline Type const *        data()                                     const {return v_.data();}
   //@}
   //*******************************************************************************************************************

   //**Euler rotations*****************************************************************************
   //! Order of the Euler rotation
   /*! This codes are needed for the EulerAngles function in order to calculate the Euler angles
       for a specific combination of rotations. */
   enum EulerRotation {
      XYZs =  0,  //!< Rotation order x, y, z in a static frame.
      ZYXr =  1,  //!< Rotation order z, y, x in a rotating frame.
      XYXs =  2,  //!< Rotation order x, y, x in a static frame.
      XYXr =  3,  //!< Rotation order x, y, z in a rotating frame.
      XZYs =  4,  //!< Rotation order x, z, y in a static frame.
      YZXr =  5,  //!< Rotation order y, z, x in a rotating frame.
      XZXs =  6,  //!< Rotation order x, z, x in a static frame.
      XZXr =  7,  //!< Rotation order x, z, x in a rotating frame.
      YZXs =  8,  //!< Rotation order y, z, x in a static frame.
      XZYr =  9,  //!< Rotation order x, z, y in a rotating frame.
      YZYs = 10,  //!< Rotation order y, z, y in a static frame.
      YZYr = 11,  //!< Rotation order y, z, y in a rotating frame.
      YXZs = 12,  //!< Rotation order y, x, z in a static frame.
      ZXYr = 13,  //!< Rotation order z, x, y in a rotating frame.
      YXYs = 14,  //!< Rotation order y, x, y in a static frame.
      YXYr = 15,  //!< Rotation order y, x, y in a rotating frame.
      ZXYs = 16,  //!< Rotation order z, x, y in a static frame.
      YXZr = 17,  //!< Rotation order y, x, z in a rotating frame.
      ZXZs = 18,  //!< Rotation order z, x, z in a static frame.
      ZXZr = 19,  //!< Rotation order z, x, z in a rotating frame.
      ZYXs = 20,  //!< Rotation order z, y, x in a static frame.
      XYZr = 21,  //!< Rotation order x, y, z in a rotating frame.
      ZYZs = 22,  //!< Rotation order z, y, z in a static frame.
      ZYZr = 23   //!< Rotation order z, y, z in a rotating frame.
   };

   /*!\name Euler rotations
   // \brief For the classification of the Euler rotation, the following characteristics are
   // \brief defined:\n
   // \brief  - Inner axis: the inner axis is the axis of the first rotation matrix multiplied
   // \brief    to a vector.
   // \brief  - Parity: the parity is even, if the inner axis X is followed by the middle axis
   // \brief    Y, or Y is followed by Z, or Z is followed by X; otherwise parity is odd.
   // \brief  - Repetition: repetition tells, if the first and last axes are the same or different.
   // \brief  - Frame: the frame refers to the frame from which the Euler angles are calculated.
   // \brief
   // \brief Altogether, there are 24 possible Euler rotations. The possibilities are consisting
   // \brief of the choice of the inner axis (X,Y or Z), the parity (even or odd), repetition
   // \brief (yes or no) and the frame (static or rotating). E.g., an Euler order of XYZs stands
   // \brief for the rotation order of x-, y- and z-axis in a static frame (inner axis: X, parity:
   // \brief even, repetition: no, frame: static), whereas YXYr stands for the rotation order y-,
   // \brief x- and y-axis in a rotating frame ( inner axis: Y, parity: odd, repetition: yes,
   // \brief frame: rotating).
   */
   //@{
          const Vector3<Type> getEulerAngles( EulerRotation order ) const;
   inline const Vector3<Type> getEulerAnglesXYZ()                   const;
   //@}
   //*******************************************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   /**
    * The nine statically allocated matrix elements.
    * Access to the matrix elements is gained via the subscript or function call
    * operator. The order of the elements is
    * \f[\left(\begin{array}{*{3}{c}}
    * 0 & 1 & 2 \\
    * 3 & 4 & 5 \\
    * 6 & 7 & 8 \\
    * \end{array}\right)\f]
   **/
   std::array<Type, 9> v_ = {{Type(1), Type(0), Type(0),
                              Type(0), Type(1), Type(0),
                              Type(0), Type(0), Type(1)}};
   //@}
   //*******************************************************************************************************************
};
static_assert( std::is_trivially_copyable<Matrix3<real_t>>::value, "Matrix3<real_t> has to be trivially copyable!");
//**********************************************************************************************************************




//======================================================================================================================
//
//  CONSTRUCTORS
//
//======================================================================================================================


//**********************************************************************************************************************
/*!\fn Matrix3<Type>::Matrix3( Type init )
// \brief Constructor for a homogeneous initialization of all elements.
//
// \param init Initial value for all matrix elements.
*/
template< typename Type >
inline constexpr Matrix3<Type>::Matrix3( Type init )
{
   v_[0] = v_[1] = v_[2] = v_[3] = v_[4] = v_[5] = v_[6] = v_[7] = v_[8] = init;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
// \brief Constructor for a direct initialization of all matrix elements
//
// \param a The first column of the matrix.
// \param b The second column of the matrix.
// \param c The third column of the matrix.
//**********************************************************************************************************************
template< typename Type >
inline constexpr Matrix3<Type>::Matrix3( const Vector3<Type>& a, const Vector3<Type>& b, const Vector3<Type>& c )
{
   v_[0] = a[0]; v_[1] = b[0]; v_[2] = c[0];
   v_[3] = a[1]; v_[4] = b[1]; v_[5] = c[1];
   v_[6] = a[2]; v_[7] = b[2]; v_[8] = c[2];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>::Matrix3( Type xx, Type xy, Type xz, Type yx, Type yy, Type yz, Type zx, Type zy, Type zz )
// \brief Constructor for a direct initialization of all matrix elements
//
// \param xx The initial value for the xx-component.
// \param xy The initial value for the xy-component.
// \param xz The initial value for the xz-component.
// \param yx The initial value for the yx-component.
// \param yy The initial value for the yy-component.
// \param yz The initial value for the yz-component.
// \param zx The initial value for the zx-component.
// \param zy The initial value for the zy-component.
// \param zz The initial value for the zz-component.
*/
template< typename Type >
inline constexpr Matrix3<Type>::Matrix3( Type xx, Type xy, Type xz,
                                         Type yx, Type yy, Type yz,
                                         Type zx, Type zy, Type zz )
{
   v_[0] = xx; v_[1] = xy; v_[2] = xz;
   v_[3] = yx; v_[4] = yy; v_[5] = yz;
   v_[6] = zx; v_[7] = zy; v_[8] = zz;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>::Matrix3( const Type* init )
// \brief Constructor for an array initializer
//
// \param init Pointer to the initialization array.
//
// The array is assumed to have at least nine valid elements.
*/
template< typename Type >
inline constexpr Matrix3<Type>::Matrix3( const Type* init )
{
   v_[0] = init[0];
   v_[1] = init[1];
   v_[2] = init[2];
   v_[3] = init[3];
   v_[4] = init[4];
   v_[5] = init[5];
   v_[6] = init[6];
   v_[7] = init[7];
   v_[8] = init[8];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>::Matrix3( Vector3<Axis> axis, Angle angle )
// \brief Rotation matrix constructor.
//
// \param axis The rotation axis.
// \param angle The rotation angle.
//
// This constructor is only defined for floating point vectors. The attempt to use this
// constructor for vectors of integral data type results in a compile time error.
*/
template< typename Type >
template< typename Axis, typename Angle >
Matrix3<Type>::Matrix3( Vector3<Axis> axis, Angle angle )
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer  );
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Axis>::is_integer  );
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Angle>::is_integer );

   const Angle sina( std::sin(angle) );
   const Angle cosa( std::cos(angle) );
   const Angle tmp( Angle(1)-cosa );

   normalize(axis);

   v_[0] = cosa + axis[0]*axis[0]*tmp;
   v_[1] = axis[0]*axis[1]*tmp - axis[2]*sina;
   v_[2] = axis[0]*axis[2]*tmp + axis[1]*sina;
   v_[3] = axis[1]*axis[0]*tmp + axis[2]*sina;
   v_[4] = cosa + axis[1]*axis[1]*tmp;
   v_[5] = axis[1]*axis[2]*tmp - axis[0]*sina;
   v_[6] = axis[2]*axis[0]*tmp - axis[1]*sina;
   v_[7] = axis[2]*axis[1]*tmp + axis[0]*sina;
   v_[8] = cosa + axis[2]*axis[2]*tmp;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>::Matrix3( const Matrix3<Other>& m )
// \brief Conversion constructor from different Matrix3 instances.
//
// \param m Matrix to be copied.
*/
template< typename Type >
template< typename Other >
inline Matrix3<Type>::Matrix3( const Matrix3<Other>& m )
{
   v_[0] = m.v_[0];
   v_[1] = m.v_[1];
   v_[2] = m.v_[2];
   v_[3] = m.v_[3];
   v_[4] = m.v_[4];
   v_[5] = m.v_[5];
   v_[6] = m.v_[6];
   v_[7] = m.v_[7];
   v_[8] = m.v_[8];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type> Matrix3<Type>::makeDiagonalMatrix( const Type xx, const Type yy, const Type zz )
// \brief Named constructor to create a diagonal matrix. All non-diagonal elements are initialized with zero.
//
// \param xx value for element (0,0).
// \param yy value for element (1,1).
// \param zz value for element (2,2).
*/
template< typename Type >
Matrix3<Type> Matrix3<Type>::makeDiagonalMatrix( const Type xx, const Type yy, const Type zz )
{
   return Matrix3<Type>( xx, Type(), Type(), Type(), yy, Type(), Type(), Type(), zz );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type> Matrix3<Type>::makeDiagonalMatrix( const Type d )
// \brief Named constructor to create a diagonal matrix. All non-diagonal elements are initialized with zero.
//
// \param d value for diagonal elements.
*/
template< typename Type >
Matrix3<Type> Matrix3<Type>::makeDiagonalMatrix( const Type d )
{
   return makeDiagonalMatrix( d, d, d );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type> Matrix3<Type>::makeIdentityMatrix()
// \brief Named constructor to create the identity matrix.
// 
// All diagonal elements are initialized to one, alls others to zero.
//
*/
template< typename Type >
Matrix3<Type> Matrix3<Type>::makeIdentityMatrix()
{
   return makeDiagonalMatrix( Type(1) );
}
//**********************************************************************************************************************


//======================================================================================================================
//
//  OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Matrix3<Type>& Matrix3<Type>::operator=( Type set )
// \brief Homogenous assignment to all matrix elements.
//
// \param set Scalar value to be assigned to all matrix elements.
// \return Reference to the assigned matrix.
*/
template< typename Type >
inline Matrix3<Type>& Matrix3<Type>::operator=( Type set )
{
   v_[0] = set;
   v_[1] = set;
   v_[2] = set;
   v_[3] = set;
   v_[4] = set;
   v_[5] = set;
   v_[6] = set;
   v_[7] = set;
   v_[8] = set;
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>& Matrix3<Type>::operator=( const Matrix3<Other>& set )
// \brief Assignment operator for different Matrix3 instances.
//
// \param set Matrix to be copied.
// \return Reference to the assigned matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix3<Type>& Matrix3<Type>::operator=( const Matrix3<Other>& set )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = set.v_[0];
   v_[1] = set.v_[1];
   v_[2] = set.v_[2];
   v_[3] = set.v_[3];
   v_[4] = set.v_[4];
   v_[5] = set.v_[5];
   v_[6] = set.v_[6];
   v_[7] = set.v_[7];
   v_[8] = set.v_[8];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Matrix3<Type>::operator==( const Matrix3<Other>& rhs ) const
// \brief Equality operator for the comparison of two matrices.
//
// \param rhs The right-hand-side matrix for the comparison.
// \return bool
*/
template< typename Type >
template< typename Other >
inline bool Matrix3<Type>::operator==( const Matrix3<Other>& rhs ) const
{
   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type.
   if( !equal( v_[0], rhs.v_[0] ) ||
       !equal( v_[1], rhs.v_[1] ) ||
       !equal( v_[2], rhs.v_[2] ) ||
       !equal( v_[3], rhs.v_[3] ) ||
       !equal( v_[4], rhs.v_[4] ) ||
       !equal( v_[5], rhs.v_[5] ) ||
       !equal( v_[6], rhs.v_[6] ) ||
       !equal( v_[7], rhs.v_[7] ) ||
       !equal( v_[8], rhs.v_[8] ) )
      return false;
   else return true;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Matrix3<Type>::operator!=( const Matrix3<Other>& rhs ) const
// \brief Inequality operator for the comparison of two matrices.
//
// \param rhs The right-hand-side matrix for the comparison.
// \return bool
*/
template< typename Type >
template< typename Other >
inline bool Matrix3<Type>::operator!=( const Matrix3<Other>& rhs ) const
{
   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type.
   if( !equal( v_[0], rhs.v_[0] ) ||
       !equal( v_[1], rhs.v_[1] ) ||
       !equal( v_[2], rhs.v_[2] ) ||
       !equal( v_[3], rhs.v_[3] ) ||
       !equal( v_[4], rhs.v_[4] ) ||
       !equal( v_[5], rhs.v_[5] ) ||
       !equal( v_[6], rhs.v_[6] ) ||
       !equal( v_[7], rhs.v_[7] ) ||
       !equal( v_[8], rhs.v_[8] ) )
      return true;
   else return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type& Matrix3<Type>::operator[]( uint_t index )
// \brief 1D-access to the matrix elements.
//
// \param index Access index. The index has to be in the range \f$[0..8]\f$.
// \return Reference to the accessed value.
*/
template< typename Type >
inline Type& Matrix3<Type>::operator[]( uint_t index )
{
   WALBERLA_ASSERT_LESS( index, 9 ,"Invalid matrix access index" );
   return v_[index];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Type& Matrix3<Type>::operator[]( uint_t index ) const
// \brief 1D-access to the matrix elements.
//
// \param index Access index. The index has to be in the range \f$[0..8]\f$.
// \return Reference-to-const to the accessed value.
*/
template< typename Type >
inline const Type& Matrix3<Type>::operator[]( uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, 9 ,"Invalid matrix access index" );
   return v_[index];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type& Matrix3<Type>::operator()( uint_t i, uint_t j )
// \brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..2].
// \param j Access index for the column. The index has to be in the range [0..2].
// \return Reference to the accessed value.
*/
template< typename Type >
inline Type& Matrix3<Type>::operator()( uint_t i, uint_t j )
{
   WALBERLA_ASSERT( i<3 && j<3,"Invalid matrix access index" );
   return v_[i*3+j];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Type& Matrix3<Type>::operator()( uint_t i, uint_t j ) const
// \brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..2].
// \param j Access index for the column. The index has to be in the range [0..2].
// \return Reference-to-const to the accessed value.
*/
template< typename Type >
inline const Type& Matrix3<Type>::operator()( uint_t i, uint_t j ) const
{
   WALBERLA_ASSERT( i<3 && j<3 ,"Invalid matrix access index" );
   return v_[i*3+j];
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  ARITHMETIC OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Matrix3<Type>& Matrix3<Type>::operator+=( const Matrix3<Other>& rhs )
// \brief Addition assignment operator for the addition of two matrices (\f$ A+=B \f$).
//
// \param rhs The right-hand-side matrix to be added to the matrix.
// \return Reference to the matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix3<Type>& Matrix3<Type>::operator+=( const Matrix3<Other>& rhs )
{
   v_[0] += rhs.v_[0];
   v_[1] += rhs.v_[1];
   v_[2] += rhs.v_[2];
   v_[3] += rhs.v_[3];
   v_[4] += rhs.v_[4];
   v_[5] += rhs.v_[5];
   v_[6] += rhs.v_[6];
   v_[7] += rhs.v_[7];
   v_[8] += rhs.v_[8];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>& Matrix3<Type>::operator-=( const Matrix3<Other>& rhs )
// \brief Subtraction assignment operator for the subtraction of two matrices (\f$ A-=B \f$).
//
// \param rhs The right-hand-side matrix to be subtracted from the matrix.
// \return Reference to the matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix3<Type>& Matrix3<Type>::operator-=( const Matrix3<Other>& rhs )
{
   v_[0] -= rhs.v_[0];
   v_[1] -= rhs.v_[1];
   v_[2] -= rhs.v_[2];
   v_[3] -= rhs.v_[3];
   v_[4] -= rhs.v_[4];
   v_[5] -= rhs.v_[5];
   v_[6] -= rhs.v_[6];
   v_[7] -= rhs.v_[7];
   v_[8] -= rhs.v_[8];
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>& Matrix3<Type>::operator*=( Other rhs )
// \brief Multiplication assignment operator for the multiplication between a matrix and
// \brief a scalar value (\f$ A*=s \f$).
//
// \param rhs The right-hand-side scalar value for the multiplication.
// \return Reference to the matrix.
*/
template< typename Type >
template< typename Other >
inline typename std::enable_if< std::is_arithmetic<Other>::value, Matrix3<Type>& >::type Matrix3<Type>::operator*=( Other rhs )
{
   v_[0] *= rhs;
   v_[1] *= rhs;
   v_[2] *= rhs;
   v_[3] *= rhs;
   v_[4] *= rhs;
   v_[5] *= rhs;
   v_[6] *= rhs;
   v_[7] *= rhs;
   v_[8] *= rhs;
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>& Matrix3<Type>::operator*=( const Matrix3<Other>& rhs )
// \brief Multiplication assignment operator for the multiplication between two matrices
// \brief (\f$ A*=B \f$).
//
// \param rhs The right-hand-side matrix for the multiplication.
// \return Reference to the matrix.
*/
template< typename Type >
template< typename Other >
inline Matrix3<Type>& Matrix3<Type>::operator*=( const Matrix3<Other>& rhs )
{
   // Creating a temporary due to data dependencies
   Matrix3 tmp( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[3] + v_[2]*rhs.v_[6],
                v_[0]*rhs.v_[1] + v_[1]*rhs.v_[4] + v_[2]*rhs.v_[7],
                v_[0]*rhs.v_[2] + v_[1]*rhs.v_[5] + v_[2]*rhs.v_[8],
                v_[3]*rhs.v_[0] + v_[4]*rhs.v_[3] + v_[5]*rhs.v_[6],
                v_[3]*rhs.v_[1] + v_[4]*rhs.v_[4] + v_[5]*rhs.v_[7],
                v_[3]*rhs.v_[2] + v_[4]*rhs.v_[5] + v_[5]*rhs.v_[8],
                v_[6]*rhs.v_[0] + v_[7]*rhs.v_[3] + v_[8]*rhs.v_[6],
                v_[6]*rhs.v_[1] + v_[7]*rhs.v_[4] + v_[8]*rhs.v_[7],
                v_[6]*rhs.v_[2] + v_[7]*rhs.v_[5] + v_[8]*rhs.v_[8] );

   return this->operator=( tmp );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<HIGH> Matrix3<Type>::operator+( const Matrix3<Other>& rhs ) const
// \brief Addition operator for the addition of two matrices (\f$ A=B+C \f$).
//
// \param rhs The right-hand-side matrix to be added to the matrix.
// \return The sum of the two matrices.
*/
template< typename Type >
template< typename Other >
inline const Matrix3<HIGH> Matrix3<Type>::operator+( const Matrix3<Other>& rhs ) const
{
   return Matrix3<HIGH>( v_[0] + rhs.v_[0],
                         v_[1] + rhs.v_[1],
                         v_[2] + rhs.v_[2],
                         v_[3] + rhs.v_[3],
                         v_[4] + rhs.v_[4],
                         v_[5] + rhs.v_[5],
                         v_[6] + rhs.v_[6],
                         v_[7] + rhs.v_[7],
                         v_[8] + rhs.v_[8] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<HIGH> Matrix3<Type>::operator-( const Matrix3<Other>& rhs ) const
// \brief Subtraction operator for the subtraction of two matrices (\f$ A=B-C \f$).
//
// \param rhs The right-hand-side matrix to be subtracted from the matrix.
// \return The difference of the two matrices.
*/
template< typename Type >
template< typename Other >
inline const Matrix3<HIGH> Matrix3<Type>::operator-( const Matrix3<Other>& rhs ) const
{
   return Matrix3<HIGH>( v_[0] - rhs.v_[0],
                         v_[1] - rhs.v_[1],
                         v_[2] - rhs.v_[2],
                         v_[3] - rhs.v_[3],
                         v_[4] - rhs.v_[4],
                         v_[5] - rhs.v_[5],
                         v_[6] - rhs.v_[6],
                         v_[7] - rhs.v_[7],
                         v_[8] - rhs.v_[8] );
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn const Matrix3<Type> operator-( const Matrix3<Type>& rhs )
// \brief Negation operator for the negation of one matrix.
//
// \param rhs The right-hand-side matrix of the operator.
// \return The negative matrix.
*/
template< typename Type >
inline const Matrix3<Type> operator-( const Matrix3<Type>& rhs )
{
   return Matrix3<Type>( -rhs[0],
                         -rhs[1],
                         -rhs[2],
                         -rhs[3],
                         -rhs[4],
                         -rhs[5],
                         -rhs[6],
                         -rhs[7],
                         -rhs[8] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<HIGH> Matrix3<Type>::operator*( Other rhs ) const
// \brief Multiplication operator for the multiplication of a matrix and a scalar value
// \brief (\f$ A=B*s \f$).
//
// \param rhs The right-hand-side scalar value for the multiplication.
// \return The scaled result matrix.
*/
template< typename Type >
template< typename Other >
inline typename std::enable_if< std::is_arithmetic<Other>::value ,const Matrix3<HIGH> >::type Matrix3<Type>::operator*( Other rhs ) const
{
   return Matrix3<HIGH>( v_[0]*rhs, v_[1]*rhs, v_[2]*rhs,
                         v_[3]*rhs, v_[4]*rhs, v_[5]*rhs,
                         v_[6]*rhs, v_[7]*rhs, v_[8]*rhs );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector3<HIGH> Matrix3<Type>::operator*( const Vector3<Other>& rhs ) const
// \brief Multiplication operator for the multiplication of a matrix and a vector
// \brief (\f$ \vec{a}=B*\vec{c} \f$).
//
// \param rhs The right-hand-side vector for the multiplication.
// \return The resulting vector.
*/
template< typename Type >
template< typename Other >
inline const Vector3<HIGH> Matrix3<Type>::operator*( const Vector3<Other>& rhs ) const
{
   return Vector3<HIGH>( v_[0]*rhs[0] + v_[1]*rhs[1] + v_[2]*rhs[2],
                         v_[3]*rhs[0] + v_[4]*rhs[1] + v_[5]*rhs[2],
                         v_[6]*rhs[0] + v_[7]*rhs[1] + v_[8]*rhs[2] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<HIGH> Matrix3<Type>::operator*( const Matrix3<Other>& rhs ) const
// \brief Multiplication operator for the multiplication of two matrices (\f$ A=B*C \f$).
//
// \param rhs The right-hand-side matrix for the multiplication.
// \return The resulting matrix.
*/
template< typename Type >
template< typename Other >
inline const Matrix3<HIGH> Matrix3<Type>::operator*( const Matrix3<Other>& rhs ) const
{
   return Matrix3<HIGH>( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[3] + v_[2]*rhs.v_[6],
                         v_[0]*rhs.v_[1] + v_[1]*rhs.v_[4] + v_[2]*rhs.v_[7],
                         v_[0]*rhs.v_[2] + v_[1]*rhs.v_[5] + v_[2]*rhs.v_[8],
                         v_[3]*rhs.v_[0] + v_[4]*rhs.v_[3] + v_[5]*rhs.v_[6],
                         v_[3]*rhs.v_[1] + v_[4]*rhs.v_[4] + v_[5]*rhs.v_[7],
                         v_[3]*rhs.v_[2] + v_[4]*rhs.v_[5] + v_[5]*rhs.v_[8],
                         v_[6]*rhs.v_[0] + v_[7]*rhs.v_[3] + v_[8]*rhs.v_[6],
                         v_[6]*rhs.v_[1] + v_[7]*rhs.v_[4] + v_[8]*rhs.v_[7],
                         v_[6]*rhs.v_[2] + v_[7]*rhs.v_[5] + v_[8]*rhs.v_[8] );
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Type Matrix3<Type>::getDeterminant() const
// \brief Calculation of the determinant of the matrix.
//
// \return The determinant of the matrix.
*/
template< typename Type >
inline Type Matrix3<Type>::getDeterminant() const
{
   return v_[0]*v_[4]*v_[8] + v_[1]*v_[5]*v_[6] + v_[2]*v_[3]*v_[7] -
          v_[6]*v_[4]*v_[2] - v_[7]*v_[5]*v_[0] - v_[8]*v_[3]*v_[1];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>& Matrix3<Type>::transpose()
// \brief Transposing the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type >
inline Matrix3<Type>& Matrix3<Type>::transpose()
{
   std::swap( v_[1], v_[3] );
   std::swap( v_[2], v_[6] );
   std::swap( v_[5], v_[7] );
   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<Type> Matrix3<Type>::getTranspose() const
// \brief Calculation of the transpose of the matrix.
//
// \return The transpose of the matrix.
*/
template< typename Type >
inline const Matrix3<Type> Matrix3<Type>::getTranspose() const
{
   return Matrix3( v_[0], v_[3], v_[6], v_[1], v_[4], v_[7], v_[2], v_[5], v_[8] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Matrix3<Type>& Matrix3<Type>::invert()
// \brief Inverting the matrix.
//
// \return Reference to the inverted matrix.
//
// The calculation is performed with the matrix inversion by Cramer. This function is only
// defined for matrices of floating point type. The attempt to use this function with matrices
// of integral data types will result in a compile time error.
*/
template< typename Type >
inline Matrix3<Type>& Matrix3<Type>::invert()
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   Type det = v_[0] * ( ( v_[4] * v_[8] ) - ( v_[7] * v_[5] ) ) +
              v_[1] * ( ( v_[5] * v_[6] ) - ( v_[8] * v_[3] ) ) +
              v_[2] * ( ( v_[3] * v_[7] ) - ( v_[4] * v_[6] ) );

   WALBERLA_ASSERT_NOT_IDENTICAL( det, Type(0), "Matrix is singular and cannot be inverted" );

   det = Type(1) / det;

   // Creating a temporary due to data dependencies
   return *this = Matrix3( det * ( ( v_[4]*v_[8] ) - ( v_[5]*v_[7] ) ),
                           det * ( ( v_[7]*v_[2] ) - ( v_[8]*v_[1] ) ),
                           det * ( ( v_[1]*v_[5] ) - ( v_[2]*v_[4] ) ),
                           det * ( ( v_[5]*v_[6] ) - ( v_[3]*v_[8] ) ),
                           det * ( ( v_[8]*v_[0] ) - ( v_[6]*v_[2] ) ),
                           det * ( ( v_[2]*v_[3] ) - ( v_[0]*v_[5] ) ),
                           det * ( ( v_[3]*v_[7] ) - ( v_[4]*v_[6] ) ),
                           det * ( ( v_[6]*v_[1] ) - ( v_[7]*v_[0] ) ),
                           det * ( ( v_[0]*v_[4] ) - ( v_[1]*v_[3] ) ) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<Type> Matrix3<Type>::getInverse() const
// \brief Calculation of the inverse of the matrix.
//
// \return The inverse of the matrix.
//
// The calculation is performed with the matrix inversion by Cramer. This function is only
// defined for matrices of floating point type. The attempt to use this function with matrices
// of integral data types will result in a compile time error.
*/
template< typename Type >
inline const Matrix3<Type> Matrix3<Type>::getInverse() const
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   Type det = v_[0] * ( ( v_[4] * v_[8] ) - ( v_[7] * v_[5] ) ) +
              v_[1] * ( ( v_[5] * v_[6] ) - ( v_[8] * v_[3] ) ) +
              v_[2] * ( ( v_[3] * v_[7] ) - ( v_[4] * v_[6] ) );

   WALBERLA_ASSERT_NOT_IDENTICAL( det, Type(0), "Matrix is singular and cannot be inverted" );

   det = Type(1) / det;

   return Matrix3( det * ( ( v_[4]*v_[8] ) - ( v_[5]*v_[7] ) ),
                   det * ( ( v_[7]*v_[2] ) - ( v_[8]*v_[1] ) ),
                   det * ( ( v_[1]*v_[5] ) - ( v_[2]*v_[4] ) ),
                   det * ( ( v_[5]*v_[6] ) - ( v_[3]*v_[8] ) ),
                   det * ( ( v_[8]*v_[0] ) - ( v_[6]*v_[2] ) ),
                   det * ( ( v_[2]*v_[3] ) - ( v_[0]*v_[5] ) ),
                   det * ( ( v_[3]*v_[7] ) - ( v_[4]*v_[6] ) ),
                   det * ( ( v_[6]*v_[1] ) - ( v_[7]*v_[0] ) ),
                   det * ( ( v_[0]*v_[4] ) - ( v_[1]*v_[3] ) ) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector3<HIGH> Matrix3<Type>::multTranspose( const Vector3<Other>& rhs ) const
// \brief Multiplication of the transpose of the matrix and a vector (\f$ \vec{a}=B^T*\vec{c} \f$).
//
// \param rhs The right-hand-side vector for the multiplication.
// \return The resulting vector.
*/
template< typename Type >
template< typename Other >
inline const Vector3<HIGH> Matrix3<Type>::multTranspose( const Vector3<Other>& rhs ) const
{
   return Vector3<HIGH>( v_[0]*rhs[0] + v_[3]*rhs[1] + v_[6]*rhs[2],
                         v_[1]*rhs[0] + v_[4]*rhs[1] + v_[7]*rhs[2],
                         v_[2]*rhs[0] + v_[5]*rhs[1] + v_[8]*rhs[2] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<HIGH> Matrix3<Type>::rotate( const Matrix3<Other>& m ) const
// \brief Rotation of a matrix M (\f$ ROT=R*M*R^{-1} \f$).
//
// \param m The matrix to be rotated.
// \return The rotated matrix.
//
// This function is only defined for matrices of floating point type. The attempt to use this
// function with matrices of integral data type will result in a compile time error.
*/
template< typename Type >
template< typename Other >
inline const Matrix3<HIGH> Matrix3<Type>::rotate( const Matrix3<Other>& m ) const
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer  );
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Other>::is_integer );

   //--Multiplication in two steps (number of FLOP = 90, 1 additional temporary matrix)------------

   // Precalculation of tmp = m * R(-1)
   const Matrix3<HIGH> tmp( m.v_[0]*v_[0] + m.v_[1]*v_[1] + m.v_[2]*v_[2],
                            m.v_[0]*v_[3] + m.v_[1]*v_[4] + m.v_[2]*v_[5],
                            m.v_[0]*v_[6] + m.v_[1]*v_[7] + m.v_[2]*v_[8],
                            m.v_[3]*v_[0] + m.v_[4]*v_[1] + m.v_[5]*v_[2],
                            m.v_[3]*v_[3] + m.v_[4]*v_[4] + m.v_[5]*v_[5],
                            m.v_[3]*v_[6] + m.v_[4]*v_[7] + m.v_[5]*v_[8],
                            m.v_[6]*v_[0] + m.v_[7]*v_[1] + m.v_[8]*v_[2],
                            m.v_[6]*v_[3] + m.v_[7]*v_[4] + m.v_[8]*v_[5],
                            m.v_[6]*v_[6] + m.v_[7]*v_[7] + m.v_[8]*v_[8] );

   // Calculating ROT = R * tmp
   return Matrix3<HIGH>( v_[0]*tmp.v_[0] + v_[1]*tmp.v_[3] + v_[2]*tmp.v_[6],
                         v_[0]*tmp.v_[1] + v_[1]*tmp.v_[4] + v_[2]*tmp.v_[7],
                         v_[0]*tmp.v_[2] + v_[1]*tmp.v_[5] + v_[2]*tmp.v_[8],
                         v_[3]*tmp.v_[0] + v_[4]*tmp.v_[3] + v_[5]*tmp.v_[6],
                         v_[3]*tmp.v_[1] + v_[4]*tmp.v_[4] + v_[5]*tmp.v_[7],
                         v_[3]*tmp.v_[2] + v_[4]*tmp.v_[5] + v_[5]*tmp.v_[8],
                         v_[6]*tmp.v_[0] + v_[7]*tmp.v_[3] + v_[8]*tmp.v_[6],
                         v_[6]*tmp.v_[1] + v_[7]*tmp.v_[4] + v_[8]*tmp.v_[7],
                         v_[6]*tmp.v_[2] + v_[7]*tmp.v_[5] + v_[8]*tmp.v_[8] );

   //--Multiplication in one step (number of FLOP = 180, no additional temporary matrix)-----------
   /*
   return Matrix3<HIGH>( m.v_[0]*v_[0]*v_[0] + m.v_[4]*v_[1]*v_[1] + m.v_[8]*v_[2]*v_[2] + v_[0]*v_[1]*( m.v_[1]+m.v_[3] ) + v_[0]*v_[2]*( m.v_[2]+m.v_[6] ) + v_[1]*v_[2]*( m.v_[5]+m.v_[7] ),
                         v_[0]*( m.v_[0]*v_[3] + m.v_[1]*v_[4] + m.v_[2]*v_[5] ) + v_[1]*( m.v_[3]*v_[3] + m.v_[4]*v_[4] + m.v_[5]*v_[5] ) + v_[2]*( m.v_[6]*v_[3] + m.v_[7]*v_[4] + m.v_[8]*v_[5] ),
                         v_[0]*( m.v_[0]*v_[6] + m.v_[1]*v_[7] + m.v_[2]*v_[8] ) + v_[1]*( m.v_[3]*v_[6] + m.v_[4]*v_[7] + m.v_[5]*v_[8] ) + v_[2]*( m.v_[6]*v_[6] + m.v_[7]*v_[7] + m.v_[8]*v_[8] ),
                         v_[3]*( m.v_[0]*v_[0] + m.v_[1]*v_[1] + m.v_[2]*v_[2] ) + v_[4]*( m.v_[3]*v_[0] + m.v_[4]*v_[1] + m.v_[5]*v_[2] ) + v_[5]*( m.v_[6]*v_[0] + m.v_[7]*v_[1] + m.v_[8]*v_[2] ),
                         m.v_[0]*v_[3]*v_[3] + m.v_[4]*v_[4]*v_[4] + m.v_[8]*v_[5]*v_[5] + v_[3]*v_[4]*( m.v_[1]+m.v_[3] ) + v_[3]*v_[5]*( m.v_[2]+m.v_[6] ) + v_[4]*v_[5]*( m.v_[5]+m.v_[7] ),
                         v_[3]*( m.v_[0]*v_[6] + m.v_[1]*v_[7] + m.v_[2]*v_[8] ) + v_[4]*( m.v_[3]*v_[6] + m.v_[4]*v_[7] + m.v_[5]*v_[8] ) + v_[5]*( m.v_[6]*v_[6] + m.v_[7]*v_[7] + m.v_[8]*v_[8] ),
                         v_[6]*( m.v_[0]*v_[0] + m.v_[1]*v_[1] + m.v_[2]*v_[2] ) + v_[7]*( m.v_[3]*v_[0] + m.v_[4]*v_[1] + m.v_[5]*v_[2] ) + v_[8]*( m.v_[6]*v_[0] + m.v_[7]*v_[1] + m.v_[8]*v_[2] ),
                         v_[6]*( m.v_[0]*v_[3] + m.v_[1]*v_[4] + m.v_[2]*v_[5] ) + v_[7]*( m.v_[3]*v_[3] + m.v_[4]*v_[4] + m.v_[5]*v_[5] ) + v_[8]*( m.v_[6]*v_[3] + m.v_[7]*v_[4] + m.v_[8]*v_[5] ),
                         m.v_[0]*v_[6]*v_[6] + m.v_[4]*v_[7]*v_[7] + m.v_[8]*v_[8]*v_[8] + v_[6]*v_[7]*( m.v_[1]+m.v_[3] ) + v_[6]*v_[8]*( m.v_[2]+m.v_[6] ) + v_[7]*v_[8]*( m.v_[5]+m.v_[7] ) );
   */
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<HIGH> Matrix3<Type>::diagRotate( const Matrix3<Other>& m ) const
// \brief Rotation of a diagonal matrix M (\f$ ROT=R*M*R^{-1} \f$).
//
// \param m The diagonal matrix to be rotated.
// \return The rotated matrix.
//
// The DiagRotate function is a special case of the Rotate function. The matrix is assumed to
// be a diagonal matrix, which reduces the number of floating point operations of the rotation.
// This function is only defined for matrices of floating point type. The attempt to use this
// function with matrices of integral data type will result in a compile time error.
*/
template< typename Type >
template< typename Other >
inline const Matrix3<HIGH> Matrix3<Type>::diagRotate( const Matrix3<Other>& m ) const
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer  );
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Other>::is_integer );

   // Precalculating tmp = m * R(-1)
   const Matrix3<HIGH> tmp( m.v_[0]*v_[0], m.v_[0]*v_[3], m.v_[0]*v_[6],
                            m.v_[4]*v_[1], m.v_[4]*v_[4], m.v_[4]*v_[7],
                            m.v_[8]*v_[2], m.v_[8]*v_[5], m.v_[8]*v_[8] );

   // Calculating ROT = R * tmp
   return Matrix3<HIGH>( v_[0]*tmp.v_[0] + v_[1]*tmp.v_[3] + v_[2]*tmp.v_[6],
                         v_[0]*tmp.v_[1] + v_[1]*tmp.v_[4] + v_[2]*tmp.v_[7],
                         v_[0]*tmp.v_[2] + v_[1]*tmp.v_[5] + v_[2]*tmp.v_[8],
                         v_[3]*tmp.v_[0] + v_[4]*tmp.v_[3] + v_[5]*tmp.v_[6],
                         v_[3]*tmp.v_[1] + v_[4]*tmp.v_[4] + v_[5]*tmp.v_[7],
                         v_[3]*tmp.v_[2] + v_[4]*tmp.v_[5] + v_[5]*tmp.v_[8],
                         v_[6]*tmp.v_[0] + v_[7]*tmp.v_[3] + v_[8]*tmp.v_[6],
                         v_[6]*tmp.v_[1] + v_[7]*tmp.v_[4] + v_[8]*tmp.v_[7],
                         v_[6]*tmp.v_[2] + v_[7]*tmp.v_[5] + v_[8]*tmp.v_[8] );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Matrix3<Type>::isSingular() const
// \brief Singularity check for the matrix (det=0).
//
// \return \a true if the matrix is singular, \a false if not.
*/
template< typename Type >
inline bool Matrix3<Type>::isSingular() const
{
   if( equal( getDeterminant(), Type(0) ) )
      return true;
   else return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Matrix3<Type>::isSymmetric() const
// \brief Checks if the matrix is symmetric.
//
// \return \a true if the matrix is symmetric, \a false if not.
*/
template< typename Type >
inline bool Matrix3<Type>::isSymmetric() const
{
   if( !equal( v_[1], v_[3] ) || !equal( v_[2], v_[6] ) || !equal( v_[5], v_[7] ) )
      return false;
   else return true;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn bool Matrix3<Type>::isZero() const
// \brief Checks if all matrix entries are zero.
//
// \return \a true if all entries are zero otherwise \a false.
*/
template< typename Type >
inline bool Matrix3<Type>::isZero() const
{
   if( equal( v_[0], Type(0) ) &&
       equal( v_[1], Type(0) ) &&
       equal( v_[2], Type(0) ) &&
       equal( v_[3], Type(0) ) &&
       equal( v_[4], Type(0) ) &&
       equal( v_[5], Type(0) ) &&
       equal( v_[6], Type(0) ) &&
       equal( v_[7], Type(0) ) &&
       equal( v_[8], Type(0) ) )
      return true;
   else return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<Type> Matrix3<Type>::getCholesky() const
// \brief Cholesky decomposition of the matrix (\f$ A = LR \f$).
//
// \return The decomposed matrix \f$ L \f$.
//
// This function is only defined for matrices of floating point type. The attempt to use this
// function with matrices of integral data type will result in a compile time error.
*/
template< typename Type >
inline const Matrix3<Type> Matrix3<Type>::getCholesky() const
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   Matrix3 tmp( Type(0) );

   WALBERLA_ASSERT( !isSingular() ,"Matrix is singular and cannot be decomposed" );

   Type sum( 0 );

   for( int k=0; k<3; ++k ) {
      for( int p=0; p<k; ++p ) {
         sum += tmp.v_[k*3+p] * tmp.v_[k*3+p];
      }
      tmp.v_[k*4] = std::sqrt( v_[k*4]-sum );
      sum = Type(0);
      for( int i=(k+1); i<3; ++i ) {
         for( int p=0; p<k; ++p ) {
            sum += tmp.v_[k*3+p] * tmp.v_[i*3+p];
         }
         tmp.v_[i*3+k] = (v_[i*3+k]-sum) / tmp.v_[k*4];
         sum = Type(0);
      }
   }

   return tmp;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector3<HIGH> Matrix3<Type>::solve( const Vector3<Other> &rhs ) const
// \brief Solving the linear system of equations \f$A*x=b\f$ with the decomposed matrix
// \brief \f$LR = A\f$.
//
// \param rhs The right-hand-side of the linear system of equations.
// \return The solution of the linear system of equations.
//
// This function is only defined for matrices of floating point type. The attempt to use this
// function with matrices of integral data type will result in a compile time error.
*/
template< typename Type >
template< typename Other >
inline const Vector3<HIGH> Matrix3<Type>::solve( const Vector3<Other> &rhs ) const
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer  );
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Other>::is_integer );

   Vector3<HIGH> tmp1, tmp2;
   HIGH sum;

   // Solving the equation L*y = b
   for( int i=0; i<3; ++i ) {
      sum = rhs[i];
      for( int j=0; j<i; ++j ) {
         sum -= v_[i*3+j] * tmp1[j];
      }
      tmp1[i] = sum / v_[i*4];
   }

   // Solving the equation R*x = y
   for( int i=2; i>=0; --i ) {
      sum = tmp1[i];
      for( int j=2; j>i ; --j ) {
         sum -= v_[j*3+i] * tmp2[j];
      }
      tmp2[i] = sum / v_[i*4];
   }

   return tmp2;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Type Matrix3<Type>::trace() const
// \brief Computes the trace of the matrix (sum of diagonal elements).
//
// \return The trace of the matrix.
*/
template< typename Type >
inline Type Matrix3<Type>::trace() const
{
   return v_[0] + v_[4] + v_[8];
}
//**********************************************************************************************************************



//======================================================================================================================
//
//  EULER ROTATIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn const Vector3<Type> Matrix3<Type>::getEulerAnglesXYZ() const
// \brief Calculation of the Euler angles (in radian measure).
//
// \return The Euler angles for a rotation order of x, y, z.
//
// The Euler angles are calculated for a rotation order of x-, y- and z-axis. This function
// is only defined for matrices of floating point type. The attempt to use this function with
// matrices of integral data type will result in a compile time error.
*/
template< typename Type >
inline const Vector3<Type> Matrix3<Type>::getEulerAnglesXYZ() const
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   const Type cy( std::sqrt( v_[0]*v_[0] + v_[3]*v_[3] ) );

   if( cy > Type(16)*std::numeric_limits<Type>::epsilon() ) {
      return Vector3<Type>( std::atan2( v_[7], v_[8] ), std::atan2( -v_[6], cy ), std::atan2( v_[3], v_[0] ) );
   }
   else {
      return Vector3<Type>( std::atan2( -v_[5], v_[4] ), std::atan2( -v_[6], cy ), Type(0) );
   }
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Vector3<Type> Matrix3<Type>::getEulerAngles( EulerRotation order ) const
// \brief Calculation of the Euler angles for a specific rotation order.
//
// \param order The specific rotation order.
// \return The specific Euler angles.
//
// This function is only defined for matrices of floating point type. The attempt to use this
// function with matrices of integral data type will result in a compile time error.
*/
template< typename Type >
const Vector3<Type> Matrix3<Type>::getEulerAngles( EulerRotation order ) const
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );

   static const uint_t eulSafe[4] = { 0, 1, 2, 0 };
   static const uint_t eulNext[4] = { 1, 2, 0, 1 };

   Vector3<Type> ea;

   // Unpacking the euler order
   const uint_t frame( order&1 );
   const uint_t repetition( (order&2)>>1 );
   const uint_t parity( (order&4)>>2 );
   const uint_t i( eulSafe[(order&24)>>3] );
   const uint_t j( eulNext[i+parity] );
   const uint_t k( eulNext[i+1-parity] );

   // Treatment of rotations with repetition
   if( repetition ) {
      const Type sy( std::sqrt( v_[i*3+j]*v_[i*3+j] + v_[i*3+k]*v_[i*3+k] ) );
      if( sy > Type(16)*std::numeric_limits<Type>::epsilon() ) {
         ea[0] = std::atan2( v_[i*3+j], v_[i*3+k] );
         ea[1] = std::atan2( sy, v_[i*3+i] );
         ea[2] = std::atan2( v_[j*3+i], -v_[k*3+i] );
      }
      else {
         ea[0] = std::atan2( -v_[j*3+k], v_[j*3+j] );
         ea[1] = std::atan2( sy, v_[i*3+i] );
         ea[2] = Type(0);
      }
   }

   // Treatment of rotations without repetition
   else {
      const Type cy( std::sqrt( v_[i*3+i]*v_[i*3+i] + v_[j*3+i]*v_[j*3+i] ) );
      if( cy > Type(16)*std::numeric_limits<Type>::epsilon() ) {
         ea[0] = std::atan2( v_[k*3+j], v_[k*3+k] );
         ea[1] = std::atan2( -v_[k*3+i], cy );
         ea[2] = std::atan2( v_[j*3+i], v_[i*3+i] );
      }
      else {
         ea[0] = std::atan2( -v_[j*3+k], v_[j*3+j] );
         ea[1] = std::atan2( -v_[k*3+i], cy );
         ea[2] = Type(0);
      }
   }

   // Treatment of an odd partity
   if( parity ) {
      ea[0] = -ea[0];
      ea[1] = -ea[1];
      ea[2] = -ea[2];
   }

   // Treatment of a rotating frame
   if( frame ) {
      real_t tmp = ea[0];
      ea[0] = ea[2];
      ea[2] = tmp;
   }

   return ea;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  GLOBAL OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\name Matrix3 operators */
//@{
//template< typename Type, typename Other >
//inline const Matrix3<HIGH> operator*( Other scalar, const Matrix3<Type>& matrix );

template< typename Type >
std::ostream& operator<<( std::ostream& os, const Matrix3<Type>& m );

template< typename Type >
inline bool isnan( const Matrix3<Type>& m );

template< typename Type >
inline const Matrix3<Type> abs( const Matrix3<Type>& m );

template< typename Type >
inline const Matrix3<Type> fabs( const Matrix3<Type>& m );
//@}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<Type> operator*( Other scalar, const Matrix3<Type>& matrix )
// \brief Multiplication operator for the multiplication of a scalar value and a matrix
// \brief (\f$ A=s*B \f$).
//
// \param scalar The left-hand-side scalar value for the multiplication.
// \param matrix The right-hand-side matrix for the multiplication.
// \return The scaled result matrix.
*/
//template< typename Type, typename Other >
//inline const Matrix3<HIGH> operator*( Other scalar, const Matrix3<Type>& matrix )
//{
//   static_assert( ! std::is_scalar<Other>::value, "Only scalar types allowed" );
//   return matrix*scalar;
//}
//**********************************************************************************************************************

//*************************************************************************************************
/*!\brief Cross product (outer product) between a matrix and a vector (\f$ R = M \cdot r^{\times} \f$).
 *
 * \param mat The left-hand side matrix for the skew-symmetric cross product.
 * \param vec The right-hand side vector for the skew-symmetric cross product.
 * \return The matrix \f$ M \cdot r^{\times} \f$.
 *
 * This operator multiplies the matrix with the skew-symmetric matrix resulting from the
 * vector, which corresponds to a cross product between the columns vectors of the matrix
 * and the right-hand side vector:

                           \f[
                           M \cdot r^{\times} =

                           \left(\begin{array}{*{3}{c}}
                           m_{00} & m_{01} & m_{02} \\
                           m_{10} & m_{11} & m_{12} \\
                           m_{20} & m_{21} & m_{22} \\
                           \end{array}\right) \cdot

                           \left(\begin{array}{*{1}{c}}
                           r_0 \\
                           r_1 \\
                           r_2 \\
                           \end{array}\right)^{\times} =

                           \left(\begin{array}{*{3}{c}}
                           m_{00} & m_{01} & m_{02} \\
                           m_{10} & m_{11} & m_{12} \\
                           m_{20} & m_{21} & m_{22} \\
                           \end{array}\right) \cdot

                           \left(\begin{array}{*{3}{c}}
                           0    & -r_2 & r_1  \\
                           r_2  & 0    & -r_0 \\
                           -r_1 & r_0  & 0    \\
                           \end{array}\right)
                           \f]

 * The operator returns a matrix of the higher-order data type of the two involved data types
 * \a T1 and \a T2. In case \a T1 and \a T2 match, the operator works for any data type as long
 * as the data type has a multiplication and subtraction operator. In case \a T1 and \a T2 differ,
 * the operator only works for data types supported by the MathTrait class template.
 */
template< typename Type >
inline const Matrix3< Type > skewSymCrossProduct( const Matrix3<Type>& mat, const Vector3<Type>& vec )
{
   return Matrix3<Type>( mat[1] * vec[2] - mat[2] * vec[1],
                         mat[2] * vec[0] - mat[0] * vec[2],
                         mat[0] * vec[1] - mat[1] * vec[0],
                         mat[4] * vec[2] - mat[5] * vec[1],
                         mat[5] * vec[0] - mat[3] * vec[2],
                         mat[3] * vec[1] - mat[4] * vec[0],
                         mat[7] * vec[2] - mat[8] * vec[1],
                         mat[8] * vec[0] - mat[6] * vec[2],
                         mat[6] * vec[1] - mat[7] * vec[0] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Cross product (outer product) between a vector and a matrix (\f$ R = r^{\times} \cdot M \f$).
 *
 * \param vec The left-hand side vector for the skew-symmetric product.
 * \param mat The right-hand side matrix for the skew-symmetric product.
 * \return The matrix \f$ r^{\times} \cdot M \f$.
 *
 * This operator multiplies the skew-symmetric matrix resulting from the vector with the
 * right-hand side matrix, which corresponds to a cross product between the vector and the
 * columns vectors of the matrix:

                           \f[
                           r^{\times} \cdot M =

                           \left(\begin{array}{*{1}{c}}
                           r_0 \\
                           r_1 \\
                           r_2 \\
                           \end{array}\right)^{\times} \cdot

                           \left(\begin{array}{*{3}{c}}
                           m_{00} & m_{01} & m_{02} \\
                           m_{10} & m_{11} & m_{12} \\
                           m_{20} & m_{21} & m_{22} \\
                           \end{array}\right) =

                           \left(\begin{array}{*{3}{c}}
                           0    & -r_2 & r_1  \\
                           r_2  & 0    & -r_0 \\
                           -r_1 & r_0  & 0    \\
                           \end{array}\right) \cdot

                           \left(\begin{array}{*{3}{c}}
                           m_{00} & m_{01} & m_{02} \\
                           m_{10} & m_{11} & m_{12} \\
                           m_{20} & m_{21} & m_{22} \\
                           \end{array}\right)
                           \f]

 * The operator returns a matrix of the higher-order data type of the two involved data types
 * \a T1 and \a T2. In case \a T1 and \a T2 match, the operator works for any data type as long
 * as the data type has a multiplication and subtraction operator. In case \a T1 and \a T2 differ,
 * the operator only works for data types supported by the MathTrait class template.
 */
template< typename Type >
inline const Matrix3<Type> skewSymCrossProduct( const Vector3<Type>& vec, const Matrix3<Type>& mat )
{
   return Matrix3<Type>( vec[1] * mat[6] - vec[2] * mat[3],
                         vec[1] * mat[7] - vec[2] * mat[4],
                         vec[1] * mat[8] - vec[2] * mat[5],
                         vec[2] * mat[0] - vec[0] * mat[6],
                         vec[2] * mat[1] - vec[0] * mat[7],
                         vec[2] * mat[2] - vec[0] * mat[8],
                         vec[0] * mat[3] - vec[1] * mat[0],
                         vec[0] * mat[4] - vec[1] * mat[1],
                         vec[0] * mat[5] - vec[1] * mat[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Dyadic product of two vectors (\f$ M = u \otimes v \f$).
 *
 * \param vec1 The first vector argument.
 * \param vec2 The second vector argument.
 * \return The matrix \f$ u \otimes v \f$.
 */
template< typename Type >
inline const Matrix3<Type> dyadicProduct( const Vector3<Type>& vec1, const Vector3<Type>& vec2 )
{
   return Matrix3<Type>( vec1[0] * vec2[0], vec1[0] * vec2[1], vec1[0] * vec2[2],
                         vec1[1] * vec2[0], vec1[1] * vec2[1], vec1[1] * vec2[2],
                         vec1[2] * vec2[0], vec1[2] * vec2[1], vec1[2] * vec2[2]);
}
//*************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::ostream& operator<<( std::ostream& os, const Matrix3<Type>& m )
// \brief Global output operator for 3x3 matrices.
//
// \param os Reference to the output stream.
// \param m Reference to a constant matrix object.
// \return Reference to the output stream.
*/
template< typename Type >
std::ostream& operator<<( std::ostream& os, const Matrix3<Type>& m )
{
   return os << " ( " << m[0] << " , " << m[1] << " , " << m[2] << " )\n"
             << " ( " << m[3] << " , " << m[4] << " , " << m[5] << " )\n"
             << " ( " << m[6] << " , " << m[7] << " , " << m[8] << " )\n";
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool math::isnan( const Matrix3<Type>& m )
// \brief Checks the given matrix for not-a-number elements.
//
// \param m The matrix to be checked for not-a-number elements.
// \return \a true if at least one element of the matrix is not-a-number, \a false otherwise.
*/
template< typename Type >
inline bool isnan( const Matrix3<Type>& m )
{
   if( math::isnan( m[0] ) || math::isnan( m[1] ) || math::isnan( m[2] ) ||
       math::isnan( m[3] ) || math::isnan( m[4] ) || math::isnan( m[5] ) ||
       math::isnan( m[6] ) || math::isnan( m[7] ) || math::isnan( m[8] ) )
      return true;
   else return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<Type> abs( const Matrix3<Type>& m )
// \brief Returns a matrix containing the absolute values of each single element of \a m.
//
// \param m The integral input matrix.
// \return The absolute value of each single element of \a m.
//
// The \a abs function calculates the absolute value of each element of the input matrix
// \a m. This function can only be applied to matrices of integral data type. For floating
// point matrices, the pe::fabs( const Matrix3& ) function can be used.
*/
template< typename Type >
inline const Matrix3<Type> abs( const Matrix3<Type>& m )
{
   WALBERLA_STATIC_ASSERT( std::numeric_limits<Type>::is_integer );
   return Matrix3<Type>( std::abs(m[0]), std::abs(m[1]), std::abs(m[2]),
                         std::abs(m[3]), std::abs(m[4]), std::abs(m[5]),
                         std::abs(m[6]), std::abs(m[7]), std::abs(m[8]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const Matrix3<Type> fabs( const Matrix3<Type>& m )
// \brief Returns a matrix containing the absolute values of each single element of \a m.
//
// \param m The floating point input matrix.
// \return The absolute value of each single element of \a m.
//
// The \a fabs function calculates the absolute value of each element of the input matrix
// \a m. This function can only be applied to  floating point matrices. For matrices of
// integral data type, the pe::abs( const Matrix3& ) function can be used.
*/
template< typename Type >
inline const Matrix3<Type> fabs( const Matrix3<Type>& m )
{
   WALBERLA_STATIC_ASSERT( !std::numeric_limits<Type>::is_integer );
   return Matrix3<Type>( std::fabs(m[0]), std::fabs(m[1]), std::fabs(m[2]),
                         std::fabs(m[3]), std::fabs(m[4]), std::fabs(m[5]),
                         std::fabs(m[6]), std::fabs(m[7]), std::fabs(m[8]) );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool math::isinf( const Matrix3<Type>& m )
// \brief Checks the given matrix for infinite elements.
//
// \param m The matrix to be checked for infinite elements.
// \return \a true if at least one element of the matrix is infinite, \a false otherwise.
*/
template< typename Type >
inline bool isinf( const Matrix3<Type>& m )
{
   if( math::isinf( m[0] ) || math::isinf( m[1] ) || math::isinf( m[2] ) ||
       math::isinf( m[3] ) || math::isinf( m[4] ) || math::isinf( m[5] ) ||
       math::isinf( m[6] ) || math::isinf( m[7] ) || math::isinf( m[8] ) )
      return true;
   else return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
// \brief Calculates the tensor product of two vectors.
//**********************************************************************************************************************
template< typename T0, typename T1>
Matrix3< typename MathTrait<T0,T1>::High > tensorProduct( Vector3<T0> v0, Vector3<T1> v1 ){
   return Matrix3< typename MathTrait<T0,T1>::High > (
      v0[0] * v1[0], v0[0] * v1[1], v0[0] * v1[2],
      v0[1] * v1[0], v0[1] * v1[1], v0[1] * v1[2],
      v0[2] * v1[0], v0[2] * v1[1], v0[2] * v1[2]);
}
//**********************************************************************************************************************

/**
 * Equivalent to R*A*R.getTranspose().
 */
template< typename Type>
inline Matrix3< Type > transformMatrixRART( const Matrix3<Type>& R, const Matrix3<Type>& A )
{
   const auto r0 = A[0]*R[0]*R[0] + A[1]*R[0]*R[1] + A[2]*R[0]*R[2] + A[3]*R[0]*R[1] + A[4]*R[1]*R[1] + A[5]*R[1]*R[2] + A[6]*R[0]*R[2] + A[7]*R[1]*R[2] + A[8]*R[2]*R[2];
   const auto r1 = A[0]*R[0]*R[3] + A[1]*R[0]*R[4] + A[2]*R[0]*R[5] + A[3]*R[1]*R[3] + A[4]*R[1]*R[4] + A[5]*R[1]*R[5] + A[6]*R[2]*R[3] + A[7]*R[2]*R[4] + A[8]*R[2]*R[5];
   const auto r2 = A[0]*R[0]*R[6] + A[1]*R[0]*R[7] + A[2]*R[0]*R[8] + A[3]*R[1]*R[6] + A[4]*R[1]*R[7] + A[5]*R[1]*R[8] + A[6]*R[2]*R[6] + A[7]*R[2]*R[7] + A[8]*R[2]*R[8];
   const auto r3 = A[0]*R[0]*R[3] + A[1]*R[1]*R[3] + A[2]*R[2]*R[3] + A[3]*R[0]*R[4] + A[4]*R[1]*R[4] + A[5]*R[2]*R[4] + A[6]*R[0]*R[5] + A[7]*R[1]*R[5] + A[8]*R[2]*R[5];
   const auto r4 = A[0]*R[3]*R[3] + A[1]*R[3]*R[4] + A[2]*R[3]*R[5] + A[3]*R[3]*R[4] + A[4]*R[4]*R[4] + A[5]*R[4]*R[5] + A[6]*R[3]*R[5] + A[7]*R[4]*R[5] + A[8]*R[5]*R[5];
   const auto r5 = A[0]*R[3]*R[6] + A[1]*R[3]*R[7] + A[2]*R[3]*R[8] + A[3]*R[4]*R[6] + A[4]*R[4]*R[7] + A[5]*R[4]*R[8] + A[6]*R[5]*R[6] + A[7]*R[5]*R[7] + A[8]*R[5]*R[8];
   const auto r6 = A[0]*R[0]*R[6] + A[1]*R[1]*R[6] + A[2]*R[2]*R[6] + A[3]*R[0]*R[7] + A[4]*R[1]*R[7] + A[5]*R[2]*R[7] + A[6]*R[0]*R[8] + A[7]*R[1]*R[8] + A[8]*R[2]*R[8];
   const auto r7 = A[0]*R[3]*R[6] + A[1]*R[4]*R[6] + A[2]*R[5]*R[6] + A[3]*R[3]*R[7] + A[4]*R[4]*R[7] + A[5]*R[5]*R[7] + A[6]*R[3]*R[8] + A[7]*R[4]*R[8] + A[8]*R[5]*R[8];
   const auto r8 = A[0]*R[6]*R[6] + A[1]*R[6]*R[7] + A[2]*R[6]*R[8] + A[3]*R[6]*R[7] + A[4]*R[7]*R[7] + A[5]*R[7]*R[8] + A[6]*R[6]*R[8] + A[7]*R[7]*R[8] + A[8]*R[8]*R[8];

   return Matrix3<Type>(r0,r1,r2,r3,r4,r5,r6,r7,r8);
}

} // namespace math

using math::Matrix3;

} // namespace walberla

//======================================================================================================================
//
//  Vector Trait Specialization
//
//======================================================================================================================

namespace walberla {

// Specialization of VectorTrait for Matrix3s
template<typename T>
struct VectorTrait< Matrix3<T> >
{
   typedef T OutputType;

   static const uint_t F_SIZE =  9u;
   static T    get( const Matrix3<T> & v, uint_t f )       { return v[f]; }
   static void set(       Matrix3<T> & v, uint_t f, T val) { v[f] = val;  }
};

} // namespace walberla

//======================================================================================================================
//
//  comparison backend for Matrix3<real_t>
//
//======================================================================================================================

namespace walberla {
namespace debug {
namespace check_functions_detail {

template< >
inline bool check_float_equal( const math::Matrix3<real_t> & lhs, const math::Matrix3<real_t> & rhs )
{
   return floatIsEqual( lhs[0], rhs[0] ) && floatIsEqual( lhs[1], rhs[1] ) && floatIsEqual( lhs[2], rhs[2] )
       && floatIsEqual( lhs[3], rhs[3] ) && floatIsEqual( lhs[4], rhs[4] ) && floatIsEqual( lhs[5], rhs[5] )
       && floatIsEqual( lhs[6], rhs[6] ) && floatIsEqual( lhs[7], rhs[7] ) && floatIsEqual( lhs[8], rhs[8] );
}

template< >
inline bool check_float_equal_eps( const math::Matrix3<real_t> & lhs, const math::Matrix3<real_t> & rhs, const real_t epsilon )
{
   return floatIsEqual( lhs[0], rhs[0], epsilon ) && floatIsEqual( lhs[1], rhs[1], epsilon ) && floatIsEqual( lhs[2], rhs[2], epsilon )
       && floatIsEqual( lhs[3], rhs[3], epsilon ) && floatIsEqual( lhs[4], rhs[4], epsilon ) && floatIsEqual( lhs[5], rhs[5], epsilon )
       && floatIsEqual( lhs[6], rhs[6], epsilon ) && floatIsEqual( lhs[7], rhs[7], epsilon ) && floatIsEqual( lhs[8], rhs[8], epsilon );
}

}
}
}


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
                typename MT >  // Element type of matrix
      mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const Matrix3<MT> & m )
      {
         buf.addDebugMarker( "m3" );
         static_assert ( std::is_trivially_copyable< Matrix3<MT> >::value,
                         "type has to be trivially copyable for the memcpy to work correctly" );
         auto pos = buf.forward(sizeof(Matrix3<MT>));
         std::memcpy(pos, &m, sizeof(Matrix3<MT>));
         return buf;
      }

      template< typename T,    // Element type  of RecvBuffer
                typename MT >  // Element type of matrix
      mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, Matrix3<MT> & m )
      {
         buf.readDebugMarker( "m3" );
         static_assert ( std::is_trivially_copyable< Matrix3<MT> >::value,
                         "type has to be trivially copyable for the memcpy to work correctly" );
         auto pos = buf.skip(sizeof(Matrix3<MT>));
         //suppress https://gcc.gnu.org/onlinedocs/gcc/C_002b_002b-Dialect-Options.html#index-Wclass-memaccess
         std::memcpy(static_cast<void*>(&m), pos, sizeof(Matrix3<MT>));
         return buf;
      }

      template<typename VT>
      struct BufferSizeTrait< walberla::math::Matrix3<VT> > {
         static const bool constantSize = true;
         static const uint_t size = 9 * BufferSizeTrait<VT>::size + mpi::BUFFER_DEBUG_OVERHEAD;
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
   struct MPITrait< Matrix3<T> >
   {
      static inline MPI_Datatype type()
      {
         // cannot use mpi::Datatype here because its destructor calls MPI_Type_free and static variables are destroyed after the MPI_Finalize
         static MPI_Datatype datatype;
         static bool initialized = false;

         if( ! initialized ) {
            MPI_Type_contiguous(9, MPITrait<T>::type(), &datatype );
            MPI_Type_commit( &datatype );
            initialized = true;
         }
         return datatype;
      }
   };
} // namespace walberla
