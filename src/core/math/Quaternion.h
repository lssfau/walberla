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
//! \file Quaternion.h
//! \ingroup core
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the implementation of a quaternion
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "core/math/MathTrait.h"
#include "core/math/Utility.h"
#include "core/debug/Debug.h"
#include "core/debug/CheckFunctions.h"
#include "core/DataTypes.h"

#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include <core/logging/Logging.h>

#include <type_traits>

#include <cmath>
#include <istream>
#include <ostream>
#include <limits>

namespace walberla {
namespace math {

//=================================================================================================
//
//  NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< typename Type >       class Matrix3;
template< typename Type >       class Vector3;




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Efficient implementation of a quaternion.
 *
 * Quaternions are a superior way to deal with rotations and orientations. This quaternion
 * consists of 4 statically allocated elements, where the first element represents the real
 * part and the three other elements represent the three imaginary parts. The naming
 * convention of the elements is as following:

                             \f[\left(\begin{array}{*{4}{c}}
                             r & i & j & k \\
                             \end{array}\right)\f]

 * These elements can be accessed directly with the subscript operator. The numbering of the
 * quaternion elements is

                             \f[\left(\begin{array}{*{4}{c}}
                             0 & 1 & 2 & 3 \\
                             \end{array}\right)\f]

 * \b Note: The Quaternion class can only be instantiated for non-cv-qualified floating point
 * types! Therefore the only possible Quaternion instantiations are
 *
 *  - Quaternion<float>
 *  - Quaternion<double>
 *  - Quaternion<long double>
 *
 * The attempt to create a quaternion with an integral data type results in a compile time
 * error.
 */
template< typename Type >  // Data type of the quaternion
class Quaternion
{
   //**Compile time checks*************************************************************************
   /*! \cond internal */
   static_assert(std::is_floating_point<Type>::value, "T has to be floating point!");
   static_assert(!std::is_const<Type>::value, "T has to be non const!");
   static_assert(!std::is_volatile<Type>::value, "T has to be non volatile!");
   /*! \endcond */
   //**********************************************************************************************

public:
   //**Type definitions****************************************************************************
   using ElementType = Type;  //!< Type of the quaternion elements.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Quaternion() = default;

   explicit inline Quaternion( Type r, Type i, Type j, Type k );

   template< typename Axis >
   explicit inline Quaternion( Vector3<Axis> axis, Type angle );

   explicit inline Quaternion( Type xangle, Type yangle, Type zangle );

   template< typename Other >
   explicit inline Quaternion( const Vector3<Other>& euler );

   inline Quaternion( const Quaternion& q ) = default;

   template< typename Other >
   inline Quaternion( const Quaternion<Other>& q );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
                              inline Quaternion& operator= ( const Quaternion& rhs ) = default;
   template< typename Other > inline Quaternion& operator= ( const Quaternion<Other>& rhs );
                              inline Type        operator[]( size_t index ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline Quaternion&                set( Type r, Type i, Type j, Type k );
                              inline void                       reset();
                              inline Quaternion&                invert();
                              inline const Quaternion           getInverse()       const;
                              inline const Matrix3<Type>        toRotationMatrix() const;
                              inline Type                       getAngle()         const;
                              inline const Vector3<Type>        getAxis()          const;
                              inline void                       rotateX( Type angle );
                              inline void                       rotateY( Type angle );
                              inline void                       rotateZ( Type angle );
                              inline void                       swap( Quaternion& q ) /* throw() */;
                              inline const Vector3<Type>        getEulerAnglesXYZ() const;
                              inline Type*                      data()                         {return v_;}
                              inline Type const *               data()                         const {return v_;}
   //@}
   //**********************************************************************************************

   //**Math functions******************************************************************************
   /*!\name Math functions
   //
   // The return type of the math functions depends on the involved data types of the
   // quaternions, matrices and vectors (for further detail see the MathTrait class
   // description).
   */
   //@{
   template< typename Other >
   inline const Vector3< typename MathTrait<Type,Other>::MultType >
      rotate( const Vector3<Other>& v ) const;

   template< typename Other >
   inline const Matrix3< typename MathTrait<Type,Other>::MultType >
      rotate( const Matrix3<Other>& m ) const;

   template< typename Other >
   inline const Matrix3< typename MathTrait<Type,Other>::MultType >
      diagRotate( const Matrix3<Other>& m ) const;

   template< typename Other >
   inline typename MathTrait<Type,Other>::HighType
      calcAngle( const Vector3<Other>& axis ) const;
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   /**
    * The four statically allocated quaternion elements.
    *
    * Access to the quaternion values is gained via the subscript operator.
    * The order of the elements is
    * \f[\left(\begin{array}{*{4}{c}}
    * 0 & 1 & 2 & 3 \\
    * \end{array}\right)\f]
   **/
   Type v_[4] = {Type(1), Type(0), Type(0), Type(0)};
   //@}
   //**********************************************************************************************
};
static_assert( std::is_trivially_copyable<Quaternion<real_t>>::value, "Quaternion<real_t> has to be trivially copyable!");
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Constructor for a direct initialization of all quaternion elements.
 *
 * \param r The initial value for the real part.
 * \param i The initial value for the first imaginary part.
 * \param j The initial value for the second imaginary part.
 * \param k The initial value for the third imaginary part.
 *
 * The initial values for the quaterion have to be chosen such that the length of the
 * quaternion is 1.
 */
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>::Quaternion( Type r, Type i, Type j, Type k )
{
   v_[0] = r; v_[1] = i; v_[2] = j; v_[3] = k;

   WALBERLA_CHECK_FLOAT_EQUAL( r*r + i*i + j*j + k*k, Type(1), "Invalid quaternion parameters:\n" <<
                               "r: " << r << "\n" <<
                               "i: " << i << "\n" <<
                               "j: " << j << "\n" <<
                               "k: " << k);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a quaternion depending on a rotation axis and angle.
 *
 * \param axis The rotation axis.
 * \param angle The rotation angle (radian measure).
 *
 * This constructor creates a quaternion from the rotation axis \a axis and the rotation angle
 * \a angle. \a axis may be an arbitrary, non-zero vector of any length. However, it is allowed
 * to use the zero vector (0,0,0) in combination with an angle of 0. This combination results
 * in a default quaternion

           \f[ \left(\begin{array}{c} 1 \\ 0 \\ 0 \\ 0 \end{array}\right) \f]
 */
template< typename Type >  // Data type of the quaternion
template< typename Axis >  // Data type of the rotation axis
inline Quaternion<Type>::Quaternion( Vector3<Axis> axis, Type angle )
{
   static_assert(std::is_floating_point<Axis>::value, "Axis has to be floating point!" );

   auto axisLength = axis.length();
   if( (floatIsEqual(axisLength, 0)) || (math::equal(std::fabs(angle), real_t(0)))  ) {
      reset();
      return;
   }

   const Type sina( std::sin( angle*Type(0.5) ) );
   const Type cosa( std::cos( angle*Type(0.5) ) );

   auto invAxisLength = real_t(1) / axisLength;
   axis *= invAxisLength;

   v_[0] = cosa;
   v_[1] = sina * axis[0];
   v_[2] = sina * axis[1];
   v_[3] = sina * axis[2];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a quaternion rotated by the Euler angles \a xangle, \a yangle and \a zangle.
 *
 * \param xangle Rotation around the x-axis (radian measure).
 * \param yangle Rotation around the y-axis (radian measure).
 * \param zangle Rotation around the z-axis (radian measure).
 *
 * This constructor creates a quaternion rotated by the Euler angles \a xangle, \a yangle and
 * \a zangle (all in radian measure). The rotations are applied in the order x, y, and z.
 */
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>::Quaternion( Type xangle, Type yangle, Type zangle )
{
   reset();
   rotateX( xangle );
   rotateY( yangle );
   rotateZ( zangle );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a quaternion rotated by the Euler angles \a euler.
 *
 * \param euler 3-dimensional vector of the three rotation angles (radian measure).
 *
 * This constructor creates a quaternion rotated by the Euler angles \a euler (all components
 * in radian measure). The rotations are applied in the order x, y, and z.
 */
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the Euler angle vector
inline Quaternion<Type>::Quaternion( const Vector3<Other>& euler )
{
   reset();
   rotateX( euler[0] );
   rotateY( euler[1] );
   rotateZ( euler[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different Quaternion instances.
 *
 * \param q Quaternion to be copied.
 */
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the foreign quaternion
inline Quaternion<Type>::Quaternion( const Quaternion<Other>& q )
{
   v_[0] = q[0];
   v_[1] = q[1];
   v_[2] = q[2];
   v_[3] = q[3];
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Assignment operator for different Quaternion instances.
 *
 * \param rhs Quaternion to be copied.
 * \return Reference to the assigned quaternion.
 */
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the foreign quaternion
inline Quaternion<Type>& Quaternion<Type>::operator=( const Quaternion<Other>& rhs )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = rhs[0];
   v_[1] = rhs[1];
   v_[2] = rhs[2];
   v_[3] = rhs[3];
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the quaternion elements.
 *
 * \param index Access index. The index has to be in the range \f$[0..3]\f$.
 * \return Copy of the accessed element.
 *
 * When compiled in Debug mode, this operator performs an index check.
 */
template< typename Type >  // Data type of the quaternion
inline Type Quaternion<Type>::operator[]( size_t index ) const
{
   WALBERLA_ASSERT( index < 4, "Invalid quaternion access index" );
   return v_[index];
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Setting the value of the quaternion elements.
 *
 * \param r The value for the real part.
 * \param i The value for the first imaginary part.
 * \param j The value for the second imaginary part.
 * \param k The value for the third imaginary part.
 */
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>& Quaternion<Type>::set( Type r, Type i, Type j, Type k )
{
   WALBERLA_CHECK_FLOAT_EQUAL( std::fabs( r*r + i*i + j*j + k*k ), Type(1),
                               "Invalid quaternion parameters: " << r << ", "<< i << ", "<< j << ", "<< k );
   v_[0] = r;
   v_[1] = i;
   v_[2] = j;
   v_[3] = k;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
 *
 * \return void
 *
 * This function resets the quaternion to the default initial values. The real part of the
 * quaternion is reset to 1, whereas the imaginary parts are reset to 0:

                             \f[\left(\begin{array}{*{4}{c}}
                             1 & 0 & 0 & 0 \\
                             \end{array}\right)\f]
 */
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::reset()
{
   v_[0] = Type(1);
   v_[1] = v_[2] = v_[3] = Type(0);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inversion of the quaternion (\f$ \hat{q} = \hat{q}^{-1} \f$).
 *
 * \return Reference to the inverted quaternion.
 */
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>& Quaternion<Type>::invert()
{
   v_[1] *= -Type(1);
   v_[2] *= -Type(1);
   v_[3] *= -Type(1);
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the inverse quaternion (\f$ \hat{q} = \hat{p}^{-1} \f$).
 *
 * \return The inverse quaternion.
 */
template< typename Type >  // Data type of the quaternion
inline const Quaternion<Type> Quaternion<Type>::getInverse() const
{
   return Quaternion( v_[0], -v_[1], -v_[2], -v_[3] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion to a rotation matrix.
 *
 * \return The resulting rotation matrix.
 */
template< typename Type >  // Data type of the quaternion
inline const Matrix3<Type> Quaternion<Type>::toRotationMatrix() const
{
   return Matrix3<Type>( Type(1) - Type(2)*v_[2]*v_[2] - Type(2)*v_[3]*v_[3],
                                Type(2)*( v_[1]*v_[2] - v_[0]*v_[3] ),
                                Type(2)*( v_[1]*v_[3] + v_[0]*v_[2] ),
                                Type(2)*( v_[1]*v_[2] + v_[0]*v_[3] ),
                                Type(1) - Type(2)*v_[1]*v_[1] - Type(2)*v_[3]*v_[3],
                                Type(2)*( v_[2]*v_[3] - v_[0]*v_[1] ),
                                Type(2)*( v_[1]*v_[3] - v_[0]*v_[2] ),
                                Type(2)*( v_[2]*v_[3] + v_[0]*v_[1] ),
                                Type(1) - Type(2)*v_[1]*v_[1] - Type(2)*v_[2]*v_[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the inherent rotation angle.
 *
 * \return The rotation angle (radian measure).
 */
template< typename Type >
inline Type Quaternion<Type>::getAngle() const
{
   //due to numerical accuracy v_[0] might be slightly larger than 1
   //which will results in nan for the acos function.
   if (v_[0]<Type(1)) return Type(2)*std::acos(v_[0]); else return Type(0);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the inherent rotation axis.
 *
 * \return The rotation axis
 *
 * If the angle is (close to) zero, the return value is set to (1,0,0)
 */
template< typename Type >
inline const Vector3<Type> Quaternion<Type>::getAxis() const
{
   // tmp might be negative due to numerical accuracy
   // check before putting it into sqrt!!!
   auto tmp = Type(1)-v_[0]*v_[0];

   if (tmp < std::numeric_limits<Type>::epsilon())
   {
      return Vector3<Type>( 1, 0, 0 );
   } else {
      Type s( std::sqrt( tmp ) );
      return Vector3<Type>( v_[1] / s, v_[2] / s, v_[3] / s );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotating the quaternion around the global x-axis by \a angle degrees (radian measure).
 *
 * \param angle The rotation angle (radian measure).
 * \return void
 */
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::rotateX( Type angle )
{
   const Type sina( std::sin( angle*Type(0.5) ) );
   const Type cosa( std::cos( angle*Type(0.5) ) );

   const Quaternion q( cosa*v_[0] - sina*v_[1],
                       cosa*v_[1] + sina*v_[0],
                       cosa*v_[2] - sina*v_[3],
                       cosa*v_[3] + sina*v_[2] );

   this->operator=( q );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotating the quaternion around the global y-axis by \a angle degrees (radian measure).
 *
 * \param angle The rotation angle (radian measure).
 * \return void
 */
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::rotateY( Type angle )
{
   const Type sina( std::sin( angle*Type(0.5) ) );
   const Type cosa( std::cos( angle*Type(0.5) ) );

   const Quaternion q( cosa*v_[0] - sina*v_[2],
                       cosa*v_[1] + sina*v_[3],
                       cosa*v_[2] + sina*v_[0],
                       cosa*v_[3] - sina*v_[1] );

   this->operator=( q );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotating the quaternion around the global z-axis by \a angle degrees (radian measure).
 *
 * \param angle The rotation angle (radian measure).
 * \return void
 */
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::rotateZ( Type angle )
{
   const Type sina( std::sin( angle*Type(0.5) ) );
   const Type cosa( std::cos( angle*Type(0.5) ) );

   const Quaternion q( cosa*v_[0] - sina*v_[3],
                       cosa*v_[1] - sina*v_[2],
                       cosa*v_[2] + sina*v_[1],
                       cosa*v_[3] + sina*v_[0] );

   this->operator=( q );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two quaternions.
 *
 * \param q The quaternion to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::swap( Quaternion& q ) /* throw() */
{
   std::swap( v_[0], q.v_[0] );
   std::swap( v_[1], q.v_[1] );
   std::swap( v_[2], q.v_[2] );
   std::swap( v_[3], q.v_[3] );
}
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Returns the euler angles in xyz order
*/
template< typename Type >
inline const Vector3<Type> Quaternion<Type>::getEulerAnglesXYZ() const
{
   WALBERLA_STATIC_ASSERT(!std::numeric_limits<Type>::is_integer);

   Vector3<Type> eulerAngles;
   
   eulerAngles[0] = std::atan2( Type(2) * ( v_[0] * v_[1] + v_[2] * v_[3] ), Type(1) - Type(2) * ( v_[1] * v_[1] + v_[2] * v_[2] ) );
   eulerAngles[1] = std::asin ( Type(2) * ( v_[0] * v_[2] - v_[3] * v_[1] ) );
   eulerAngles[2] = std::atan2( Type(2) * ( v_[0] * v_[3] + v_[1] * v_[2] ), Type(1) - Type(2) * ( v_[2] * v_[2] + v_[3] * v_[3] ) );

   return eulerAngles;
}
//*************************************************************************************************


//=================================================================================================
//
//  MATH FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Rotation of a vector v (\f$ \vec{rot} = \hat{q} \cdot \vec{v} \cdot \hat{q}^{-1} \f$).
 *
 * \param v The vector to be rotated.
 * \return The rotated vector.
 *
 * The function is selected for vectors of different data type (in case \a Type and \a Other
 * are supported by the MathTrait class). The function returns a vector of the higher-order
 * data type of the two involved data types.
 */
template< typename Type >  // Data type of the quaternion
template< typename Other>   // Data type of the vector
inline const Vector3< typename MathTrait<Type,Other>::MultType >
   Quaternion<Type>::rotate( const Vector3<Other>& v ) const
{
   using MT = typename MathTrait<Type, Other>::MultType;

   // Multiplication in two steps
   const MT w( v_[1]*v[0] + v_[2]*v[1] + v_[3]*v[2] );
   const MT x( v_[0]*v[0] - v_[3]*v[1] + v_[2]*v[2] );
   const MT y( v_[0]*v[1] - v_[1]*v[2] + v_[3]*v[0] );
   const MT z( v_[0]*v[2] - v_[2]*v[0] + v_[1]*v[1] );

   return Vector3<MT>( v_[0]*x + v_[1]*w + v_[2]*z - v_[3]*y,
                          v_[0]*y + v_[2]*w + v_[3]*x - v_[1]*z,
                          v_[0]*z + v_[3]*w + v_[1]*y - v_[2]*x );

   // Multiplication in one step
   /*
   const MT v0( v_[0] * v_[0] );
   const MT v1( v_[1] * v_[1] );
   const MT v2( v_[2] * v_[2] );
   const MT v3( v_[3] * v_[3] );

   return Vector3<MT>( ( v[0]*( v0 + v1 - v2 - v3 ) + 2.0*v_[0]*( v_[2]*v[2] - v_[3]*v[1] ) + 2.0*v_[1]*( v_[2]*v[1] + v_[3]*v[2] ) ),
                          ( v[1]*( v0 - v1 + v2 - v3 ) + 2.0*v_[0]*( v_[3]*v[0] - v_[1]*v[2] ) + 2.0*v_[2]*( v_[1]*v[0] + v_[3]*v[2] ) ),
                          ( v[2]*( v0 - v1 - v2 + v3 ) + 2.0*v_[0]*( v_[1]*v[1] - v_[2]*v[0] ) + 2.0*v_[3]*( v_[1]*v[0] + v_[2]*v[1] ) ) );
   */
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of a matrix.
 *
 * \param m The matrix to be rotated.
 * \return The rotated matrix.
 *
 * The function is selected for matrices of different data type (in case \a Type and \a Other
 * are supported by the MathTrait class). The function returns a matrix of the higher-order
 * data type of the two involved data types.
 */
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the matrix
inline const Matrix3< typename MathTrait<Type,Other>::MultType >
   Quaternion<Type>::rotate( const Matrix3<Other>& m ) const
{
   using MT = typename MathTrait<Type, Other>::MultType;
   const Matrix3<MT> R( this->toRotationMatrix() );
   return R.rotate( m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of a diagonal matrix.
 *
 * \param m The diagonal matrix to be rotated.
 * \return The rotated matrix.
 *
 * The DiagRotate function is a special case of the rotate function. The matrix is assumed to
 * be a diagonal matrix, which reduces the number of floating point operations of the rotation.
 * The function is selected for matrices of different data type (in case \a Type and \a Other
 * are supported by the MathTrait class). The function returns a matrix of the higher-order
 * data type of the two involved data types.
 */
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the diagonal matrix
inline const Matrix3< typename MathTrait<Type,Other>::MultType >
   Quaternion<Type>::diagRotate( const Matrix3<Other>& m ) const
{
   using MT = typename MathTrait<Type, Other>::MultType;
   const Matrix3<MT> R( this->toRotationMatrix() );
   return R.diagRotate( m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the angle in radian measure between the quaterion and a given axis.
 *
 * \param axis The given axis.
 * \return The angle in radian measure.
 */
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the axis
inline typename MathTrait<Type,Other>::HighType
   Quaternion<Type>::calcAngle( const Vector3<Other>& axis ) const
{
   using High = typename MathTrait<Type, Other>::HighType;

   const Vector3<High> u( v_[1], v_[2], v_[3] );
   const High y  ( u.length() );
   const High x  ( v_[0] );
   const High dot( u * axis );

   return High(2) * std::atan2( y, ( dot < real_c(0) )?( -x ):( x ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Quaternion operators */
//@{
template< typename T1, typename T2 >
inline bool operator==( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs );

template< typename T1, typename T2 >
inline bool operator!=( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs );

template< typename Type >
std::ostream& operator<<( std::ostream& os, const Quaternion<Type>& q );

template< typename Type >
std::istream& operator>>( std::istream& is, Quaternion<Type>& q );

template< typename Type >
inline bool isnan( const Quaternion<Type>& q );

template< typename Type >
inline void reset( Quaternion<Type>& q );

template< typename Type >
inline void clear( Quaternion<Type>& q );

template< typename Type >
inline bool isDefault( const Quaternion<Type>& q );

template< typename Type >
inline const Quaternion<Type> inv( const Quaternion<Type>& m );

template< typename Type >
inline const Quaternion<Type> sq( const Quaternion<Type>& m );

template< typename Type >
inline void swap( Quaternion<Type>& a, Quaternion<Type>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two quaternions.
 * \ingroup quaternion
 *
 * \param lhs The left-hand side quaternion for the comparison.
 * \param rhs The right-hand side quaternion for the comparison.
 * \return \a true if the two quaternions are equal, \a false if not.
 */
template< typename T1    // Data type of the left-hand side quaternion
        , typename T2 >  // Data type of the right-hand side quaternion
inline bool operator==( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs )
{
   // In order to compare the two quaternions, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   return equal( lhs[0], rhs[0] ) &&
          equal( lhs[1], rhs[1] ) &&
          equal( lhs[2], rhs[2] ) &&
          equal( lhs[2], rhs[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two quaternions.
 * \ingroup quaternion
 *
 * \param lhs The left-hand side quaternion for the comparison.
 * \param rhs The right-hand side quaternion for the comparison.
 * \return \a true if the two quaternions are not equal, \a false if they are equal.
 */
template< typename T1    // Data type of the left-hand side quaternion
        , typename T2 >  // Data type of the right-hand side quaternion
inline bool operator!=( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for quaternions.
 * \ingroup quaternion
 *
 * \param os Reference to the output stream.
 * \param q Reference to a constant quaternion object.
 * \return Reference to the output stream.
 */
template< typename Type >  // Data type of the quaternion
std::ostream& operator<<( std::ostream& os, const Quaternion<Type>& q )
{
   return os << "<" << q[0] << "," << q[1] << "," << q[2] << "," << q[3] << ">";
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global input operator for quaternions.
 * \ingroup quaternion
 *
 * \param is Reference to the input stream.
 * \param q Reference to a quaternion object.
 * \return The input stream.
 */
template< typename Type >  // Data type of the quaternion
std::istream& operator>>( std::istream& is, Quaternion<Type>& q )
{
   if( !is ) return is;

   char bracket1, bracket2, comma1, comma2, comma3;
   Type r(0), i(0), j(0), k(0);
   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );

   // Setting the 'skip whitespaces' flag
   is >> std::skipws;

   // Extracting the quaternion
   if( !(is >> bracket1 >> r >> comma1 >> i >> comma2 >> j >> comma3 >> k >> bracket2) ||
       bracket1 != '<' || comma1 != ',' || comma2 != ',' || comma3 != ',' || bracket2 != '>' ) {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      is.flags( oldFlags );
      return is;
   }

   // Transfering the input to the quaternion values
   q.set( r, i, j, k );

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given quaternion for not-a-number elements.
 * \ingroup quaternion
 *
 * \param q The quaternion to be checked for not-a-number elements.
 * \return \a true if at least one element of the quaternion is not-a-number, \a false otherwise.
 */
template< typename Type >  // Data type of the quaternion
inline bool isnan( const Quaternion<Type>& q )
{
   return isnan( q[0] ) || isnan( q[1] ) || isnan( q[2] ) || isnan( q[3] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given quaternion.
 * \ingroup quaternion
 *
 * \param q The quaternion to be resetted.
 * \return void
 */
template< typename Type >  // Data type of the quaternion
inline void reset( Quaternion<Type>& q )
{
   q.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given quaternion.
 * \ingroup quaternion
 *
 * \param q The quaternion to be cleared.
 * \return void
 *
 * Clearing a quaternion is equivalent to resetting it via the reset() function.
 */
template< typename Type >  // Data type of the quaternion
inline void clear( Quaternion<Type>& q )
{
   q.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given quaternion is in default state.
 * \ingroup quaternion
 *
 * \param q The quaternion to be tested for its default state.
 * \return \a true in case the given quaternion is in default state, \a false otherwise.
 *
 * The function returns \a true in case the real part of the quaternion is 1 and the imaginary
 * parts are 0, otherwise it returns \a false.

                             \f[\left(\begin{array}{*{4}{c}}
                             1 & 0 & 0 & 0 \\
                             \end{array}\right)\f]
 */
template< typename Type >  // Data type of the quaternion
inline bool isDefault( const Quaternion<Type>& q )
{
   return ( q[0] == Type(1) ) && ( q[1] == Type(0) ) && ( q[2] == Type(0) ) && ( q[3] == Type(0) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given quaternion.
 * \ingroup quaternion
 *
 * \param q The quaternion to be inverted.
 * \return The inverse quaternion.
 *
 * This function returns the inverse of the given quaternion. It has the same effect as
 * calling the getInverse() member function of the quaternion.
 */
template< typename Type >  // Data type of the quaternion
inline const Quaternion<Type> inv( const Quaternion<Type>& q )
{
   return q.getInverse();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Squaring the given quaternion.
 * \ingroup quaternion
 *
 * \param q The quaternion to be squared.
 * \return The result of the square operation.
 *
 * This function squares the given quaternion \a q. This function has the same effect as
 * multiplying the quaternion with itself (\f$ q * q \f$).
 */
template< typename Type >  // Data type of the quaternion
inline const Quaternion<Type> sq( const Quaternion<Type>& q )
{
   return q * q;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two quaternions.
 * \ingroup quaternion
 *
 * \param a The first quaternion to be swapped.
 * \param b The second quaternion to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
template< typename Type >  // Data type of the quaternions
inline void swap( Quaternion<Type>& a, Quaternion<Type>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Quaternion arithmetic operators
 *
 * These operators support operations between quaternions of different element types. They work
 * for all element types supported by the MathTrait class template.
 */
//@{
template< typename T1, typename T2 >
inline const Quaternion< typename MathTrait<T1,T2>::MultType >
   operator*( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of two quaternions
 *        (\f$ \hat{q}=\hat{p}*\hat{r} \f$).
 * \ingroup quaternion
 *
 * \param lhs The left-hand side quaternion for the multiplication.
 * \param rhs The right-hand side quaternion for the multiplication.
 * \return The resulting quaternion.
 */
template< typename T1    // Data type of the left-hand side quaternion
        , typename T2 >  // Data type of the right-hand side quaternion
inline const Quaternion< typename MathTrait<T1,T2>::MultType >
   operator*( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs )
{
   using MT = typename MathTrait<T1, T2>::MultType;

   const MT r( lhs[0]*rhs[0] - lhs[1]*rhs[1] - lhs[2]*rhs[2] - lhs[3]*rhs[3] );
   const MT i( lhs[0]*rhs[1] + lhs[1]*rhs[0] + lhs[2]*rhs[3] - lhs[3]*rhs[2] );
   const MT j( lhs[0]*rhs[2] + lhs[2]*rhs[0] + lhs[3]*rhs[1] - lhs[1]*rhs[3] );
   const MT k( lhs[0]*rhs[3] + lhs[3]*rhs[0] + lhs[1]*rhs[2] - lhs[2]*rhs[1] );

   const MT len2( r*r + i*i + j*j + k*k );
   WALBERLA_ASSERT(!std::isnan(len2), lhs << "\n" << rhs);

   if( std::fabs( len2 - MT(1) ) < MT(1E-8) ) {
      return Quaternion<MT>( r, i, j, k );
   }
   else {
      const MT ilen( MT(1) / std::sqrt( len2 ) );
      return Quaternion<MT>( r*ilen, i*ilen, j*ilen, k*ilen );
   }
}
//*************************************************************************************************

} // namespace math
}

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

   template< typename T,    // Element type of SendBuffer
             typename G,    // Growth policy of SendBuffer
             typename VT >  // Element type of Quaternion
   mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const math::Quaternion<VT> & quat )
   {
      buf.addDebugMarker( "q4" );
      static_assert ( std::is_trivially_copyable< math::Quaternion<VT> >::value,
                      "type has to be trivially copyable for the memcpy to work correctly" );
      auto pos = buf.forward(sizeof(math::Quaternion<VT>));
      std::memcpy(pos, &quat, sizeof(math::Quaternion<VT>));
      return buf;
   }

   template< typename T,    // Element type  of RecvBuffer
             typename VT >  // Element type of vector
   mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, math::Quaternion<VT> & quat )
   {
      buf.readDebugMarker( "q4" );
      static_assert ( std::is_trivially_copyable< math::Quaternion<VT> >::value,
                      "type has to be trivially copyable for the memcpy to work correctly" );
      auto pos = buf.skip(sizeof(math::Quaternion<VT>));
      //suppress https://gcc.gnu.org/onlinedocs/gcc/C_002b_002b-Dialect-Options.html#index-Wclass-memaccess
      std::memcpy(static_cast<void*>(&quat), pos, sizeof(math::Quaternion<VT>));
      return buf;
   }

   template<typename VT>
   struct BufferSizeTrait< walberla::math::Quaternion<VT> > {
      static const bool constantSize = true;
      static const uint_t size = 4 * BufferSizeTrait<VT>::size + mpi::BUFFER_DEBUG_OVERHEAD;
   };
}
}

//======================================================================================================================
//
//  Vector Trait Specialization
//
//======================================================================================================================

namespace walberla {

using math::Quaternion;

// Specialization of VectorTrait for Quaternions
template<typename T>
struct VectorTrait< Quaternion<T> >
{
   using OutputType = T;

   static const uint_t F_SIZE =  4u;
   static T    get( const Quaternion<T> & v, uint_t f )       { return v[f]; }
   static void set(       Quaternion<T> & v, uint_t f, T val) { v[f] = val;  }
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
inline bool check_float_equal( const math::Quaternion<real_t> & lhs, const math::Quaternion<real_t> & rhs )
{
   return floatIsEqual( lhs[0], rhs[0] ) && floatIsEqual( lhs[1], rhs[1] ) && floatIsEqual( lhs[2], rhs[2] ) && floatIsEqual( lhs[3], rhs[3] );
}

template< >
inline bool check_float_equal_eps( const math::Quaternion<real_t> & lhs, const math::Quaternion<real_t> & rhs, const real_t epsilon )
{
   return floatIsEqual( lhs[0], rhs[0], epsilon ) && floatIsEqual( lhs[1], rhs[1], epsilon ) && floatIsEqual( lhs[2], rhs[2], epsilon ) && floatIsEqual( lhs[3], rhs[3], epsilon );
}

}
}
}
