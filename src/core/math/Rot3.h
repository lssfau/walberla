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
//! \file Rot3.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/math/Matrix3.h>
#include <core/math/Quaternion.h>

#include <type_traits>

namespace walberla {
namespace math {

/**
 * Rotation class which merges quaternion and matrix representation.
 *
 * For numerical reasons the representation of a rotation via a quaternion
 * is favourable. However application of the rotation to vectors and matrices
 * is numerical more efficient with matrices. Therefore this class combines both
 * representations and takes care that both are in sync.
 */
template <typename Type> // floating point type
class Rot3
{
   //**Compile time checks*************************************************************************
   /*! \cond internal */
   static_assert(std::is_floating_point<Type>::value, "T has to be floating point!");
   static_assert(!std::is_const<Type>::value, "T has to be non const!");
   static_assert(!std::is_volatile<Type>::value, "T has to be non volatile!");
   /*! \endcond */
   //**********************************************************************************************
public:
   Rot3();
   Rot3(const Vector3<Type>& rot);
   Rot3(const Vector3<Type>& axis, const real_t& angle);
   Rot3(const Quaternion<Type> q);

   const Quaternion<Type>& getQuaternion() const { return quat_; }
   const Matrix3<Type>&    getMatrix() const { return mat_; }

   void rotate(const Vector3<Type>& rot);
   void rotate(const Vector3<Type>& axis, const real_t& angle);
private:
   Quaternion<Type> quat_;
   Matrix3<Type>    mat_;
};

template< typename Type >  // floating point type
inline Rot3<Type>::Rot3()
   : quat_(Vector3<Type>(Type(1),Type(0),Type(0)), Type(0))
   , mat_(quat_.toRotationMatrix())
{
   WALBERLA_ASSERT_FLOAT_EQUAL( mat_.getDeterminant(), real_t(1), "Corrupted rotation matrix determinant" );
}

template< typename Type >  // floating point type
inline Rot3<Type>::Rot3(const Vector3<Type>& phi)
   : Rot3<Type>()
{
   rotate(phi);
}

template< typename Type >  // floating point type
inline Rot3<Type>::Rot3(const Vector3<Type>& axis, const real_t& angle)
   : Rot3<Type>()
{
   rotate(axis, angle);
}

template< typename Type >  // floating point type
inline Rot3<Type>::Rot3(const Quaternion<Type> q)
   : quat_(q)
   , mat_(quat_.toRotationMatrix())
{}

template< typename Type >  // floating point type
inline void Rot3<Type>::rotate(const Vector3<Type>& phi)
{
   auto len = phi.length();
   if (!floatIsEqual(len, 0))
   {
      auto q = Quaternion<Type>( phi, len );
      quat_ = q * quat_;
      mat_  = quat_.toRotationMatrix();
      WALBERLA_ASSERT_FLOAT_EQUAL( mat_.getDeterminant(), real_t(1), "Corrupted rotation matrix determinant" );
   }
}

template< typename Type >  // floating point type
inline void Rot3<Type>::rotate(const Vector3<Type>& axis, const real_t& angle)
{
   if (!floatIsEqual(angle, 0))
   {
      auto q = Quaternion<Type>( axis, angle );
      quat_ = q * quat_;
      mat_  = quat_.toRotationMatrix();
      WALBERLA_ASSERT_FLOAT_EQUAL( mat_.getDeterminant(), real_t(1), "Corrupted rotation matrix determinant" );
   }
}

template< typename Type >  // floating point type
inline std::ostream& operator<<( std::ostream& os, const Rot3<Type>& r )
{
   os << r.getMatrix();
   return os;
}

} // math
} // walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename V >   // value type
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const math::Rot3<V>& obj )
{
   buf.addDebugMarker( "ro" );
   buf << obj.getQuaternion();
   return buf;
}

template< typename T,    // Element type  of RecvBuffer
          typename V >   // value type
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, math::Rot3<V>& objparam )
{
   math::Quaternion<V> q;
   buf.readDebugMarker( "ro" );
   buf >> q;
   objparam = math::Rot3<V>(q);
   return buf;
}

template< typename V > // value type
struct BufferSizeTrait< math::Rot3<V> > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait< math::Quaternion<V> >::size + mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla
