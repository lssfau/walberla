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
//! \file GenericAABB.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "core/math/Shims.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <random>
#include <type_traits>


namespace walberla {
namespace math {

/**
 * \brief GenericAABB represents an axis-aligned bounding box
 *
 * \tparam T The data type used to store the minimal and maximal corner points. May be float, double or long double.
 *
 * \invariant minCorner.x <= maxCorner.x && minCorner.y <= maxCorner.y && minCorner.z <= maxCorner.z
 */
template< typename T >
class GenericAABB
{
   static_assert( std::is_floating_point< T >::value, "GenericAABB only works with floating point types for T!" );

public:
   // Typedefs
   using value_type = T;  /// scalar data type
   using vector_type = Vector3<T>; /// data type for three dimensional vectors

   // Constructors
   inline GenericAABB();
   template< typename U >
   inline GenericAABB( const GenericAABB< U > & other );
   inline GenericAABB( const vector_type & corner0, const vector_type & corner1 );
   inline GenericAABB( const value_type x0, const value_type y0, const value_type z0,
                       const value_type x1, const value_type y1, const value_type z1 );
   template< typename InputIterator >
   inline GenericAABB( InputIterator first, InputIterator last );

   // Static Factory Methods
   static GenericAABB createFromMinMaxCorner( const vector_type & theMinCorner, const vector_type & theMaxCorner );
   static GenericAABB createFromMinMaxCorner( const value_type minX, const value_type minY, const value_type minZ,
                                              const value_type maxX, const value_type maxY, const value_type maxZ );

   // Operators ( also see free functions below )
   template< typename U >
   inline GenericAABB & operator=( const GenericAABB< U > & other );

   // Observers
   inline const vector_type & minCorner() const;
   inline const vector_type & maxCorner() const;

   inline const vector_type & min() const { return minCorner(); }
   inline const vector_type & max() const { return maxCorner(); }

   inline value_type xMin() const;
   inline value_type yMin() const;
   inline value_type zMin() const;

   inline value_type min( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, uint_t(3) ); return minCorner_[index]; }

   inline value_type xMax() const;
   inline value_type yMax() const;
   inline value_type zMax() const;

   inline value_type max( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, uint_t(3) ); return maxCorner_[index]; }

   inline bool empty() const;

   inline vector_type sizes() const;
   inline value_type size( const uint_t index ) const;
   inline value_type xSize() const;
   inline value_type ySize() const;
   inline value_type zSize() const;

   inline value_type volume() const;

   inline vector_type center() const;

   inline bool contains( const value_type x, const value_type y, const value_type z ) const;
   inline bool contains( const value_type x, const value_type y, const value_type z, const value_type dx ) const;
   inline bool contains( const vector_type & point ) const;
   inline bool contains( const vector_type & point, const value_type dx ) const;
   inline bool contains( const GenericAABB & other ) const;
   inline bool contains( const GenericAABB & other, const value_type dx ) const;

   inline bool containsClosedInterval( const vector_type & point ) const;
   inline bool containsClosedInterval( const vector_type & point, const value_type dx ) const;

   inline GenericAABB getExtended( const value_type eps ) const;
   inline GenericAABB getExtended( const vector_type & eps ) const;

   inline GenericAABB getTranslated( const vector_type & translation ) const;

   inline GenericAABB getScaled( const value_type factor ) const;
   inline GenericAABB getScaled( const vector_type & factors ) const;

   inline GenericAABB getMerged( const vector_type & point ) const;
   inline GenericAABB getMerged( const GenericAABB & other ) const;

   template< typename InputIterator >
   inline GenericAABB getMerged( InputIterator first, InputIterator last ) const;

   inline bool intersects( const GenericAABB & other ) const;
   inline bool intersects( const GenericAABB & other, const value_type dx ) const;
   inline bool intersectsClosedInterval( const GenericAABB & other ) const;
   inline bool intersectsClosedInterval( const GenericAABB & other, const value_type dx ) const;
   inline value_type intersectionVolume( const GenericAABB & other ) const;
   inline GenericAABB getIntersection( const GenericAABB & other ) const;

   inline bool isIdentical( const GenericAABB & rhs ) const;

   inline bool isEqual( const GenericAABB & rhs ) const;

   inline value_type sqDistance( const vector_type & point ) const;
   inline value_type sqSignedDistance( const vector_type & point ) const;
   inline value_type sqMaxDistance( const vector_type & point ) const;
   inline value_type distance( const vector_type & point ) const;
   inline value_type signedDistance( const vector_type & point ) const;
   inline value_type maxDistance( const vector_type & point ) const;

   inline value_type sqDistance( const GenericAABB & other ) const;
   inline value_type sqMaxDistance( const GenericAABB & other ) const;

   inline std::array< vector_type, 8 > corners() const;

   // Modifiers
   inline void init();
   inline void init( const vector_type & corner0, const vector_type & corner1 );
   inline void init( const value_type x0, const value_type y0, const value_type z0,
                     const value_type x1, const value_type y1, const value_type z1 );
   template< typename InputIterator >
   inline void init( InputIterator first, InputIterator last );

   inline void initMinMaxCorner( const vector_type & theMinCorner, const vector_type & theMaxCorner );
   inline void initMinMaxCorner( const value_type minX, const value_type minY, const value_type minZ,
                                 const value_type maxX, const value_type maxY, const value_type maxZ );

   inline void setAxisBounds( const uint_t index, const value_type value1, const value_type value2 );

   inline void extend( const value_type eps );
   inline void extend( const vector_type & eps );

   inline void setCenter( const vector_type & center );
   inline void translate( const vector_type & translation );

   inline void scale( const value_type factor );
   inline void scale( const vector_type & factors );

   inline void merge( const vector_type & point );
   inline void merge( const GenericAABB & other );

   template< typename InputIterator >
   inline void merge( InputIterator first, InputIterator last );

   inline void intersect( const GenericAABB & other );

   template< typename Engine >
   inline vector_type randomPoint( Engine & engine ) const;

   bool checkInvariant() const;

   /**
    * \brief Deserializes an GenericAABB from a mpi::GenericRecvBuffer
    *
    * \param buf  The mpi::GenericRecvBuffer written to
    * \param aabb The GenericAABB the be deserialized
    *
    * \returns A reference to buf
    */
   template< typename ET >    // Element type  of RecvBuffer
   inline friend mpi::GenericRecvBuffer< ET > & operator>>( mpi::GenericRecvBuffer< ET > & buf, GenericAABB< T > & aabb )
   {
      buf.readDebugMarker( "bb" );
      static_assert ( std::is_trivially_copyable< GenericAABB< T > >::value,
                      "type has to be trivially copyable for the memcpy to work correctly" );
      auto pos = buf.skip(sizeof(GenericAABB< T >));
      std::memcpy(&aabb, pos, sizeof(GenericAABB< T >));
      WALBERLA_ASSERT( aabb.checkInvariant() );
      return buf;
   }

protected:

   struct MinMaxCornerGivenT { }; /// Helper type to distinguish ctors with equal signature.

   inline GenericAABB( const vector_type & theMinCorner, const vector_type & theMaxCorner, MinMaxCornerGivenT );
   inline GenericAABB( const value_type minX, const value_type minY, const value_type minZ,
                       const value_type maxX, const value_type maxY, const value_type maxZ, MinMaxCornerGivenT );
   inline GenericAABB( const GenericAABB & lhs, const GenericAABB & rhs );

private:
   vector_type minCorner_; /// minimal values
   vector_type maxCorner_; /// maximal values
};
static_assert( std::is_trivially_copyable<GenericAABB<real_t>>::value, "GenericAABB<real_t> has to be trivially copyable!");



template< typename T, typename U >
inline bool operator==( const GenericAABB< T > & lhs, const GenericAABB< U > & rhs );

template< typename T, typename U >
inline bool operator!=( const GenericAABB< T > & lhs, const GenericAABB< U > & rhs );

template< typename T >
inline std::ostream& operator<<( std::ostream& os, const GenericAABB< T > & aabb );

template< typename T >
inline std::istream& operator>>( std::istream& is, GenericAABB< T > & aabb );


template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename VT >  // Element type of GenericAABB
inline mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const GenericAABB< VT > & aabb );


} // namespace math


namespace mpi {

template< typename T >
struct BufferSizeTrait< walberla::math::GenericAABB< T > > {
   static const bool constantSize = true;
   static const uint_t size = 2u * BufferSizeTrait< typename walberla::math::GenericAABB< T >::vector_type >::size + mpi::BUFFER_DEBUG_OVERHEAD;
};

} // namespace mpi

} // namespace walberla

#include "GenericAABB.impl.h"
