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
//! \file GenericAABB.impl.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================



namespace walberla {
namespace math {

/**
 * \brief Builds a default GenericAABB
 *
 * \post minCorner = maxCorner = (0,0,0)
 */
template< typename T >
GenericAABB< T >::GenericAABB() : minCorner_( value_type(0), value_type(0), value_type(0) ),
                                  maxCorner_( value_type(0), value_type(0), value_type(0) )
{
   WALBERLA_ASSERT( checkInvariant(), *this );
   WALBERLA_ASSERT( empty() );
}


/**
 * \brief Copy constructor for GenericAABBs of different value_type
 *
 * \tparam U value_type of other GenericAABB
 *
 * \param other other GenericAABB
 */
template< typename T >
template< typename U >
inline GenericAABB< T >::GenericAABB( const GenericAABB< U > & other )
   : minCorner_( other.minCorner() ),
     maxCorner_( other.maxCorner() )
{
   WALBERLA_ASSERT( checkInvariant(), *this );
}


/**
 * \brief Builds a GenericAABB from two arbitrary points
 *
 * \param corner0 first corner spanning the AABB
 * \param corner1 second corner spanning the AABB
 */
template< typename T >
GenericAABB< T >::GenericAABB( const vector_type & corner0, const vector_type & corner1 )
{
   for( uint_t i = 0; i < 3; ++i )
   {
      if( corner0[i] < corner1[i] )
      {
         minCorner_[i] = corner0[i];
         maxCorner_[i] = corner1[i];
      }
      else
      {
         minCorner_[i] = corner1[i];
         maxCorner_[i] = corner0[i];
      }
   }

   WALBERLA_ASSERT( checkInvariant(), *this );
}



/**
 * \brief Builds a GenericAABB from two arbitrary points given with their coordinates
 *
 * \param x0  x-coordinate of the first corner spanning the AABB
 * \param y0  y-coordinate of the first corner spanning the AABB
 * \param z0  z-coordinate of the first corner spanning the AABB
 * \param x1  x-coordinate of the second corner spanning the AABB
 * \param y1  y-coordinate of the second corner spanning the AABB
 * \param z1  z-coordinate of the second corner spanning the AABB
 */
template< typename T >
GenericAABB< T >::GenericAABB( const value_type x0, const value_type y0, const value_type z0,
                               const value_type x1, const value_type y1, const value_type z1 )
{
   if( x0 < x1 )
   {
      minCorner_[0] = x0;
      maxCorner_[0] = x1;
   }
   else
   {
      minCorner_[0] = x1;
      maxCorner_[0] = x0;
   }

   if( y0 < y1 )
   {
      minCorner_[1] = y0;
      maxCorner_[1] = y1;
   }
   else
   {
      minCorner_[1] = y1;
      maxCorner_[1] = y0;
   }

   if( z0 < z1 )
   {
      minCorner_[2] = z0;
      maxCorner_[2] = z1;
   }
   else
   {
      minCorner_[2] = z1;
      maxCorner_[2] = z0;
   }

   WALBERLA_ASSERT( checkInvariant(), *this );
}



/**
 * \brief Builds a default GenericAABB by merging a Sequence of points or GenericAABBs
 *
 * \tparam InputIterator ForwardIterator. Must dereference to Vector3< T > or GenericAABB< T >.
 *
 * \param first first element of sequence [first, last)
 * \param last  final element of sequence [first, last)
 *
 * \post first == last => minCorner = maxCorner = (0,0,0)
 */
template< typename T >
template< typename InputIterator >
GenericAABB< T >::GenericAABB( InputIterator first, InputIterator last )
   : minCorner_(  std::numeric_limits<value_type>::max(),  std::numeric_limits<value_type>::max(),  std::numeric_limits<value_type>::max() ),
     maxCorner_( -std::numeric_limits<value_type>::max(), -std::numeric_limits<value_type>::max(), -std::numeric_limits<value_type>::max() )
{
   if( first != last )
   {
      merge( first, last );
   }
   else
   {
      init();
   }

   WALBERLA_ASSERT( checkInvariant(), *this );
}



/**
 * \brief Builds a GenericAABB from a min- and a maxCorner
 *
 * \param theMinCorner corner with minimal values
 * \param theMaxCorner corner with maximal values
 *
 * \pre theMinCorner.x <= theMaxCorner.x && theMinCorner.y <= theMaxCorner.y && theMinCorner.z <= theMaxCorner.z
 */
template< typename T >
GenericAABB< T >::GenericAABB( const vector_type & theMinCorner, const vector_type & theMaxCorner, MinMaxCornerGivenT )
    : minCorner_( theMinCorner ), maxCorner_( theMaxCorner )
{
   WALBERLA_ASSERT( checkInvariant(), *this );
}



/**
 * \brief Builds a GenericAABB as an intersection of two other GenericAABBs
 *
 * \param lhs an arbitrary GenericAABB
 * \param rhs an arbitrary GenericAABB
 *
 * \post lhs does not intersect rhs => minCorner = maxCorner )
 */
template< typename T >
GenericAABB< T >::GenericAABB( const GenericAABB & lhs, const GenericAABB & rhs )
   : minCorner_( std::max( lhs.minCorner_[0], rhs.minCorner_[0] ),
                 std::max( lhs.minCorner_[1], rhs.minCorner_[1] ),
                 std::max( lhs.minCorner_[2], rhs.minCorner_[2] ) ),

     maxCorner_( std::max( minCorner_[0], std::min( lhs.maxCorner_[0], rhs.maxCorner_[0] ) ),
                 std::max( minCorner_[1], std::min( lhs.maxCorner_[1], rhs.maxCorner_[1] ) ),
                 std::max( minCorner_[2], std::min( lhs.maxCorner_[2], rhs.maxCorner_[2] ) ) )
{
   WALBERLA_ASSERT( checkInvariant(), *this );
}



/**
 * \brief Builds a GenericAABB from a min- and a maxCorner
 *
 * \param minX  x-coordinate of the minCorner
 * \param minY  y-coordinate of the minCorner
 * \param minZ  z-coordinate of the minCorner
 * \param maxX  x-coordinate of the maxCorner
 * \param maxY  y-coordinate of the maxCorner
 * \param maxZ  z-coordinate of the maxCorner
 *
 * \pre x0 <= x1 && y0 <= y1 && z0 <= z1
 */
template< typename T >
GenericAABB< T >::GenericAABB( const value_type minX, const value_type minY, const value_type minZ,
                               const value_type maxX, const value_type maxY, const value_type maxZ, MinMaxCornerGivenT )
   : minCorner_( minX, minY, minZ ), maxCorner_( maxX, maxY, maxZ )
{
   WALBERLA_ASSERT( checkInvariant(), *this );
}



/**
 * \brief Builds a GenericAABB from a min- and a maxCorner
 *
 * \param minCorner corner with minimal values
 * \param maxCorner corner with maximal values
 *
 * \pre minCorner.x <= maxCorner.x && minCorner.y <= maxCorner.y && minCorner.z <= maxCorner.z
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::createFromMinMaxCorner( const vector_type & theMinCorner, const vector_type & theMaxCorner )
{
   return GenericAABB< T >( theMinCorner, theMaxCorner, MinMaxCornerGivenT() );
}



/**
 * \brief Builds a GenericAABB from a min- and a maxCorner
 *
 * \param minX  x-coordinate of the minCorner
 * \param minY  y-coordinate of the minCorner
 * \param minZ  z-coordinate of the minCorner
 * \param maxX  x-coordinate of the maxCorner
 * \param maxY  y-coordinate of the maxCorner
 * \param maxZ  z-coordinate of the maxCorner
 *
 * \pre x0 <= x1 && y0 <= y1 && z0 <= z1
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::createFromMinMaxCorner(
   const value_type minX, const value_type minY, const value_type minZ,
   const value_type maxX, const value_type maxY, const value_type maxZ )
{
   return GenericAABB< T >( minX, minY, minZ, maxX, maxY, maxZ, MinMaxCornerGivenT() );
}



/**
 * \brief Assignment operator for GenericAABBs of different value_type
 *
 * \tparam U value_type of other GenericAABB
 *
 * \param other other GenericAABB
 */
template< typename T >
template< typename U >
inline GenericAABB< T > & GenericAABB< T >::operator=( const GenericAABB< U > & other )
{
   minCorner_ = other.minCorner_();
   maxCorner_ = other.maxCorner_();

   return *this;
}



/**
 * \brief Provides access to the corner with minimal values
 *
 * \returns the corner with minimal values
 */
template< typename T >
const typename GenericAABB< T >::vector_type & GenericAABB< T >::minCorner() const
{
   return minCorner_;
}



/**
 * \brief Provides access to the corner with maximal values
 *
 * \returns the corner with maximal values
 */
template< typename T >
const typename GenericAABB< T >::vector_type & GenericAABB< T >::maxCorner() const
{
   return maxCorner_;
}



/**
 * \brief Provides access to the x-component of the corner with minimal values
 *
 * \returns the x-component of the corner with minimal values
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::xMin() const
{
   return minCorner_[0];
}



/**
 * \brief Provides access to the y-component of the corner with minimal values
 *
 * \returns the y-component of the corner with minimal values
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::yMin() const
{
   return minCorner_[1];
}



/**
 * \brief Provides access to the z-component of the corner with minimal values
 *
 * \returns the z-component of the corner with minimal values
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::zMin() const
{
   return minCorner_[2];
}



/**
 * \brief Provides access to the x-component of the corner with maximal values
 *
 * \returns the x-component of the corner with maximal values
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::xMax() const
{
   return maxCorner_[0];
}



/**
 * \brief Provides access to the y-component of the corner with maximal values
 *
 * \returns the y-component of the corner with maximal values
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::yMax() const
{
   return maxCorner_[1];
}



/**
 * \brief Provides access to the z-component of the corner with maximal values
 *
 * \returns the z-component of the corner with maximal values
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::zMax() const
{
   return maxCorner_[2];
}



/**
 * \brief Checks whether the GenericAABB is empty
 *
 * \returns true, if the GenericAABB is empty(), false else
 */
template< typename T >
bool GenericAABB< T >::empty() const
{
   return walberla::isIdentical( minCorner_[0], maxCorner_[0] ) ||
          walberla::isIdentical( minCorner_[1], maxCorner_[1] ) ||
          walberla::isIdentical( minCorner_[2], maxCorner_[2] );
}



/**
 * \brief Provides access to a vector of the genericAABBs extends
 *
 * \returns A vector with the GenericAABBs extends in x-, y-, and z-Direction
 */
template< typename T >
typename GenericAABB< T >::vector_type GenericAABB< T >::sizes() const
{
   WALBERLA_ASSERT( checkInvariant() );

   return vector_type( maxCorner_[0] - minCorner_[0], maxCorner_[1] - minCorner_[1], maxCorner_[2] - minCorner_[2] );
}



/**
 * \brief Provides access to one of the extends of the genericAABB
 *
 * \param index Determines the axis of which the extend is given: 0 -> x, 1 -> y, 2 -> z
 *
 * \returns The extend of the GenericAABB on the requested axis
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::size( const uint_t index ) const
{
   WALBERLA_ASSERT( checkInvariant() );

   return maxCorner_[index] - minCorner_[index];
}



/**
 * \brief Provides access to the extend of the genericAABB in x direction
 *
 * \returns The extend of the GenericAABB on the x axis
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::xSize() const
{
   WALBERLA_ASSERT( checkInvariant() );

   return maxCorner_[0] - minCorner_[0];
}



/**
 * \brief Provides access to the extend of the genericAABB in y direction
 *
 * \returns The extend of the GenericAABB on the y axis
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::ySize() const
{
   WALBERLA_ASSERT( checkInvariant() );

   return maxCorner_[1] - minCorner_[1];
}



/**
 * \brief Provides access to the extend of the genericAABB in z direction
 *
 * \returns The extend of the GenericAABB on the z axis
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::zSize() const
{
   WALBERLA_ASSERT( checkInvariant() );

   return maxCorner_[2] - minCorner_[2];
}



/**
 * \brief Calculates the volume of the genericAABB
 *
 * \returns The calculated volume
 *
 * \invariant empty() => 0
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::volume() const
{
   WALBERLA_ASSERT( checkInvariant() );
   return ( maxCorner_[0] - minCorner_[0] ) * ( maxCorner_[1] - minCorner_[1] ) * ( maxCorner_[2] - minCorner_[2] );
}



/**
 * \brief Calculates the center of the genericAABB
 *
 * \returns The calculated center
 */
template< typename T >
typename GenericAABB< T >::vector_type GenericAABB< T >::center() const
{
   WALBERLA_ASSERT( checkInvariant() );

   return vector_type( ( maxCorner_[0] + minCorner_[0] ) * value_type(0.5),
                       ( maxCorner_[1] + minCorner_[1] ) * value_type(0.5),
                       ( maxCorner_[2] + minCorner_[2] ) * value_type(0.5) );
}



/**
 * \brief Tests whether a point lies within the GenericAABB
 *
 * This test interprets the boundaries of the GenericAABB as an half-open interval. The minimal values stored in minCorner
 * are considered to belong to the bounding box while the maximal values stored in maxCorner are not.
 *
 * \param x x coordinate of the point to be tested for containment
 * \param y y coordinate of the point to be tested for containment
 * \param z z coordinate of the point to be tested for containment
 *
 * \returns if point lies in the GenericAABB true, else false
 */
template< typename T >
bool GenericAABB< T >::contains( const value_type x, const value_type y, const value_type z ) const
{
   return x >= minCorner_[0] && x < maxCorner_[0] &&
          y >= minCorner_[1] && y < maxCorner_[1] &&
          z >= minCorner_[2] && z < maxCorner_[2];
}



/**
 * \brief Tests whether a point lies within the extended GenericAABB.
 *
 * This test interprets the boundaries of the GenericAABB as an half-open interval. The minimal values stored in minCorner
 * are considered to belong to the bounding box while the maximal values stored in maxCorner are not.
 *
 * \param x  x coordinate of the point to be tested for containment
 * \param y  y coordinate of the point to be tested for containment
 * \param z  z coordinate of the point to be tested for containment
 * \param dx An epsilon the box is extended by in each direction before the test
 *
 * \returns if point lies in the extended GenericAABB true, else false
 */
template< typename T >
bool GenericAABB< T >::contains( const value_type x, const value_type y, const value_type z, const value_type dx ) const
{
   return x >= ( minCorner_[0] - dx ) && x < ( maxCorner_[0] + dx ) &&
          y >= ( minCorner_[1] - dx ) && y < ( maxCorner_[1] + dx ) &&
          z >= ( minCorner_[2] - dx ) && z < ( maxCorner_[2] + dx );
}



/**
 * \brief Tests whether a point lies within the GenericAABB
 *
 * This test interprets the boundaries of the GenericAABB as an half-open interval. The minimal values stored in minCorner
 * are considered to belong to the bounding box while the maximal values stored in maxCorner are not.
 *
 * \param point The point to be tested for containment
 *
 * \returns if point lies in the GenericAABB true, else false
 */
template< typename T >
bool GenericAABB< T >::contains( const vector_type & point ) const
{
   return contains( point[0], point[1], point[2] );
}



/**
 * \brief Tests whether a point lies within the extended GenericAABB.
 *
 * This test interprets the boundaries of the GenericAABB as an half-open interval. The minimal values stored in minCorner
 * are considered to belong to the bounding box while the maximal values stored in maxCorner are not.
 *
 * \param point The point to be tested for containment
 * \param dx    An epsilon the box is extended by in each direction before the test
 *
 * \returns if point lies in the extended GenericAABB true, else false
 */
template< typename T >
bool GenericAABB< T >::contains( const vector_type & point, const value_type dx ) const
{
   return contains( point[0], point[1], point[2], dx );
}



/**
 * \brief Tests whether an other GenericAABB lies within the GenericAABB
 *
 * \param other An other GenericAABB
 *
 * \returns if other lies in the GenericAABB true, else false
 */
template< typename T >
bool GenericAABB< T >::contains( const GenericAABB & other ) const
{
   return other.minCorner_[0] >= minCorner_[0] && other.maxCorner_[0] <= maxCorner_[0] &&
          other.minCorner_[1] >= minCorner_[1] && other.maxCorner_[1] <= maxCorner_[1] &&
          other.minCorner_[2] >= minCorner_[2] && other.maxCorner_[2] <= maxCorner_[2];
}



/**
 * \brief Tests whether an other GenericAABB lies within the extended GenericAABB
 *
 * \param other An other GenericAABB
 * \param dx    An epsilon the box is extended by in each direction before the test
 *
 * \returns if other lies in the extended GenericAABB true, else false
 */
template< typename T >
bool GenericAABB< T >::contains( const GenericAABB & other, const value_type dx ) const
{
   return other.minCorner_[0] >= ( minCorner_[0] - dx ) && other.maxCorner_[0] <= ( maxCorner_[0] + dx ) &&
          other.minCorner_[1] >= ( minCorner_[1] - dx ) && other.maxCorner_[1] <= ( maxCorner_[1] + dx ) &&
          other.minCorner_[2] >= ( minCorner_[2] - dx ) && other.maxCorner_[2] <= ( maxCorner_[2] + dx );
}



/**
 * \brief Tests whether a point lies within the GenericAABB
 *
 * This test interprets the boundaries of the GenericAABB as a closed interval. The minimal values stored in minCorner
 * are considered to belong to the bounding box as well as the maximal values stored in maxCorner.
 *
 * \param point The point to be tested for containment
 *
 * \returns if point lies in the GenericAABB true, else false
 */
template< typename T >
bool GenericAABB< T >::containsClosedInterval( const vector_type & point ) const
{
   return point[0] >= minCorner_[0] && point[0] <= maxCorner_[0] &&
          point[1] >= minCorner_[1] && point[1] <= maxCorner_[1] &&
          point[2] >= minCorner_[2] && point[2] <= maxCorner_[2];
}



/**
 * \brief Tests whether a point lies within the extended GenericAABB.
 *
 * This test interprets the boundaries of the GenericAABB as a closed interval. The minimal values stored in minCorner
 * are considered to belong to the bounding box as well as the maximal values stored in maxCorner.
 *
 * \param point The point to be tested for containment
 * \param dx    An epsilon the box is extended by in each direction before the test
 *
 * \returns if point lies in the extended GenericAABB true, else false
 */
template< typename T >
bool GenericAABB< T >::containsClosedInterval( const vector_type & point, const value_type dx ) const
{
   return point[0] >= ( minCorner_[0] - dx ) && point[0] <= ( maxCorner_[0] + dx ) &&
          point[1] >= ( minCorner_[1] - dx ) && point[1] <= ( maxCorner_[1] + dx ) &&
          point[2] >= ( minCorner_[2] - dx ) && point[2] <= ( maxCorner_[2] + dx );
}



/**
 * \brief Creates a new GenericAABB by extending this one
 *
 * \param e epsilon by which the bounding box is extended in each direction
 *
 * \returns The extended GenericAABB
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::getExtended( const value_type e ) const
{
   vector_type newMinCorner( minCorner_[0] - e, minCorner_[1] - e, minCorner_[2] - e );

   return createFromMinMaxCorner( newMinCorner[0], newMinCorner[1], newMinCorner[2],
                                  std::max( newMinCorner[0], maxCorner_[0] + e ),
                                  std::max( newMinCorner[1], maxCorner_[1] + e ),
                                  std::max( newMinCorner[2], maxCorner_[2] + e ) );
}



/**
 * \brief Creates a new GenericAABB by extending this one
 *
 * \param e epsilon vector by which the bounding box is extended. The box is extended in each direction by
 *          the corresponding vector component.
 *
 * \returns The extended GenericAABB
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::getExtended( const vector_type & e ) const
{
   vector_type newMinCorner( minCorner_ - e );
   return createFromMinMaxCorner( newMinCorner[0], newMinCorner[1], newMinCorner[2],
                                  std::max( newMinCorner[0], maxCorner_[0] + e[0] ),
                                  std::max( newMinCorner[1], maxCorner_[1] + e[0] ),
                                  std::max( newMinCorner[2], maxCorner_[2] + e[0] ) );
}



/**
 * \brief Creates a new GenericAABB by translating this one
 *
 * \param translation Vector by which the GenericAABB is translated.
 *
 * \returns The translated GenericAABB
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::getTranslated( const vector_type & translation ) const
{
   return createFromMinMaxCorner( minCorner_ + translation, maxCorner_ + translation );
}



/**
 * \brief Creates a new GenericAABB by scaling this one
 *
 * \param factor Factor by which the bounding box gets scaled.
 *
 * \returns The scaled GenericAABB
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::getScaled( const value_type factor ) const
{
   vector_type theCenter = center();

   GenericAABB result = getTranslated( -theCenter );
   result.minCorner_ *= factor;
   result.maxCorner_ *= factor;
   result.translate( theCenter );

   return result;
}



/**
 * \brief Creates a new GenericAABB by scaling this one
 *
 * \param factor Vector of scaling factors by which the bounding box gets scaled on the
 *               respective axises.
 *
 * \returns The scaled GenericAABB
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::getScaled( const vector_type & factors ) const
{
   vector_type theCenter = center();

   GenericAABB result = getTranslated( -theCenter );

   for( uint_t i = 0; i < 3; ++i )
   {
      result.minCorner_[i] *= factors[i];
      result.maxCorner_[i] *= factors[i];
   }

   result.translate( theCenter );

   return result;
}



/**
 * \brief Creates a new GenericAABB by merging this one with an additional point
 *
 * Note that for the resulting GenericAABB containsClosedInterval( point ) will be true
 * but contains( point ) may be false!
 *
 * \param point The point that will be covered by the resulting GenericAABB
 *
 * \returns A GenericAABB covering point
 *
 * \post containsClosedInterval( point )
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::getMerged( const vector_type & point ) const
{
   return createFromMinMaxCorner( std::min( minCorner_[0], point[0] ),
                                  std::min( minCorner_[1], point[1] ),
                                  std::min( minCorner_[2], point[2] ),
                                  std::max( maxCorner_[0], point[0] ),
                                  std::max( maxCorner_[1], point[1] ),
                                  std::max( maxCorner_[2], point[2] ) );
}



/**
 * \brief Creates a new GenericAABB by merging this one with an other GenericAABB
 *
 * \param other The GenericAABB that will be covered by the resulting GenericAABB
 *
 * \returns A GenericAABB covering other
 *
 * \post contains( other )
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::getMerged( const GenericAABB & other ) const
{
   return createFromMinMaxCorner( std::min( minCorner_[0], other.minCorner_[0] ),
                                  std::min( minCorner_[1], other.minCorner_[1] ),
                                  std::min( minCorner_[2], other.minCorner_[2] ),
                                  std::max( maxCorner_[0], other.maxCorner_[0] ),
                                  std::max( maxCorner_[1], other.maxCorner_[1] ),
                                  std::max( maxCorner_[2], other.maxCorner_[2] ) );
}



/**
 * \brief Creates a new GenericAABB by merging this one to a Sequence of points or other GenericAABBs.
 *
 * \tparam InputIterator A ForwardIterator. Must dereference to Vector3< T > or GenericAABB< T >.
 *
 * \param first First element of sequence [first, last)
 * \param last  Final element of sequence [first, last)
 *
 * \returns A GenericAABB covering *this and all the Points/GenericAABBs in the sequence
 */
template< typename T >
template< typename InputIterator >
GenericAABB< T > GenericAABB< T >::getMerged( InputIterator first, InputIterator last ) const
{
   if( first == last  )
   {
      return *this;
   }
   else
   {
      GenericAABB< T > result( first, last );
      result.merge( *this );
      return result;
   }
}



/**
 * \brief Tests whether this genericAABB intersects another GenericAABB
 *
 * This test interprets the boundaries of the GenericAABB as an half-open interval. The minimal values stored in minCorner
 * are considered to belong to the bounding boxes while the maximal values stored in maxCorner are not.
 *
 * \param other Other GenericAABB to be tested for intersection
 *
 * \returns true, if *this intersects other, false else
 */
template< typename T >
bool GenericAABB< T >::intersects( const GenericAABB & other ) const
{
   return other.maxCorner_[0] > minCorner_[0] && other.minCorner_[0] < maxCorner_[0] &&
          other.maxCorner_[1] > minCorner_[1] && other.minCorner_[1] < maxCorner_[1] &&
          other.maxCorner_[2] > minCorner_[2] && other.minCorner_[2] < maxCorner_[2] && 
          !empty() && !other.empty();
}



/**
 * \brief Tests whether this extended GenericAABB intersects another GenericAABB
 *
 * This test interprets the boundaries of the GenericAABB as an half-open interval. The minimal values stored in minCorner
 * are considered to belong to the bounding boxes while the maximal values stored in maxCorner are not.
 *
 * \param other Other GenericAABB to be tested for intersection
 * \param dx    An epsilon the box is extended by in each direction before the test
 *
 * \returns true, if *this extended by dx intersects other, false else
 */
template< typename T >
bool GenericAABB< T >::intersects( const GenericAABB & other, const value_type dx ) const
{
   WALBERLA_ASSERT_LESS_EQUAL( minCorner_[0] - dx, maxCorner_[0] + dx );
   WALBERLA_ASSERT_LESS_EQUAL( minCorner_[1] - dx, maxCorner_[1] + dx );
   WALBERLA_ASSERT_LESS_EQUAL( minCorner_[2] - dx, maxCorner_[2] + dx );

   return other.maxCorner_[0] > ( minCorner_[0] - dx ) && other.minCorner_[0] < ( maxCorner_[0] + dx ) &&
          other.maxCorner_[1] > ( minCorner_[1] - dx ) && other.minCorner_[1] < ( maxCorner_[1] + dx ) &&
          other.maxCorner_[2] > ( minCorner_[2] - dx ) && other.minCorner_[2] < ( maxCorner_[2] + dx ) &&
          !other.empty() &&
          minCorner_[0] - dx < maxCorner_[0] + dx &&
          minCorner_[1] - dx < maxCorner_[1] + dx &&
          minCorner_[2] - dx < maxCorner_[2] + dx;
}



/**
 * \brief Tests whether this genericAABB intersects another GenericAABB
 *
 * This test interprets the boundaries of the GenericAABB as a closed interval. The minimal values stored in minCorner
 * are considered to belong to the bounding boxes as well as the maximal values stored in maxCorner.
 *
 * \param other Other GenericAABB to be tested for intersection
 *
 * \returns true, if *this intersects other, false else
 */
template< typename T >
bool GenericAABB< T >::intersectsClosedInterval( const GenericAABB & other ) const
{
   return other.maxCorner_[0] >= minCorner_[0] && other.minCorner_[0] <= maxCorner_[0] &&
          other.maxCorner_[1] >= minCorner_[1] && other.minCorner_[1] <= maxCorner_[1] &&
          other.maxCorner_[2] >= minCorner_[2] && other.minCorner_[2] <= maxCorner_[2];
}



/**
 * \brief Tests whether this extended GenericAABB intersects another GenericAABB
 *
 * This test interprets the boundaries of the GenericAABB as a closed interval. The minimal values stored in minCorner
 * are considered to belong to the bounding boxes as well as the maximal values stored in maxCorner.
 *
 * \param other Other GenericAABB to be tested for intersection
 * \param dx    An epsilon the box is extended by in each direction before the test
 *
 * \returns true, if *this extended by dx intersects other, false else
 */
template< typename T >
bool GenericAABB< T >::intersectsClosedInterval( const GenericAABB & other, const value_type dx ) const
{
   return other.maxCorner_[0] >= ( minCorner_[0] - dx ) && other.minCorner_[0] <= ( maxCorner_[0] + dx ) &&
          other.maxCorner_[1] >= ( minCorner_[1] - dx ) && other.minCorner_[1] <= ( maxCorner_[1] + dx ) &&
          other.maxCorner_[2] >= ( minCorner_[2] - dx ) && other.minCorner_[2] <= ( maxCorner_[2] + dx );
}



/**
 * \brief Computes the volume of the intersection of this and another GenericAABB
 *
 * \param other Other GenericAABB to be intersected with *this
 *
 * \returns The volume of *this intersected with other
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::intersectionVolume( const GenericAABB & other ) const
{
   return getIntersection( other ).volume();
}



/**
 * \brief Computes the intersection of this and another GenericAABB
 *
 * \param other Other GenericAABB to be intersected with *this
 *
 * \returns The intersection of *this and other
 */
template< typename T >
GenericAABB< T > GenericAABB< T >::getIntersection( const GenericAABB & other ) const
{
   return GenericAABB( *this, other );
}



/**
 * \brief Test whether this and another GenericAABB are exactly identical
 *
 * \param other Other GenericAABB to be tested for identity
 *
 * \returns true, if *this is exactly identical to other, false else
 */
template< typename T >
bool GenericAABB< T >::isIdentical( const GenericAABB< T > & other ) const
{
   return walberla::isIdentical( minCorner_[0], other.minCorner_[0] ) &&
          walberla::isIdentical( minCorner_[1], other.minCorner_[1] ) &&
          walberla::isIdentical( minCorner_[2], other.minCorner_[2] ) &&
          walberla::isIdentical( maxCorner_[0], other.maxCorner_[0] ) &&
          walberla::isIdentical( maxCorner_[1], other.maxCorner_[1] ) &&
          walberla::isIdentical( maxCorner_[2], other.maxCorner_[2] );
}



/**
 * \brief Test whether this and another GenericAABB are equal
 *
 * Uses walberla::floatIsEqual method for determination of equality. Small roundoff errors are neglected in this check.
 *
 * \param other Other GenericAABB to be tested for equality
 *
 * \returns true, if *this is equal to other, false else
 */
template< typename T >
bool GenericAABB< T >::isEqual( const GenericAABB< T > & other ) const
{
   return walberla::floatIsEqual( minCorner_[0], other.minCorner_[0] ) &&
          walberla::floatIsEqual( minCorner_[1], other.minCorner_[1] ) &&
          walberla::floatIsEqual( minCorner_[2], other.minCorner_[2] ) &&
          walberla::floatIsEqual( maxCorner_[0], other.maxCorner_[0] ) &&
          walberla::floatIsEqual( maxCorner_[1], other.maxCorner_[1] ) &&
          walberla::floatIsEqual( maxCorner_[2], other.maxCorner_[2] );
}



/**
 * \brief Computes the squared distance of point to this GenericAABB
 *
 * \param point The point of which the distance to *this should be computed
 *
 * \returns 0, if point lies inside or on the surface of *this
 * \returns The squared (positive) distance of point to the surface of *this, if point lies outside of *this
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::sqDistance( const vector_type & point ) const
{
   const vector_type d = point - minCorner_;
   const vector_type theSizes = sizes();
   value_type sqDist(0);

   for( uint_t i = 0; i < 3; i++ )
   {
      if( d[i] < 0 )
         sqDist += d[i] * d[i];
      else if( d[i] > theSizes[i] )
      {
         const value_type axisDist = d[i] - theSizes[i];
         sqDist += axisDist * axisDist;
      }
   }

   WALBERLA_ASSERT_GREATER_EQUAL( sqDist, value_type( 0 ) );

   return sqDist;
}



/**
 * \brief Computes the squared signed distance of point to this GenericAABB
 *
 * \param point The point of which the distance to *this should be computed
 *
 * \returns The negative squared distance of point to the surface of *this, if point lies inside of *this
 * \returns 0, if point lies on the surface of *this
 * \returns The positive squared distance of point to the surface of *this, if point lies outside of *this
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::sqSignedDistance( const vector_type & point ) const
{
   const vector_type d        = point - minCorner_;
   const vector_type theSizes = sizes();

   bool inside = d[0] >= value_type( 0 ) && d[0] < theSizes[0] &&
                 d[1] >= value_type( 0 ) && d[1] < theSizes[1] && 
                 d[2] >= value_type( 0 ) && d[2] < theSizes[2];

   if( !inside )
      return sqDistance( point );

   value_type sqAxisDist[3];
   
   for( uint_t i = 0; i < 3; ++i )
   {
      WALBERLA_ASSERT_GREATER_EQUAL( d[i], 0 );
      WALBERLA_ASSERT_LESS( d[i], theSizes[i] );
      
      if( d[i] < theSizes[i] * value_type(0.5) )
      {
         sqAxisDist[i] = d[i] * d[i];
      }
      else
      {
         const value_type axisDist = d[i] - theSizes[i];
         sqAxisDist[i] = axisDist * axisDist;
      }
   }

   return -std::min( sqAxisDist[0], std::min( sqAxisDist[1], sqAxisDist[2] ) );
}



/**
 * \brief Computes the squared maximal distance of point to the far side of this GenericAABB
 *
 * \param point The point of which the maximal distance to *this should be computed
 *
 * \returns The maximal squared (positive) distance of point to the surface of *this.
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::sqMaxDistance( const vector_type & point ) const
{
   const vector_type d0 = point - minCorner_;
   const vector_type d1 = point - maxCorner_;
   value_type sqDist(0);

   for( uint_t i = 0; i < 3; i++ )
   {
      if( std::fabs( d0[i] ) > std::fabs( d1[i] ) )
         sqDist += d0[i] * d0[i];
      else
         sqDist += d1[i] * d1[i];
   }

   WALBERLA_ASSERT_GREATER_EQUAL( sqDist, value_type( 0 ) );

   return sqDist;
}



/**
 * \brief Computes the distance of point to this GenericAABB
 *
 * \param point The point of which the distance to *this should be computed
 *
 * \returns 0, if point lies inside or on the surface of *this
 * \returns The (positive) distance of point to the surface of *this, if point lies outside of *this
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::distance( const vector_type & point ) const
{
   return std::sqrt( sqDistance( point ) );
}



/**
 * \brief Computes the signed distance of point to this GenericAABB
 *
 * \param point The point of which the distance to *this should be computed
 *
 * \returns The negative distance of point to the surface of *this, if point lies inside of *this
 * \returns 0, if point lies on the surface of *this
 * \returns The positive distance of point to the surface of *this, if point lies outside of *this
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::signedDistance( const vector_type & point ) const
{
   const vector_type d        = point - minCorner_;
   const vector_type theSizes = sizes();

   bool inside = d[0] >= value_type( 0 ) && d[0] < theSizes[0] &&
                 d[1] >= value_type( 0 ) && d[1] < theSizes[1] &&
                 d[2] >= value_type( 0 ) && d[2] < theSizes[2];

   if( !inside )
      return distance( point );

   value_type axisDist[3];

   for( uint_t i = 0; i < 3; ++i )
   {
      WALBERLA_ASSERT_GREATER_EQUAL( d[i], 0 );
      WALBERLA_ASSERT_LESS( d[i], theSizes[i] );

      if( d[i] < theSizes[i] * value_type( 0.5 ) )
      {
         axisDist[i] = d[i];
      }
      else
      {
         axisDist[i] = theSizes[i] - d[i];
      }
   }

   return -std::min( axisDist[0], std::min( axisDist[1], axisDist[2] ) );
}



/**
 * \brief Computes the maximal distance of point to the far side of this GenericAABB
 *
 * \param point The point of which the maximal distance to *this should be computed
 *
 * \returns The maximal (positive) distance of point to the surface of *this.
 */
template< typename T >
typename GenericAABB< T >::value_type GenericAABB< T >::maxDistance( const vector_type & point ) const
{
   return std::sqrt( sqMaxDistance( point ) );
}


/**
* \brief Computes the distance between two GenericAABBs
*
* \param other The other AABB to which the distance to *this should be computed
*
* \returns The (positive) distance of the surface of other to the surface of *this. 0 if *this and other are intersecting.
*/
template< typename T >
inline typename GenericAABB< T >::value_type GenericAABB< T >::sqDistance( const GenericAABB & other ) const
{
   value_type theSqDistance = value_type(0);

   if( other.maxCorner_[0] < minCorner_[0] )
   {
      theSqDistance += sq( minCorner_[0] - other.maxCorner_[0] );
   }
   else if( other.minCorner_[0] > maxCorner_[0] )
   {
      theSqDistance += sq( other.minCorner_[0] - maxCorner_[0] );
   }

   if( other.maxCorner_[1] < minCorner_[1] )
   {
      theSqDistance += sq( minCorner_[1] - other.maxCorner_[1] );
   }
   else if( other.minCorner_[1] > maxCorner_[1] )
   {
      theSqDistance += sq( other.minCorner_[1] - maxCorner_[1] );
   }

   if( other.maxCorner_[2] < minCorner_[2] )
   {
      theSqDistance += sq( minCorner_[2] - other.maxCorner_[2] );
   }
   else if( other.minCorner_[2] > maxCorner_[2] )
   {
      theSqDistance += sq( other.minCorner_[2] - maxCorner_[2] );
   }

   WALBERLA_ASSERT_GREATER_EQUAL( theSqDistance, value_type(0) );
   WALBERLA_ASSERT( !intersects( other ) || walberla::isIdentical( theSqDistance, value_type(0) ) ); // intersect => distance == 0

   return theSqDistance;
}


/**
* \brief Computes the maximal distance of any two points from two GenericAABBs
*
* \param other The other AABB to which the maximal distance to *this should be computed
*
* \returns The maximal (positive) distance of any point in other to any point in this.
*/
template< typename T >
inline typename GenericAABB< T >::value_type GenericAABB< T >::sqMaxDistance( const GenericAABB & other ) const
{
   value_type theSqMaxDistance = value_type(0);

   theSqMaxDistance += sq( std::max( maxCorner_[0] - other.minCorner_[0], other.maxCorner_[0] - minCorner_[0] ) );
   theSqMaxDistance += sq( std::max( maxCorner_[1] - other.minCorner_[1], other.maxCorner_[1] - minCorner_[1] ) );
   theSqMaxDistance += sq( std::max( maxCorner_[2] - other.minCorner_[2], other.maxCorner_[2] - minCorner_[2] ) );

   WALBERLA_ASSERT_GREATER_EQUAL( theSqMaxDistance, value_type(0) );

   return theSqMaxDistance;
}



/**
 * \brief Computes the eight corners of this GenericAABB
 *
 * \returns An array of size eight containing the corner points. Indices match those of the stencil::D3CornerStencil
 */
template< typename T >
inline std::array< typename GenericAABB< T >::vector_type, 8 > GenericAABB< T >::corners() const
{
   std::array< vector_type, 8 > cornerArray;

   // TNE
   cornerArray[0][0] = maxCorner_[0];
   cornerArray[0][1] = maxCorner_[1];
   cornerArray[0][2] = maxCorner_[2];

   // TNW
   cornerArray[1][0] = minCorner_[0];
   cornerArray[1][1] = maxCorner_[1];
   cornerArray[1][2] = maxCorner_[2];

   // TSE
   cornerArray[2][0] = maxCorner_[0];
   cornerArray[2][1] = minCorner_[1];
   cornerArray[2][2] = maxCorner_[2];

   // TSW
   cornerArray[3][0] = minCorner_[0];
   cornerArray[3][1] = minCorner_[1];
   cornerArray[3][2] = maxCorner_[2];

   // BNE
   cornerArray[4][0] = maxCorner_[0];
   cornerArray[4][1] = maxCorner_[1];
   cornerArray[4][2] = minCorner_[2];

   // BNW
   cornerArray[5][0] = minCorner_[0];
   cornerArray[5][1] = maxCorner_[1];
   cornerArray[5][2] = minCorner_[2];

   // BSE
   cornerArray[6][0] = maxCorner_[0];
   cornerArray[6][1] = minCorner_[1];
   cornerArray[6][2] = minCorner_[2];

   // BSW
   cornerArray[7][0] = minCorner_[0];
   cornerArray[7][1] = minCorner_[1];
   cornerArray[7][2] = minCorner_[2];

   return cornerArray;
}



/**
 * \brief Reinitializes this GenericAABB with default values
 *
 * \post minCorner = maxCorner = (0,0,0)
 */
template< typename T >
void GenericAABB< T >::init()
{
   minCorner_.set( value_type( 0 ), value_type( 0 ), value_type( 0 ) );
   maxCorner_.set( value_type( 0 ), value_type( 0 ), value_type( 0 ) );

   WALBERLA_ASSERT( checkInvariant() );
   WALBERLA_ASSERT( empty() );
}



/**
 * \brief Reinitializes this GenericAABB from two arbitrary points
 *
 * \param corner0 first corner spanning the AABB
 * \param corner1 second corner spanning the AABB
 */
template< typename T >
void GenericAABB< T >::init( const vector_type & corner0, const vector_type & corner1 )
{
   for( uint_t i = 0; i < 3; ++i )
   {
      if( corner0[i] < corner1[i] )
      {
         minCorner_[i] = corner0[i];
         maxCorner_[i] = corner1[i];
      }
      else
      {
         minCorner_[i] = corner1[i];
         maxCorner_[i] = corner0[i];
      }
   }

   WALBERLA_ASSERT( checkInvariant() );
}



/**
 * \brief Reinitializes this GenericAABB from two arbitrary points
 *
 * \param x0  x-coordinate of the first corner spanning the AABB
 * \param y0  y-coordinate of the first corner spanning the AABB
 * \param z0  z-coordinate of the first corner spanning the AABB
 * \param x1  x-coordinate of the second corner spanning the AABB
 * \param y1  y-coordinate of the second corner spanning the AABB
 * \param z1  z-coordinate of the second corner spanning the AABB
 */
template< typename T >
void GenericAABB< T >::init( const value_type x0, const value_type y0, const value_type z0,
                             const value_type x1, const value_type y1, const value_type z1 )
{
   if( x0 < x1 )
   {
      minCorner_[0] = x0;
      maxCorner_[0] = x1;
   }
   else
   {
      minCorner_[0] = x1;
      maxCorner_[0] = x0;
   }

   if( y0 < y1 )
   {
      minCorner_[1] = y0;
      maxCorner_[1] = y1;
   }
   else
   {
      minCorner_[1] = y1;
      maxCorner_[1] = y0;
   }

   if( z0 < z1 )
   {
      minCorner_[2] = z0;
      maxCorner_[2] = z1;
   }
   else
   {
      minCorner_[2] = z1;
      maxCorner_[2] = z0;
   }

   WALBERLA_ASSERT( checkInvariant() );
}



/**
 * \brief Reinitializes this GenericAABB by merging a Sequence of points or GenericAABBs
 *
 * \tparam InputIterator ForwardIterator. Must dereference to Vector3< T > or GenericAABB< T >.
 *
 * \param first first element of sequence [first, last)
 * \param last  final element of sequence [first, last)
 *
 * \post first == last => minCorner = maxCorner = (0,0,0)
 */
template< typename T >
template< typename InputIterator >
void GenericAABB< T >::init( InputIterator first, InputIterator last )
{
   if( first != last )
   {
      minCorner_.set(  std::numeric_limits<value_type>::max(),  std::numeric_limits<value_type>::max(),  std::numeric_limits<value_type>::max() );
      maxCorner_.set( -std::numeric_limits<value_type>::max(), -std::numeric_limits<value_type>::max(), -std::numeric_limits<value_type>::max() );

      merge( first, last );
   }
   else
   {
      init();
   }

   WALBERLA_ASSERT( checkInvariant() );
}



/**
 * \brief Reinitializes this GenericAABB from a min- and a maxCorner
 *
 * \param minCorner corner with minimal values
 * \param maxCorner corner with maximal values
 *
 * \pre minCorner.x <= maxCorner.x && minCorner.y <= maxCorner.y && minCorner.z <= maxCorner.z
 */
template< typename T >
void GenericAABB< T >::initMinMaxCorner( const value_type minX, const value_type minY, const value_type minZ,
                                         const value_type maxX, const value_type maxY, const value_type maxZ )
{
   minCorner_[0] = minX;
   minCorner_[1] = minY;
   minCorner_[2] = minZ;

   maxCorner_[0] = maxX;
   maxCorner_[1] = maxY;
   maxCorner_[2] = maxZ;

   WALBERLA_ASSERT( checkInvariant() );
}



/**
 * \brief Reinitializes this GenericAABB from a min- and a maxCorner
 *
 * \param minX  x-coordinate of the minCorner
 * \param minY  y-coordinate of the minCorner
 * \param minZ  z-coordinate of the minCorner
 * \param maxX  x-coordinate of the maxCorner
 * \param maxY  y-coordinate of the maxCorner
 * \param maxZ  z-coordinate of the maxCorner
 *
 * \pre x0 <= x1 && y0 <= y1 && z0 <= z1
 */
template< typename T >
void GenericAABB< T >::initMinMaxCorner( const vector_type & theMinCorner, const vector_type & theMaxCorner )
{
   minCorner_ = theMinCorner;
   maxCorner_ = theMaxCorner;

   WALBERLA_ASSERT( checkInvariant() );
}


/**
 * \brief Sets the minimum and maximum for one axis
 *
 * \param index   0 for x, 1 for y, 2 for z axis
 * \param value*  the smaller of the two values is taken as minimum, the other as maximum
 */
template< typename T >
void GenericAABB< T >::setAxisBounds( const uint_t index, const value_type value1, const value_type value2 )
{
   WALBERLA_ASSERT_LESS( index, 3 );

   if ( value1 < value2 ) {
      minCorner_[index] = value1;
      maxCorner_[index] = value2;
   } else {
      minCorner_[index] = value2;
      maxCorner_[index] = value1;
   }
   WALBERLA_ASSERT( checkInvariant() );
}


/**
 * \brief Extends this GenericAABB
 *
 * \param e epsilon by which the bounding box is extended in each direction
 */
template< typename T >
void GenericAABB< T >::extend( const value_type e )
{
   minCorner_[0] -= e;
   minCorner_[1] -= e;
   minCorner_[2] -= e;

   maxCorner_[0] += e;
   maxCorner_[1] += e;
   maxCorner_[2] += e;

   for( uint_t i = 0; i < 3; ++i )
      if( minCorner_[i] > maxCorner_[i] )
         maxCorner_[i] = minCorner_[i];

   WALBERLA_ASSERT( checkInvariant() );
}



/**
 * \brief Extends this GenericAABB
 *
 * \param e epsilon vector by which the bounding box is extended. The box is extended in each direction by
 *          the corresponding vector component.
 */
template< typename T >
void GenericAABB< T >::extend( const vector_type & e )
{
   minCorner_ -= e;
   maxCorner_ += e;

   for( uint_t i = 0; i < 3; ++i )
      if( minCorner_[i] > maxCorner_[i] )
         maxCorner_[i] = minCorner_[i];

   WALBERLA_ASSERT( checkInvariant() );
}



/**
 * \brief Sets center of this GenericAABB
 *
 * AABB gets translated such that its center matches the given center.
 *
 * \param new center location
 */
template< typename T >
void GenericAABB< T >::setCenter( const vector_type & center )
{
   translate(center - this->center());
}



/**
 * \brief Translates this GenericAABB
 *
 * \param translation Vector by which the GenericAABB is translated.
 */
template< typename T >
void GenericAABB< T >::translate( const vector_type & translation )
{
   minCorner_ += translation;
   maxCorner_ += translation;
}



/**
 * \brief Scales this GenericAABB
 *
 * \param factor Factor by which the bounding box gets scaled.
 */
template< typename T >
void GenericAABB< T >::scale( const value_type factor )
{
   WALBERLA_ASSERT( factor > value_type( 0 ) );

   const vector_type theCenter = center();

   translate( -theCenter );

   minCorner_ *= factor;
   maxCorner_ *= factor;

   translate( theCenter );
}



/**
 * \brief Scales this GenericAABB
 *
 * \param factor Vector of scaling factors by which the bounding box gets scaled on the
 *               respective axises.
 */
template< typename T >
void GenericAABB< T >::scale( const vector_type & factors )
{
   WALBERLA_ASSERT( factors[0] > value_type( 0 ) );
   WALBERLA_ASSERT( factors[1] > value_type( 0 ) );
   WALBERLA_ASSERT( factors[2] > value_type( 0 ) );

   const vector_type theCenter = center();

   translate( -theCenter );

   for( uint_t i = 0; i < 3; ++i )
   {
      minCorner_[i] *= factors[i];
      maxCorner_[i] *= factors[i];
   }

   translate( theCenter );
}



/**
 * \brief Merges this GenericAABB with an additional point
 *
 * Note that for the resulting GenericAABB containsClosedInterval( point ) will be true
 * but contains( point ) may be false!
 *
 * \param point The point that will be covered by the resulting GenericAABB
 *
 * \post containsClosedInterval( point )
 */
template< typename T >
void GenericAABB< T >::merge( const vector_type & point )
{
   for( uint_t i = 0; i < 3; ++i )
   {
      if( point[i] < minCorner_[i] )
         minCorner_[i] = point[i];
      if( point[i] > maxCorner_[i] )
         maxCorner_[i] = point[i];
   }
}



/**
 * \brief Merges this GenericAABB with an other GenericAABB
 *
 * \param other The GenericAABB that will be covered by the resulting GenericAABB
 *
 * \post contains( other )
 */
template< typename T >
void GenericAABB< T >::merge( const GenericAABB & other )
{
   for( uint_t i = 0; i < 3; ++i )
   {
      if( other.minCorner_[i] < minCorner_[i] )
         minCorner_[i] = other.minCorner_[i];
      if( other.maxCorner_[i] > maxCorner_[i] )
         maxCorner_[i] = other.maxCorner_[i];
   }
}



/**
 * \brief Merges this GenericAABB with a Sequence of points or other GenericAABBs.
 *
 * \tparam InputIterator A ForwardIterator. Must dereference to Vector3< T > or GenericAABB< T >.
 *
 * \param first First element of sequence [first, last)
 * \param last  Final element of sequence [first, last)
 */
template< typename T >
template< typename InputIterator >
void GenericAABB< T >::merge( InputIterator first, InputIterator last )
{
   while( first != last )
      merge( *first++ );
}



/**
 * \brief Intersects this with another GenericAABB
 *
 * \param other Other GenericAABB to be intersected with *this
 */
template< typename T >
void GenericAABB< T >::intersect( const GenericAABB & other )
{
   minCorner_[0] = std::max( minCorner_[0], other.minCorner_[0] );
   minCorner_[1] = std::max( minCorner_[1], other.minCorner_[1] );
   minCorner_[2] = std::max( minCorner_[2], other.minCorner_[2] );

   maxCorner_[0] = std::max( minCorner_[0], std::min( maxCorner_[0], other.maxCorner_[0] ) );
   maxCorner_[1] = std::max( minCorner_[1], std::min( maxCorner_[1], other.maxCorner_[1] ) );
   maxCorner_[2] = std::max( minCorner_[2], std::min( maxCorner_[2], other.maxCorner_[2] ) );

   WALBERLA_ASSERT( checkInvariant() );
}



/**
 * \brief Generates a random point uniformly distributed within the AABB
 *
 * The point is in ( [ xMin(), xMax() ), [ yMin(), yMax() ), [ zMin(), zMax() ) )
 *
 * \pre !empty()
 * \param engine  An Uniform Random Number Generator (e.g. std::mt19937)
 * \returns Random point within *this
 */
template< typename T >
template< typename Engine >
typename GenericAABB< T >::vector_type GenericAABB< T >::randomPoint( Engine & engine ) const
{
   WALBERLA_ASSERT( !empty() );
   std::uniform_real_distribution< T > randX( xMin(), xMax() );
   std::uniform_real_distribution< T > randY( yMin(), yMax() );
   std::uniform_real_distribution< T > randZ( zMin(), zMax() );

   return vector_type( randX( engine ), randY( engine ), randZ( engine ) );
}


/**
 * \brief Tests whether the class invariant is satisfied
 *
 * The invariant is minCorner.x <= maxCorner.x && minCorner.y <= maxCorner.y && minCorner.z <= maxCorner.z
 *
 * \returns true if the invariant is satisfied, false else
 */
template< typename T >
bool GenericAABB< T >::checkInvariant() const
{
   return minCorner_[0] <= maxCorner_[0] && minCorner_[1] <= maxCorner_[1] && minCorner_[2] <= maxCorner_[2];
}



/**
 * \brief Compares two GenericeAABBs for equality
 *
 * see isEqual
 *
 * ATTENTION: Does not check for identity! If you need to check for identity,
 *            use member function isIdentical of class GenericAABB
 *
 * \param lhs An arbitrary GenericAABB
 * \param rhs An arbitrary GenericAABB
 *
 * \returns true, if lhs is equal to rhs, false else
 */
template< typename T, typename U >
bool operator==( const GenericAABB< T > & lhs, const GenericAABB< U > & rhs )
{
   return lhs.isEqual( rhs );
}



/**
 * \brief Compares two GenericeAABBs for inequality
 *
 * see isEqual
 *
 * ATTENTION: Does not check for non-identity! If you need to check for non-identity,
 *            use member function isIdentical of class GenericAABB
 *
 * \param lhs An arbitrary GenericAABB
 * \param rhs An arbitrary GenericAABB
 *
 * \returns false, if lhs is equal to rhs, true else
 */
template< typename T, typename U >
bool operator!=( const GenericAABB< T > & lhs, const GenericAABB< U > & rhs )
{
   return !lhs.isEqual( rhs );
}



/**
 * \brief Writes a GenericAABB to a std::ostream
 *
 * The format is [ <0,0,0>, <1,1,1> ]
 *
 * \param os   The std::ostream written to
 * \param aabb The GenericAABB written
 *
 * \returns A reference to os
 */
template< typename T >
std::ostream& operator<<( std::ostream& os, const GenericAABB< T > & aabb )
{
   return os << "[ " << aabb.minCorner() << ", " << aabb.maxCorner() << " ]";
}



/**
 * \brief Reads an GenericAABB from a std::istream
 *
 * The format is [ <0,0,0>, <1,1,1> ]
 *
 * \param is   The std::istream read from
 * \param aabb The GenericAABB the read values are to be stored in
 *
 * \returns A reference to is
 */
template< typename T >
std::istream& operator>>( std::istream& is, GenericAABB< T > & aabb )
{
   if( !is ) return is;

   char bracket0, bracket1;
   char comma;
   typename GenericAABB< T >::vector_type corner0, corner1;

   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );

   // Setting the 'skip whitespaces' flag
   is >> std::skipws;

   // Extracting the vector
   if( !( is >> bracket0 >> corner0 >> comma >> corner1 >> bracket1 ) ||
       comma != ',' || bracket0 != '[' || bracket1 != ']' )
   {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      is.flags( oldFlags );
      return is;
   }

   // Transferring the input to the aabb values
   aabb.init( corner0, corner1 );

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}



/**
 * \brief Serializes an GenericAABB into a mpi::GenericSendBuffer
 *
 * \param buf  The mpi::GenericSendBuffer written to
 * \param aabb The GenericAABB the be serialized
 *
 * \returns A reference to buf
 */
template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename VT >  // Element type of GenericAABB
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const GenericAABB< VT > & aabb )
{
   buf.addDebugMarker( "bb" );
   static_assert ( std::is_trivially_copyable< GenericAABB< VT > >::value,
                   "type has to be trivially copyable for the memcpy to work correctly" );
   auto pos = buf.forward(sizeof(GenericAABB< VT >));
   std::memcpy(pos, &aabb, sizeof(GenericAABB< VT >));
   return buf;
}



} // namespace math
} // namespace walberla
