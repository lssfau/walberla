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
//! \file GenericAABBTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/GenericAABB.h"
#include "core/math/Shims.h"
#include "core/math/Utility.h"

#include "stencil/D3CornerStencil.h"

#include <random>
#include <array>


using namespace walberla;
using math::GenericAABB;

template< typename T >
void testEmptyAABB( const GenericAABB< T > & aabb )
{
   using VecT = Vector3<T>;

   WALBERLA_CHECK_EQUAL( aabb.minCorner(), VecT( T(0), T(0), T(0) ) );
   WALBERLA_CHECK_EQUAL( aabb.maxCorner(), VecT( T(0), T(0), T(0) ) );
   WALBERLA_CHECK( aabb.empty() );

   WALBERLA_CHECK_EQUAL( aabb.sizes(), VecT( T(0), T(0), T(0) ) );
   for( uint_t i = 0; i < 3; ++i )
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.size( i ), T(0) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.xSize(), T(0) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.ySize(), T(0) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.zSize(), T(0) );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb.volume(), T(0) );
   WALBERLA_CHECK_EQUAL( aabb.center(), VecT( T(0), T(0), T(0) ) );

   WALBERLA_CHECK( !aabb.contains( VecT( T(0), T(0), T(0) ) ) );
   WALBERLA_CHECK( aabb.contains( VecT( T(0), T(0), T(0) ), std::numeric_limits<T>::epsilon() ) );
   WALBERLA_CHECK( aabb.contains( VecT( T(-1), T(-1), T(-1) ), T(1) ) );
   WALBERLA_CHECK( !aabb.contains( VecT( T(1), T(1), T(1) ), T(1) ) );

   auto corners = aabb.corners();
   using namespace stencil;
   for( auto dir = D3CornerStencil::begin(); dir != D3CornerStencil::end(); ++dir )
   {
      WALBERLA_CHECK( !aabb.contains( corners[ dir.toIdx() ] ) );
      WALBERLA_CHECK( aabb.containsClosedInterval( corners[ dir.toIdx() ] ) );
   }
}

template< typename T >
void testNonEmptyAABB( const GenericAABB< T > & aabb )
{
   WALBERLA_CHECK( aabb.contains( aabb.minCorner() ) );
   WALBERLA_CHECK( aabb.contains( aabb.center() ) );
   WALBERLA_CHECK( !aabb.contains( aabb.maxCorner() ) );

   auto corners = aabb.corners();
   using namespace stencil;
   for( auto dir = D3CornerStencil::begin(); dir != D3CornerStencil::end(); ++dir )
   {
      if( *dir == BSW )
      {
         WALBERLA_CHECK( aabb.contains( corners[ dir.toIdx() ] ) );
      }
      else
      {
         WALBERLA_CHECK( !aabb.contains( corners[ dir.toIdx() ] ) );
      }

      WALBERLA_CHECK( aabb.containsClosedInterval( corners[ dir.toIdx() ] ) );
   }

   WALBERLA_CHECK( aabb.intersects( aabb ) );

   GenericAABB< T > intersectingBox( aabb.minCorner() + aabb.sizes() * T(0.5), aabb.maxCorner() + aabb.sizes() * T(0.5) );
   WALBERLA_CHECK( aabb.intersects( intersectingBox ) );
   WALBERLA_CHECK( intersectingBox.intersects( aabb ) );
   WALBERLA_CHECK( aabb.intersectsClosedInterval( intersectingBox ) );
   WALBERLA_CHECK( intersectingBox.intersectsClosedInterval( aabb ) );
   GenericAABB< T > tmpAABB = aabb;
   tmpAABB.intersect( intersectingBox );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getIntersection( intersectingBox ) );
   WALBERLA_CHECK_EQUAL( aabb.getIntersection( intersectingBox ), tmpAABB );
   WALBERLA_CHECK_FLOAT_EQUAL( tmpAABB.volume(), aabb.volume() / T(8) );
   WALBERLA_CHECK_FLOAT_EQUAL( tmpAABB.volume(), aabb.intersectionVolume( intersectingBox ) );

   intersectingBox.init( aabb.minCorner() + aabb.sizes(), aabb.maxCorner() + aabb.sizes() );
   tmpAABB = aabb;
   tmpAABB.intersect( intersectingBox );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getIntersection( intersectingBox ) );
   WALBERLA_CHECK_EQUAL( intersectingBox.getIntersection( aabb ), aabb.getIntersection( intersectingBox ) );
   WALBERLA_CHECK_FLOAT_EQUAL( tmpAABB.volume(), T(0) );
   WALBERLA_CHECK_IDENTICAL( tmpAABB.volume(), aabb.intersectionVolume( intersectingBox ) );

   std::mt19937 urng;
   for( int i = 0; i < 100; ++i )
   {
      auto p = aabb.randomPoint( urng );
      WALBERLA_CHECK( aabb.contains( p ) );
   }
}

template< typename T >
void testAnyAABB( const GenericAABB< T > & aabb )
{
   using VecT = Vector3<T>;


   WALBERLA_CHECK_IDENTICAL( aabb.minCorner()[0], aabb.xMin() );
   WALBERLA_CHECK_IDENTICAL( aabb.minCorner()[1], aabb.yMin() );
   WALBERLA_CHECK_IDENTICAL( aabb.minCorner()[2], aabb.zMin() );
   WALBERLA_CHECK_IDENTICAL( aabb.maxCorner()[0], aabb.xMax() );
   WALBERLA_CHECK_IDENTICAL( aabb.maxCorner()[1], aabb.yMax() );
   WALBERLA_CHECK_IDENTICAL( aabb.maxCorner()[2], aabb.zMax() );

   mpi::SendBuffer sendBuffer;
   sendBuffer << aabb;
   mpi::RecvBuffer receiveBuffer( sendBuffer );
   GenericAABB< T > aabb0;
   receiveBuffer >> aabb0;
   WALBERLA_CHECK_EQUAL( aabb, aabb0 );

   std::stringstream ss;
   ss << aabb;
   GenericAABB< T > aabb1;
   ss >> aabb1;
   WALBERLA_CHECK( !ss.fail() );

   WALBERLA_CHECK( aabb.contains( aabb ) );
   WALBERLA_CHECK( !aabb.contains( aabb.getExtended( T(1) ) ) );
   WALBERLA_CHECK( !aabb.contains( aabb.getTranslated( VecT(T(1), T(-1), T(1) ) ) ) );

   WALBERLA_CHECK_EQUAL( aabb, GenericAABB< T >::createFromMinMaxCorner( aabb.minCorner(), aabb.maxCorner() ) );

   WALBERLA_CHECK_EQUAL( aabb, aabb.getIntersection( aabb ) );

   GenericAABB< T > tmpAABB = aabb;
   tmpAABB.extend( T(1) );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getExtended( T(1) ) );

   tmpAABB = aabb;
   tmpAABB.extend( T(-1) );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getExtended( T(-1) ) );

   tmpAABB = aabb;
   tmpAABB.translate( VecT( T(-1), T(0), T(1) ) );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getTranslated( VecT( T(-1), T(0), T(1) ) ) );
   WALBERLA_CHECK( !aabb.isIdentical( tmpAABB ) );
   WALBERLA_CHECK( !aabb.isEqual( tmpAABB ) );

   tmpAABB = aabb;
   tmpAABB.scale( T(2) );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getScaled( VecT( T(2) ) ) );

   tmpAABB = aabb;
   tmpAABB.scale( VecT( T(1), T(2), T(3) ) );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getScaled( VecT( T(1), T(2), T(3) ) ) );

   tmpAABB = aabb;
   VecT newPoint = tmpAABB.minCorner() - VecT( T(1), T(2), T(3) );
   tmpAABB.merge( newPoint );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getMerged( newPoint ) );
   WALBERLA_CHECK( tmpAABB.contains( aabb ) );
   WALBERLA_CHECK( tmpAABB.contains( newPoint ) );
   WALBERLA_CHECK( tmpAABB.containsClosedInterval( newPoint ) );
   WALBERLA_CHECK( !aabb.isIdentical( tmpAABB ) );
   WALBERLA_CHECK( !aabb.isEqual( tmpAABB ) );

   tmpAABB = aabb;
   newPoint = tmpAABB.maxCorner() + VecT( T(1), T(2), T(3) );
   tmpAABB.merge( newPoint );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getMerged( newPoint ) );
   WALBERLA_CHECK( tmpAABB.contains( aabb ) );
   WALBERLA_CHECK( !tmpAABB.contains( newPoint ) );
   WALBERLA_CHECK( tmpAABB.containsClosedInterval( newPoint ) );
   WALBERLA_CHECK( !aabb.isIdentical( tmpAABB ) );
   WALBERLA_CHECK( !aabb.isEqual( tmpAABB ) );

   tmpAABB = aabb;
   GenericAABB< T > newBox( tmpAABB.maxCorner(), tmpAABB.maxCorner() + VecT( T(1), T(2), T(3) ) );
   tmpAABB.merge( newBox );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getMerged( newBox ) );
   WALBERLA_CHECK( tmpAABB.contains( aabb ) );
   WALBERLA_CHECK( tmpAABB.contains( newBox ) );
   WALBERLA_CHECK( !aabb.isIdentical( tmpAABB ) );
   WALBERLA_CHECK( !aabb.isEqual( tmpAABB ) );

   std::vector< VecT > points;
   points.push_back( aabb.maxCorner() + VecT( T(1), T(2), T(3) ) );
   points.push_back( aabb.minCorner() + VecT( T(1), T(2), T(3) ) );
   tmpAABB = aabb;
   tmpAABB.merge( points.begin(), points.end() );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getMerged( points.begin(), points.end() ) );
   WALBERLA_CHECK( tmpAABB.contains( aabb ) );
   WALBERLA_CHECK( tmpAABB.containsClosedInterval( points[0] ) );
   WALBERLA_CHECK( tmpAABB.containsClosedInterval( points[1] ) );
   WALBERLA_CHECK( !aabb.isIdentical( tmpAABB ) );
   WALBERLA_CHECK( !aabb.isEqual( tmpAABB ) );

   std::vector< GenericAABB< T > > boxes;
   boxes.push_back( GenericAABB< T >( aabb.minCorner() + VecT( T(1), T(2), T(3) ), aabb.maxCorner() + VecT( T(1), T(2), T(3) ) ) );
   boxes.push_back( GenericAABB< T >( aabb.minCorner() - VecT( T(1), T(2), T(3) ), aabb.maxCorner() - VecT( T(1), T(2), T(3) ) ) );
   tmpAABB = aabb;
   tmpAABB.merge( boxes.begin(), boxes.end() );
   WALBERLA_CHECK_EQUAL( tmpAABB, aabb.getMerged( boxes.begin(), boxes.end() ) );
   WALBERLA_CHECK( tmpAABB.contains( aabb ) );
   WALBERLA_CHECK( tmpAABB.contains( boxes[0] ) );
   WALBERLA_CHECK( tmpAABB.contains( boxes[1] ) );
   WALBERLA_CHECK( !aabb.isIdentical( tmpAABB ) );
   WALBERLA_CHECK( !aabb.isEqual( tmpAABB ) );

   WALBERLA_CHECK( aabb.intersectsClosedInterval( aabb ) );
   WALBERLA_CHECK_EQUAL( aabb.getIntersection( aabb ), aabb );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.intersectionVolume( aabb ), aabb.volume() );

   WALBERLA_CHECK( aabb.isIdentical( aabb ) );
   WALBERLA_CHECK( aabb.isEqual( aabb ) );
   GenericAABB< T > oneBox( VecT( T(1), T(1), T(1) ), VecT( T(1), T(1), T(1) ) );
   GenericAABB< T > epsilonShiftedBox = oneBox.getTranslated( VecT( std::numeric_limits<T>::epsilon() ) );
   WALBERLA_CHECK( !oneBox.isIdentical( epsilonShiftedBox ) );
   WALBERLA_CHECK( oneBox.isEqual( epsilonShiftedBox ) );

   using math::sqr;

   WALBERLA_CHECK_IDENTICAL( aabb.sqDistance( aabb.center() ), T(0) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.sqMaxDistance( aabb.center() ), ( aabb.minCorner() - aabb.center() ).sqrLength() );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.sqSignedDistance( aabb.center() ), -math::sqr( std::min( aabb.xSize(), std::min( aabb.ySize(), aabb.zSize() ) ) * T(0.5) ) );
   WALBERLA_CHECK_IDENTICAL( aabb.distance( aabb.center() ), T(0) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.maxDistance( aabb.center() ), ( aabb.minCorner() - aabb.center() ).length() );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.signedDistance( aabb.center() ), -std::min( aabb.xSize(), std::min( aabb.ySize(), aabb.zSize() ) ) * T( 0.5 ) );

   auto corners = aabb.corners();
   for( auto it = corners.begin(); it != corners.end(); ++it )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.sqDistance( *it ), T(0) );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.sqMaxDistance( *it ), sqr( aabb.xSize() ) + sqr( aabb.ySize() ) + sqr( aabb.zSize() ) );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.sqSignedDistance( *it ), T(0) );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.distance( *it ), T(0) );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.maxDistance( *it ), std::sqrt( sqr( aabb.xSize() ) + sqr( aabb.ySize() ) + sqr( aabb.zSize() ) ) );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.signedDistance( *it ), T(0) );

      VecT ep = *it + ( *it - aabb.center() );
      T sqDist = ( *it - aabb.center() ).sqrLength();
      T dist = std::sqrt( sqDist );

      WALBERLA_CHECK_FLOAT_EQUAL( aabb.sqDistance( ep ), sqDist );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.sqMaxDistance( ep ), sqr( T(3) * dist ) );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.sqSignedDistance( ep ), sqDist );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.distance( ep ), dist );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.maxDistance( ep ), T(3) * dist );
      WALBERLA_CHECK_FLOAT_EQUAL( aabb.signedDistance( ep ), dist );
   }
}

template< typename T >
void testFixedAABB()
{
   using VecT = Vector3<T>;

   GenericAABB<T> aabb0;
   testEmptyAABB( aabb0 );
   testAnyAABB( aabb0 );

   {
      std::stringstream ss;
      ss << aabb0;
      GenericAABB< T > tmpAABB;
      ss >> tmpAABB;
      WALBERLA_CHECK_EQUAL( aabb0, tmpAABB );
   }


   GenericAABB<T> aabb1( VecT( T(1), T(1), T(5) ), VecT( T(0), T(3), T(2) ) );

   testNonEmptyAABB( aabb1 );
   testAnyAABB( aabb1 );

   {
      std::stringstream ss;
      ss << aabb1;
      GenericAABB< T > tmpAABB;
      ss >> tmpAABB;
      WALBERLA_CHECK_EQUAL( aabb1, tmpAABB );
   }

   WALBERLA_CHECK_EQUAL( aabb1.minCorner(), VecT( T(0), T(1), T(2) ) );
   WALBERLA_CHECK_EQUAL( aabb1.maxCorner(), VecT( T(1), T(3), T(5) ) );
   WALBERLA_CHECK( !aabb1.empty() );

   WALBERLA_CHECK_EQUAL( aabb1.sizes(), VecT( T(1), T(2), T(3) ) );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb1.size( 0 ), T(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb1.size( 1 ), T(2) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb1.size( 2 ), T(3) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb1.xSize(), T(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb1.ySize(), T(2) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb1.zSize(), T(3) );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb1.volume(), T(6) );

   WALBERLA_CHECK_EQUAL( aabb1.center(), VecT( T(0.5), T(2), T(3.5) ) );

   WALBERLA_CHECK( aabb1.contains( aabb1.maxCorner(), T(5) * std::numeric_limits<T>::epsilon() ) );
   WALBERLA_CHECK( aabb1.contains( VecT( T(-1), T(0), T(1) ), T(1) ) );
   WALBERLA_CHECK( !aabb1.contains( VecT( T(2), T(4), T(6) ), T(1) ) );

   aabb1.init();
   testEmptyAABB( aabb1 );
}

template< typename T >
void testConstructors( const T x0, const T y0, const T z0, const T x1, const T y1, const T z1 )
{
   using VecT = Vector3<T>;

   {
      GenericAABB< T > refAABB;
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], T(0) );
   }

   {
      GenericAABB< T > refAABB( x0, y0, z0, x1, y1, z1 );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > refAABB( VecT( x0, y0, z0 ), VecT( x1, y1, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > toBeCopied( x0, y0, z0, x1, y1, z1 );
      const GenericAABB< T > refAABB( toBeCopied ); // NOLINT
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > toBeCopied( x0, y0, z0, x1, y1, z1 );
      GenericAABB< T > refAABB;
      refAABB = toBeCopied;
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      std::vector< VecT > points;
      points.push_back( VecT( x0, y0, z0 ) );
      points.push_back( VecT( x1, y1, z1 ) );
      GenericAABB< T > refAABB( points.begin(), points.end() );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > toBeCopied( x0, y0, z0, x1, y1, z1 );
      GenericAABB< T > refAABB( &toBeCopied, &toBeCopied + 1u );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > refAABB = GenericAABB< T >::createFromMinMaxCorner(
         std::min( x0, x1 ), std::min( y0, y1 ), std::min( z0, z1 ),
         std::max( x0, x1 ), std::max( y0, y1 ), std::max( z0, z1 ) );

      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > refAABB = GenericAABB< T >::createFromMinMaxCorner(
         VecT( std::min( x0, x1 ), std::min( y0, y1 ), std::min( z0, z1 ) ),
         VecT( std::max( x0, x1 ), std::max( y0, y1 ), std::max( z0, z1 ) ) );

      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > refAABB( x0, y0, z0, x1, y1, z1 );
      refAABB.init();

      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], T(0) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], T(0) );
   }

   {
      GenericAABB< T > refAABB;
      refAABB.init( x0, y0, z0, x1, y1, z1 );

      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > refAABB;
      refAABB.init( VecT( x0, y0, z0 ), VecT( x1, y1, z1 ) );

      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      std::vector< VecT > points;
      points.push_back( VecT( x0, y0, z0 ) );
      points.push_back( VecT( x1, y1, z1 ) );
      GenericAABB< T > refAABB;
      refAABB.init( points.begin(), points.end() );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > toBeCopied( x0, y0, z0, x1, y1, z1 );
      GenericAABB< T > refAABB;
      refAABB.init( &toBeCopied, &toBeCopied + 1u );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > refAABB;
      refAABB.initMinMaxCorner(
         std::min( x0, x1 ), std::min( y0, y1 ), std::min( z0, z1 ),
         std::max( x0, x1 ), std::max( y0, y1 ), std::max( z0, z1 ) );

      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }

   {
      GenericAABB< T > refAABB;
      refAABB.initMinMaxCorner(
         VecT( std::min( x0, x1 ), std::min( y0, y1 ), std::min( z0, z1 ) ),
         VecT( std::max( x0, x1 ), std::max( y0, y1 ), std::max( z0, z1 ) ) );

      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[0], std::min( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[1], std::min( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.minCorner()[2], std::min( z0, z1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[0], std::max( x0, x1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[1], std::max( y0, y1 ) );
      WALBERLA_CHECK_IDENTICAL( refAABB.maxCorner()[2], std::max( z0, z1 ) );
   }
}

template< typename T >
void randomTest()
{
   typedef std::mersenne_twister_engine< walberla::uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 0xa37d3c92 > mt11213b;
   mt11213b rng;

   std::uniform_real_distribution<T> dist( -T(10), T(10) );

   for( int i = 0; i < 1000; ++i )
   {
      T minX = dist(rng);
      T minY = dist(rng);
      T minZ = dist(rng);
      T maxX = dist(rng);
      T maxY = dist(rng);
      T maxZ = dist(rng);

      testConstructors( minX, minY, minZ, maxX, maxY, maxZ );

      GenericAABB< T > aabb( minX, minY, minZ, maxX, maxY, maxZ );

      testAnyAABB( aabb );
      if( aabb.empty() )
         testEmptyAABB( aabb );
      else
         testNonEmptyAABB( aabb );
   }
}


template< typename T >
void testAABBDistancesFixed()
{
   GenericAABB< T > aabb0( T(-1), T(-1), T(-1), T(1), T(1), T(1) );
   GenericAABB< T > aabb1( T(-0.5), T(-0.5), T(-0.5), T(0.5), T(0.5), T(0.5) );
      
   WALBERLA_CHECK_IDENTICAL( aabb0.sqDistance( aabb1 ), T(0) );
   WALBERLA_CHECK_IDENTICAL( aabb1.sqDistance( aabb0 ), T(0) );

   WALBERLA_CHECK_IDENTICAL( aabb0.sqMaxDistance( aabb1 ), T(3) * math::sq( T(1.5) ) );
   WALBERLA_CHECK_IDENTICAL( aabb1.sqMaxDistance( aabb0 ), T(3) * math::sq( T(1.5) ) );

   aabb1.init( T(1.5), T(1.5), T(1.5), T(2), T(2), T(2) );

   WALBERLA_CHECK_IDENTICAL( aabb0.sqDistance( aabb1 ), T(3) * math::sq( T(0.5) ) );
   WALBERLA_CHECK_IDENTICAL( aabb1.sqDistance( aabb0 ), T(3) * math::sq( T(0.5) ) );

   WALBERLA_CHECK_IDENTICAL( aabb0.sqMaxDistance( aabb1 ), T(3) * math::sq( T(3) ) );
   WALBERLA_CHECK_IDENTICAL( aabb1.sqMaxDistance( aabb0 ), T(3) * math::sq( T(3) ) );


   aabb1.init( T(0), T(0), T(0), T(2), T(2), T(2) );

   WALBERLA_CHECK_IDENTICAL( aabb0.sqDistance( aabb1 ), T(0) );
   WALBERLA_CHECK_IDENTICAL( aabb1.sqDistance( aabb0 ), T(0) );

   WALBERLA_CHECK_IDENTICAL( aabb0.sqMaxDistance( aabb1 ), T(3) * math::sq( T(3) ) );
   WALBERLA_CHECK_IDENTICAL( aabb1.sqMaxDistance( aabb0 ), T(3) * math::sq( T(3) ) );

   aabb1.init( T(-2), T(-2), T(-2), T(0), T(0), T(0) );

   WALBERLA_CHECK_IDENTICAL( aabb0.sqDistance( aabb1 ), T(0) );
   WALBERLA_CHECK_IDENTICAL( aabb1.sqDistance( aabb0 ), T(0) );

   WALBERLA_CHECK_IDENTICAL( aabb0.sqMaxDistance( aabb1 ), T(3) * math::sq( T(3) ) );
   WALBERLA_CHECK_IDENTICAL( aabb1.sqMaxDistance( aabb0 ), T(3) * math::sq( T(3) ) );
}


template< typename T >
void testAABBDistancesRandom( const GenericAABB< T > & baseAABB )
{
   static const uint_t NUM_BOXES  = 100;
   static const uint_t NUM_POINTS = 1000;
   std::mt19937 rng;

   for( uint_t i = 0; i < NUM_BOXES; ++i )
   {
      math::GenericAABB< T > subAabb0(  baseAABB.randomPoint( rng ), baseAABB.randomPoint( rng )  );
      math::GenericAABB< T > subAabb1(  baseAABB.randomPoint( rng ), baseAABB.randomPoint( rng )  );

      WALBERLA_CHECK_IDENTICAL( subAabb0.sqDistance( subAabb1 ), subAabb1.sqDistance( subAabb0 ) );
      WALBERLA_CHECK_IDENTICAL( subAabb0.sqMaxDistance( subAabb1 ), subAabb1.sqMaxDistance( subAabb0 ) );

      const T minSqDistance = subAabb0.sqDistance( subAabb1 );
      const T maxSqDistance = subAabb0.sqMaxDistance( subAabb1 );

      WALBERLA_CHECK_GREATER_EQUAL( maxSqDistance, minSqDistance );

      for( uint_t j = 0; j < NUM_POINTS; ++j )
      {
         auto p0 = subAabb0.randomPoint( rng );
         auto p1 = subAabb1.randomPoint( rng );

         const auto sqPointDistance = (p0 - p1).sqrLength();

         WALBERLA_CHECK_GREATER_EQUAL( sqPointDistance, minSqDistance );
         WALBERLA_CHECK_LESS_EQUAL( sqPointDistance, maxSqDistance );
      }
   }
}


int main(int argc, char**argv)
{
   walberla::debug::enterTestMode();

   walberla::Environment env( argc, argv );

   testConstructors<float>( 1.0f, 2.0f, 3.0f, -1.0f, 3.0f, 2.0f );
   testConstructors<double>( 1.0, 2.0, 3.0, -1.0, 3.0, 2.0 );
   
   testFixedAABB<float>();
   testFixedAABB<double>();
   
   GenericAABB< float > floatAABB( 1.0f, 1.0f, 1.0f, 2.0f, 2.0f, 2.0f );
   GenericAABB< double > doubleAABB( 1.0, 1.0, 1.0, 2.0, 2.0 ,2.0 );
   
   GenericAABB< double > copiedAABB0( floatAABB );
   const GenericAABB< double > copiedAABB1( doubleAABB );

   WALBERLA_UNUSED(copiedAABB0);
   WALBERLA_UNUSED(copiedAABB1);
   
   randomTest<float>();
   randomTest<double>();

   testAABBDistancesFixed<float>();
   testAABBDistancesFixed<double>();

   testAABBDistancesRandom<float>( GenericAABB<float>( -1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f ) );
   testAABBDistancesRandom<double>( GenericAABB<double>( -1.0, -1.0, -1.0, 1.0, 1.0, 1.0 ) );

   return 0;
}
