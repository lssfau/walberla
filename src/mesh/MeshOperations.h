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
//! \file MeshOperations.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "MatrixVectorOperations.h"
#include "TriangleMeshes.h"

#include "core/math/GenericAABB.h"
#include "core/Optional.h"

#include <set>
#include <iterator>
#include <queue>

namespace walberla {
namespace mesh {

template< typename MeshType >
math::GenericAABB< typename MeshType::Scalar > computeAABB( const MeshType & mesh );

template< typename MeshType, typename InputIterator >
math::GenericAABB< typename MeshType::Scalar > computeAABBForVertices( const MeshType & mesh, InputIterator beginVh, InputIterator endVh );

template< typename MeshType, typename InputIterator >
math::GenericAABB< typename MeshType::Scalar > computeAABBForFaces( const MeshType & mesh, InputIterator beginFh, InputIterator endFh );

template< typename MeshType >
void translate( MeshType & mesh, const Vector3< typename MeshType::Scalar > & offset );

template< typename MeshType >
void scale( MeshType & mesh, const Vector3< typename MeshType::Scalar > & scaleFactors );

template< typename MeshType >
typename MeshType::Scalar computeVolume( const MeshType & mesh );

template< typename MeshType >
typename MeshType::Scalar computeSurfaceArea( const MeshType & mesh );

template< typename MeshType >
typename MeshType::Point computeCentroid( const MeshType & mesh );

template< typename MeshType >
Matrix3<typename MeshType::Scalar> computeInertiaTensor( const MeshType & mesh );

template< typename MeshType >
typename MeshType::Point computeCentroid( const MeshType & mesh, const typename MeshType::FaceHandle fh );

template< typename MeshType, typename InputIterator >
std::vector< typename MeshType::VertexHandle > findConnectedVertices( const MeshType & mesh, InputIterator facesBegin, InputIterator facesEnd );

template< typename MeshType >
void findConnectedVertices( const MeshType & mesh, const typename MeshType::FaceHandle & face, std::vector< typename MeshType::VertexHandle > & outVertices );

template< typename DistanceObject, typename T, typename U >
walberla::optional< bool > isIntersecting( const DistanceObject & distanceObject, const math::GenericAABB< T > & aabb, const T & maxError );

template< typename DistanceObject, typename T, typename U >
walberla::optional< bool > fullyCoversAABB( const DistanceObject & distanceObject, const math::GenericAABB< T > & aabb, const T & maxError );

template< typename DistanceObject, typename T, typename U, typename V >
walberla::optional< bool > intersectsSurface( const DistanceObject & distanceObject, const math::GenericAABB< T > & aabb, const U & maxError, const V & surfaceDistance );

template< typename MeshType, typename InputIterator >
typename MeshType::Point principalComponent( const MeshType & mesh, InputIterator beginFh, InputIterator endFh, const uint_t iterations = uint_t(10) );

template< typename MeshType >
typename MeshType::Point principalComponent( const MeshType & mesh, const uint_t iterations = uint_t( 10 ) );


template< typename MeshType >
math::GenericAABB< typename MeshType::Scalar > computeAABB( const MeshType & mesh )
{
   return computeAABBForVertices( mesh, mesh.vertices_begin(), mesh.vertices_end() );
}


template< typename MeshType, typename InputIterator >
math::GenericAABB< typename MeshType::Scalar > computeAABBForVertices( const MeshType & mesh, InputIterator beginVh, InputIterator endVh )
{
   if( beginVh == endVh )
      return math::GenericAABB< typename MeshType::Scalar >();

   auto vIt = beginVh;

   auto & pFirst = mesh.point( *vIt++ );

   Vector3< typename MeshType::Scalar > min, max;
   min[0] = max[0] = pFirst[0];
   min[1] = max[1] = pFirst[1];
   min[2] = max[2] = pFirst[2];

   while( vIt != endVh )
   {
      auto & p = mesh.point( *vIt++ );

      min[0] = std::min( min[0], p[0] );
      min[1] = std::min( min[1], p[1] );
      min[2] = std::min( min[2], p[2] );

      max[0] = std::max( max[0], p[0] );
      max[1] = std::max( max[1], p[1] );
      max[2] = std::max( max[2], p[2] );
   }

   return math::GenericAABB< typename MeshType::Scalar >::createFromMinMaxCorner( min, max );
}

template< typename MeshType, typename InputIterator >
math::GenericAABB< typename MeshType::Scalar > computeAABBForFaces( const MeshType & mesh, InputIterator beginFh, InputIterator endFh )
{
   if( beginFh == endFh )
      return math::GenericAABB< typename MeshType::Scalar >();

   auto fIt = beginFh;

   const auto & pFirst = mesh.point( *(mesh.cfv_iter( *fIt )) );

   Vector3< typename MeshType::Scalar > min, max;
   min[0] = max[0] = pFirst[0];
   min[1] = max[1] = pFirst[1];
   min[2] = max[2] = pFirst[2];

   while( fIt != endFh )
   {
      for( auto vIt = mesh.cfv_begin( *fIt ); vIt != mesh.cfv_end( *fIt ); ++vIt )
      {
         const auto & p = mesh.point( *vIt );

         min[0] = std::min( min[0], p[0] );
         min[1] = std::min( min[1], p[1] );
         min[2] = std::min( min[2], p[2] );

         max[0] = std::max( max[0], p[0] );
         max[1] = std::max( max[1], p[1] );
         max[2] = std::max( max[2], p[2] );
      }
      ++fIt;
   }

   return math::GenericAABB< typename MeshType::Scalar >::createFromMinMaxCorner( min, max );
}


template< typename MeshType >
void translate( MeshType & mesh, const Vector3< typename MeshType::Scalar > & offset )
{
   const auto o = toOpenMesh( offset );

   for( auto vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt )
      mesh.point( *vIt ) += o;
}


template< typename MeshType >
void scale( MeshType & mesh, const Vector3< typename MeshType::Scalar > & scaleFactors )
{
   for( auto vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt )
   {
      auto & p = mesh.point( *vIt );
      p[0] *= scaleFactors[0];
      p[1] *= scaleFactors[1];
      p[2] *= scaleFactors[2];
   }
}


template< typename MeshType >
typename MeshType::Scalar computeVolume( const MeshType & mesh )
{
   static_assert( MeshType::IsTriMesh == 1, "computeVolume only works with triangular meshes!" );

   typedef typename MeshType::Point  Point;
   typedef typename MeshType::Scalar Scalar;

   Point v0, v1, v2;
   Scalar result(0);

   for( auto fIt = mesh.faces_begin(); fIt != mesh.faces_end(); ++fIt )
   {
      getVertexPositions( mesh, *fIt, v0, v1, v2 );
      result += ( v0 | ( v1 % v2 ) ) / Scalar(6);
   }

   return std::fabs( result );
}


template< typename MeshType >
typename MeshType::Scalar computeSurfaceArea( const MeshType & mesh )
{
   static_assert( MeshType::IsTriMesh == 1, "computeSurfaceArea only works with triangular meshes!" );

   typedef typename MeshType::Point  Point;
   typedef typename MeshType::Scalar Scalar;

   Point v0, v1, v2;
   Scalar result(0);

   for( auto fIt = mesh.faces_begin(); fIt != mesh.faces_end(); ++fIt )
   {
      getVertexPositions( mesh, *fIt, v0, v1, v2 );
      result += ( ( v1 - v0 ) % ( v2 - v0 ) ).length();
   }

   return result * Scalar( 0.5 );
}


template< typename MeshType >
typename MeshType::Point computeCentroid( const MeshType & mesh )
{
   static_assert( MeshType::IsTriMesh == 1, "computeCentroid only works with triangular meshes!" );

   typedef typename MeshType::Point  Point;
   typedef typename MeshType::Scalar Scalar;

   Point v0, v1, v2;
   Point centroid(Scalar(0), Scalar(0), Scalar(0));
   Scalar volume(0);

   for( auto fIt = mesh.faces_begin(); fIt != mesh.faces_end(); ++fIt )
   {
      getVertexPositions( mesh, *fIt, v0, v1, v2 );
      Scalar tetraVolume = ( v0 | ( v1 % v2 ) ) / Scalar(6);
      // v0 + v1 + v2 + 0 / 4 is the tetraeder centroid
      centroid += tetraVolume * (v0 + v1 + v2);
      volume += tetraVolume;
   }

   return centroid / ( Scalar(4) * volume );
}


template< typename MeshType >
Matrix3<typename MeshType::Scalar> computeInertiaTensor( const MeshType & mesh )
{
   static_assert( MeshType::IsTriMesh == 1, "computeInertiaTensor only works with triangular meshes!" );

   typedef typename MeshType::Point  Point;
   typedef typename MeshType::Scalar Scalar;

   Point v0, v1, v2;
   Scalar p00, p01, p02, p11, p12, p22;
   p00 = p01 = p02 = p11 = p12 = p22 = Scalar(0);

   auto deltaP_jk = [&]( uint_t j, uint_t k ) {
      return Scalar(2) * v0[j]*v0[k]
           + Scalar(2) * v1[j]*v1[k]
           + Scalar(2) * v2[j]*v2[k]
           + v0[j] * v1[k] + v0[k] * v1[j]
           + v0[j] * v2[k] + v0[k] * v2[j]
           + v1[j] * v2[k] + v1[k] * v2[j];
   };

   for( auto fIt = mesh.faces_begin(); fIt != mesh.faces_end(); ++fIt )
   {
      getVertexPositions( mesh, *fIt, v0, v1, v2 );
      Scalar tetraVolume = ( v0 | ( v1 % v2 ) );
      p00 += tetraVolume * deltaP_jk(0,0);
      p01 += tetraVolume * deltaP_jk(0,1);
      p02 += tetraVolume * deltaP_jk(0,2);
      p11 += tetraVolume * deltaP_jk(1,1);
      p12 += tetraVolume * deltaP_jk(1,2);
      p22 += tetraVolume * deltaP_jk(2,2);
   }

   p00 /= Scalar(20 * 6);
   p01 /= Scalar(20 * 6);
   p02 /= Scalar(20 * 6);
   p11 /= Scalar(20 * 6);
   p12 /= Scalar(20 * 6);
   p22 /= Scalar(20 * 6);

   return Matrix3<Scalar>(  p11 + p22, -p01,       -p02,
                           -p01,        p00 + p22, -p12,
                           -p02,       -p12,        p00 + p11 );
}


template< typename MeshType >
typename MeshType::Point computeCentroid( const MeshType & mesh, const typename MeshType::FaceHandle fh )
{
   typedef typename MeshType::Point Point;
   typedef typename MeshType::Scalar Scalar;

   Point centroid = Point( Scalar(0), Scalar(0), Scalar(0) );
   uint_t numVertices(0);
   for(auto vh : mesh.fv_range(fh))
   {
      centroid += mesh.point(vh);
      ++numVertices;
   }

   centroid /= numeric_cast<Scalar>(numVertices);

   return centroid;
}


template< typename MeshType, typename InputIterator >
std::vector< typename MeshType::VertexHandle > findConnectedVertices( const MeshType & mesh, const InputIterator facesBegin, const InputIterator facesEnd )
{
   std::vector< typename MeshType::VertexHandle > vertexHandles;

   for( auto it = facesBegin; it != facesEnd; ++it )
      findConnectedVertices( mesh, *it, vertexHandles );

   std::sort( vertexHandles.begin(), vertexHandles.end() );
   vertexHandles.erase( std::unique( vertexHandles.begin(), vertexHandles.end() ), vertexHandles.end() );

   return vertexHandles;
}


template< typename MeshType >
void findConnectedVertices( const MeshType & mesh, const typename MeshType::FaceHandle & face, std::vector< typename MeshType::VertexHandle > & outVertices )
{
   std::copy( mesh.cfv_begin( face ), mesh.cfv_end( face ), std::back_inserter( outVertices ) );
}



template< typename DistanceObject, typename T, typename U >
walberla::optional< bool > isIntersecting( const DistanceObject & distanceObject, const math::GenericAABB< T > & aabb, const U & maxError )
{
   typedef typename DistanceObject::Scalar Scalar;

   if( aabb.empty() )
      return false;

   Scalar maxErrorScalar = numeric_cast<Scalar>( maxError      );

   typedef math::GenericAABB< Scalar > Box;

   Box rootAABB( numeric_cast<Scalar>( aabb.xMin() ),
                 numeric_cast<Scalar>( aabb.yMin() ),
                 numeric_cast<Scalar>( aabb.zMin() ),
                 numeric_cast<Scalar>( aabb.xMax() ),
                 numeric_cast<Scalar>( aabb.yMax() ),
                 numeric_cast<Scalar>( aabb.zMax() ) );

   std::queue< Box > boxQueue;
   boxQueue.push( rootAABB );

   while( !boxQueue.empty() )
   {
      const Box & curAabb = boxQueue.front();
      const Scalar minEdgeLength = std::min( curAabb.xSize(), std::min( curAabb.ySize(), curAabb.zSize() ) );
      const Scalar circumRadius = curAabb.sizes().length() * Scalar(0.5);
      const Scalar sqCircumRadius = circumRadius * circumRadius;
      const Scalar inRadius = minEdgeLength * Scalar(0.5);
      const Scalar sqInRadius = inRadius * inRadius;

      const Scalar sqSignedDistance = distanceObject.sqSignedDistance( toOpenMesh( curAabb.center() ) );


      if( sqSignedDistance > sqCircumRadius )
      {
         boxQueue.pop();
         continue; // clearly the mesh does not intersect this box
      }

      if( sqSignedDistance < sqInRadius )
         return true; // clearly the box must intersect the mesh

      // The point closest to the box center is located in the spherical shell between the box's insphere and ciricumsphere
      const Scalar error = circumRadius - inRadius;

      if( error < maxErrorScalar )
      {
         return walberla::nullopt; // we still don't know if there is an intersection but the error margin is already small enough
      }

      const auto &    min = curAabb.minCorner();
      const auto &    max = curAabb.maxCorner();
      const auto   center = curAabb.center();

      boxQueue.push( Box::createFromMinMaxCorner(    min[0],    min[1],    min[2], center[0], center[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0],    min[1], center[2], center[0], center[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0], center[1],    min[2], center[0],    max[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0], center[1], center[2], center[0],    max[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0],    min[1],    min[2],    max[0], center[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0],    min[1], center[2],    max[0], center[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], center[1],    min[2],    max[0],    max[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], center[1], center[2],    max[0],    max[1],    max[2] ) );

      boxQueue.pop();
   }

   return false;
}



template< typename DistanceObject, typename T, typename U >
walberla::optional< bool > fullyCoversAABB( const DistanceObject & distanceObject, const math::GenericAABB< T > & aabb, const U & maxError )
{
   typedef typename DistanceObject::Scalar Scalar;

   if( aabb.empty() )
      return false;

   Scalar maxErrorScalar = numeric_cast<Scalar>( maxError );

   typedef math::GenericAABB< Scalar > Box;

   Box rootAABB( numeric_cast<Scalar>( aabb.xMin() ),
                 numeric_cast<Scalar>( aabb.yMin() ),
                 numeric_cast<Scalar>( aabb.zMin() ),
                 numeric_cast<Scalar>( aabb.xMax() ),
                 numeric_cast<Scalar>( aabb.yMax() ),
                 numeric_cast<Scalar>( aabb.zMax() ) );

   std::queue< Box > boxQueue;
   boxQueue.push( rootAABB );

   while( !boxQueue.empty() )
   {
      const Box & curAabb = boxQueue.front();
      const Scalar minEdgeLength = std::min( curAabb.xSize(), std::min( curAabb.ySize(), curAabb.zSize() ) );
      const Scalar circumRadius = curAabb.sizes().length() * Scalar(0.5);
      const Scalar sqCircumRadius = circumRadius * circumRadius;
      const Scalar inRadius = minEdgeLength * Scalar(0.5);
      const Scalar sqInRadius = inRadius * inRadius;

      const Scalar sqSignedDistance = distanceObject.sqSignedDistance( toOpenMesh( curAabb.center() ) );

      if( sqSignedDistance < -sqCircumRadius )
      {
         boxQueue.pop();
         continue; // clearly the box is fully covered by the mesh
      }

      if( sqSignedDistance > -sqInRadius )
         return false; // clearly the box must partially be outside of the mesh

      // The point closest to the box center is located in the spherical shell between the box's insphere and ciricumsphere
      const Scalar error = circumRadius - inRadius;

      if( error < maxErrorScalar )
      {
         return walberla::nullopt; // we still don't know if there is an intersection but the error margin is already small enough
      }

      const auto &    min = curAabb.minCorner();
      const auto &    max = curAabb.maxCorner();
      const auto   center = curAabb.center();

      boxQueue.push( Box::createFromMinMaxCorner(    min[0],    min[1],    min[2], center[0], center[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0],    min[1], center[2], center[0], center[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0], center[1],    min[2], center[0],    max[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0], center[1], center[2], center[0],    max[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0],    min[1],    min[2],    max[0], center[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0],    min[1], center[2],    max[0], center[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], center[1],    min[2],    max[0],    max[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], center[1], center[2],    max[0],    max[1],    max[2] ) );

      boxQueue.pop();
   }

   return true;
}



template< typename DistanceObject, typename T, typename U, typename V >
walberla::optional< bool > intersectsSurface( const DistanceObject & distanceObject, const math::GenericAABB< T > & aabb, const U & maxError, const V & surfaceDistance )
{
   typedef typename DistanceObject::Scalar Scalar;

   if(aabb.empty())
      return false;

   Scalar maxErrorScalar        = numeric_cast<Scalar>( maxError         );
   Scalar surfaceDistanceScalar = numeric_cast<Scalar>( surfaceDistance );

   typedef math::GenericAABB< Scalar > Box;

   Box rootAABB( numeric_cast<Scalar>( aabb.xMin() ),
      numeric_cast<Scalar>( aabb.yMin() ),
      numeric_cast<Scalar>( aabb.zMin() ),
      numeric_cast<Scalar>( aabb.xMax() ),
      numeric_cast<Scalar>( aabb.yMax() ),
      numeric_cast<Scalar>( aabb.zMax() ) );

   std::queue< Box > boxQueue;
   boxQueue.push( rootAABB );

   while(!boxQueue.empty())
   {
      const Box & curAabb = boxQueue.front();
      const Scalar minEdgeLength = std::min( curAabb.xSize(), std::min( curAabb.ySize(), curAabb.zSize() ) );
      const Scalar circumRadius = curAabb.sizes().length() * Scalar( 0.5 ) + surfaceDistanceScalar;
      const Scalar sqCircumRadius = circumRadius * circumRadius;
      const Scalar inRadius = minEdgeLength * Scalar( 0.5 ) + surfaceDistanceScalar;
      const Scalar sqInRadius = inRadius * inRadius;

      const Scalar sqDistance = distanceObject.sqDistance( toOpenMesh( curAabb.center() ) );


      if( sqDistance > sqCircumRadius )
      {
         boxQueue.pop();
         continue; // clearly the mesh does not intersect this box
      }

      if( sqDistance < sqInRadius)
         return true; // clearly the box must intersect the mesh

                      // The point closest to the box center is located in the spherical shell between the box's insphere and ciricumsphere
      const Scalar error = circumRadius - inRadius;

      if(error < maxErrorScalar)
      {
         return walberla::nullopt; // we still don't know if there is an intersection but the error margin is already small enough
      }

      const auto &    min = curAabb.minCorner();
      const auto &    max = curAabb.maxCorner();
      const auto   center = curAabb.center();

      boxQueue.push( Box::createFromMinMaxCorner( min[0], min[1], min[2], center[0], center[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( min[0], min[1], center[2], center[0], center[1], max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( min[0], center[1], min[2], center[0], max[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( min[0], center[1], center[2], center[0], max[1], max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], min[1], min[2], max[0], center[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], min[1], center[2], max[0], center[1], max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], center[1], min[2], max[0], max[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], center[1], center[2], max[0], max[1], max[2] ) );

      boxQueue.pop();
   }

   return false;
}



template< typename MeshType, typename InputIterator >
typename MeshType::Point principalComponent( const MeshType & mesh, InputIterator beginFh, InputIterator endFh, const uint_t iterations )
{
   typedef typename MeshType::Point Point;
   typedef typename MeshType::Scalar Scalar;
   typedef typename MeshType::VertexHandle VertexHandle;

   Point r( numeric_cast<Scalar>(1), numeric_cast<Scalar>(-1), numeric_cast<Scalar>(2) );

   std::vector< VertexHandle > vertices = findConnectedVertices( mesh, beginFh, endFh );
   Point centroid( Scalar(0), Scalar(0), Scalar(0) );
   uint_t ctr = 0;
   for( auto vh : vertices )
   {
      centroid += mesh.point( vh );
      ++ctr;
   }
   centroid /= numeric_cast<Scalar>( ctr );

   std::vector< Point > decenteredPoints( vertices.size() );
   auto vIt = vertices.cbegin();
   auto pIt = decenteredPoints.begin();
   while(vIt != vertices.end())
   {
      *pIt = mesh.point( *vIt ) - centroid;
      ++pIt;
      ++vIt;
   }

   for(uint_t i = 0; i < iterations; ++i)
   {
      Point s( 0, 0, 0 );

      for( const auto & x : decenteredPoints )
      {
         s += ( x | r ) * x;
      }
      const Scalar sLength = s.length();
      r = s / sLength;
   }

   return r;
}

template< typename MeshType >
typename MeshType::Point principalComponent( const MeshType & mesh, const uint_t iterations )
{
   return principalComponent( mesh, mesh.faces_begin(), mesh.faces_end(), iterations );
}



} // namespace mesh
} // namespace walberla
