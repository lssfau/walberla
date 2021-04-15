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
#include "mesh_common/TriangleMeshes.h"

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
void rotate( MeshType& mesh, Vector3<typename  MeshType::Scalar > axis, typename MeshType::Scalar angle, Vector3< typename MeshType::scalar> axis_foot);

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
void rotate( MeshType& mesh, Vector3<typename  MeshType::Scalar > axis, typename MeshType::Scalar angle, Vector3< typename MeshType::Scalar> axis_foot)
{
    Matrix3< typename MeshType::Scalar > mat(axis, angle);

    for( auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        auto &p = mesh.point(*v_it);
        p -= mesh::toOpenMesh(axis_foot);
        p = mat*p;
        p += mesh::toOpenMesh(axis_foot);
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
   // Inertia tensor is calculated relative to the origin of the meshes coordinate system!
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

/**
 * \brief Computes all mass properties (mass, centroid, inertia tensor at once).
 *
 * This function computes the mass, centroid and the inertia tensor relative to the calculated centroid.
 * Source: https://www.cs.upc.edu/~virtual/SGI/docs/3.%20Further%20Reading/Polyhedral%20Mass%20Properties%20Revisited.pdf
 *
 * \tparam MeshType
 * \param mesh      Triangular input mesh.
 * \param density   Density of the mesh.
 * \param centroid  Output centroid point.
 * \param inertiaTensor Output inertia matrix.
 * \param mass      Output mass.
 * \attention The inertia tensor is computed relative to the centroid.
 */
template< typename MeshType >
void computeMassProperties(const MeshType & mesh, typename MeshType::Scalar density,
      Vector3<typename MeshType::Scalar>& centroid, Matrix3<typename MeshType::Scalar>& inertiaTensor,
      typename MeshType::Scalar& mass) {
   static_assert( MeshType::IsTriMesh == 1, "computeMassProperties only works with triangular meshes!" );

   typedef typename MeshType::Point Point;
   typedef typename MeshType::Scalar Scalar;

   const Scalar mult[10] = {Scalar(1)/Scalar(6),
                            Scalar(1)/Scalar(24), Scalar(1)/Scalar(24), Scalar(1)/Scalar(24),
                            Scalar(1)/Scalar(60), Scalar(1)/Scalar(60), Scalar(1)/Scalar(60),
                            Scalar(1)/Scalar(120), Scalar(1)/Scalar(120), Scalar(1)/Scalar(120)};

   Scalar intg[10] = {0,0,0,0,0,0,0,0,0,0};

   auto subExpr = [](Scalar& w0, Scalar& w1, Scalar& w2,
         Scalar& f1, Scalar& f2, Scalar& f3,
         Scalar& g0, Scalar& g1, Scalar& g2){
      Scalar temp0 = w0+w1;
      f1 = temp0 + w2;
      Scalar temp1 = w0*w0;
      Scalar temp2 = temp1 + w1*temp0;
      f2 = temp2 + w2*f1;
      f3 = w0*temp1 + w1*temp2 + w2*f2;
      g0 = f2 + w0*(f1+w0);
      g1 = f2 + w1*(f1+w1);
      g2 = f2 + w2*(f1+w2);
   };

   Point v0, v1, v2;
   for (auto fIt = mesh.faces_begin(); fIt != mesh.faces_end(); ++fIt) {
      getVertexPositions(mesh, *fIt, v0, v1, v2);

      Scalar a1 = v1[0]-v0[0];
      Scalar b1 = v1[1]-v0[1];
      Scalar c1 = v1[2]-v0[2];
      Scalar a2 = v2[0]-v0[0];
      Scalar b2 = v2[1]-v0[1];
      Scalar c2 = v2[2]-v0[2];

      Scalar d0 = b1*c2 - b2*c1;
      Scalar d1 = a2*c1 - a1*c2;
      Scalar d2 = a1*b2 - a2*b1;

      Scalar f1x, f2x, f3x, g0x, g1x, g2x;
      subExpr(v0[0], v1[0], v2[0], f1x, f2x, f3x, g0x, g1x, g2x);
      Scalar f1y, f2y, f3y, g0y, g1y, g2y;
      subExpr(v0[1], v1[1], v2[1], f1y, f2y, f3y, g0y, g1y, g2y);
      Scalar f1z, f2z, f3z, g0z, g1z, g2z;
      subExpr(v0[2], v1[2], v2[2], f1z, f2z, f3z, g0z, g1z, g2z);

      intg[0] += d0*f1x;
      intg[1] += d0*f2x; intg[2] += d1*f2y; intg[3] += d2*f2z;
      intg[4] += d0*f3x; intg[5] += d1*f3y; intg[6] += d2*f3z;
      intg[7] += d0*(v0[1]*g0x + v1[1]*g1x + v2[1]*g2x);
      intg[8] += d1*(v0[2]*g0y + v1[2]*g1y + v2[2]*g2y);
      intg[9] += d2*(v0[0]*g0z + v1[0]*g1z + v2[0]*g2z);
   }

   for (uint_t i = 0; i < 10; ++i) {
      intg[i] *= mult[i];
   }

   mass = intg[0];

   centroid[0] = intg[1] / mass;
   centroid[1] = intg[2] / mass;
   centroid[2] = intg[3] / mass;

   inertiaTensor[0] = density * (intg[5] + intg[6] - mass*(centroid[1]*centroid[1] + centroid[2]*centroid[2])); //xx
   inertiaTensor[4] = density * (intg[4] + intg[6] - mass*(centroid[2]*centroid[2] + centroid[0]*centroid[0])); // yy
   inertiaTensor[8] = density * (intg[4] + intg[5] - mass*(centroid[0]*centroid[0] + centroid[1]*centroid[1])); // zz
   inertiaTensor[1] = density * (-(intg[7] - mass * centroid[0]*centroid[1]));
   inertiaTensor[5] = density * (-(intg[8] - mass * centroid[1]*centroid[2]));
   inertiaTensor[2] = density * (-(intg[9] - mass * centroid[2]*centroid[0]));

   mass *= density;
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

      // The point closest to the box center is located in the spherical shell between the box's insphere and circumsphere
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

      // The point closest to the box center is located in the spherical shell between the box's insphere and circumsphere
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

                      // The point closest to the box center is located in the spherical shell between the box's insphere and circumsphere
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
