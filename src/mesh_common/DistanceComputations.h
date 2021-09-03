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
//! \file DistanceComputations.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"
#include "core/debug/Debug.h"
#include "core/math/GenericAABB.h"
#include "core/math/Matrix3.h"
#include "core/logging/Logging.h"

#include "mesh_common/MatrixVectorOperations.h"
#include "mesh_common/TriangleMeshes.h"

#include <OpenMesh/Core/Utils/PropertyManager.hh>

namespace walberla {
namespace mesh {

template< typename MeshType >
struct DistanceProperties
{
   typedef OpenMesh::VectorT< typename MeshType::Scalar, 2 > Vec2;
   typedef OpenMesh::VectorT< typename MeshType::Scalar, 3 > Vec3;
   typedef typename MeshType::Scalar Scalar;
   typedef math::Matrix3<Scalar> Matrix;

   // Dummy constructor to suppress GCC 7 warnings
   DistanceProperties() : e0(real_t(0)), e1(real_t(0)), e2(real_t(0)),
                          e1_normal(real_t(0)), e2_normal(real_t(0)),
                          e1_normalized(real_t(0)), e2_normalized(real_t(0)), e0_normalized(real_t(0)),
                          e0l(real_t(0)), e1l(real_t(0)), e2l(real_t(0)),
                          translation(real_t(0))
   {}

   Vec2 e0, e1, e2;
   Vec2 e1_normal, e2_normal;
   Vec2 e1_normalized, e2_normalized, e0_normalized;
   Scalar e0l, e1l, e2l;

   Vec3 translation;
   Matrix rotation;

   math::GenericAABB< Scalar > aabb;
};


/**
* \brief Adds information required to compute signed distances from a point to a triangle
* 
* The class adds vertex, edge and face normals to mesh. The normals are calculated according to \cite Baerentzen2005
* to allow for a numerically stable sign computation. To compute the point<->triangle distance, we use the rotation
* method described in \cite Jones1995. Some data is precomputed to allow for faster computations and is stored
* at each face in an object of class DistanceProperties. The numbering of vertices, edges and voronoi regions in
* the rotated triangles are shown here:
*
* \image html mesh/triangle_topo.svg "A rotated triangle in the planar coordinate system. The vertex numbering is shown in red, the edge numbering in blue and the numbering of the voronoi regions in green."
*
* The class offers multiple methods to get the signed squared distance from a point to a single triangle or the whole mesh.
* Please note that the distance computation for whole meshes is rather inefficient. Instead you should an object of this
* class into a \ref mesh::distance_octree::DistanceOctree "mesh::DistanceOctree" reduce the computational complexity
* from \f$\mathcal O(n)\f$ to \f$\mathcal O(\log n)\f$, where \f$\mathcal O(n)\f$ where \f$n\f$ is the number triangles.
*
* Additionally to the signed squared distance you can also retrieve the closest point on the triangle or mesh to
* your point of inquiry. You may also retrieve the corresponding normal and the closest voronoi region.
*/
template< typename MeshType >
class TriangleDistance
{
public:
   typedef typename MeshType::Scalar     Scalar;
   typedef typename MeshType::Point      Point;
   typedef typename MeshType::Normal     Normal;
   typedef typename MeshType::FaceHandle FaceHandle;
   typedef math::Vector3<Scalar>         Vec3;
   typedef math::GenericAABB<Scalar>     BoundingBox;

   TriangleDistance( const shared_ptr<MeshType> & mesh ) : mesh_(mesh), distanceProperties_( *mesh, "DistanceProperties" ) { computeNormals(); computeDistanceProperties(); }

   Scalar sqDistance      ( const FaceHandle fh, const Point & p ) const;
   Scalar sqSignedDistance( const FaceHandle fh, const Point & p ) const;
   Scalar distance        ( const FaceHandle fh, const Point & p ) const;
   Scalar signedDistance  ( const FaceHandle fh, const Point & p ) const;

   Scalar sqDistance      ( const FaceHandle fh, const Point & p, Point & closestPoint ) const;
   Scalar sqSignedDistance( const FaceHandle fh, const Point & p, Point & closestPoint ) const;
   Scalar distance        ( const FaceHandle fh, const Point & p, Point & closestPoint ) const;
   Scalar signedDistance  ( const FaceHandle fh, const Point & p, Point & closestPoint ) const;

   Scalar sqDistance      ( const FaceHandle fh, const Point & p, Point & closestPoint, Normal & normal ) const;
   Scalar sqSignedDistance( const FaceHandle fh, const Point & p, Point & closestPoint, Normal & normal ) const;
   Scalar distance        ( const FaceHandle fh, const Point & p, Point & closestPoint, Normal & normal ) const;
   Scalar signedDistance  ( const FaceHandle fh, const Point & p, Point & closestPoint, Normal & normal ) const;

   Scalar sqDistance      ( const Point & p ) const;
   Scalar sqSignedDistance( const Point & p ) const;
   Scalar distance        ( const Point & p ) const;
   Scalar signedDistance  ( const Point & p ) const;

   Scalar sqDistance      ( const Point & p, FaceHandle & closestTriangle ) const;
   Scalar sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const;
   Scalar distance        ( const Point & p, FaceHandle & closestTriangle ) const;
   Scalar signedDistance  ( const Point & p, FaceHandle & closestTriangle ) const;

   Scalar sqDistance      ( const Point & p, Point & closestPoint ) const;
   Scalar sqSignedDistance( const Point & p, Point & closestPoint ) const;
   Scalar distance        ( const Point & p, Point & closestPoint ) const;
   Scalar signedDistance  ( const Point & p, Point & closestPoint ) const;

   Scalar sqDistance      ( const Point & p, Point & closestPoint, Normal & normal ) const;
   Scalar sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const;
   Scalar distance        ( const Point & p, Point & closestPoint, Normal & normal ) const;
   Scalar signedDistance  ( const Point & p, Point & closestPoint, Normal & normal ) const;

   template< typename InputIterator >
   Scalar sqDistance      ( InputIterator fhBegin, InputIterator fhEnd, const Point & p ) const;
   template< typename InputIterator >
   Scalar sqSignedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p ) const;
   template< typename InputIterator >
   Scalar distance        ( InputIterator fhBegin, InputIterator fhEnd, const Point & p ) const;
   template< typename InputIterator >
   Scalar signedDistance  ( InputIterator fhBegin, InputIterator fhEnd, const Point & p ) const;

   template< typename InputIterator >
   Scalar sqDistance      ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, FaceHandle & closestTriangle ) const;
   template< typename InputIterator >
   Scalar sqSignedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, FaceHandle & closestTriangle ) const;
   template< typename InputIterator >
   Scalar distance        ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, FaceHandle & closestTriangle ) const;
   template< typename InputIterator >
   Scalar signedDistance  ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, FaceHandle & closestTriangle ) const;

   template< typename InputIterator >
   Scalar sqDistance      ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint ) const;
   template< typename InputIterator >
   Scalar sqSignedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint ) const;
   template< typename InputIterator >
   Scalar distance        ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint ) const;
   template< typename InputIterator >
   Scalar signedDistance  ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint ) const;

   template< typename InputIterator >
   Scalar sqDistance      ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint, Normal & normal ) const;
   template< typename InputIterator >
   Scalar sqSignedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint, Normal & normal ) const;
   template< typename InputIterator >
   Scalar distance        ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint, Normal & normal ) const;
   template< typename InputIterator >
   Scalar signedDistance  ( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint, Normal & normal ) const;

   template< typename InputIterator, typename OutputIterator >
   void filterTrianglesForAABB( const BoundingBox & aabb, const InputIterator fhInBegin, const InputIterator fhInEnd, OutputIterator fhOutBegin ) const;

   const MeshType & getMesh() const { return *mesh_; }
   shared_ptr<MeshType> getMeshPtr() { return mesh_; }

   const BoundingBox & getAabb( FaceHandle fh ) const { return distanceProperties_[ fh ].aabb; }

   void triangleToStream( const FaceHandle fh, std::ostream & os ) const;

protected:
   typedef typename OpenMesh::FPropHandleT< DistanceProperties<MeshType> > DistancePropertyHandle;

   void computeNormals();
   void computeDistanceProperties();

   Scalar sqDistance( const FaceHandle fh, const Point & p, int & region, Point & closestPoint ) const;
   template< typename InputIterator >
   Scalar sqDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, int & region, Point & closestPoint, FaceHandle & fh ) const;

   inline Scalar toSignedDistance( const Scalar & d, const FaceHandle fh, const Point & p, const int region, Point & closestPoint ) const;
   inline Scalar toSignedDistance( const Scalar & d, const Point & p, Point & closestPoint, const Normal & normal ) const;

   inline const Normal & getNormal( const FaceHandle fh, const int region ) const;

   shared_ptr<MeshType> mesh_;
   OpenMesh::PropertyManager< DistancePropertyHandle, MeshType > distanceProperties_;
};



template< typename MeshType >
void TriangleDistance<MeshType>::computeNormals()
{
   // J. B�rentzen and H. Aan�s. Signed distance computation using the angle weighted pseudonormal.
   // Visualization and Computer Graphics, IEEE Transactions on, 11(3):243�253, 2005.

   static_assert( MeshType::IsTriMesh == 1, "computeNormals only works with triangular meshes!" );

   mesh_->request_face_normals();
   mesh_->request_halfedge_normals();
   mesh_->request_vertex_normals();

   WALBERLA_ASSERT( mesh_->has_vertex_normals() );
   WALBERLA_ASSERT( mesh_->has_halfedge_normals() );
   WALBERLA_ASSERT( mesh_->has_face_normals() );

   mesh_->update_face_normals();

   bool watertightWarningLogged = false;

   // Update half edge normals
   for(auto he_it = mesh_->halfedges_begin(); he_it != mesh_->halfedges_end(); ++he_it)
   {
      auto face0 = mesh_->face_handle( *he_it );
      auto face1 = mesh_->face_handle( mesh_->opposite_halfedge_handle( *he_it ) );

      Normal n;
      if( !face0.is_valid() || !face1.is_valid() )
      {
         if(!watertightWarningLogged)
         {
            WALBERLA_LOG_WARNING( "Mesh is not watertight and has boundary edges! Make sure the boundary edges are located outside of the domain or at least are located at its edge" );
            watertightWarningLogged = true;
         }

         if( face0.is_valid() )
         {
            WALBERLA_ASSERT_FLOAT_EQUAL( mesh_->normal( face0 ).sqrnorm(), typename MeshType::Scalar(1) );
            n  = mesh_->normal( face0 );
         }
         else if( face1.is_valid() )
         {
            WALBERLA_ASSERT_FLOAT_EQUAL( mesh_->normal( face1 ).sqrnorm(), typename MeshType::Scalar(1) );
            n  = mesh_->normal( face1 );
         }
         else
         {
            WALBERLA_ABORT("Found loose halfedge!");
         }
      }
      else 
      {
         WALBERLA_ASSERT_FLOAT_EQUAL( mesh_->normal( face0 ).sqrnorm(), typename MeshType::Scalar(1) );
         WALBERLA_ASSERT_FLOAT_EQUAL( mesh_->normal( face1 ).sqrnorm(), typename MeshType::Scalar(1) );

         n  = ( mesh_->normal( face0 ) + mesh_->normal( face1 ) ).normalized();
         
      }
      mesh_->set_normal( *he_it, n );
      mesh_->set_normal( mesh_->opposite_halfedge_handle( *he_it ), n );
   }

   // Update vertex normals
   for(auto v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
   {
      if( mesh_->vih_cwbegin( *v_it ) == mesh_->vih_cwend( *v_it ) ) // unreferenced vertex which we can ignore....
         continue;

      typename MeshType::Normal normal( typename MeshType::Scalar(0) );
      auto vih_it     = mesh_->vih_ccwbegin( *v_it );
      auto prev_edge = mesh_->point( mesh_->from_vertex_handle( *vih_it ) ) - mesh_->point( *v_it );
      ++vih_it;
      while( vih_it != mesh_->vih_ccwend( *v_it ) )
      {
         auto edge = mesh_->point( mesh_->from_vertex_handle( *vih_it ) ) - mesh_->point( *v_it );
         auto angle = std::acos( dot(edge, prev_edge) / ( edge.length() * prev_edge.length() ) );
         normal += angle * mesh_->normal( mesh_->face_handle( *vih_it ) );
         ++vih_it;
         prev_edge = edge;
      }
      vih_it = mesh_->vih_ccwbegin( *v_it );
      auto edge = mesh_->point( mesh_->from_vertex_handle( *vih_it ) ) - mesh_->point( *v_it );
      auto angle = std::acos( dot(edge, prev_edge) / ( edge.length() * prev_edge.length() ) );
      normal += angle * mesh_->normal( mesh_->face_handle( *vih_it ) );
      normal.normalize();

      mesh_->set_normal( *v_it, normal );
   }
}


template< typename MeshType >
void TriangleDistance<MeshType>::computeDistanceProperties()
{
   typedef DistanceProperties<MeshType> DP;

   WALBERLA_ASSERT( mesh_->has_vertex_normals() );
   WALBERLA_ASSERT( mesh_->has_halfedge_normals() );
   WALBERLA_ASSERT( mesh_->has_face_normals() );

   //Precompute properties for distance computation

   for(auto f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
   {
      DP & dp = distanceProperties_[ *f_it ];

      const typename DP::Vec3 & v0 = mesh_->point( getVertexHandle( *mesh_, *f_it, 0U ) );
      const typename DP::Vec3 & v1 = mesh_->point( getVertexHandle( *mesh_, *f_it, 1U ) );
      const typename DP::Vec3 & v2 = mesh_->point( getVertexHandle( *mesh_, *f_it, 2U ) );

      dp.translation = -v0;

      typename DP::Vec3 e0 = v1 - v0;
      typename DP::Vec3 e1 = v2 - v0;

      typename DP::Vec3 newYAxis = e0;
      newYAxis.normalize();
      typename DP::Vec3 newZAxis = e0 % e1;
      newZAxis.normalize();
      typename DP::Vec3 newXAxis = newYAxis % newZAxis;
      newXAxis.normalize();

      if( ( newXAxis | e1 ) < typename DP::Scalar(0) )
      {
         newZAxis = -newZAxis;
         newXAxis = -newXAxis;
      }

      WALBERLA_ASSERT_FLOAT_EQUAL( newXAxis | newYAxis, typename DP::Scalar(0) );
      WALBERLA_ASSERT_FLOAT_EQUAL( newXAxis | newZAxis, typename DP::Scalar(0) );
      WALBERLA_ASSERT_FLOAT_EQUAL( newYAxis | newZAxis, typename DP::Scalar(0) );

      for( uint_t i = 0; i < 3; ++i )
      {
         dp.rotation(0, i) = newXAxis[i];
         dp.rotation(1, i) = newYAxis[i];
         dp.rotation(2, i) = newZAxis[i];
      }

      WALBERLA_ASSERT_FLOAT_EQUAL( ( dp.rotation * dp.rotation.getTranspose() ).trace(), typename DP::Scalar(3) ); // Matrix is orthogonal
      WALBERLA_ASSERT_FLOAT_EQUAL( dp.rotation.getDeterminant(), typename DP::Scalar(1) );

      e0 = dp.rotation * e0;
      e1 = dp.rotation * e1;

      typename DP::Vec3 e2 = e1 - e0;

      dp.e0 = typename DP::Vec2( e0[0], e0[1] );
      dp.e1 = typename DP::Vec2( e1[0], e1[1] );
      dp.e2 = typename DP::Vec2( e2[0], e2[1] );

      dp.e0l = dp.e0.length();
      dp.e1l = dp.e1.length();
      dp.e2l = dp.e2.length();

      dp.e0_normalized = dp.e0 / dp.e0l;
      dp.e1_normalized = dp.e1 / dp.e1l;
      dp.e2_normalized = dp.e2 / dp.e2l;

      dp.e1_normal = typename DP::Vec2(  dp.e1[1], -dp.e1[0] );
      dp.e2_normal = typename DP::Vec2( -dp.e2[1],  dp.e2[0] );

      dp.e1_normal.normalize();
      dp.e2_normal.normalize();

      dp.aabb.init( toWalberla( v0 ), toWalberla( v1 ) );
      dp.aabb.merge( toWalberla( v2 ) );
   }
}



template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( const FaceHandle fh, const Point & p, int & region, Point & closestPoint ) const
{
   // Jones, Mark W. "3D distance from a point to a triangle."
   // Department of Computer Science, University of Wales Swansea Technical Report CSR-5 (1995).

   WALBERLA_ASSERT( mesh_->is_valid_handle( fh ), "Given face handle is invalid!" );

   typedef DistanceProperties<MeshType> DP;

   const DP & dp = distanceProperties_[ fh ];

   typename DP::Vec3 pt3 = dp.rotation * typename DP::Vec3( p + dp.translation );

   closestPoint[2] = 0.0;
   typename DP::Scalar result = pt3[2] * pt3[2];

   typename DP::Vec2 pt( pt3[0], pt3[1] );

   typename DP::Scalar e0p = dp.e0_normalized | pt;
   typename DP::Scalar e1p = dp.e1_normalized | pt;
   typename DP::Scalar e2p = dp.e2_normalized | ( pt - dp.e0 );

   typename DP::Scalar e0d = typename DP::Vec2( -1, 0 ) | pt;
   typename DP::Scalar e1d = dp.e1_normal | pt;
   typename DP::Scalar e2d = dp.e2_normal | ( pt - dp.e0 );

   if( e0p <= 0 && e1p <= 0  )
   {
      // Voronoi area of vertex 0
      region = 1;
      result += pt.sqrnorm(); // distance from v0
      closestPoint[0] = closestPoint[1] = 0;
   }
   else if( e0p >= dp.e0l && e2p <= 0 )
   {
      // Voronoi area of vertex 1
      region = 2;
      result += (pt - dp.e0).sqrnorm(); // distance from v1
      closestPoint[0] = dp.e0[0];
      closestPoint[1] = dp.e0[1];
   }
   else if( e1p >= dp.e1l && e2p >= dp.e2l )
   {
      // Voronoi area of vertex 2
      region = 3;
      result += (pt - dp.e1).sqrnorm(); // distance from v2
      closestPoint[0] = dp.e1[0];
      closestPoint[1] = dp.e1[1];
   }
   else if( e0d <= 0 && e1d <= 0 && e2d <= 0 )
   {
      // Voronoi area of face
      region = 0;
      // result += 0;
      closestPoint[0] = pt[0];
      closestPoint[1] = pt[1];
   }
   else if( e0d >= 0 && e0p > 0 && e0p < dp.e0l )
   {
      // Voronoi area of edge 0
      region = 4;
      result += pt[0] * pt[0];
      closestPoint[0] = 0;
      closestPoint[1] = pt[1];
   }
   else if( e1d >= 0 && e1p > 0 && e1p < dp.e1l )
   {
      // Voronoi area of edge 1
      region = 5;
      result += e1d * e1d;
      closestPoint[0] = e1p * dp.e1_normalized[0];
      closestPoint[1] = e1p * dp.e1_normalized[1];
   }
   else if( e2d >= 0 && e2p > 0 && e2p < dp.e2l )
   {
      // Voronoi area of edge 2
      region = 6;
      result += e2d * e2d;
      closestPoint[0] = dp.e0[0] + e2p * dp.e2_normalized[0];
      closestPoint[1] = dp.e0[1] + e2p * dp.e2_normalized[1];
   }
   else
   {
      region = -1; // Silence compiler warning
      WALBERLA_ASSERT(false);
   }

   closestPoint = dp.rotation.getTranspose() * closestPoint - dp.translation;

   return result;
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::toSignedDistance( const Scalar & d, const FaceHandle fh, const Point & p, const int region, Point & closestPoint ) const
{
   return toSignedDistance( d, p, closestPoint, getNormal( fh, region ) );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::toSignedDistance( const Scalar & d, const Point & p, Point & closestPoint, const Normal & normal ) const
{
   WALBERLA_ASSERT_GREATER_EQUAL( d, Scalar(0) );

   Scalar dot = ( p - closestPoint ) | normal;

   return dot >= Scalar(0) ? d : -d;
}


template< typename MeshType >
void TriangleDistance<MeshType>::triangleToStream( const FaceHandle fh, std::ostream & os ) const
{
   os << "V0: " << mesh_->point(  getVertexHandle( *mesh_, fh, 0U ) ) << "\n"
      << "V1: " << mesh_->point(  getVertexHandle( *mesh_, fh, 1U ) ) << "\n"
      << "V2: " << mesh_->point(  getVertexHandle( *mesh_, fh, 2U ) ) << "\n"
      << "\n"
      << "FN:  " << mesh_->normal( fh ) << "\n"
      << "V0N: " << mesh_->normal( getVertexHandle( *mesh_, fh, 0U ) ) << "\n"
      << "V1N: " << mesh_->normal( getVertexHandle( *mesh_, fh, 1U ) ) << "\n"
      << "V2N: " << mesh_->normal( getVertexHandle( *mesh_, fh, 2U ) ) << "\n"
      << "E0N: " << mesh_->normal( getHalfedgeHandle( *mesh_, fh, 0U, 1U ) ) << "\n"
      << "E1N: " << mesh_->normal( getHalfedgeHandle( *mesh_, fh, 0U, 2U ) ) << "\n"
      << "E2N: " << mesh_->normal( getHalfedgeHandle( *mesh_, fh, 1U, 2U ) );
}

template< typename MeshType >
const typename MeshType::Normal & TriangleDistance<MeshType>::getNormal( const FaceHandle fh, const int region ) const
{
   switch(region)
   {
   case 1: return mesh_->normal( getVertexHandle( *mesh_, fh, 0U ) );
   case 2: return mesh_->normal( getVertexHandle( *mesh_, fh, 1U ) );
   case 3: return mesh_->normal( getVertexHandle( *mesh_, fh, 2U ) );
   case 4: return mesh_->normal( getHalfedgeHandle( *mesh_, fh, 0U, 1U ) );
   case 5: return mesh_->normal( getHalfedgeHandle( *mesh_, fh, 0U, 2U ) );
   case 6: return mesh_->normal( getHalfedgeHandle( *mesh_, fh, 1U, 2U ) );
   default:
   case 0: return mesh_->normal( fh );
   }
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( const FaceHandle fh, const Point & p ) const
{
   int region;
   Point closestPoint;
   return sqDistance( fh, p, region, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( const FaceHandle fh, const Point & p ) const
{
   int region;
   Point closestPoint;
   Scalar d = sqDistance( fh, p, region, closestPoint );
   return toSignedDistance( d, fh, p, region, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( const FaceHandle fh, const Point & p ) const
{
   int region;
   Point closestPoint;
   return std::sqrt( sqDistance( fh, p, region, closestPoint ) );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( const FaceHandle fh, const Point & p ) const
{
   int region;
   Point closestPoint;
   Scalar d = std::sqrt( sqDistance( fh, p, region, closestPoint ) );
   return toSignedDistance( d, fh, p, region, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( const FaceHandle fh, const Point & p, Point & closestPoint ) const
{
   int region;
   return sqDistance( fh, p, region, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( const FaceHandle fh, const Point & p, Point & closestPoint ) const
{
   int region;
   Scalar d = sqDistance( fh, p, region, closestPoint );
   return toSignedDistance( d, fh, p, region, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( const FaceHandle fh, const Point & p, Point & closestPoint ) const
{
   int region;
   return std::sqrt( sqDistance( fh, p, region, closestPoint ) );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( const FaceHandle fh, const Point & p, Point & closestPoint ) const
{
   int region;
   Scalar d = std::sqrt( sqDistance( fh, p, region, closestPoint ) );
   return toSignedDistance( d, fh, p, region, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( const FaceHandle fh, const Point & p, Point & closestPoint, Normal & normal ) const
{
   int region;
   auto d = sqDistance( fh, p, region, closestPoint );
   normal = getNormal( fh, region );
   return d;
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( const FaceHandle fh, const Point & p, Point & closestPoint, Normal & normal ) const
{
   int region;
   Scalar d = sqDistance( fh, p, region, closestPoint );
   normal = getNormal( fh, region );
   return toSignedDistance( d, p, closestPoint, normal );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( const FaceHandle fh, const Point & p, Point & closestPoint, Normal & normal ) const
{
   int region;
   auto d = sqDistance( fh, p, region, closestPoint );
   normal = getNormal( fh, region );
   return std::sqrt( d );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( const FaceHandle fh, const Point & p, Point & closestPoint, Normal & normal ) const
{
   int region;
   Scalar d = std::sqrt( sqDistance( fh, p, region, closestPoint ) );
   normal = getNormal( fh, region );
   return toSignedDistance( d, p, closestPoint, normal );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, int & region, Point & closestPoint, FaceHandle & fh ) const
{
   WALBERLA_ASSERT_GREATER( std::distance(fhBegin, fhEnd), 0, "Empty face list!" );
   Scalar minDistance = std::numeric_limits<Scalar>::max();

   int tmpRegion;
   Point tmpClosestPoint;
   for( auto fIt = fhBegin; fIt != fhEnd; ++fIt )
   {
      Scalar d = sqDistance( *fIt, p, tmpRegion, tmpClosestPoint );
      if( d < minDistance )
      {
         minDistance = d;
         region = tmpRegion;
         closestPoint = tmpClosestPoint;
         fh = *fIt;
      }
   }

   return minDistance;
}

template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( const Point & p ) const
{
   return sqDistance( mesh_->faces_begin(), mesh_->faces_end(), p );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( const Point & p ) const
{
   return sqSignedDistance( mesh_->faces_begin(), mesh_->faces_end(), p );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( const Point & p ) const
{
   return distance( mesh_->faces_begin(), mesh_->faces_end(), p );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( const Point & p ) const
{
   return signedDistance( mesh_->faces_begin(), mesh_->faces_end(), p );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( const Point & p, FaceHandle & closestTriangle ) const
{
   return sqDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestTriangle );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const
{
   return sqSignedDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestTriangle );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( const Point & p, FaceHandle & closestTriangle ) const
{
   return distance( mesh_->faces_begin(), mesh_->faces_end(), p, closestTriangle );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( const Point & p, FaceHandle & closestTriangle ) const
{
   return signedDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestTriangle );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( const Point & p, Point & closestPoint ) const
{
   return sqDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( const Point & p, Point & closestPoint ) const
{
   return sqSignedDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( const Point & p, Point & closestPoint ) const
{
   return distance( mesh_->faces_begin(), mesh_->faces_end(), p, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( const Point & p, Point & closestPoint ) const
{
   return signedDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestPoint );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( const Point & p, Point & closestPoint, Normal & normal ) const
{
   return sqDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestPoint, normal );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const
{
   return sqSignedDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestPoint, normal );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( const Point & p, Point & closestPoint, Normal & normal ) const
{
   return distance( mesh_->faces_begin(), mesh_->faces_end(), p, closestPoint, normal );
}


template< typename MeshType >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( const Point & p, Point & closestPoint, Normal & normal ) const
{
   return signedDistance( mesh_->faces_begin(), mesh_->faces_end(), p, closestPoint, normal );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p ) const
{
   int region = 0;
   Point closestPoint;
   FaceHandle fh;
   return sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p ) const
{
   int region = 0;
   Point closestPoint;
   FaceHandle fh;
   Scalar d = sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh );
   return toSignedDistance( d, fh, p, region, closestPoint );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( InputIterator fhBegin, InputIterator fhEnd, const Point & p ) const
{
   int region = 0;
   Point closestPoint;
   FaceHandle fh;
   return std::sqrt( sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh ) );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p ) const
{
   int region = 0;
   Point closestPoint;
   FaceHandle fh;
   Scalar d = std::sqrt( sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh ) );
   return toSignedDistance( d, fh, p, region, closestPoint );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, FaceHandle & closestTriangle ) const
{
   int region = 0;
   Point closestPoint;
   return sqDistance( fhBegin, fhEnd, p, region, closestPoint, closestTriangle );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, FaceHandle & closestTriangle ) const
{
   int region = 0;
   Point closestPoint;
   Scalar d = sqDistance( fhBegin, fhEnd, p, region, closestPoint, closestTriangle );
   return toSignedDistance( d, closestTriangle, p, region, closestPoint );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, FaceHandle & closestTriangle ) const
{
   int region = 0;
   Point closestPoint;
   return std::sqrt( sqDistance( fhBegin, fhEnd, p, region, closestPoint, closestTriangle ) );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, FaceHandle & closestTriangle ) const
{
   int region = 0;
   Point closestPoint;
   Scalar d = std::sqrt( sqDistance( fhBegin, fhEnd, p, region, closestPoint, closestTriangle ) );
   return toSignedDistance( d, closestTriangle, p, region, closestPoint );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint ) const
{
   int region = 0;
   FaceHandle fh;
   return sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint ) const
{
   int region = 0;
   FaceHandle fh;
   Scalar d = sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh );
   return toSignedDistance( d, fh, p, region, closestPoint );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint ) const
{
   int region = 0;
   FaceHandle fh;
   return std::sqrt( sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh ) );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint ) const
{
   int region = 0;
   FaceHandle fh;
   Scalar d = std::sqrt( sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh ) );
   return toSignedDistance( d, fh, p, region, closestPoint );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint, Normal & normal ) const
{
   int region = 0;
   FaceHandle fh;
   auto d = sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh );
   normal = getNormal( fh, region );
   return d;
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::sqSignedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint, Normal & normal ) const
{
   int region = 0;
   FaceHandle fh;
   Scalar d = sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh );
   normal = getNormal( fh, region );
   return toSignedDistance( d, p, closestPoint, normal );
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::distance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint, Normal & normal ) const
{
   int region = 0;
   FaceHandle fh;
   auto d = std::sqrt( sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh ) );
   normal = getNormal( fh, region );
   return d;
}


template< typename MeshType >
template< typename InputIterator >
typename MeshType::Scalar TriangleDistance<MeshType>::signedDistance( InputIterator fhBegin, InputIterator fhEnd, const Point & p, Point & closestPoint, Normal & normal ) const
{
   int region = 0;
   FaceHandle fh;
   Scalar d = std::sqrt( sqDistance( fhBegin, fhEnd, p, region, closestPoint, fh ) );
   normal = getNormal( fh, region );
   return toSignedDistance( d, p, closestPoint, normal );
}


template< typename MeshType >
template< typename InputIterator, typename OutputIterator >
void TriangleDistance<MeshType>::filterTrianglesForAABB( const BoundingBox & aabb, const InputIterator fhInBegin, const InputIterator fhInEnd, OutputIterator fhOut ) const
{
   struct TringleDistance
   {
      TringleDistance( FaceHandle _fh, const BoundingBox & _aabb, const BoundingBox & triAabb )
         : fh( _fh )
      {
         minSqDistance = _aabb.sqDistance( triAabb );
         maxSqDistance = _aabb.sqMaxDistance( triAabb );
      }

      bool operator<( const TringleDistance & other ) const { return this->minSqDistance < other.minSqDistance; }

      FaceHandle fh;
      Scalar     minSqDistance;
      Scalar     maxSqDistance;
   };

   if( fhInBegin == fhInEnd )
      return;

   std::vector<TringleDistance> td;
   for( auto it = fhInBegin; it != fhInEnd; ++it )
      td.push_back( TringleDistance( *it, aabb, getAabb( *it ) ) );

   std::sort( td.begin(), td.end() );

   Scalar sqMaxDist = std::numeric_limits<Scalar>::max();

   for( auto it = td.begin(); it != td.end(); ++it )
      sqMaxDist = std::min( sqMaxDist, it->maxSqDistance );


   for( auto it = td.begin(); it != td.end(); ++it )
      if( it->minSqDistance <= sqMaxDist )
      {
         *fhOut++ = it->fh;
      }
      else
         break;

}

} // namespace mesh
} // namespace walberla
