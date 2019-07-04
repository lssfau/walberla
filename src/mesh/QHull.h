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
//! \file QHull.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include "mesh/TriangleMeshes.h"
#include "mesh/MatrixVectorOperations.h"
#include "mesh/vtk/VTKMeshWriter.h"
#include "mesh/vtk/CommonDataSources.h"

#include "vtk/VTKOutput.h"
#include "vtk/PointDataSource.h"

#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include <algorithm>
#include <iterator>
#include <vector>
#include <queue>

namespace walberla {
namespace mesh {

template< typename MeshType >
inline shared_ptr<MeshType> convexHull( const Vector3<real_t> & points );

template< typename MeshType >
inline shared_ptr<MeshType> convexHull( const MeshType & mesh );


// Forward declarations
template< typename MeshType >
class QHullFaceSorter;

template< typename MeshType >
class QHullPointDataSource;
  
template< typename MeshType >
class QHull
{
public:
   static_assert( MeshType::IsTriMesh == 1, "QHull only works on triangular meshes!" );

   typedef typename MeshType::Point Point;
   typedef typename MeshType::Scalar Scalar;
   typedef typename MeshType::FaceHandle FaceHandle;
   typedef typename MeshType::VertexHandle VertexHandle;
   typedef typename MeshType::HalfedgeHandle HalfedgeHandle;

   QHull( const std::vector< Vector3<real_t> > & pointCloud, const shared_ptr<MeshType> & mesh = make_shared< MeshType >());
   QHull( const std::vector< Point > & pointCloud, const shared_ptr<MeshType> & mesh = make_shared< MeshType >() );

   uint_t run();

   inline const shared_ptr< MeshType > & meshPtr() const { return mesh_; }
   inline const MeshType & mesh() const { return *mesh_; }
   inline const std::vector<Point> & getVisiblePoints( const FaceHandle fh ) const { return visiblePoints_[fh]; }

   void enableDebugVTKOutput( const std::string & identifierPrefix = "QHull", const std::string & baseFolder = "vtk_out" );
   void disableDebugVTKOutput();
   
private:
   inline void initMesh();

   inline bool pointIsVisible( const FaceHandle fh, const Point & p ) const
   { 
      return ( ( p - mesh_->point( *( mesh_->cfv_begin(fh) ) ) ) | mesh_->normal(fh) ) > Scalar(0);
   }

   inline Scalar pointDistance( const FaceHandle fh, const Point & p ) const
   { 
      return std::fabs( ( p - mesh_->point( *( mesh_->cfv_begin(fh) ) ) ) | mesh_->normal(fh) );
   }

   inline void createInitialSimplex();
   
   inline void iteration();
   
   inline void deleteVisibleFaces( const FaceHandle startFaceHandle, const Point & p );
   
   inline void addNewFaces( const Point & p );

   template< typename InputIterator >
   inline void assignPointsToFaces( const std::vector<Point> & points, InputIterator facesBegin, InputIterator facesEnd );


   
   std::vector< Point > pointCloud_; /// The initial point cloud
   shared_ptr< MeshType > mesh_; /// The genereated convex mesh
   typedef typename OpenMesh::FPropHandleT< std::vector<Point> > VisiblePointsPropertyHandle;
   OpenMesh::PropertyManager< VisiblePointsPropertyHandle, MeshType > visiblePoints_; /// Property storing the points of a certain face
   std::priority_queue<FaceHandle, std::vector<FaceHandle>, QHullFaceSorter<MeshType> > queue_; /// queue to proptize faces

   // Vectors to be reused in between iterations
   std::vector<Point> orphanPoints_; /// Points getting orphaned during face removal
   std::vector<HalfedgeHandle> horizon_; /// the new horizon
   std::vector<FaceHandle> newFaces_; /// new faces created in an iteration replacing the removed ones

   // Debug VTK Output
   bool writeDebugOutput_; /// Should debug ouput be written?
   shared_ptr<VTKMeshWriter< MeshType >> meshWriter_; /// Writes the current state of the mesh
   shared_ptr<vtk::VTKOutput> remainingPointsWriter_; /// Writes the remaining points belonging to the faces
};


template< typename MeshType >
class QHullFaceSorter
{
public:
   QHullFaceSorter( const QHull<MeshType> & qhull ) : qhull_(qhull) {}
   bool operator()( const typename QHull<MeshType>::FaceHandle lhs, const typename QHull<MeshType>::FaceHandle rhs )
   {
      return qhull_.getVisiblePoints(lhs).size() < qhull_.getVisiblePoints(rhs).size();
   }

private:
   const QHull<MeshType> & qhull_;
};


template< typename MeshType >
shared_ptr<MeshType> convexHull( const Vector3<real_t> & points )
{
   QHull<MeshType> qhull( points );
   qhull.run();
   return qhull.meshPtr();
}

template< typename MeshType >
inline shared_ptr<MeshType> convexHull( const MeshType & mesh )
{
   QHull<MeshType> qhull( std::vector<typename MeshType::Point>( mesh.vertices_begin(), mesh.vertices_.end() ) );
   qhull.run();
   return qhull.meshPtr();
}

template< typename MeshType >
class QHullPointDataSource : public vtk::PointDataSource
{
public:
   typedef typename MeshType::Point Point;
   typedef typename MeshType::Scalar Scalar;

   QHullPointDataSource( const QHull<MeshType> & qhull ) : qhull_( qhull ) {}

   virtual std::vector< Attributes > getAttributes() const 
   {
      std::vector< Attributes > attributes; 
      attributes.push_back( Attributes("Int32", "index", uint_t(1)) );
      return attributes;
   }

   virtual std::vector< Vector3< real_t > > getPoints() { return points_; }
   virtual void configure()
   {
      points_.clear();
      indices_.clear();
      for(auto fhIt = qhull_.mesh().faces_sbegin(); fhIt != qhull_.mesh().faces_end(); ++fhIt )
      {
         for(const auto & p : qhull_.getVisiblePoints( *fhIt ))
            points_.push_back( Vector3<real_t>( real_c(p[0]), real_c(p[1]), real_c(p[2]) ) );

         indices_.insert( indices_.end(), qhull_.getVisiblePoints(*fhIt).size(), fhIt->idx() );
      }
   };

   virtual void push( std::ostream& os,  const uint_t /*data*/, const uint_t point, const uint_t /*component*/ )
   {
      os << indices_[point];
   };

   virtual void push( vtk::Base64Writer& b64, const uint_t /*data*/, const uint_t point, const uint_t /*component*/ )
   {
      b64 << indices_[point];
   };

private:
   const QHull<MeshType> & qhull_;
   std::vector<Vector3<real_t>> points_;
   std::vector<int> indices_;
};


template< typename MeshType >
QHull<MeshType>::QHull( const std::vector< Vector3<real_t> > & pointCloud, const shared_ptr<MeshType> & _mesh ) 
   : pointCloud_(pointCloud.size()) , mesh_(_mesh), visiblePoints_( *mesh_, "VisiblePoints"),
   queue_(*this), writeDebugOutput_(false)
{
   std::transform( pointCloud.begin(), pointCloud.end(), pointCloud_.begin(), toOpenMesh<real_t> );
   initMesh();
}


template< typename MeshType >
QHull<MeshType>::QHull( const std::vector< Point > & pointCloud, const shared_ptr<MeshType> & _mesh ) 
   : pointCloud_(pointCloud) , mesh_( _mesh ), visiblePoints_( *mesh_, "VisiblePoints"), 
   queue_(*this), writeDebugOutput_(false)
{
   initMesh();
}


template< typename MeshType >
uint_t QHull<MeshType>::run()
{
   createInitialSimplex();

   uint_t iterations(0);

   while( !queue_.empty() )
   {
      iteration();
      ++iterations;
   }

   mesh_->garbage_collection();

   return iterations;
}

template< typename MeshType >
void QHull<MeshType>::iteration()
{
   FaceHandle fh = queue_.top();
   queue_.pop();

   WALBERLA_ASSERT( !visiblePoints_[fh].empty() );

   // Find the point farthest away from the current face
   auto it = std::max_element( visiblePoints_[fh].begin(), visiblePoints_[fh].end(), [&](const Point & lhs, const Point & rhs){
      return pointDistance(fh, lhs) < pointDistance(fh, rhs);
   });

   // Remove the point from the remaining points
   Point p = *it;
   visiblePoints_[fh].erase(it);

   // Delete all faces visible by p
   deleteVisibleFaces( fh, p );

   // Fill the hole by connecting the horizon edges to p
   addNewFaces( p );

   // reassign the orphaned points from the deleted faces to the newly created faces
   assignPointsToFaces( orphanPoints_, newFaces_.begin(), newFaces_.end() );

   // Add new faces to the queue
   for( auto nfh : newFaces_)
   {
      if( !visiblePoints_[nfh].empty() )
         queue_.push(nfh);
   }

   // Remove deleted faces from the top of the queue
   while(!queue_.empty() && mesh_->status( queue_.top() ).deleted())
   {
      queue_.pop();
   }

   if(writeDebugOutput_)
   {
      (*meshWriter_)();
      remainingPointsWriter_->write();
   }
}

template< typename MeshType >
void QHull<MeshType>::createInitialSimplex()
{
   // Find extreme points on coordinate axis
   std::vector< Point > ep( size_t(6),  pointCloud_.front() );
   for(const auto & p : pointCloud_)
   {
      if(p[0] < ep[0][0])
         ep[0] = p;
      if(p[1] < ep[1][1])
         ep[1] = p;
      if(p[2] < ep[2][2])
         ep[2] = p;
      if(p[0] > ep[3][0])
         ep[3] = p;
      if(p[1] > ep[4][1])
         ep[4] = p;
      if(p[2] > ep[5][2])
         ep[5] = p;
   }

   // Find the most distant pair of points out of ep
   Point base0, base1;
   Scalar maxSqDist = Scalar(0);

   for( auto p1It = ep.begin(); p1It != ep.end(); ++p1It )
      for(auto p2It = std::next( p1It ); p2It != ep.end(); ++p2It)
      {
         const Scalar d = (*p1It - *p2It).sqrnorm() ;
         if( d > maxSqDist )
         {
            base0 = *p1It;
            base1 = *p2It;
            maxSqDist = d;
         }
      }

   // find the most distant point in ep to the previously found line
   Point dir = base0 - base1;
   Scalar maxDist = Scalar(0);
   Point base2;

   for(const auto & p : ep)
   {
      const Scalar d = ((p - base0) % dir).sqrnorm();
      if(d > maxDist)
      {
         base2 = p;
         maxDist = d;
      }
   }

   // find the most distant point in the cloud to the previously found triangle
   Point basePlaneNormal = (base0 - base2) % (base1 - base2);
   Scalar basePlaneOffset = base0 | basePlaneNormal;

   Point tip;
   Scalar maxDistanceToBasePlane =  Scalar(0);
   for(const auto & p : pointCloud_)
   {
      const Scalar d = std::fabs( ( basePlaneNormal | p ) - basePlaneOffset);
      if(d > maxDistanceToBasePlane)
      {
         tip = p;
         maxDistanceToBasePlane = d;
      }
   }

   const VertexHandle base0h = mesh_->add_vertex( base0 );
   const VertexHandle base1h = mesh_->add_vertex( base1 );
   const VertexHandle base2h = mesh_->add_vertex( base2 );
   const VertexHandle tiph   = mesh_->add_vertex( tip   );

   if( ( basePlaneNormal | tip ) - basePlaneOffset < Scalar(0) )
   {
      mesh_->add_face( base0h, base1h, base2h );
      mesh_->add_face( base0h, base2h, tiph );
      mesh_->add_face( base2h, base1h, tiph );
      mesh_->add_face( base1h, base0h, tiph );
   }
   else
   {
      mesh_->add_face( base2h, base1h, base0h );
      mesh_->add_face( base2h, base0h, tiph );
      mesh_->add_face( base1h, base2h, tiph );
      mesh_->add_face( base0h, base1h, tiph );
   }

   mesh_->update_face_normals();

   pointCloud_.erase( std::remove_if( pointCloud_.begin(), pointCloud_.end(), [&](const Point & p) {
      return p == base0 || p == base1 || p == base2 || p == tip;
   }), pointCloud_.end() );

   assignPointsToFaces( pointCloud_, mesh_->faces_begin(), mesh_->faces_end() );

   for(auto fh : mesh_->faces())
   {
      if( !visiblePoints_[fh].empty() )
         queue_.push(fh);
   }

   if(writeDebugOutput_)
   {
      (*meshWriter_)();
      remainingPointsWriter_->write();
   }
}


template< typename MeshType >
void QHull<MeshType>::deleteVisibleFaces( const FaceHandle startFaceHandle, const Point & p )
{
   horizon_.clear();
   orphanPoints_.clear();

   std::queue<HalfedgeHandle> q;

   // Push all opposing halfedges from the start face to the queue and delete the start face
   // The half edges mark an initial horizon that gets pushed further outwards later
   for(auto heh : mesh_->fh_range( startFaceHandle ))
   {
      q.push( mesh_->opposite_halfedge_handle(heh) );
   }

   orphanPoints_.insert( orphanPoints_.end(), visiblePoints_[startFaceHandle].begin(), visiblePoints_[startFaceHandle].end() );

   mesh_->delete_face( startFaceHandle, true );

   // Move the horizon farther outwards
   while(!q.empty())
   {
      HalfedgeHandle heh = q.front();
      q.pop();

      // If the face belonging to the half edge  has already been removed just go on with the next half edge
      if( mesh_->status( mesh_->edge_handle(heh) ).deleted() )
         continue;

      // Test if p is visible from the face that lies beyond the current horizon half edge
      auto fh = mesh_->face_handle(heh);
      if(pointIsVisible( fh, p ))
      {
         // the face will be deleted, but first push its opposing half edges to the stack
         for(auto neighbor_heh : mesh_->fh_range( fh ))
         {
            if( !mesh_->status(neighbor_heh).deleted() )
               q.push( mesh_->opposite_halfedge_handle(neighbor_heh) );
         }
         orphanPoints_.insert( orphanPoints_.end(), visiblePoints_[fh].begin(), visiblePoints_[fh].end() );
         mesh_->delete_face( fh, true );
      }
      else // p not visible from face fh
      {
         // The half edge is part of our final horizon
         horizon_.push_back( heh );
      }
   }
}


template< typename MeshType >
void QHull<MeshType>::addNewFaces( const Point & p )
{
   newFaces_.clear();

   VertexHandle v0h = mesh_->add_vertex(p);
   for(const auto heh : horizon_)
   {
      VertexHandle v1h = mesh_->to_vertex_handle(heh);
      VertexHandle v2h = mesh_->from_vertex_handle(heh);

      FaceHandle fh = mesh_->add_face( v0h, v1h, v2h );
      mesh_->update_normal(fh);
      newFaces_.push_back(fh);         
   }

}


template< typename MeshType >
template< typename InputIterator >
void QHull<MeshType>::assignPointsToFaces( const std::vector<Point> & points, InputIterator facesBegin, InputIterator facesEnd )
{
   for(const auto & p : points)
   {
      for(auto fhIt = facesBegin; fhIt != facesEnd; ++fhIt)
      {
         if(pointIsVisible( *fhIt, p ))
         {
            visiblePoints_[*fhIt].push_back(p);
            break;
         }
      }
   }
}


template< typename MeshType >
void QHull<MeshType>::enableDebugVTKOutput( const std::string & identifierPrefix, const std::string & baseFolder ) 
{
   meshWriter_ = walberla::make_shared< VTKMeshWriter< MeshType > >( meshPtr(), identifierPrefix + "Mesh", 1, baseFolder );
   meshWriter_->addDataSource( make_shared<IndexFaceDataSource<MeshType>>() );

   auto pointDataSource = make_shared<QHullPointDataSource<MeshType>>( *this );
   remainingPointsWriter_ = vtk::createVTKOutput_PointData( pointDataSource, identifierPrefix + "RemainingPoints", 1, baseFolder);

   writeDebugOutput_ = true;
}


template< typename MeshType >
void QHull<MeshType>::disableDebugVTKOutput() 
{
   writeDebugOutput_ = false;
   meshWriter_.reset();
   remainingPointsWriter_.reset();
}


template< typename MeshType >
void QHull<MeshType>::initMesh() 
{
   mesh_->clear();
   mesh_->request_face_normals();
   mesh_->request_face_status();
   mesh_->request_edge_status();
   mesh_->request_vertex_status();
   mesh_->request_halfedge_status();
}

} // namespace mesh
} // namespace walberla
