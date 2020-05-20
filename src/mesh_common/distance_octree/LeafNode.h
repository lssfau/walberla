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
//! \file LeafNode.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Node.h"

#include "mesh_common/DistanceComputations.h"

namespace walberla {
namespace mesh {
namespace distance_octree {

template <typename MeshType>
class LeafNode : public Node<MeshType>
{
public:
   typedef typename Node<MeshType>::Point      Point;
   typedef typename Node<MeshType>::Normal     Normal;
   typedef typename Node<MeshType>::Scalar     Scalar;
   typedef typename Node<MeshType>::FaceHandle FaceHandle; 
   typedef typename Node<MeshType>::AABB       AABB;
   
   LeafNode( const shared_ptr< TriangleDistance<MeshType> > & triDistance, const std::vector<FaceHandle> & triangles )
      : Node<MeshType>( triDistance->getMesh(), triangles.begin(), triangles.end() ), triangles_( triangles ), triDistance_( triDistance ) { }

   virtual Scalar sqSignedDistance( const Point & p ) const;
   virtual Scalar sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const;
   virtual Scalar sqSignedDistance( const Point & p, Point & closestPoint ) const;
   virtual Scalar sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const;

   virtual Scalar sqDistance( const Point & p ) const;
   virtual Scalar sqDistance( const Point & p, FaceHandle & closestTriangle ) const;
   virtual Scalar sqDistance( const Point & p, Point & closestPoint ) const;
   virtual Scalar sqDistance( const Point & p, Point & closestPoint, Normal & normal ) const;

   uint_t numTriangles() const { return uint_c( triangles_.size() ); }
   void numTrianglesToStream( std::ostream & os, const uint_t level ) const;
   virtual uint_t height() const { return 0; }
   virtual uint_t numChildren() const { return 0; };
   virtual const Node<MeshType> * getChild( const uint_t /*idx*/ ) const { WALBERLA_ABORT("DistanceOctree: You are requesting access to children of a Leaf Node!"); return 0; }
   
protected:
   std::vector<FaceHandle> triangles_;
   shared_ptr< TriangleDistance<MeshType> > triDistance_;
};


template <typename MeshType>
void LeafNode<MeshType>::numTrianglesToStream( std::ostream & os, const uint_t level ) const
{
   for( uint_t i = 0; i < level; ++i )
      os << "   ";
   os << numTriangles() << "\n";
}


template <typename MeshType>
typename LeafNode<MeshType>::Scalar LeafNode<MeshType>::sqSignedDistance( const Point & p ) const
{
   if( triangles_.empty() )
      return std::numeric_limits<Scalar>::max();

   return triDistance_->sqSignedDistance( triangles_.begin(), triangles_.end(), p ) ;
}


template <typename MeshType>
typename LeafNode<MeshType>::Scalar LeafNode<MeshType>::sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const
{
   if( triangles_.empty() )
      return std::numeric_limits<Scalar>::max();

   return triDistance_->sqSignedDistance( triangles_.begin(), triangles_.end(), p, closestTriangle ) ;
}


template <typename MeshType>
typename LeafNode<MeshType>::Scalar LeafNode<MeshType>::sqSignedDistance( const Point & p, Point & closestPoint ) const
{
   if( triangles_.empty() )
      return std::numeric_limits<Scalar>::max();

   return triDistance_->sqSignedDistance( triangles_.begin(), triangles_.end(), p, closestPoint ) ;
}


template <typename MeshType>
typename LeafNode<MeshType>::Scalar LeafNode<MeshType>::sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const
{
   if( triangles_.empty() )
      return std::numeric_limits<Scalar>::max();

   return triDistance_->sqSignedDistance( triangles_.begin(), triangles_.end(), p, closestPoint, normal ) ;
}


template <typename MeshType>
typename LeafNode<MeshType>::Scalar LeafNode<MeshType>::sqDistance( const Point & p ) const
{
   if(triangles_.empty())
      return std::numeric_limits<Scalar>::max();

   return triDistance_->sqDistance( triangles_.begin(), triangles_.end(), p );
}


template <typename MeshType>
typename LeafNode<MeshType>::Scalar LeafNode<MeshType>::sqDistance( const Point & p, FaceHandle & closestTriangle ) const
{
   if(triangles_.empty())
      return std::numeric_limits<Scalar>::max();

   return triDistance_->sqDistance( triangles_.begin(), triangles_.end(), p, closestTriangle );
}


template <typename MeshType>
typename LeafNode<MeshType>::Scalar LeafNode<MeshType>::sqDistance( const Point & p, Point & closestPoint ) const
{
   if(triangles_.empty())
      return std::numeric_limits<Scalar>::max();

   return triDistance_->sqDistance( triangles_.begin(), triangles_.end(), p, closestPoint );
}


template <typename MeshType>
typename LeafNode<MeshType>::Scalar LeafNode<MeshType>::sqDistance( const Point & p, Point & closestPoint, Normal & normal ) const
{
   if(triangles_.empty())
      return std::numeric_limits<Scalar>::max();

   return triDistance_->sqDistance( triangles_.begin(), triangles_.end(), p, closestPoint, normal );
}


} // namespace distance_octree
} // namespace mesh
} // namespace walberla