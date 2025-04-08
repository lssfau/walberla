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
   using Point = typename Node<MeshType>::Point;
   using Normal = typename Node<MeshType>::Normal;
   using Scalar = typename Node<MeshType>::Scalar;
   using FaceHandle = typename Node<MeshType>::FaceHandle; 
   using AABB = typename Node<MeshType>::AABB;
   
   LeafNode( const shared_ptr< TriangleDistance<MeshType> > & triDistance, const std::vector<FaceHandle> & triangles )
      : Node<MeshType>( triDistance->getMesh(), triangles.begin(), triangles.end() ), triangles_( triangles ), triDistance_( triDistance ) { }

   Scalar sqSignedDistance( const Point & p ) const override;
   Scalar sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const override;
   Scalar sqSignedDistance( const Point & p, Point & closestPoint ) const override;
   Scalar sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const override;

   Scalar sqDistance( const Point & p ) const override;
   Scalar sqDistance( const Point & p, FaceHandle & closestTriangle ) const override;
   Scalar sqDistance( const Point & p, Point & closestPoint ) const override;
   Scalar sqDistance( const Point & p, Point & closestPoint, Normal & normal ) const override;

   Scalar getRayDistanceToMeshObject(const Point & ray_origin, const Point & normalised_ray_direction) const override;
   
   uint_t numTriangles() const override { return uint_c( triangles_.size() ); }
   void numTrianglesToStream( std::ostream & os, const uint_t level ) const override;
   uint_t height() const override { return 0; }
   uint_t numChildren() const override { return 0; };
   const Node<MeshType> * getChild( const uint_t /*idx*/ ) const override { WALBERLA_ABORT("DistanceOctree: You are requesting access to children of a Leaf Node!"); return nullptr; }
   
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
   if(triangles_.empty()){
      return std::numeric_limits<Scalar>::max();
   }


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

template< typename MeshType >
typename LeafNode<MeshType>::Scalar LeafNode< MeshType >::getRayDistanceToMeshObject(const Point& ray_origin, const Point& normalised_ray_direction) const
{
   if(triangles_.empty())
      return std::numeric_limits<Scalar>::max();

   return triDistance_->getRayDistanceToMeshObject( triangles_.begin(), triangles_.end(), ray_origin, normalised_ray_direction );
}

} // namespace distance_octree
} // namespace mesh
} // namespace walberla