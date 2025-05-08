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
//! \file BranchNode.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Node.h"

#include "core/debug/Debug.h"

#include "mesh_common/DistanceComputations.h"
#include "mesh_common/MatrixVectorOperations.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/distance_octree/Ray.h"

#include <algorithm>
#include <iterator>
#include <optional>

namespace walberla {
namespace mesh {
namespace distance_octree {

template <typename MeshType>
class BranchNode : public Node<MeshType>
{
public:
   using Point = typename Node<MeshType>::Point;
   using Normal = typename Node<MeshType>::Normal;
   using Scalar = typename Node<MeshType>::Scalar;
   using FaceHandle = typename Node<MeshType>::FaceHandle;
   using AABB = typename Node<MeshType>::AABB;

   template< typename InputIterator >
   BranchNode( const shared_ptr< TriangleDistance<MeshType> > & triDistance, InputIterator beginFh, InputIterator endFh,
               uint_t maxDepth, uint_t minNumTriangles );

   ~BranchNode() override { for( uint_t i = 0; i < 8; ++i ) delete children_[i]; }


   Scalar sqSignedDistance( const Point & p ) const override;
   Scalar sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const override;
   Scalar sqSignedDistance( const Point & p, Point & closestPoint ) const override ;
   Scalar sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const override;

   Scalar sqDistance( const Point & p ) const override;
   Scalar sqDistance( const Point & p, FaceHandle & closestTriangle ) const override;
   Scalar sqDistance( const Point & p, Point & closestPoint ) const override;
   Scalar sqDistance( const Point & p, Point & closestPoint, Normal & normal ) const override;

   Scalar getRayDistanceToMeshObject(const Point & ray_origin, const Point & normalised_ray_direction) const override;
   
   inline uint_t numTriangles() const override;
   void numTrianglesToStream( std::ostream & os, const uint_t level ) const override;
   inline uint_t height() const override;
   uint_t numChildren() const override { return uint_t(8); };
   const Node<MeshType> * getChild( const uint_t idx ) const override { WALBERLA_ASSERT_LESS( idx, 8 ); return children_[idx]; };

   BranchNode( const BranchNode & other ) = delete;
   BranchNode & operator=( const BranchNode & other ) = delete;
private:

   struct ChildInfo
   {
      ChildInfo( const Node<MeshType> * _child, const Point & p )
         : child( _child ), minSqBoxDist( child->getAABB().sqDistance( toWalberla( p ) ) )
      { }

      bool operator<( const ChildInfo & other ) const { return minSqBoxDist < other.minSqBoxDist; }

      const Node<MeshType> * child;
      Scalar minSqBoxDist;
   };
/*
   struct ChildInfoIntersects
   {
      ChildInfoIntersects(const Node<MeshType> * _child, const Point & ray_origin, const Point & ray_direction)
         :child( _child ), intersectsAabb( pe::raytracing::intersects( child->getAABB(), 
                                                                       pe::raytracing::Ray(toWalberla( ray_origin ), (toWalberla( ray_direction )).getNormalized()),
                                                                       parametricDistance,
                                                                       real_t(0.0), 
                                                                       &normal
                                                                     )) {}

      bool operator<( const ChildInfoIntersects & other ) const 
         { return parametricDistance < other.parametricDistance; }
   const Node<MeshType> * child;
      real_t parametricDistance;
      pe::Vec3 normal;
      bool intersectsAabb;
   };
*/
   struct ChildInfoIntersects
   {
      ChildInfoIntersects(const Node<MeshType> * _child, const real_t & _parametricDistance, 
                          const Vector3<real_t> & _normal, const bool & _intersects)
         :child( _child ), parametricDistance(_parametricDistance), normal(_normal), intersectsAabb(_intersects){}

      static ChildInfoIntersects fromRay(const Node<MeshType> * child, const Point & ray_origin, const Point & ray_direction) {
         real_t distance;
         Vector3<real_t> ray_normal;

         bool intersects { mesh::rayAABBIntersection(   child->getAABB(),
                                                        mesh::Ray(
                                                               toWalberla( ray_origin ), 
                                                               (toWalberla( ray_direction )).getNormalized()),
                                                         distance, real_t(0.0), &ray_normal
                                                      )
                        };
                  
         return ChildInfoIntersects(child, distance, ray_normal, intersects);
      }

      bool operator<(const ChildInfoIntersects & other) const {
         return parametricDistance < other.parametricDistance;
      }

      const Node<MeshType> * child;
      const real_t parametricDistance;
      const Vector3<real_t> normal;
      const bool intersectsAabb;
   };

protected:
   std::array<const Node<MeshType> *, 8> children_;
};


template <typename MeshType>
uint_t BranchNode<MeshType>::numTriangles() const
{
   return children_[0]->numTriangles() + children_[1]->numTriangles()
        + children_[2]->numTriangles() + children_[3]->numTriangles()
        + children_[4]->numTriangles() + children_[5]->numTriangles()
        + children_[6]->numTriangles() + children_[7]->numTriangles();
}


template <typename MeshType>
uint_t BranchNode<MeshType>::height() const
{
   uint_t maxChildHeight = children_[0]->height();
   for( uint_t i = 1; i < 8; ++i )
   {
      uint_t childHeight = children_[i]->height();
      if( childHeight > maxChildHeight )
         maxChildHeight = childHeight;
   }

   return maxChildHeight + 1;
}


template <typename MeshType>
template< typename InputIterator >
BranchNode<MeshType>::BranchNode( const shared_ptr< TriangleDistance<MeshType> > & triDistance, InputIterator beginFh, InputIterator endFh,
                                  uint_t maxDepth, uint_t minNumTriangles )
   : Node<MeshType>( triDistance->getMesh(), beginFh, endFh )
{
   for( uint_t i = 0; i < 8; ++i )
      children_[i] = nullptr;

   const auto &    min = this->aabb_.minCorner();
   const auto &    max = this->aabb_.maxCorner();
   const auto   center = this->aabb_.center();

   std::array<AABB, 8> childAABBs = {
      AABB::createFromMinMaxCorner(    min[0],    min[1],    min[2], center[0], center[1], center[2] ),
      AABB::createFromMinMaxCorner(    min[0],    min[1], center[2], center[0], center[1],    max[2] ),
      AABB::createFromMinMaxCorner(    min[0], center[1],    min[2], center[0],    max[1], center[2] ),
      AABB::createFromMinMaxCorner(    min[0], center[1], center[2], center[0],    max[1],    max[2] ),
      AABB::createFromMinMaxCorner( center[0],    min[1],    min[2],    max[0], center[1], center[2] ),
      AABB::createFromMinMaxCorner( center[0],    min[1], center[2],    max[0], center[1],    max[2] ),
      AABB::createFromMinMaxCorner( center[0], center[1],    min[2],    max[0],    max[1], center[2] ),
      AABB::createFromMinMaxCorner( center[0], center[1], center[2],    max[0],    max[1],    max[2] )
   };

   uint_t theNumTriangles = uint_c( std::distance( beginFh, endFh ) );

   std::vector<bool> const triangleUsed( theNumTriangles, false );
   std::array<std::vector<FaceHandle>, 8> childTriangles;

   for( auto fhIt = beginFh; fhIt != endFh; ++fhIt )
   {
      auto centroid = computeCentroid( triDistance->getMesh(), *fhIt );

      uint_t minIdx = 0;
      Scalar minSqDistance = ( childAABBs[0].center() - toWalberla( centroid ) ).sqrLength();
      for(uint_t i = 1; i < uint_t( 8 ); ++i)
      {
         Scalar theSqDistance = ( childAABBs[i].center() - toWalberla( centroid ) ).sqrLength();
         if(theSqDistance < minSqDistance)
         {
            minSqDistance = theSqDistance;
            minIdx = i;
         }
      }
      childTriangles[minIdx].push_back( *fhIt );
   }

   for( uint_t i = 0; i < uint_t(8); ++i )
   {
      if( maxDepth == 0 || childTriangles[i].size() < minNumTriangles || theNumTriangles == childTriangles[i].size() )
         children_[i] = new LeafNode<MeshType>( triDistance, childTriangles[i] );
      else
         children_[i] = new BranchNode<MeshType>( triDistance, childTriangles[i].begin(), childTriangles[i].end(), maxDepth - 1, minNumTriangles );
   }
}


template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::sqSignedDistance( const Point & p ) const
{
   std::array<ChildInfo, 8> childinfos = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( std::begin(childinfos), std::end(childinfos) );

   Scalar absMinSqSignedDistance = childinfos[0].child->sqSignedDistance( p );
   for( uint_t i = 1; i < 8; ++i )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( childinfos[i].child );
      if( std::fabs( absMinSqSignedDistance ) < childinfos[i].minSqBoxDist )
         continue;

      Scalar theSqSignedDistance = childinfos[i].child->sqSignedDistance( p );
      if( std::fabs( theSqSignedDistance )  <  std::fabs( absMinSqSignedDistance ) )
         absMinSqSignedDistance = theSqSignedDistance;
   }

   WALBERLA_ASSERT_LESS( absMinSqSignedDistance, std::numeric_limits<Scalar>::max() );
   return absMinSqSignedDistance;

}


template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const
{
   std::array<ChildInfo, 8> childinfos = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( std::begin(childinfos), std::end(childinfos) );

   Scalar absMinSqSignedDistance = childinfos[0].child->sqSignedDistance( p, closestTriangle );

   for( uint_t i = 1; i < 8; ++i )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( childinfos[i].child );
      if( std::fabs( absMinSqSignedDistance ) < childinfos[i].minSqBoxDist )
         continue;

      FaceHandle triangle;
      Scalar theSqSignedDistance = childinfos[i].child->sqSignedDistance( p, triangle );
      if( std::fabs( theSqSignedDistance )  <  std::fabs( absMinSqSignedDistance ) )
      {
         absMinSqSignedDistance = theSqSignedDistance;
         closestTriangle = triangle;
      }
   }

   WALBERLA_ASSERT_LESS( absMinSqSignedDistance, std::numeric_limits<Scalar>::max() );
   return absMinSqSignedDistance;
}


template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::sqSignedDistance( const Point & p, Point & closestPoint ) const
{
   std::array<ChildInfo, 8> childinfos = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( std::begin(childinfos), std::end(childinfos) );

   Scalar absMinSqSignedDistance = childinfos[0].child->sqSignedDistance( p, closestPoint );

   for( uint_t i = 1; i < 8; ++i )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( childinfos[i].child );
      if( std::fabs( absMinSqSignedDistance ) < childinfos[i].minSqBoxDist )
         continue;

      Point point;
      Scalar theSqSignedDistance = childinfos[i].child->sqSignedDistance( p, point );
      if( std::fabs( theSqSignedDistance )  <  std::fabs( absMinSqSignedDistance ) )
      {
         absMinSqSignedDistance = theSqSignedDistance;
         closestPoint = point;
      }
   }

   WALBERLA_ASSERT_LESS( absMinSqSignedDistance, std::numeric_limits<Scalar>::max() );
   return absMinSqSignedDistance;
}


template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const
{
   std::array<ChildInfo, 8> childinfos = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( std::begin(childinfos), std::end(childinfos) );

   Scalar absMinSqSignedDistance = childinfos[0].child->sqSignedDistance( p, closestPoint, normal );

   for( uint_t i = 1; i < 8; ++i )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( childinfos[i].child );
      if( std::fabs( absMinSqSignedDistance ) < childinfos[i].minSqBoxDist )
         continue;

      Point point;
      Normal tmpNormal;
      Scalar theSqSignedDistance = childinfos[i].child->sqSignedDistance( p, point, tmpNormal );
      if( std::fabs( theSqSignedDistance )  <  std::fabs( absMinSqSignedDistance ) )
      {
         absMinSqSignedDistance = theSqSignedDistance;
         closestPoint = point;
         normal = tmpNormal;
      }
   }

   WALBERLA_ASSERT_LESS( absMinSqSignedDistance, std::numeric_limits<Scalar>::max() );
   return absMinSqSignedDistance;
}



template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::sqDistance( const Point & p ) const
{
   //WALBERLA_ASSERT( getAABB().contains( toWalberla( p ) ) );

   std::array<ChildInfo, 8> childinfos = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( std::begin(childinfos), std::end(childinfos) );

   Scalar absMinSqDistance = childinfos[0].child->sqDistance( p );
   for(uint_t i = 1; i < 8; ++i)
   {
      WALBERLA_ASSERT_NOT_NULLPTR( childinfos[i].child );
      if( absMinSqDistance < childinfos[i].minSqBoxDist)
         continue;

      Scalar theSqDistance = childinfos[i].child->sqDistance( p );
      if( theSqDistance < absMinSqDistance )
         absMinSqDistance = theSqDistance;
   }

   WALBERLA_ASSERT_LESS( absMinSqDistance, std::numeric_limits<Scalar>::max() );
   return absMinSqDistance;

}


template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::sqDistance( const Point & p, FaceHandle & closestTriangle ) const
{
   std::array<ChildInfo, 8> childinfos = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( std::begin(childinfos), std::end(childinfos) );

   Scalar absMinSqDistance = childinfos[0].child->sqDistance( p, closestTriangle );

   for(uint_t i = 1; i < 8; ++i)
   {
      WALBERLA_ASSERT_NOT_NULLPTR( childinfos[i].child );
      if( absMinSqDistance < childinfos[i].minSqBoxDist)
         continue;

      FaceHandle triangle;
      Scalar theSqDistance = childinfos[i].child->sqDistance( p, triangle );
      if( theSqDistance < absMinSqDistance )
      {
         absMinSqDistance = theSqDistance;
         closestTriangle = triangle;
      }
   }

   WALBERLA_ASSERT_LESS( absMinSqDistance, std::numeric_limits<Scalar>::max() );
   return absMinSqDistance;
}


template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::sqDistance( const Point & p, Point & closestPoint ) const
{
   std::array<ChildInfo, 8> childinfos = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( std::begin(childinfos), std::end(childinfos) );

   Scalar absMinSqDistance = childinfos[0].child->sqDistance( p, closestPoint );

   for(uint_t i = 1; i < 8; ++i)
   {
      WALBERLA_ASSERT_NOT_NULLPTR( childinfos[i].child );
      if( absMinSqDistance < childinfos[i].minSqBoxDist)
         continue;

      Point point;
      Scalar theSqDistance = childinfos[i].child->sqDistance( p, point );
      if(theSqDistance < absMinSqDistance )
      {
         absMinSqDistance = theSqDistance;
         closestPoint = point;
      }
   }

   WALBERLA_ASSERT_LESS( absMinSqDistance, std::numeric_limits<Scalar>::max() );
   return absMinSqDistance;
}


template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::sqDistance( const Point & p, Point & closestPoint, Normal & normal ) const
{
   std::array<ChildInfo, 8> childinfos = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( std::begin(childinfos), std::end(childinfos) );

   Scalar absMinSqDistance = childinfos[0].child->sqDistance( p, closestPoint, normal );

   for(uint_t i = 1; i < 8; ++i)
   {
      WALBERLA_ASSERT_NOT_NULLPTR( childinfos[i].child );
      if(absMinSqDistance < childinfos[i].minSqBoxDist)
         continue;

      Point point;
      Normal tmpNormal;
      Scalar theSqDistance = childinfos[i].child->sqDistance( p, point, tmpNormal );
      if(theSqDistance < absMinSqDistance)
      {
         absMinSqDistance = theSqDistance;
         closestPoint = point;
         normal = tmpNormal;
      }
   }

   WALBERLA_ASSERT_LESS( absMinSqDistance, std::numeric_limits<Scalar>::max() );
   return absMinSqDistance;
}


template <typename MeshType>
typename BranchNode<MeshType>::Scalar BranchNode<MeshType>::getRayDistanceToMeshObject(const Point & ray_origin, const Point & normalised_ray_direction) const
{
   ChildInfoIntersects childinfos[8] = {
      ChildInfoIntersects::fromRay( children_[0], ray_origin, normalised_ray_direction ), 
      ChildInfoIntersects::fromRay( children_[1], ray_origin, normalised_ray_direction ),
      ChildInfoIntersects::fromRay( children_[2], ray_origin, normalised_ray_direction ), 
      ChildInfoIntersects::fromRay( children_[3], ray_origin, normalised_ray_direction ),
      ChildInfoIntersects::fromRay( children_[4], ray_origin, normalised_ray_direction ), 
      ChildInfoIntersects::fromRay( children_[5], ray_origin, normalised_ray_direction ),
      ChildInfoIntersects::fromRay( children_[6], ray_origin, normalised_ray_direction ), 
      ChildInfoIntersects::fromRay( children_[7], ray_origin, normalised_ray_direction )
   };

   Scalar distance( std::numeric_limits<Scalar>::max() );

   for (const auto& childinfo : childinfos) {
      if (childinfo.intersectsAabb) 
      {
         WALBERLA_ASSERT_NOT_NULLPTR(childinfo.child);

         Scalar newDistance = childinfo.child->getRayDistanceToMeshObject(ray_origin, normalised_ray_direction);
         
         if (newDistance < distance)
            distance = newDistance;
      }
   }

   return distance;
}

template <typename MeshType>
void BranchNode<MeshType>::numTrianglesToStream( std::ostream & os, const uint_t level ) const
{
   for( uint_t i = 0; i < level; ++i )
      os << "   ";
   os << numTriangles() << "\n";
   for( uint_t i = 0; i < 8; ++i )
      children_[i]->numTrianglesToStream(os, level + 1);

}


} // namespace distance_octree
} // namespace mesh
} // namespace walberla
