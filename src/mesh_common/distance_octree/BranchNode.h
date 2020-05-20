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

#include <algorithm>
#include <iterator>

namespace walberla {
namespace mesh {
namespace distance_octree {

template <typename MeshType>
class BranchNode : public Node<MeshType>
{
public:
   typedef typename Node<MeshType>::Point      Point;
   typedef typename Node<MeshType>::Normal     Normal;
   typedef typename Node<MeshType>::Scalar     Scalar;  
   typedef typename Node<MeshType>::FaceHandle FaceHandle; 
   typedef typename Node<MeshType>::AABB       AABB;
   
   template< typename InputIterator >
   BranchNode( const shared_ptr< TriangleDistance<MeshType> > & triDistance, InputIterator beginFh, InputIterator endFh,
               uint_t maxDepth, uint_t minNumTriangles );

   virtual ~BranchNode() { for( int i = 0; i < 8; ++i ) delete children_[i]; }


   virtual Scalar sqSignedDistance( const Point & p ) const;
   virtual Scalar sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const;
   virtual Scalar sqSignedDistance( const Point & p, Point & closestPoint ) const ;
   virtual Scalar sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const;

   virtual Scalar sqDistance( const Point & p ) const;
   virtual Scalar sqDistance( const Point & p, FaceHandle & closestTriangle ) const;
   virtual Scalar sqDistance( const Point & p, Point & closestPoint ) const;
   virtual Scalar sqDistance( const Point & p, Point & closestPoint, Normal & normal ) const;

   inline uint_t numTriangles() const;
   void numTrianglesToStream( std::ostream & os, const uint_t level ) const;
   inline virtual uint_t height() const;
   virtual uint_t numChildren() const { return uint_t(8); };
   virtual const Node<MeshType> * getChild( const uint_t idx ) const { WALBERLA_ASSERT_LESS( idx, 8 ); return children_[idx]; };

private:
   BranchNode( const BranchNode & other );
   BranchNode & operator=( const BranchNode & other );

   struct ChildInfo
   {
      ChildInfo( const Node<MeshType> * _child, const Point & p )
         : child( _child ), minSqBoxDist( child->getAABB().sqDistance( toWalberla( p ) ) )
      { }

      bool operator<( const ChildInfo & other ) const { return minSqBoxDist < other.minSqBoxDist; }

      const Node<MeshType> * child;
      Scalar minSqBoxDist;
   };

protected:
   const Node<MeshType> * children_[8];
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
   for( int i = 1; i < 8; ++i )
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
   for( int i = 0; i < 8; ++i )
      children_[i] = NULL;

   const auto &    min = this->aabb_.minCorner();
   const auto &    max = this->aabb_.maxCorner();
   const auto   center = this->aabb_.center();

   AABB childAABBs[8] = {
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

   std::vector<bool> triangleUsed( theNumTriangles, false );
   std::vector<FaceHandle> childTriangles[8];

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
   ChildInfo childinfos[8] = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( childinfos, childinfos + 8 );

   Scalar absMinSqSignedDistance = childinfos[0].child->sqSignedDistance( p );
   for( int i = 1; i < 8; ++i )
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
   ChildInfo childinfos[8] = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( childinfos, childinfos + 8 );

   Scalar absMinSqSignedDistance = childinfos[0].child->sqSignedDistance( p, closestTriangle );

   for( int i = 1; i < 8; ++i )
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
   ChildInfo childinfos[8] = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( childinfos, childinfos + 8 );

   Scalar absMinSqSignedDistance = childinfos[0].child->sqSignedDistance( p, closestPoint );

   for( int i = 1; i < 8; ++i )
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
   ChildInfo childinfos[8] = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( childinfos, childinfos + 8 );

   Scalar absMinSqSignedDistance = childinfos[0].child->sqSignedDistance( p, closestPoint, normal );

   for( int i = 1; i < 8; ++i )
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

   ChildInfo childinfos[8] = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( childinfos, childinfos + 8 );

   Scalar absMinSqDistance = childinfos[0].child->sqDistance( p );
   for(int i = 1; i < 8; ++i)
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
   ChildInfo childinfos[8] = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( childinfos, childinfos + 8 );

   Scalar absMinSqDistance = childinfos[0].child->sqDistance( p, closestTriangle );

   for(int i = 1; i < 8; ++i)
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
   ChildInfo childinfos[8] = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( childinfos, childinfos + 8 );

   Scalar absMinSqDistance = childinfos[0].child->sqDistance( p, closestPoint );

   for(int i = 1; i < 8; ++i)
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
   ChildInfo childinfos[8] = {
      ChildInfo( children_[0], p ), ChildInfo( children_[1], p ),
      ChildInfo( children_[2], p ), ChildInfo( children_[3], p ),
      ChildInfo( children_[4], p ), ChildInfo( children_[5], p ),
      ChildInfo( children_[6], p ), ChildInfo( children_[7], p )
   };
   std::sort( childinfos, childinfos + 8 );

   Scalar absMinSqDistance = childinfos[0].child->sqDistance( p, closestPoint, normal );

   for(int i = 1; i < 8; ++i)
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
void BranchNode<MeshType>::numTrianglesToStream( std::ostream & os, const uint_t level ) const
{
   for( uint_t i = 0; i < level; ++i )
      os << "   ";
   os << numTriangles() << "\n";
   for( int i = 0; i < 8; ++i )
      children_[i]->numTrianglesToStream(os, level + 1);

}


} // namespace distance_octree
} // namespace mesh
} // namespace walberla