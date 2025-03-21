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
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Node.h"
#include "IndeterminateLeafNode.h"
#include "InsideLeafNode.h"
#include "OutsideLeafNode.h"

#include "core/DataTypes.h"

namespace walberla {
namespace geometry {
namespace containment_octree {

template< typename ContainmentOctreeT >
class BranchNode : public Node<ContainmentOctreeT>
{
public:
   using DistanceObject = typename Node<ContainmentOctreeT>::DistanceObject;
   using Scalar = typename Node<ContainmentOctreeT>::Scalar;
   using Point = typename Node<ContainmentOctreeT>::Point;
   using AABB = typename Node<ContainmentOctreeT>::AABB;

   using KahanAccumulator = typename Node<ContainmentOctreeT>::KahanAccumulator;

   inline BranchNode( const shared_ptr<const DistanceObject> & distanceObject, const AABB & aabb, const Scalar epsilon,
                      const uint_t maxDepth, const Scalar minAABBVolume );

   ~BranchNode() override { for(auto & i : children_) delete i; }

   inline bool contains( const Point & p ) const override;

   inline uint_t height() const override;
   inline uint_t numNodes() const override;
   inline void numNodes( uint_t & numInside, uint_t & numOutside, uint_t & numIndeterminate, uint_t & numBranch ) const override;
   void volumes( KahanAccumulator & insideVolume, KahanAccumulator & outsideVolume, KahanAccumulator & indeterminateVolume, Scalar volume ) const override;
   uint_t numChildren() const override { return 8; }
   const Point & center() const { return center_; }

   const Node<ContainmentOctreeT> * getChild( const uint_t idx ) const override { WALBERLA_ASSERT_LESS( idx, 8 ); return children_[idx]; }

   BranchNode( const BranchNode & other ) = delete;
   BranchNode & operator=( const BranchNode & other ) = delete;

protected:
   std::array< const Node<ContainmentOctreeT> *, 8 >  children_;
   Point center_;
};


template< typename ContainmentOctreeT >
BranchNode<ContainmentOctreeT>::BranchNode( const shared_ptr<const DistanceObject> & distanceObject, const AABB & aabb, const Scalar epsilon,
                                            const uint_t maxDepth, const Scalar minAABBVolume ) : center_( this->toPoint( aabb.center() ) )
{
   for( uint_t i = 0; i < uint_t(8); ++i )
      children_[i] = nullptr;

   const auto & min = aabb.minCorner();
   const auto & max = aabb.maxCorner();
   const auto & ctr = center_;

   std::array<AABB, 8> childAABBs = {
      AABB::createFromMinMaxCorner( min[0], min[1], min[2], ctr[0], ctr[1], ctr[2] ),
      AABB::createFromMinMaxCorner( min[0], min[1], ctr[2], ctr[0], ctr[1], max[2] ),
      AABB::createFromMinMaxCorner( min[0], ctr[1], min[2], ctr[0], max[1], ctr[2] ),
      AABB::createFromMinMaxCorner( min[0], ctr[1], ctr[2], ctr[0], max[1], max[2] ),
      AABB::createFromMinMaxCorner( ctr[0], min[1], min[2], max[0], ctr[1], ctr[2] ),
      AABB::createFromMinMaxCorner( ctr[0], min[1], ctr[2], max[0], ctr[1], max[2] ),
      AABB::createFromMinMaxCorner( ctr[0], ctr[1], min[2], max[0], max[1], ctr[2] ),
      AABB::createFromMinMaxCorner( ctr[0], ctr[1], ctr[2], max[0], max[1], max[2] )
   };

   auto halfAABBDimensions = childAABBs[0].sizes() * Scalar(0.5);
   halfAABBDimensions[0] += epsilon;
   halfAABBDimensions[1] += epsilon;
   halfAABBDimensions[2] += epsilon;
   Scalar maxSqDistanceCenterAABB = halfAABBDimensions.sqrLength();

   auto childAABBIt = std::begin(childAABBs);
   for( auto it = std::begin(children_); it != std::end(children_); ++it, ++childAABBIt )
   {
      const auto childAABBCenter = childAABBIt->center();

      const Scalar sqSignedDist = distanceObject->sqSignedDistance( this->toPoint( childAABBCenter ) );

      if( std::fabs( sqSignedDist ) > maxSqDistanceCenterAABB )
      {
         if( sqSignedDist <= Scalar(0) )
         {
            *it = new InsideLeafNode< ContainmentOctreeT >();
         }
         else
         {
            *it = new OutsideLeafNode< ContainmentOctreeT >();
         }
      }
      else
      {
         if( maxDepth == 0 || childAABBIt->volume() < minAABBVolume )
         {
            *it = new IndeterminateLeafNode< ContainmentOctreeT >( distanceObject, epsilon );
         }
         else
         {
            *it = new BranchNode< ContainmentOctreeT >( distanceObject, *childAABBIt, epsilon, maxDepth - 1, minAABBVolume );
         }
      }
   }
}


template< typename ContainmentOctreeT >
bool BranchNode<ContainmentOctreeT>::contains( const Point & p ) const
{
   if( p[0] < center_[0] )
   {
      if( p[1] < center_[1] )
      {
         return p[2] < center_[2] ? children_[0]->contains(p) : children_[1]->contains(p);
      }
      else // p[1] >= center_[1]
      {
         return p[2] < center_[2] ? children_[2]->contains(p) : children_[3]->contains(p);
      }
   }
   else // p[0] >= center_[0]
   {
      if( p[1] < center_[1] )
      {
         return p[2] < center_[2] ? children_[4]->contains(p) : children_[5]->contains(p);
      }
      else // p[1] >= center_[1]
      {
         return p[2] < center_[2] ? children_[6]->contains(p) : children_[7]->contains(p);
      }
   }
}


template< typename ContainmentOctreeT >
uint_t BranchNode<ContainmentOctreeT>::height() const
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


template< typename ContainmentOctreeT >
uint_t BranchNode<ContainmentOctreeT>::numNodes() const
{
   uint_t nodes = 1;
   for(auto & i : children_)
      nodes += i->numNodes();

   return nodes;
}


template< typename ContainmentOctreeT >
void BranchNode<ContainmentOctreeT>::numNodes( uint_t & numInside, uint_t & numOutside, uint_t & numIndeterminate, uint_t & numBranch ) const
{
   ++numBranch;
   for(auto & i : children_)
      i->numNodes( numInside, numOutside, numIndeterminate, numBranch );
}


template< typename ContainmentOctreeT >
void BranchNode<ContainmentOctreeT>::volumes( KahanAccumulator & insideVolume, KahanAccumulator & outsideVolume, KahanAccumulator & indeterminateVolume, Scalar volume ) const
{
   static const Scalar ONE_OVER_EIGHT = Scalar(1) / Scalar(8);
   for(auto const & i : children_)
      i->volumes( insideVolume, outsideVolume, indeterminateVolume, volume * ONE_OVER_EIGHT );
}





} // namespace containment_octree
} // namespace geometry
} // namespace walberla