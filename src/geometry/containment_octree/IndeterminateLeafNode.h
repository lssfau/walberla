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
//! \file IndeterminateLeafNode.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "LeafNode.h"

#include "core/DataTypes.h"

namespace walberla {
namespace geometry {
namespace containment_octree {

template< typename ContainmentOctreeT >
class IndeterminateLeafNode : public LeafNode<ContainmentOctreeT>
{
public:
   using LeafNode<ContainmentOctreeT>::numNodes;
   typedef typename LeafNode<ContainmentOctreeT>::DistanceObject DistanceObject;
   typedef typename LeafNode<ContainmentOctreeT>::Scalar Scalar;
   typedef typename LeafNode<ContainmentOctreeT>::Point Point;
   typedef typename LeafNode<ContainmentOctreeT>::AABB AABB;
   
   typedef typename LeafNode<ContainmentOctreeT>::KahanAccumulator KahanAccumulator;
   
   IndeterminateLeafNode( const shared_ptr<const DistanceObject> & distanceObject, const Scalar epsilon )
      :  distanceObject_( distanceObject ), sqEpsilon_( epsilon * epsilon ) { }

   virtual ~IndeterminateLeafNode() = default;

   virtual bool contains( const Point & p ) const { return distanceObject_->sqSignedDistance(p) <= sqEpsilon_; }

   virtual void numNodes( uint_t & /*numInside*/, uint_t & /*numOutside*/, uint_t & numIndeterminate, uint_t & /*numBranch*/ ) const { ++numIndeterminate; }
   virtual void volumes( KahanAccumulator & /*insideVolume*/, KahanAccumulator & /*outsideVolume*/, KahanAccumulator & indeterminateVolume, Scalar volume ) const { indeterminateVolume += volume; }

protected:
   shared_ptr<const DistanceObject> distanceObject_;
   Scalar sqEpsilon_;
};


   
} // namespace containment_octree
} // namespace geometry
} // namespace walberla
