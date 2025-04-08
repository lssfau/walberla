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
//! \file InsideLeafNode.h
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
class InsideLeafNode : public LeafNode<ContainmentOctreeT>
{
public:
   using LeafNode<ContainmentOctreeT>::numNodes;
   using DistanceObject = typename LeafNode<ContainmentOctreeT>::DistanceObject;
   using Scalar = typename LeafNode<ContainmentOctreeT>::Scalar;
   using Point = typename LeafNode<ContainmentOctreeT>::Point;
   using AABB = typename LeafNode<ContainmentOctreeT>::AABB;
    
   using KahanAccumulator = typename LeafNode<ContainmentOctreeT>::KahanAccumulator;

   ~InsideLeafNode() override = default;

   bool contains( const Point & /*p*/ ) const override { return true; }

   void numNodes( uint_t & numInside, uint_t & /*numOutside*/, uint_t & /*numIndeterminate*/, uint_t & /*numBranch*/ ) const override { ++numInside; }
   void volumes( KahanAccumulator & insideVolume, KahanAccumulator & /*outsideVolume*/, KahanAccumulator & /*indeterminateVolume*/, Scalar volume ) const override { insideVolume += volume; }
};

   
} // namespace containment_octree
} // namespace geometry
} // namespace walberla
