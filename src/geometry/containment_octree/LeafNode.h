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
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Node.h"

#include "core/Abort.h"
#include "core/DataTypes.h"

namespace walberla {
namespace geometry {
namespace containment_octree {

template< typename ContainmentOctreeT >
class LeafNode : public Node<ContainmentOctreeT>
{
public:
   using Node<ContainmentOctreeT>::numNodes;
   using DistanceObject = typename Node<ContainmentOctreeT>::DistanceObject;
   using Scalar = typename Node<ContainmentOctreeT>::Scalar;
   using Point = typename Node<ContainmentOctreeT>::Point;
   using AABB = typename Node<ContainmentOctreeT>::AABB;
   
   using KahanAccumulator = typename Node<ContainmentOctreeT>::KahanAccumulator;

   virtual ~LeafNode() = default;

   virtual uint_t height() const { return uint_t(0); }
   virtual uint_t numNodes() const { return uint_t(0); }
   virtual uint_t numChildren() const { return uint_t(0); }

   virtual const Node<ContainmentOctreeT> * getChild( const uint_t ) const { WALBERLA_ABORT("ContainmentOctree: You are requesting access to children of a Leaf Node!"); return 0; }
};


} // namespace containment_octree
} // namespace geometry
} // namespace walberla
