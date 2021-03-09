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
//! \file Node.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/GenericAABB.h"
#include "core/math/KahanSummation.h"

namespace walberla {
namespace geometry {
namespace containment_octree {

template< typename ContainmentOctreeT >
class Node
{
public:
   using Point = typename ContainmentOctreeT::Point;
   using Scalar = typename ContainmentOctreeT::Scalar;
   using AABB = typename ContainmentOctreeT::AABB;

   using DistanceObject = typename ContainmentOctreeT::DistanceObject;
   using KahanAccumulator = typename ContainmentOctreeT::KahanAccumulator;

   virtual ~Node() = default;
   virtual bool contains( const Point & p ) const = 0;
   virtual uint_t height() const = 0;
   virtual uint_t numNodes() const = 0;
   virtual void numNodes( uint_t & numInside, uint_t & numOutside, uint_t & numIndeterminate, uint_t & numBranch ) const = 0;
   virtual void volumes( KahanAccumulator & insideVolume, KahanAccumulator & outsideVolume, KahanAccumulator & indeterminateVolume, Scalar volume ) const = 0;
   virtual uint_t numChildren() const = 0;
   virtual const Node * getChild( const uint_t idx ) const = 0;

   static inline Point             toPoint( const Vector3<real_t> & p ) { return ContainmentOctreeT::toPoint( p ); }
   static inline Vector3<real_t> fromPoint( const Point & p )           { return ContainmentOctreeT::fromPoint( p ); }

   static inline Scalar   toScalar( const real_t & x ) { return ContainmentOctreeT::toScalar( x ); }
   static inline real_t fromScalar( const Scalar & x ) { return ContainmentOctreeT::fromScalar( x ); }
};

} // namespace containment_octree
} // namespace geometry
} // namespace walberla
