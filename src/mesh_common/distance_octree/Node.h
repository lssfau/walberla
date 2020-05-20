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
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/GenericAABB.h"

#include "mesh_common/MeshOperations.h"

namespace walberla {
namespace mesh {
namespace distance_octree {

template <typename MeshType>
class Node
{
public:
   typedef typename MeshType::Point                    Point;
   typedef typename MeshType::Normal                   Normal;
   typedef typename MeshType::Scalar                   Scalar;  
   typedef typename MeshType::FaceHandle               FaceHandle; 
   typedef typename math::GenericAABB<Scalar>          AABB;

   template< typename InputIterator >
   Node( const MeshType & mesh, InputIterator beginFh, InputIterator endFh ) : aabb_( computeAABBForFaces( mesh, beginFh, endFh ) ) {}
   virtual ~Node() { }

   const AABB & getAABB() const { return aabb_; }

   virtual Scalar sqSignedDistance( const Point & p ) const = 0;
   virtual Scalar sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const = 0;
   virtual Scalar sqSignedDistance( const Point & p, Point & closestPoint ) const = 0;
   virtual Scalar sqSignedDistance( const Point & p, Point & closestPoint, Point & normal ) const = 0;

   virtual Scalar sqDistance( const Point & p ) const = 0;
   virtual Scalar sqDistance( const Point & p, FaceHandle & closestTriangle ) const = 0;
   virtual Scalar sqDistance( const Point & p, Point & closestPoint ) const = 0;
   virtual Scalar sqDistance( const Point & p, Point & closestPoint, Point & normal ) const = 0;

   virtual uint_t numTriangles() const = 0;
   virtual void numTrianglesToStream( std::ostream & os, const uint_t level ) const = 0;
   virtual uint_t height() const = 0;
   virtual uint_t numChildren() const = 0;
   virtual const Node * getChild( const uint_t idx ) const = 0;

protected:
   AABB aabb_;
};




} // namespace distance_octree
} // namespace mesh
} // namespace walberla