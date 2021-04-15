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
//! \file DefaultTesselation.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "tesselation/Box.h"
#include "tesselation/ConvexPolyhedron.h"

#include "pe/rigidbody/Box.h"

namespace walberla {
namespace mesh {
namespace pe {

template< typename MeshType >
class DefaultTesselation
{
public:
   typedef typename MeshType::VertexHandle VertexHandle;
   typedef typename MeshType::FaceHandle   FaceHandle;

   void operator()( const walberla::pe::RigidBody & body, MeshType & mesh, std::vector<VertexHandle> & newVertices, std::vector<FaceHandle> & newFaces )
   {
      const id_t typeId = body.getTypeID();
      if(typeId == walberla::pe::Box::getStaticTypeID())
      {
         (*this)( dynamic_cast<const walberla::pe::Box &>(body), mesh, newVertices, newFaces );
      }
      else if(typeId == ConvexPolyhedron::getStaticTypeID())
      {
         (*this)( dynamic_cast<const ConvexPolyhedron &>(body), mesh, newVertices, newFaces );
      }
      else
      {
         WALBERLA_ABORT( "Tessellation not implemented for your body!" );
      }
   }

   void operator()( const walberla::pe::Box & body, MeshType & mesh, std::vector<VertexHandle> & newVertices, std::vector<FaceHandle> & newFaces )
   {
      tesselateBox( body, mesh, newVertices, newFaces );
   }

   void operator()( const ConvexPolyhedron & body, MeshType & mesh, std::vector<VertexHandle> & newVertices, std::vector<FaceHandle> & newFaces )
   {
      tesselateConvexPolyhedron( body, mesh, newVertices, newFaces );
   }

};

} // namespace pe
} // namespace mesh
} // namespace walberla