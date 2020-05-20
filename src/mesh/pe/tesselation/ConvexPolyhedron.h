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
//! \file ConvexPolyhedron.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesh/pe/rigid_body/ConvexPolyhedron.h"
#include "mesh_common/MatrixVectorOperations.h"
#include "mesh_common/TriangleMeshes.h"

namespace walberla {
namespace mesh {
namespace pe {

template<typename MeshType>
void tesselateConvexPolyhedron( const ConvexPolyhedron & poly, MeshType & mesh, std::vector<typename MeshType::VertexHandle> & newVertices, std::vector<typename MeshType::FaceHandle> & newFaces )
{
   typedef typename MeshType::Scalar Scalar;
   typedef typename MeshType::Point Point;
   typedef typename MeshType::VertexHandle VertexHandle;

   const TriangleMesh & polyhedronMesh = poly.getMesh();

   newVertices.reserve( polyhedronMesh.n_vertices() );
   for(auto vh : polyhedronMesh.vertices())
   {
      const TriangleMesh::Point & p0 = polyhedronMesh.point( vh );
      Vec3 p1 = poly.pointFromBFtoWF( p0[0], p0[1], p0[2] );
      Point p2( numeric_cast<Scalar>( p1[0]), numeric_cast<Scalar>( p1[1]), numeric_cast<Scalar>( p1[2]) );
      newVertices.push_back( mesh.add_vertex(p2) );
   }

   newFaces.reserve( polyhedronMesh.n_faces() );
   for(auto fh : polyhedronMesh.faces())
   {
      auto it = polyhedronMesh.cfv_ccwbegin(fh);
      VertexHandle v0 = newVertices[ size_t( (it++)->idx() ) ];
      VertexHandle v1 = newVertices[ size_t( (it++)->idx() ) ];
      VertexHandle v2 = newVertices[ size_t( (it++)->idx() ) ];
      WALBERLA_ASSERT_EQUAL(it, polyhedronMesh.cfv_ccwend(fh) );
      newFaces.push_back( mesh.add_face(v0, v1, v2) );
   }

   if(mesh.has_face_normals())
      for( auto fh: newVertices )
         mesh.update_normal(fh);

   if(mesh.has_vertex_normals())
      for( auto vh: newVertices )
         mesh.update_normal(vh);
}

} // namespace pe
} // namespace mesh
} // namespace walberla
