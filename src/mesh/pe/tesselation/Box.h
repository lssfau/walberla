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
//! \file Box.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesh_common/MatrixVectorOperations.h"

#include "pe/rigidbody/Box.h"

namespace walberla {
namespace mesh {
namespace pe {

template<typename MeshType>
void tesselateBox( const walberla::pe::Box & box, MeshType & mesh, std::vector<typename MeshType::VertexHandle> & newVertices, std::vector<typename MeshType::FaceHandle> & newFaces )
{
   typedef typename MeshType::Scalar Scalar;
   typedef typename MeshType::VertexHandle VertexHandle;
   typedef typename MeshType::FaceHandle FaceHandle;

   const Vector3<real_t> d = real_t(0.5) * box.getLengths();

   const VertexHandle nnn = mesh.add_vertex( toOpenMeshNumericCast<Scalar>( box.pointFromBFtoWF( Vector3<real_t>( -d[0], -d[1], -d[2] ) ) ) );
   const VertexHandle nnp = mesh.add_vertex( toOpenMeshNumericCast<Scalar>( box.pointFromBFtoWF( Vector3<real_t>( -d[0], -d[1], +d[2] ) ) ) );
   const VertexHandle npn = mesh.add_vertex( toOpenMeshNumericCast<Scalar>( box.pointFromBFtoWF( Vector3<real_t>( -d[0], +d[1], -d[2] ) ) ) );
   const VertexHandle npp = mesh.add_vertex( toOpenMeshNumericCast<Scalar>( box.pointFromBFtoWF( Vector3<real_t>( -d[0], +d[1], +d[2] ) ) ) );
   const VertexHandle pnn = mesh.add_vertex( toOpenMeshNumericCast<Scalar>( box.pointFromBFtoWF( Vector3<real_t>( +d[0], -d[1], -d[2] ) ) ) );
   const VertexHandle pnp = mesh.add_vertex( toOpenMeshNumericCast<Scalar>( box.pointFromBFtoWF( Vector3<real_t>( +d[0], -d[1], +d[2] ) ) ) );
   const VertexHandle ppn = mesh.add_vertex( toOpenMeshNumericCast<Scalar>( box.pointFromBFtoWF( Vector3<real_t>( +d[0], +d[1], -d[2] ) ) ) );
   const VertexHandle ppp = mesh.add_vertex( toOpenMeshNumericCast<Scalar>( box.pointFromBFtoWF( Vector3<real_t>( +d[0], +d[1], +d[2] ) ) ) );

   newVertices.push_back( nnn );
   newVertices.push_back( nnp );
   newVertices.push_back( npn );
   newVertices.push_back( npp );
   newVertices.push_back( pnn );
   newVertices.push_back( pnp );
   newVertices.push_back( ppn );
   newVertices.push_back( ppp );

   if(mesh.has_vertex_normals())
   {
      mesh.set_normal( nnn, toOpenMeshNumericCast<Scalar>( Vector3<real_t>( -d[0], -d[1], -d[2] ).getNormalized() ) );
      mesh.set_normal( nnp, toOpenMeshNumericCast<Scalar>( Vector3<real_t>( -d[0], -d[1], +d[2] ).getNormalized() ) );
      mesh.set_normal( npn, toOpenMeshNumericCast<Scalar>( Vector3<real_t>( -d[0], +d[1], -d[2] ).getNormalized() ) );
      mesh.set_normal( npp, toOpenMeshNumericCast<Scalar>( Vector3<real_t>( -d[0], +d[1], +d[2] ).getNormalized() ) );
      mesh.set_normal( pnn, toOpenMeshNumericCast<Scalar>( Vector3<real_t>( +d[0], -d[1], -d[2] ).getNormalized() ) );
      mesh.set_normal( pnp, toOpenMeshNumericCast<Scalar>( Vector3<real_t>( +d[0], -d[1], +d[2] ).getNormalized() ) );
      mesh.set_normal( ppn, toOpenMeshNumericCast<Scalar>( Vector3<real_t>( +d[0], +d[1], -d[2] ).getNormalized() ) );
      mesh.set_normal( ppp, toOpenMeshNumericCast<Scalar>( Vector3<real_t>( +d[0], +d[1], +d[2] ).getNormalized() ) );
   }

   const FaceHandle f0 = mesh.add_face( nnn, pnn, pnp, nnp );
   const FaceHandle f1 = mesh.add_face( npn, npp, ppp, ppn );
   const FaceHandle f2 = mesh.add_face( npp, nnp, pnp, ppp );
   const FaceHandle f3 = mesh.add_face( npn, ppn, pnn, nnn );
   const FaceHandle f4 = mesh.add_face( npp, npn, nnn, nnp );
   const FaceHandle f5 = mesh.add_face( ppp, pnp, pnn, ppn );

   if(mesh.has_face_normals())
   {
      mesh.update_normal( f0 );
      mesh.update_normal( f1 );
      mesh.update_normal( f2 );
      mesh.update_normal( f3 );
      mesh.update_normal( f4 );
      mesh.update_normal( f5 );
   }

   newFaces.push_back( f0 );
   newFaces.push_back( f1 );
   newFaces.push_back( f2 );
   newFaces.push_back( f3 );
   newFaces.push_back( f4 );
   newFaces.push_back( f5 );
}

} // namespace pe
} // namespace mesh
} // namespace walberla
