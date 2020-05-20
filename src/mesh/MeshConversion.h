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
//! \file MeshConversion.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================
#pragma once

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MatrixVectorOperations.h"

#include "geometry/mesh/TriangleMesh.h"

namespace walberla {
namespace mesh {

/**
* \brief Converts an OpenMesh to a geometry::TriangleMesh
*
* Only the basic topology is converted. Additional mesh attributes like normals are ignored.
* 
* \tparam MeshType The type of the OpenMesh
*
* \param omMesh The OpenMesh source mesh
* \param wbMesh The waLBerla target mesh
* \param clear  If true, the waLBerla mesh is cleared before data from omMesh is added
*/
template< typename OpenMeshType >
void convertOpenMeshToWalberla( const OpenMeshType & omMesh, geometry::TriangleMesh & wbMesh, bool clear = false )
{
   static_assert( OpenMeshType::IsTriMesh == 1, "convertOpenMeshToWalberla only works with triangular meshes!" );

   if( clear )
      wbMesh.clear();

   std::map< typename OpenMeshType::VertexHandle, uint32_t > omVertexHandlesToIndex;

   for( auto v_it = omMesh.vertices_begin(); v_it != omMesh.vertices_end(); ++v_it )
   {
      omVertexHandlesToIndex[ *v_it ] = wbMesh.addVertex( toWalberla( omMesh.point( *v_it ) ) );
   }

   for(auto f_it = omMesh.faces_begin(); f_it != omMesh.faces_end(); ++f_it)
   {
      typename OpenMeshType::VertexHandle vh0, vh1, vh2;
      getVertexHandles( omMesh, *f_it, vh0, vh1, vh2 );
      wbMesh.addTriangle( omVertexHandlesToIndex[vh0], omVertexHandlesToIndex[vh1], omVertexHandlesToIndex[vh2] );
   }

   WALBERLA_ASSERT_EQUAL( omMesh.n_faces(), wbMesh.getNumTriangles() );
   WALBERLA_ASSERT_EQUAL( omMesh.n_vertices(), wbMesh.getNumVertices() );
}


/**
* \brief Converts a geometry::TriangleMesh to an OpenMesh
*
* Only the basic topology is converted. Additional mesh attributes like normals are ignored.
* 
* \tparam MeshType The type of the OpenMesh
*
* \param wbMesh The waLBerla source mesh
* \param omMesh The OpenMesh target mesh
* \param clear  If true, the OpenMesh is cleared before data from wbMesh is added
*/
template< typename OpenMeshType >
void convertWalberlaToOpenMesh( const geometry::TriangleMesh & wbMesh, OpenMeshType & omMesh, bool clear = false )
{
   static_assert( OpenMeshType::IsTriMesh == 1, "convertOpenMeshToWalberla only works with triangular meshes!" );

   if( clear )
      omMesh.clear();

   std::map< uint32_t, typename OpenMeshType::VertexHandle > indexToOmVertexHandles;

   uint32_t ctr = 0;
   for( auto v_it = wbMesh.getVertices().begin(); v_it != wbMesh.getVertices().end(); ++v_it )
   {
      indexToOmVertexHandles[ ctr++ ] = omMesh.add_vertex( typename OpenMeshType::Point( toOpenMesh( *v_it ) ) );
   }

   WALBERLA_ASSERT_EQUAL( wbMesh.getVertexIndices().size() % size_t(3), size_t(0) );
   auto vi_it = wbMesh.getVertexIndices().begin();
   while( vi_it != wbMesh.getVertexIndices().end() )
   {
      uint32_t vi0 = *vi_it++;
      uint32_t vi1 = *vi_it++;
      uint32_t vi2 = *vi_it++;
      omMesh.add_face( indexToOmVertexHandles[ vi0 ], indexToOmVertexHandles[ vi1 ], indexToOmVertexHandles[ vi2 ] );
   }

   WALBERLA_ASSERT_EQUAL( omMesh.n_faces(), wbMesh.getNumTriangles() );
   WALBERLA_ASSERT_EQUAL( omMesh.n_vertices(), wbMesh.getNumVertices() );
}



} // namespace mesh
} // namespace walberla