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
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Marshalling of objects for data transmission or storage.
//
//======================================================================================================================

#include "ConvexPolyhedron.h"

#include "core/mpi/BufferDataTypeExtensions.h"

#include "mesh_common/OpenMeshBufferTypeExtensions.h"

#include "pe/communication/Marshalling.h"

namespace walberla {
namespace mesh {
namespace pe {

void marshal( mpi::SendBuffer& buffer, const ConvexPolyhedron& obj ) {
   walberla::pe::communication::marshal( buffer, static_cast<const GeomPrimitive&>( obj ) );
   
   const TriangleMesh & m = obj.getMesh();

   WALBERLA_DEBUG_SECTION()
   {
      buffer << std::string("poly_begin");
   }

   buffer << m.n_vertices();

   int dbgIndex = 0;
   WALBERLA_UNUSED(dbgIndex);
   for(auto vh : m.vertices())
   {
      WALBERLA_ASSERT_EQUAL( vh.idx(), dbgIndex++ ); // assume vertices are compactly stored
      buffer << m.point(vh);
   }

   buffer << m.n_faces();

   for( auto fh: m.faces() )
      for(auto vhIt = m.cfv_ccwbegin(fh); vhIt != m.cfv_ccwend(fh); ++vhIt)
      {
         WALBERLA_ASSERT_GREATER_EQUAL( vhIt->idx(), 0 );
         WALBERLA_ASSERT_LESS( vhIt->idx(), m.n_vertices() );
         buffer << vhIt->idx();
      }

   WALBERLA_DEBUG_SECTION()
   {
      buffer << std::string("poly_end");
   }
}
//*************************************************************************************************


void unmarshal( mpi::RecvBuffer& buffer, ConvexPolyhedronParameters& objparam ) {
   walberla::pe::communication::unmarshal( buffer, static_cast<GeomPrimitiveParameters&>( objparam ) );
   
   WALBERLA_DEBUG_SECTION()
   {
      std::string s;
      buffer >> s;
      WALBERLA_ASSERT_EQUAL( s, std::string("poly_begin") );
   }

   TriangleMesh & m = objparam.mesh_;

   size_t numVertices;
   buffer >> numVertices;

   std::vector<TriangleMesh::VertexHandle> vertexHandles(numVertices);
   for(size_t i = 0; i < numVertices; ++i)
   {
      TriangleMesh::Point p;
      buffer >> p;
      vertexHandles[size_t(i)] = m.add_vertex( p );
   }

   size_t numFaces;
   buffer >> numFaces;
   for(size_t i = 0; i < numFaces; ++i)
   {
      int v0;
      int v1;
      int v2;
      buffer >> v0 >> v1 >> v2;
      WALBERLA_ASSERT_GREATER_EQUAL( v0, 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( v1, 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( v2, 0 );
      WALBERLA_ASSERT_LESS( v0, numVertices );
      WALBERLA_ASSERT_LESS( v1, numVertices );
      WALBERLA_ASSERT_LESS( v2, numVertices );

      m.add_face( vertexHandles[size_t(v0)], vertexHandles[size_t(v1)], vertexHandles[size_t(v2)] );
   }

   WALBERLA_DEBUG_SECTION()
   {
      std::string s;
      buffer >> s;
      WALBERLA_ASSERT_EQUAL( s, std::string("poly_end") );
   }

}
//*************************************************************************************************

}  // namespace pe
}  // namespace mesh
}  // namespace walberla
