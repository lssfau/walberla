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
//! \file TriangleMeshes.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace OpenMesh {
namespace Python {

struct MeshTraits : public OpenMesh::DefaultTraits {
   /** Use double precision points */
   using Point = OpenMesh::Vec3d;

   /** Use double precision normals */
   using Normal = OpenMesh::Vec3d;

   /** Use RGBA colors */
   using Color = OpenMesh::Vec4f;
};

}
}


namespace walberla {
namespace mesh {

using PythonTriangleMesh = OpenMesh::TriMesh_ArrayKernelT<OpenMesh::Python::MeshTraits>;

struct RealTraits : public OpenMesh::DefaultTraits
{
   using Point = OpenMesh::VectorT<real_t, 3>;
   using Normal = OpenMesh::VectorT<real_t, 3>;
};
using TriangleMesh = OpenMesh::TriMesh_ArrayKernelT<RealTraits>;


struct FloatTraits : public OpenMesh::DefaultTraits
{
   using Point = OpenMesh::Vec3f;
   using Normal = OpenMesh::Vec3f;
};
using FloatTriangleMesh = OpenMesh::TriMesh_ArrayKernelT<FloatTraits>;


template< typename MeshType >
inline void getVertexHandles( const MeshType & mesh, const typename MeshType::FaceHandle fh,
                              typename MeshType::VertexHandle & vh0, typename MeshType::VertexHandle & vh1, typename MeshType::VertexHandle & vh2 )
{
   static_assert( MeshType::IsTriMesh == 1, "getVertexHandles only works with triangular meshes!" );

   auto v_it = mesh.cfv_ccwbegin( fh );
   vh0 = *v_it++;
   vh1 = *v_it++;
   vh2 = *v_it;
   WALBERLA_ASSERT_EQUAL( ++v_it, mesh.cfv_ccwend( fh ) );
}


template< typename MeshType >
inline void getVertexPositions( const MeshType & mesh, const typename MeshType::FaceHandle fh,
                                typename MeshType::Point & v0, typename MeshType::Point & v1, typename MeshType::Point & v2 )
{
   static_assert( MeshType::IsTriMesh == 1, "getVertexPositions only works with triangular meshes!" );

   auto v_it = mesh.cfv_ccwbegin( fh );
   v0 = mesh.point( *v_it++ );
   v1 = mesh.point( *v_it++ );
   v2 = mesh.point( *v_it );
   WALBERLA_ASSERT_EQUAL( ++v_it, mesh.cfv_ccwend( fh ) );
}


template< typename MeshType >
inline typename MeshType::VertexHandle getVertexHandle( const MeshType & mesh, const typename MeshType::FaceHandle fh, const unsigned int vertexIdx )
{
   static_assert( MeshType::IsTriMesh == 1, "getVertexHandle only works with triangular meshes!" );
   WALBERLA_ASSERT_LESS( vertexIdx, 3U );

   auto vIt = mesh.cfv_ccwbegin( fh );
   std::advance( vIt, vertexIdx );
   return *vIt;
}

template< typename MeshType >
inline typename MeshType::HalfedgeHandle getHalfedgeHandle( const MeshType & mesh, const typename MeshType::FaceHandle fh,
                                                            const unsigned int fromVertexIdx, const unsigned int toVertexIdx )
{
   static_assert( MeshType::IsTriMesh == 1, "getHalfEdgeHandle only works with triangular meshes!" );

   auto heh = mesh.find_halfedge( getVertexHandle( mesh, fh, fromVertexIdx ), getVertexHandle( mesh, fh, toVertexIdx ) );

   WALBERLA_ASSERT( heh.is_valid() );

   return heh;
}


} // namespace mesh
} // namespace walberla