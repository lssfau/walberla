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
//! \file MeshDistanceOctreeTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

#include "mesh_common/OpenMeshBufferTypeExtensions.h"

#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>

namespace walberla {
namespace mesh {

template< typename T >
void test( const T & in )
{
   mpi::SendBuffer sendBuffer;
   sendBuffer << in;

   T out;
   mpi::RecvBuffer recvBuffer( sendBuffer );
   recvBuffer >> out;

   WALBERLA_CHECK_EQUAL( in, out );
}

int main( int /*argc*/, char * /*argv*/[] )
{
   debug::enterTestMode();
   
   test( OpenMesh::VertexHandle( 42 ) );
   test( OpenMesh::FaceHandle( 42 ) );
   test( OpenMesh::EdgeHandle( 42 ) );
   test( OpenMesh::HalfedgeHandle( 42 ) );

   test( OpenMesh::Vec3d( 1.0, 2.0, 3.0) );
   test( OpenMesh::Vec3f( 1.0f, 2.0f, 3.0f) );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}