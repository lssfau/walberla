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
//! \file MeshConversionTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh/MeshConversion.h"
#include "mesh_common/MeshIO.h"

namespace walberla {
namespace mesh {

template< typename MeshType >
void test()
{
   MeshType omMesh;
   geometry::TriangleMesh wbMesh;

   mesh::readAndBroadcast("cube.obj", omMesh);
   geometry::readAndBroadcastMesh( "cube.obj", wbMesh );

   geometry::TriangleMesh convWbMesh;
   mesh::convertOpenMeshToWalberla( omMesh, convWbMesh );

   MeshType convOmMesh;
   mesh::convertWalberlaToOpenMesh( wbMesh, convOmMesh );

   WALBERLA_CHECK_EQUAL( omMesh.n_faces(), wbMesh.getNumTriangles() );
   WALBERLA_CHECK_EQUAL( wbMesh.getNumTriangles(), convWbMesh.getNumTriangles() );
   WALBERLA_CHECK_EQUAL( convWbMesh.getNumTriangles(), convOmMesh.n_faces() );

   WALBERLA_CHECK_EQUAL( omMesh.n_vertices(), wbMesh.getNumVertices() );
   WALBERLA_CHECK_EQUAL( wbMesh.getNumVertices(), convWbMesh.getNumVertices() );
   WALBERLA_CHECK_EQUAL( convWbMesh.getNumVertices(), convOmMesh.n_vertices() );
}

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   test< mesh::TriangleMesh >();
   test< mesh::FloatTriangleMesh >();
   test< mesh::PythonTriangleMesh >();

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}