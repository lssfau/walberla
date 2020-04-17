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
//! \file MeshVTKTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/vtk/VTKMeshWriter.h"
#include "mesh_common/vtk/CommonDataSources.h"
#include "mesh_common/vtk/CommonFilters.h"

#include <vector>
#include <string>

namespace walberla {
namespace mesh {

template< typename MeshType >
void test( const std::string & meshFile )
{
   auto mesh = make_shared< MeshType >();
   mesh::readAndBroadcast( meshFile, *mesh);

   auto aabb = computeAABB( *mesh );

   aabb.scale( typename MeshType::Scalar( 0.9 ) );

   mesh->request_face_normals();
   mesh->request_vertex_normals();
   mesh->update_normals();

   mesh->request_face_status();

   bool b = true;
   for( auto it = mesh->faces_begin(); it != mesh->faces_end(); ++it )
   {
      mesh->status( *it ).set_tagged( b );
      b = !b;
   }

   mesh->request_vertex_status();

   b = true;
   for( auto it = mesh->vertices_begin(); it != mesh->vertices_end(); ++it )
   {
      mesh->status( *it ).set_tagged( b );
      b = !b;
   }

   std::vector< typename MeshType::Color > colors;
   colors.push_back( typename MeshType::Color( 255,0,0 ) );
   colors.push_back( typename MeshType::Color( 0,255,0 ) );
   colors.push_back( typename MeshType::Color( 0,0,255 ) );
   auto colorIt = colors.begin();

   mesh->request_vertex_colors();
   for( auto it = mesh->vertices_begin(); it != mesh->vertices_end(); ++it )
   {
      mesh->set_color( *it, *colorIt );

      ++colorIt;
      if( colorIt == colors.end() )
         colorIt = colors.begin();
   }

   mesh->request_face_colors();
   for( auto it = mesh->faces_begin(); it != mesh->faces_end(); ++it )
   {
      mesh->set_color( *it, *colorIt );

      ++colorIt;
      if( colorIt == colors.end() )
         colorIt = colors.begin();
   }

   VTKMeshWriter< MeshType > meshWriter( mesh, "mesh_vtk_test_unfiltered", 1 );
   meshWriter.addDataSource( make_shared< NormalsVertexDataSource< MeshType > >() );
   meshWriter.addDataSource( make_shared< NormalsFaceDataSource  < MeshType > >() );

   meshWriter.addDataSource( make_shared< AreaFaceDataSource     < MeshType > >() );

   meshWriter.addDataSource( make_shared< StatusBitFaceDataSource  < MeshType > >( OpenMesh::Attributes::TAGGED, "tagged" ) );
   meshWriter.addDataSource( make_shared< StatusBitVertexDataSource< MeshType > >( OpenMesh::Attributes::TAGGED, "tagged" ) );

   meshWriter.addDataSource( make_shared< ColorFaceDataSource  < MeshType > >() );
   meshWriter.addDataSource( make_shared< ColorVertexDataSource< MeshType > >() );

   meshWriter.addDataSource( make_shared< IndexFaceDataSource  < MeshType > >() );
   meshWriter.addDataSource( make_shared< IndexVertexDataSource< MeshType > >() );

   meshWriter();


   VTKMeshWriter< MeshType > meshWriterAABBfiltered( mesh, "mesh_vtk_test_aabb_filter", 1 );
   meshWriterAABBfiltered.addDataSource( make_shared< NormalsVertexDataSource< MeshType > >() );
   meshWriterAABBfiltered.addDataSource( make_shared< NormalsFaceDataSource  < MeshType > >() );

   meshWriterAABBfiltered.addDataSource( make_shared< AreaFaceDataSource     < MeshType > >() );

   meshWriterAABBfiltered.addDataSource( make_shared< StatusBitFaceDataSource  < MeshType > >( OpenMesh::Attributes::TAGGED, "tagged" ) );
   meshWriterAABBfiltered.addDataSource( make_shared< StatusBitVertexDataSource< MeshType > >( OpenMesh::Attributes::TAGGED, "tagged" ) );

   meshWriterAABBfiltered.addDataSource( make_shared< ColorFaceDataSource  < MeshType > >() );
   meshWriterAABBfiltered.addDataSource( make_shared< ColorVertexDataSource< MeshType > >() );

   meshWriter.addDataSource( make_shared< IndexFaceDataSource  < MeshType > >() );
   meshWriter.addDataSource( make_shared< IndexVertexDataSource< MeshType > >() );

   meshWriterAABBfiltered.setFaceFilter( mesh::AABBFaceFilter< MeshType >(aabb) );

   meshWriterAABBfiltered();


   VTKMeshWriter< MeshType > meshWriterStatusfiltered( mesh, "mesh_vtk_test_status_filter", 1 );
   meshWriterStatusfiltered.addDataSource( make_shared< NormalsVertexDataSource< MeshType > >() );
   meshWriterStatusfiltered.addDataSource( make_shared< NormalsFaceDataSource  < MeshType > >() );

   meshWriterStatusfiltered.addDataSource( make_shared< AreaFaceDataSource     < MeshType > >() );

   meshWriterStatusfiltered.addDataSource( make_shared< StatusBitFaceDataSource  < MeshType > >( OpenMesh::Attributes::TAGGED, "tagged" ) );
   meshWriterStatusfiltered.addDataSource( make_shared< StatusBitVertexDataSource< MeshType > >( OpenMesh::Attributes::TAGGED, "tagged" ) );

   meshWriterStatusfiltered.addDataSource( make_shared< ColorFaceDataSource  < MeshType > >() );
   meshWriterStatusfiltered.addDataSource( make_shared< ColorVertexDataSource< MeshType > >() );

   meshWriter.addDataSource( make_shared< IndexFaceDataSource  < MeshType > >() );
   meshWriter.addDataSource( make_shared< IndexVertexDataSource< MeshType > >() );

   meshWriterStatusfiltered.setFaceFilter( mesh::StatusFaceFilter< MeshType >( OpenMesh::Attributes::TAGGED ) );

   meshWriterStatusfiltered();
}


int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );
   if( args.size() != 2 )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " MESH_FILE" );

   const std::string & meshFile = args[1];

   test< mesh::TriangleMesh >( meshFile );
   test< mesh::FloatTriangleMesh >( meshFile );
   // test< mesh::PythonTriangleMesh >( meshFile ); // problematic due to different color type

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}