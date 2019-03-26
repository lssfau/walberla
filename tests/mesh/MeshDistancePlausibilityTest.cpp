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

#include "blockforest/BlockForest.h"
#include "blockforest/Initialization.h"
#include "blockforest/StructuredBlockForest.h"

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "field/AddToStorage.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh/TriangleMeshes.h"
#include "mesh/MeshOperations.h"
#include "mesh/DistanceComputations.h"
#include "mesh/distance_octree/DistanceOctree.h"
#include "mesh/MeshIO.h"
#include "mesh/vtk/VTKMeshWriter.h"
#include "mesh/vtk/CommonDataSources.h"
#include "mesh/vtk/CommonFilters.h"

#include "stencil/D3Q27.h"

#include <vector>
#include <string>

namespace walberla {
namespace mesh {

typedef field::GhostLayerField< real_t                        , 1 > DistanceField;
typedef field::GhostLayerField< uint8_t                       , 1 > ErrorMarkerField;
typedef field::GhostLayerField< mesh::TriangleMesh::FaceHandle, 1 > FaceHandleField;

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   bool writeVtk = false;
   std::vector<std::string> args( argv, argv + argc );
   if( args.size() < 3LU || args.size() > 4LU )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " [--vtk] MESH_FILE RESOLUTION" );

   auto vtkArgIt = std::find( args.begin(), args.end(), "--vtk" );
   if( vtkArgIt != args.end() )
   {
      writeVtk = true;
      args.erase( vtkArgIt );
   }
   const std::string & meshFile = args[1];
   real_t dx = string_to_num<real_t>( args[2] );

   auto mesh = make_shared<mesh::TriangleMesh>();
   mesh::readAndBroadcast( meshFile, *mesh);
   mesh->request_face_status();

   auto meshAabb = computeAABB( *mesh );

   auto domainAABB = meshAabb.getScaled( real_t(1.5) ); // AABB containing the test points

   WALBERLA_LOG_INFO_ON_ROOT( "Preparing distance Information..." );
   auto triDist = make_shared< mesh::TriangleDistance<mesh::TriangleMesh> >( mesh );



   DistanceOctree<mesh::TriangleMesh> distanceOctree( triDist );

   if( writeVtk )
      distanceOctree.writeVTKOutput( "vtk_out/distance_octree" );

   Vector3<uint_t> cells( numeric_cast<uint_t>( std::ceil( domainAABB.xSize() / dx ) ),
                          numeric_cast<uint_t>( std::ceil( domainAABB.ySize() / dx ) ),
                          numeric_cast<uint_t>( std::ceil( domainAABB.zSize() / dx ) ) );

   domainAABB = AABB( real_t(0), real_t(0), real_t(0), real_c( cells[0] ) * dx, real_c( cells[1] ) * dx, real_c( cells[2] ) * dx );
   domainAABB.translate( meshAabb.center() - domainAABB.center() );

   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid( domainAABB,
                                                                                     uint_t(1), uint_t(1), uint_t(1),
                                                                                     cells[0], cells[1], cells[2],
                                                                                     uint_t(1), uint_t(1), uint_t(1) );

   WALBERLA_CHECK_FLOAT_EQUAL( dx,           blocks->dx() );
   WALBERLA_CHECK_FLOAT_EQUAL( blocks->dx(), blocks->dy() );
   WALBERLA_CHECK_FLOAT_EQUAL( blocks->dx(), blocks->dz() );

   BlockDataID distanceFieldId    = field::addToStorage< DistanceField    >( blocks, "DistanceField"   , real_t(0)  );
   BlockDataID errorMarkerFieldId = field::addToStorage< ErrorMarkerField >( blocks, "ErrorMarkerField", uint8_t(0) );
   BlockDataID faceHandleFieldId  = field::addToStorage< FaceHandleField  >( blocks, "FaceHandleField" );



   WALBERLA_LOG_INFO_ON_ROOT( "Computing distance field of size " << cells << " (" << cells[0] * cells[1] * cells[2] << " cells)..." );

   for( IBlock & block : *blocks )
   {
      auto * distanceField   = block.getData<DistanceField  >( distanceFieldId   );
      auto * faceHandleField = block.getData<FaceHandleField>( faceHandleFieldId );
      for(cell_idx_t z = -1; z < cell_idx_c( distanceField->zSize() ) + cell_idx_t(1); ++z )
      {
         for( cell_idx_t y = -1; y < cell_idx_c( distanceField->ySize() ) + cell_idx_t(1); ++y )
         {
            for( cell_idx_t x = -1; x < cell_idx_c( distanceField->xSize() ) + cell_idx_t(1); ++x )
            {
               Cell c( x, y, z );
               Vector3<real_t> p = blocks->getCellCenter( c );
               real_t d = distanceOctree.sqSignedDistance( toOpenMesh( p ), faceHandleField->get( c ) );
               distanceField->get( c ) = d < real_t(0) ? -std::sqrt( -d ) : std::sqrt(d);
            }
         }
      WALBERLA_LOG_INFO( "Slice " << z + 1 << "/" << cell_idx_c( distanceField->zSize() ) + cell_idx_t(2) );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Checking distance field..." );

   bool errorFound = false;

   // Use triangle inequality to perform a plausibility check on the distance field
   for( auto blockIt = blocks->begin(); blockIt != blocks->end() && !errorFound; ++blockIt )
   {
      IBlock & block = *blockIt;

      auto * distanceField    = block.getData<DistanceField   >( distanceFieldId    );
      auto * errorMarkerField = block.getData<ErrorMarkerField>( errorMarkerFieldId );
      auto * faceHandleField  = block.getData<FaceHandleField >( faceHandleFieldId  );

      for(uint_t z = 0; z < distanceField->zSize() && !errorFound; ++z)
         for( uint_t y = 0; y < distanceField->ySize() && !errorFound; ++y )
            for( uint_t x = 0; x < distanceField->xSize() && !errorFound; ++x )
            {
               const Cell c( x, y, z );
               const real_t d = distanceField->get( c );

               for( auto dirIt = stencil::D3Q27::beginNoCenter(); dirIt != stencil::D3Q27::end() && !errorFound; ++dirIt )
               {
                  const Cell nc = c + *dirIt;
                  const real_t upperLimit = d + dirIt.length() * dx + real_comparison::Epsilon<real_t>::value;
                  const real_t lowerLimit = d - dirIt.length() * dx - real_comparison::Epsilon<real_t>::value;

                  const real_t nd = distanceField->get( nc );

                  if( nd > upperLimit || nd < lowerLimit )
                  {
                     errorMarkerField->get( c ) = 1;
                     errorMarkerField->get( nc ) = 1;

                     mesh->status( faceHandleField->get( c ) ).set_tagged( true );
                     mesh->status( faceHandleField->get( nc ) ).set_tagged( true );

                     std::ostringstream oss;
                     oss << "Distance at cell " << c  << " is " <<  d << "\n"
                         << "Distance at cell " << nc << " is " << nd << "\n"
                         << "Difference is " << std::fabs( d - nd ) << " but it should be <= " << dirIt.length() * dx << "\n"
                         << "Closest Face to " << c << ":\n";
                     triDist->triangleToStream( faceHandleField->get( c ), oss );
                     oss << "\n\n"
                         << "Closest Face to " << nc << ":\n";
                     triDist->triangleToStream( faceHandleField->get( nc ), oss );
                     oss << "\n\n"
                         << "If you enabled vtk output (--vtk) the written files will contain markers at the affected cells and faces!";
                     WALBERLA_LOG_WARNING( oss.str() );
                     errorFound = true;
                  }
               }
            }
   }

   if( writeVtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing mesh VTK output..." );

      VTKMeshWriter< mesh::TriangleMesh > meshWriter( mesh, "mesh", 1 );
      meshWriter.addDataSource( make_shared< NormalsVertexDataSource< mesh::TriangleMesh > >() );
      meshWriter.addDataSource( make_shared< NormalsFaceDataSource< mesh::TriangleMesh > >() );
      meshWriter.addDataSource( make_shared< StatusBitFaceDataSource< mesh::TriangleMesh > >( OpenMesh::Attributes::TAGGED, "error marker" ) );
      meshWriter();

      WALBERLA_LOG_INFO_ON_ROOT( "Writing field VTK output..." );

      auto vtkOutput = vtk::createVTKOutput_BlockData( blocks, "plausibility test" );

      vtkOutput->addCellDataWriter( make_shared< field::VTKWriter< DistanceField, float  > >( distanceFieldId   , "distance field" ) );
      vtkOutput->addCellDataWriter( make_shared< field::VTKWriter< ErrorMarkerField      > >( errorMarkerFieldId, "error marker"   ) );

      writeFiles( vtkOutput, true )();
   }


   return errorFound ? EXIT_FAILURE : EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}