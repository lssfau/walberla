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
#include "core/logging/Logging.h"
#include "core/stringToNum.h"

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"

#include <random>

#include <algorithm>
#include <vector>

namespace walberla {
namespace mesh {


int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );
   if( args.size() != 4 )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " MESH_FILE NUM_BOXES NUM_POINTS_TESTED_PER_BOX" );

   const std::string & meshFile = args[1];
   const uint_t numBoxes = stringToNum<uint_t>( args[2] );
   const uint_t numPointsTestedPerBox = stringToNum<uint_t>( args[3] );

   auto mesh = make_shared<TriangleMesh>();
   mesh::readAndBroadcast( meshFile, *mesh );
   mesh::TriangleDistance<TriangleMesh> triDist( mesh );

   auto aabb = computeAABB( *mesh );

   std::mt19937 rng;

   for( uint_t i = 0; i < numBoxes; ++i )
   {
      math::GenericAABB< mesh::TriangleDistance<TriangleMesh>::Scalar > subAabb(  aabb.randomPoint( rng ), aabb.randomPoint( rng )  );

      std::vector<TriangleMesh::FaceHandle> subAabbFaces;

      triDist.filterTrianglesForAABB( subAabb, mesh->faces_begin(), mesh->faces_end(), std::back_inserter( subAabbFaces ) );

      WALBERLA_CHECK( !subAabbFaces.empty() );
      WALBERLA_CHECK_EQUAL( std::adjacent_find( subAabbFaces.begin(), subAabbFaces.end() ), subAabbFaces.end(), "Duplicate faces found!" );

      WALBERLA_LOG_INFO( "Box " << i + 1 << "/" << numBoxes << ": Volume ratio: " << subAabb.volume() / aabb.volume() << " Faces ratio: " << real_c( subAabbFaces.size() ) / real_c( mesh->n_faces() ) );

      for(uint_t j = 0; j < numPointsTestedPerBox; ++j)
      {
         auto p = toOpenMesh( subAabb.randomPoint( rng ) );
         
         TriangleMesh::FaceHandle fh0, fh1;
         auto d0 = triDist.sqSignedDistance( p, fh0 );
         auto d1 = triDist.sqSignedDistance( subAabbFaces.begin(), subAabbFaces.end(), p, fh1 );

         if( !floatIsEqual( d0, d1 ) )
         {
            const bool fh0InSubAabbFaces = std::find( subAabbFaces.begin(), subAabbFaces.end(), fh0 ) != subAabbFaces.end();
            std::ostringstream msg;
            msg << "Mesh AABB:           " << aabb                << "\n";
            msg << "Sub AABB:            " << subAabb             << "\n";
            msg << "#Faces:              " << mesh->n_faces()     << "\n";
            msg << "#Faces Sub AABB:     " << subAabbFaces.size() << "\n";
            msg << "p:                   " << p                   << "\n";
            msg << "d0:                  " << d0                  << "\n";
            msg << "d1:                  " << d1                  << "\n";
            msg << "fh0:                 " << fh0                 << "\n";
            msg << "fh1:                 " << fh1                 << "\n";
            msg << "fh0 in subAABBFaces? " << fh0InSubAabbFaces   << "\n";

            WALBERLA_ABORT( msg.str() );
         }
      }

   }

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}