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
#include "core/mpi/Environment.h"

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/MeshIO.h"

#include <random>

#include <vector>
#include <string>

namespace walberla {
namespace mesh {

template<typename MeshType>
void test( const std::string & meshFile )
{
   auto mesh = make_shared<MeshType>();
   mesh::readAndBroadcast( meshFile, *mesh);

   auto aabb = computeAABB( *mesh );

   auto testVolume = aabb.getScaled( typename MeshType::Scalar(1.2) ); // AABB containing the test points

   auto triDist = make_shared< mesh::TriangleDistance<MeshType> >( mesh );
   DistanceOctree<MeshType> distanceOctree( triDist );

   //distanceOctree.writeVTKOutput( "distance_octree" );

   std::mt19937 rng;

   for( int i = 0; i < 1000; ++i )
   {
      auto p = testVolume.randomPoint( rng );

      auto d0 = triDist->sqSignedDistance( toOpenMesh(p) );
      auto d1 = distanceOctree.sqSignedDistance( toOpenMesh(p) );
      auto d2 = distanceOctree.sqDistance( toOpenMesh(p) );

      WALBERLA_CHECK_FLOAT_EQUAL( d0, d1 );
      WALBERLA_CHECK_FLOAT_EQUAL( std::fabs( d1 ), d2 );
   }
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

   test<mesh::TriangleMesh>( meshFile );
   test<mesh::FloatTriangleMesh>( meshFile );
   test<mesh::PythonTriangleMesh>( meshFile );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}