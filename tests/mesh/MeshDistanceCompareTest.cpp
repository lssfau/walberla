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
//! \file MeshDistanceCompareTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "geometry/containment_octree/ContainmentOctree.h"
#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/MeshIO.h"

#include <random>

#include <vector>
#include <string>

namespace walberla {
namespace mesh {

template< typename MeshType >
void testAABBDistance( const Vector3<real_t> & translationVector = Vector3<real_t>() )
{
   auto mesh = make_shared<MeshType>();
   mesh::readAndBroadcast( "cube.obj", *mesh);

   translate( *mesh, translationVector );

   auto aabb = computeAABB( *mesh ); // works since the mesh is a cube

   auto testVolume = aabb.getScaled( real_t(2) ); // AABB containing the test points

   TriangleDistance<MeshType> triDist( mesh );

   std::mt19937 rng;

   for( int i = 0; i < 10000; ++i )
   {
      auto p = testVolume.randomPoint( rng );
      WALBERLA_CHECK_FLOAT_EQUAL( triDist.sqSignedDistance( toOpenMesh(p) ), aabb.sqSignedDistance( p ) );
   }

}

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   testAABBDistance<mesh::TriangleMesh>();
   testAABBDistance<mesh::TriangleMesh>( Vector3<real_t>( real_t(-0.5), real_t(-0.5), real_t(-0.5) ) );

   testAABBDistance<mesh::FloatTriangleMesh>();
   testAABBDistance<mesh::FloatTriangleMesh>( Vector3<real_t>( real_t(-0.5), real_t(-0.5), real_t(-0.5) ) );

   testAABBDistance<mesh::PythonTriangleMesh>();
   testAABBDistance<mesh::PythonTriangleMesh>( Vector3<real_t>( real_t(-0.5), real_t(-0.5), real_t(-0.5) ) );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}