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
//! \file MeshOperationsTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/grid_generator/SCIterator.h"
#include "core/logging/Logging.h"
#include "core/math/AABB.h"
#include "core/math/KahanSummation.h"
#include "core/mpi/Environment.h"


#include "geometry/containment_octree/ContainmentOctree.h"
#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/MeshOperations.h"
#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/distance_octree/DistanceOctree.h"


namespace walberla {
namespace mesh {


template< typename MeshType >
void testCube()
{
   MeshType mesh;

   using Scalar = typename MeshType::Scalar;

   readAndBroadcast("cube.obj", mesh);

   auto aabb = computeAABB( mesh );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb.xMin(), Scalar(0) );
   WALBERLA_CHECK_IDENTICAL( aabb.yMin(), Scalar(0) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.zMin(), Scalar(0) );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb.xMax(), Scalar(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.yMax(), Scalar(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.zMax(), Scalar(1) );


   WALBERLA_CHECK_FLOAT_EQUAL( computeSurfaceArea( mesh ), Scalar(6) );
   WALBERLA_CHECK_FLOAT_EQUAL( computeVolume( mesh ), Scalar(1) );

   auto aabb_vertices = computeAABBForVertices( mesh, mesh.vertices_begin(), mesh.vertices_end() );
   auto aabb_faces    = computeAABBForFaces( mesh, mesh.faces_begin(), mesh.faces_end() );

   WALBERLA_CHECK_EQUAL( aabb, aabb_vertices );
   WALBERLA_CHECK_EQUAL( aabb, aabb_faces );

   translate( mesh, Vector3<Scalar>( Scalar(-0.5), Scalar(-0.5), Scalar(-0.5) ) );
   aabb = computeAABB( mesh );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb.xMin(), Scalar(-0.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.yMin(), Scalar(-0.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.zMin(), Scalar(-0.5) );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb.xMax(), Scalar(0.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.yMax(), Scalar(0.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.zMax(), Scalar(0.5) );


   WALBERLA_CHECK_FLOAT_EQUAL( computeSurfaceArea( mesh ), Scalar(6) );
   WALBERLA_CHECK_FLOAT_EQUAL( computeVolume( mesh ), Scalar(1) );
   auto centroid = computeCentroid( mesh );
   auto aabbCenter = aabb.center();
   WALBERLA_CHECK_FLOAT_EQUAL( centroid[0], aabbCenter[0] );
   WALBERLA_CHECK_FLOAT_EQUAL( centroid[1], aabbCenter[1] );
   WALBERLA_CHECK_FLOAT_EQUAL( centroid[2], aabbCenter[2] );

   Matrix3<Scalar> inertiaTensor = computeInertiaTensor(mesh);
   WALBERLA_CHECK_FLOAT_EQUAL( inertiaTensor(0,0), ( aabb.ySize() * aabb.ySize() + aabb.zSize() * aabb.zSize() ) / ( real_t(12) * aabb.volume() ) );
   WALBERLA_CHECK_FLOAT_EQUAL( inertiaTensor(1,1), ( aabb.xSize() * aabb.xSize() + aabb.zSize() * aabb.zSize() ) / ( real_t(12) * aabb.volume() ) );
   WALBERLA_CHECK_FLOAT_EQUAL( inertiaTensor(2,2), ( aabb.xSize() * aabb.xSize() + aabb.ySize() * aabb.ySize() ) / ( real_t(12) * aabb.volume() ) );

   scale( mesh, Vector3<Scalar>( Scalar(2), Scalar(3), Scalar(0.5) ) );
   aabb = computeAABB( mesh );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb.xMin(), Scalar(-1) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.yMin(), Scalar(-1.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.zMin(), Scalar(-0.25) );

   WALBERLA_CHECK_FLOAT_EQUAL( aabb.xMax(), Scalar(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.yMax(), Scalar(1.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( aabb.zMax(), Scalar(0.25) );


   WALBERLA_CHECK_FLOAT_EQUAL( computeSurfaceArea( mesh ), Scalar(2) * ( aabb.xSize() * aabb.ySize() + aabb.xSize() * aabb.zSize() + aabb.ySize() * aabb.zSize() ) );
   WALBERLA_CHECK_FLOAT_EQUAL( computeVolume( mesh ),  aabb.volume() );
   centroid = computeCentroid( mesh );
   aabbCenter = aabb.center();
   WALBERLA_CHECK_FLOAT_EQUAL( centroid[0], aabbCenter[0] );
   WALBERLA_CHECK_FLOAT_EQUAL( centroid[1], aabbCenter[1] );
   WALBERLA_CHECK_FLOAT_EQUAL( centroid[2], aabbCenter[2] );

   inertiaTensor = computeInertiaTensor(mesh);
   WALBERLA_CHECK_FLOAT_EQUAL( inertiaTensor(0,0), aabb.volume() * ( aabb.ySize() * aabb.ySize() + aabb.zSize() * aabb.zSize() ) / real_t(12) );
   WALBERLA_CHECK_FLOAT_EQUAL( inertiaTensor(1,1), aabb.volume() * ( aabb.xSize() * aabb.xSize() + aabb.zSize() * aabb.zSize() ) / real_t(12) );
   WALBERLA_CHECK_FLOAT_EQUAL( inertiaTensor(2,2), aabb.volume() * ( aabb.xSize() * aabb.xSize() + aabb.ySize() * aabb.ySize() ) / real_t(12) );
}

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   testCube< mesh::TriangleMesh >();
   testCube< mesh::FloatTriangleMesh >();
   testCube< mesh::PythonTriangleMesh >();

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}