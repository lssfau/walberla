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

#include "geometry/containment_octree/ContainmentOctree.h"
#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/vtk/VTKMeshWriter.h"
#include "mesh_common/vtk/CommonDataSources.h"

#include <random>

#include <vector>
#include <string>

namespace walberla {
namespace mesh {

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );

   bool writeVtk = false;
   auto vtkArgIt = std::find( args.begin(), args.end(), "--vtk" );
   if(vtkArgIt != args.end())
   {
      writeVtk = true;
      args.erase( vtkArgIt );
   }

   if( args.size() != 2 )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " [--vtk] MESH_FILE" );

   const std::string & meshFile = args[1];

   auto mesh = make_shared<mesh::TriangleMesh>();
   mesh::readAndBroadcast( meshFile, *mesh);

   auto aabb = computeAABB( *mesh );

  //static const mesh::TriangleMesh::Point xAxis( 1, 0, 0 );
  //static const mesh::TriangleMesh::Point yAxis( 0, 1, 0 );
  //static const mesh::TriangleMesh::Point zAxis( 0, 0, 1 );
  
  mesh::TriangleMesh::Point r = mesh::toOpenMesh( ( aabb.minCorner() - aabb.maxCorner() ).getNormalized() );
  
  WALBERLA_LOG_DEVEL( "Start: r = " << r );
  
  mesh::TriangleMesh::Scalar sLengthOld = std::numeric_limits< mesh::TriangleMesh::Scalar >::max();
  for(int i = 0; i < 100; ++i)
  {
     mesh::TriangleMesh::Point s(0,0,0);
  
     for(auto vh : mesh->vertices())
     {
        const mesh::TriangleMesh::Point & x = mesh->point(vh);
        s = s + ( x | r ) * x;
     }
     const mesh::TriangleMesh::Scalar sLength = s.length();
     r = s / sLength;
  
     const mesh::TriangleMesh::Scalar eps = sLength - sLengthOld;
     WALBERLA_LOG_DEVEL( "Iteration:" << i << " r = " << r << " eps = " << eps );
     sLengthOld = sLength;
  }

  auto testVolume = aabb.getScaled( real_t(1.5) ); // AABB containing the test points
  
  auto triDist = make_shared< mesh::TriangleDistance<mesh::TriangleMesh> >( mesh );
  auto distanceOctree = make_shared< DistanceOctree<mesh::TriangleMesh> >( triDist );
  geometry::ContainmentOctree< DistanceOctree<mesh::TriangleMesh> > containmentOctree( distanceOctree );
  
  if( writeVtk )
     containmentOctree.writeVTKOutput( "containment_octree" );
  
  std::mt19937 rng;
  
  for( int i = 0; i < 10000; ++i )
  {
     const auto p = toOpenMesh( testVolume.randomPoint( rng ) );
     const bool isContained = containmentOctree.contains( p );
     const real_t distance = distanceOctree->sqSignedDistance( p );
     WALBERLA_CHECK_EQUAL( distance <= DistanceOctree<mesh::TriangleMesh>::Scalar(0), isContained, "Point " << std::setprecision(16) << p << " is " << ( isContained ? "inside" : "outside" ) << " but has distance " << distance );
  }
  
  return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}
