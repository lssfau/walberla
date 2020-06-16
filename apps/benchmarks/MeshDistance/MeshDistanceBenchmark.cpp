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
#include "core/stringToNum.h"
#include "core/timing/Timer.h"

#include "geometry/containment_octree/ContainmentOctree.h"
#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/MeshIO.h"

#include <random>

#include <algorithm>
#include <vector>
#include <string>

namespace mesh_distance_benchmark {

using namespace walberla;

template< typename MeshType >
void runBenchmark( const std::string & meshFile, const uint_t numPoints, const uint_t numRepetitions, const bool testBruteForce )
{
   auto mesh = make_shared<MeshType>();;
   mesh::readAndBroadcast( meshFile, *mesh);

   auto aabb = computeAABB( *mesh );

   WALBERLA_LOG_INFO( "Mesh has " << mesh->n_faces() << " triangles" );
   WALBERLA_LOG_INFO( "Testing " << numPoints << " Points with " << numRepetitions << " repetitions" )

   auto testVolume = aabb.getScaled( typename MeshType::Point::value_type(1.2) ); // AABB containing the test points

   WcTimer timer;
   WALBERLA_LOG_INFO( "Preparing mesh for distance computations..." );
   timer.start();
   auto triDist = make_shared< mesh::TriangleDistance<MeshType> >( mesh );
   timer.end();
   WALBERLA_LOG_INFO( "Mesh preparation took " << timer.last() << "s" );

   std::mt19937 rng;

   std::vector< typename MeshType::Point > points( numPoints );
   for(auto it = points.begin(); it != points.end(); ++it)
   {
      *it = mesh::toOpenMesh( testVolume.randomPoint( rng ) );
   }

   if( testBruteForce )
   {
      WALBERLA_LOG_INFO( "Testing brute force distance computation..." );
      timer.reset();
      for( uint_t i = 0; i < numRepetitions; ++i )
      {
         timer.start();
         for(auto it = points.begin(); it != points.end(); ++it)
         {
            triDist->sqSignedDistance( *it );
         }
         timer.end();
      }
      WALBERLA_LOG_INFO( "Brute force distance computation took " << timer.min() << "s" );
      WALBERLA_LOG_INFO( real_c(points.size()) / timer.min() << " points / s" );
   }


   WALBERLA_LOG_INFO( "Building distance octree..." );
   timer.reset();
   timer.start();
   mesh::DistanceOctree<MeshType> distanceOctree( triDist );
   timer.end();
   WALBERLA_LOG_INFO( "Building octree took " << timer.last() << "s" );
   WALBERLA_LOG_INFO( "Distance octree has height " << distanceOctree.height() );

   WALBERLA_LOG_INFO( "Testing octree distance computation..." );
   timer.reset();
   for( uint_t i = 0; i < numRepetitions; ++i )
   {
      timer.start();
      for(auto it = points.begin(); it != points.end(); ++it)
      {
         distanceOctree.sqSignedDistance( *it );
      }
      timer.end();
   }
   WALBERLA_LOG_INFO( "Octree distance computation took " << timer.min() << "s" );
   WALBERLA_LOG_INFO( real_c( points.size() ) / timer.min() << " points / s" );


   WALBERLA_LOG_INFO( "Building containment octree..." );
   timer.reset();
   timer.start();
   geometry::ContainmentOctree< mesh::DistanceOctree<MeshType> > containmentOctree( make_shared<mesh::DistanceOctree<MeshType>>( distanceOctree ) );
   timer.end();
   WALBERLA_LOG_INFO( "Building octree took " << timer.last() << "s" );
   WALBERLA_LOG_INFO( "Containment octree has height " << containmentOctree.height() );

   WALBERLA_LOG_INFO( "Testing octree containment computation..." );
   timer.reset();
   for( uint_t i = 0; i < numRepetitions; ++i )
   {
      timer.start();
      for(auto it = points.begin(); it != points.end(); ++it)
      {
         containmentOctree.contains( *it );
      }
      timer.end();
   }
   WALBERLA_LOG_INFO( "Octree containment computation took " << timer.min() << "s" );
   WALBERLA_LOG_INFO( real_c( points.size() ) / timer.min() << " points / s" );
}

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );

   bool testBruteForce = true;
   auto argIt = std::find( args.begin(), args.end(), std::string("--no-brute-force") );
   if( argIt != args.end() )
   {
      testBruteForce = false;
      args.erase( argIt );
   }

   bool forceFloat = false;
   argIt = std::find( args.begin(), args.end(), std::string("--force-float") );
   if( argIt != args.end() )
   {
      forceFloat = true;
      args.erase( argIt );
   }
   
   if( args.size() != 4 )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " [--no-brute-force] [--force-float] MESH_FILE NUM_POINTS NUM_REPETITIONS" );

   const std::string & meshFile = args[1];
   const uint_t numPoints       = stringToNum<uint_t>( args[2] );
   const uint_t numRepetitions  = stringToNum<uint_t>( args[3] );

   if(forceFloat)
   {
      runBenchmark< mesh::FloatTriangleMesh >( meshFile, numPoints, numRepetitions, testBruteForce );
   }
   else
   {
      runBenchmark< mesh::TriangleMesh >( meshFile, numPoints, numRepetitions, testBruteForce );
   }

   return EXIT_SUCCESS;
}


} // namespace mesh_distance_benchmark

int main( int argc, char * argv[] )
{
   return mesh_distance_benchmark::main( argc, argv );
}