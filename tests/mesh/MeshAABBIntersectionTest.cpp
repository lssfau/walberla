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
//! \file MeshAABBIntersectionTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/AABB.h"
#include "core/mpi/Environment.h"
#include "core/Optional.h"

#include "mesh/MeshIO.h"
#include "mesh/MeshOperations.h"
#include "mesh/DistanceComputations.h"
#include "mesh/TriangleMeshes.h"

#include <random>

#include <vector>
#include <string>

namespace walberla {
namespace mesh {


template< typename MeshType >
void runTests( const uint_t numAABBs )
{
   auto mesh = make_shared<MeshType>();
   mesh::readAndBroadcast( "cube.obj", *mesh);

   auto meshAABB = computeAABB( *mesh ); // works since the mesh is a cube

   auto testVolume = meshAABB.getScaled( real_t(3) ); // AABB containing the test points

   TriangleDistance<MeshType> triDist( mesh );

   WALBERLA_CHECK( isIntersecting( triDist, meshAABB, real_t(0) ).value_or( false ) );

   std::mt19937 rng( uint32_t(42) );

   for(uint_t i = 0; i < numAABBs; ++i)
   {
      math::GenericAABB< typename MeshType::Scalar > testAABB( testVolume.randomPoint( rng ), testVolume.randomPoint( rng ) );

      const real_t maxErr = real_t(1e-2);

      walberla::optional< bool > result = isIntersecting( triDist, testAABB, maxErr );

      if ( result )
      {
         if(result.value())
         {
            WALBERLA_CHECK( meshAABB.intersects( testAABB ), "Box#: " << i );
         }
         else if(!result.value())
         {
            WALBERLA_CHECK( !meshAABB.intersects( testAABB ), "Box#: " << i );
         }
      }
   }
}

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );
   if( args.size() != 2 )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " NUM_AABBS" );

   const uint_t numAABBs = string_to_num< uint_t >( args[1] );

   runTests< mesh::TriangleMesh >( numAABBs );
   runTests< mesh::FloatTriangleMesh >( numAABBs );
   runTests< mesh::PythonTriangleMesh >( numAABBs );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}
