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
//! \file MeshInitializationHelperTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/loadbalancing/StaticParMetis.h"

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/IntegerFactorization.h"
#include "core/mpi/Environment.h"
#include "core/stringToNum.h"

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/MeshIO.h"
#include "mesh/blockforest/BlockForestInitialization.h"
#include "mesh/blockforest/BlockWorkloadMemory.h"
#include "mesh/blockforest/BlockExclusion.h"
#include "mesh/blockforest/RefinementSelection.h"


#include "mesh_common/distance_octree/DistanceOctree.h"

#include <cmath>
#include <vector>
#include <string>

namespace walberla {
namespace mesh {

template< typename MeshType >
void testHelperFunctions( const std::string & meshFile, const uint_t numTotalBlocks )
{
   auto mesh = make_shared<MeshType>();
   mesh::readAndBroadcast( meshFile, *mesh);
   auto triDist = make_shared< mesh::TriangleDistance<MeshType> >( mesh );
   auto distanceOctree = make_shared< DistanceOctree< MeshType > >( triDist );

   const real_t meshVolume  = real_c( computeVolume( *mesh ) );
   const real_t blockVolume = meshVolume / real_c( numTotalBlocks );
   static const real_t cellsPersBlock = real_t(1000);
   const real_t cellVolume = blockVolume / cellsPersBlock;
   const real_t dx = std::pow( cellVolume, real_t(1) / real_t(3) );

   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with createStructuredBlockStorageInsideMesh with block size" );
   auto sbf0 = mesh::createStructuredBlockStorageInsideMesh( distanceOctree, dx, numTotalBlocks );
   
   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with createStructuredBlockStorageInsideMesh with block size" );
   Vector3<uint_t> blockSize( sbf0->getNumberOfXCells(), sbf0->getNumberOfYCells(), sbf0->getNumberOfZCells() );
   auto sbf1 = mesh::createStructuredBlockStorageInsideMesh( distanceOctree, dx, blockSize );

   auto exteriorAabb = computeAABB( *mesh ).getScaled( real_t(2) );

   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with createStructuredBlockStorageInsideMesh with block size" );
   auto sbf2 = mesh::createStructuredBlockStorageOutsideMesh( exteriorAabb, distanceOctree, dx, numTotalBlocks );

   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with createStructuredBlockStorageInsideMesh with block size" );
   blockSize = Vector3<uint_t>( sbf2->getNumberOfXCells(), sbf2->getNumberOfYCells(), sbf2->getNumberOfZCells() );
   auto sbf3 = mesh::createStructuredBlockStorageOutsideMesh( exteriorAabb, distanceOctree, dx, blockSize );
}

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );
   if( args.size() != 4 )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " MESH_FILE NUM_PROCESSES NUM_BLOCKS" );

   const std::string & meshFile       = args[1];
   const uint_t        numTotalBlocks = stringToNum< uint_t >( args[3] );

   testHelperFunctions< mesh::TriangleMesh >( meshFile, numTotalBlocks );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}