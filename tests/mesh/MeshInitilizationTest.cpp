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

void commInXDirection( std::vector< std::pair< const SetupBlock*, const SetupBlock* > > blocks, std::vector< int64_t > & weights )
{
   WALBERLA_ASSERT_EQUAL( blocks.size(), weights.size() );

   for(size_t i = 0; i < blocks.size(); ++i)
   {
      auto aabb1 = blocks[i].first->getAABB();
      auto aabb2 = blocks[i].second->getAABB();

      auto centerdiff = aabb2.center() - aabb1.center();

      static const real_t eps = real_t(1e-8);

      if( std::fabs( centerdiff[0] ) > eps && std::fabs( centerdiff[1] ) < eps && std::fabs( centerdiff[2] ) < eps  )
         weights[i] = int64_t(1000);
      else
         weights[i] = int64_t(1);
   }
}

template< typename MeshType >
void test( const std::string & meshFile, const uint_t numProcesses, const uint_t numTotalBlocks )
{
   auto mesh = make_shared<MeshType>();
   mesh::readAndBroadcast( meshFile, *mesh);

   auto aabb = computeAABB( *mesh );

   auto domainAABB = aabb.getScaled( typename MeshType::Scalar(3) );

   auto triDist = make_shared< mesh::TriangleDistance<MeshType> >( mesh );
   auto distanceOctree = make_shared< DistanceOctree< MeshType > >( triDist );

   const real_t meshVolume  = real_c( computeVolume( *mesh ) );
   const real_t blockVolume = meshVolume / real_c( numTotalBlocks );
   static const real_t cellsPersBlock = real_t(1000);
   const real_t cellVolume = blockVolume / cellsPersBlock;
   const Vector3<real_t> cellSize( std::pow( cellVolume, real_t(1) / real_t(3) ) );

   ComplexGeometryStructuredBlockforestCreator bfc( domainAABB, cellSize, makeExcludeMeshInterior( distanceOctree, cellSize.min() ) );
   auto wl = mesh::makeMeshWorkloadMemory( distanceOctree, cellSize );
   wl.setInsideCellWorkload(0);
   wl.setOutsideCellWorkload(1);
   wl.setForceZeroMemoryOnZeroWorkload(true);
   bfc.setWorkloadMemorySUIDAssignmentFunction( wl );
   bfc.setRefinementSelectionFunction( makeRefinementSelection( distanceOctree, 5, cellSize[0], cellSize[0] * real_t(5) ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with StaticLevelwiseCurveBalanceWeighted Partitioner" );
   bfc.setTargetProcessAssignmentFunction( blockforest::StaticLevelwiseCurveBalanceWeighted() );
   auto sbf_default = bfc.createSetupBlockForest( Vector3<uint_t>(64,64,64), numProcesses );
   //sbf_default->writeVTKOutput("sbf_default");
   WALBERLA_LOG_INFO_ON_ROOT( sbf_default->toString() );

   return;
#ifdef WALBERLA_BUILD_WITH_PARMETIS

   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with ParMetis (PART_KWAY, no commweights)" );
   bfc.setTargetProcessAssignmentFunction( blockforest::StaticLevelwiseParMetis( blockforest::StaticLevelwiseParMetis::PARMETIS_PART_KWAY ) );
   auto sbf = bfc.createSetupBlockForest( numTotalBlocks, numProcesses );
   //sbf->writeVTKOutput("sbf");
   WALBERLA_LOG_INFO_ON_ROOT( sbf->toString() );


   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with ParMetis (PART_KWAY, commweights)" );
   bfc.setTargetProcessAssignmentFunction( blockforest::StaticLevelwiseParMetis( commInXDirection, blockforest::StaticLevelwiseParMetis::PARMETIS_PART_KWAY ) );
   auto sbf_edge = bfc.createSetupBlockForest( numTotalBlocks, numProcesses );
   //sbf_edge->writeVTKOutput("sbf_edge");
   WALBERLA_LOG_INFO_ON_ROOT( sbf_edge->toString() );

   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with ParMetis (PART_GEOM_KWAY, no commweights)" );
   bfc.setTargetProcessAssignmentFunction( blockforest::StaticLevelwiseParMetis( blockforest::StaticLevelwiseParMetis::PARMETIS_PART_GEOM_KWAY ) );
   auto sbf_geom = bfc.createSetupBlockForest( numTotalBlocks, numProcesses );
   //sbf_geom->writeVTKOutput("sbf_geom");
   WALBERLA_LOG_INFO_ON_ROOT( sbf_geom->toString() );


   WALBERLA_LOG_INFO_ON_ROOT( "Creating SBF with ParMetis (PART_GEOM_KWAY, commweights)" );
   bfc.setTargetProcessAssignmentFunction( blockforest::StaticLevelwiseParMetis( commInXDirection, blockforest::StaticLevelwiseParMetis::PARMETIS_PART_GEOM_KWAY ) );
   auto sbf_geom_edge = bfc.createSetupBlockForest( numTotalBlocks, numProcesses );
   //sbf_geom_edge->writeVTKOutput("sbf_geom_edge");
   WALBERLA_LOG_INFO_ON_ROOT( sbf_geom_edge->toString() );

#endif
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
   const uint_t        numProcesses   = stringToNum< uint_t >( args[2] );
   const uint_t        numTotalBlocks = stringToNum< uint_t >( args[3] );

   test< mesh::TriangleMesh >( meshFile, numProcesses, numTotalBlocks );
   //test< mesh::FloatTriangleMesh >( meshFile, numProcesses, numTotalBlocks );
   //test< mesh::PythonTriangleMesh >( meshFile, numProcesses, numTotalBlocks );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}