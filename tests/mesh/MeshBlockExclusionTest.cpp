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

#include "mesh/blockforest/BlockExclusion.h"

#include "mesh_common/distance_octree/DistanceOctree.h"

#include <vector>
#include <string>

namespace walberla {
namespace mesh {

struct PointInAABB
{
   PointInAABB( const AABB & aabb ) : aabb_(aabb) {}

   inline bool operator()( const Vector3<real_t> & p ) { return aabb_.contains( p ); }

   AABB aabb_;
};

struct AnyPointInAABB
{
   using Points = std::vector<Vector3<real_t> >;

   AnyPointInAABB( const Points & points ) : points_(points) {}

   inline bool operator()( const AABB & aabb )
   {
     for( auto it = points_.begin(); it != points_.end(); ++it )
        if( aabb.contains( *it ) )
           return true;
     return false;
   }

   Points points_;
};


template< typename F, typename MeshType >
void test(const shared_ptr< DistanceOctree< MeshType > > & distanceOctree, const MeshType & mesh, const AABB & domainAABB, Vector3<uint_t> numBlocks)
{
   Vector3<real_t> blockSize(domainAABB.xSize() / real_c(numBlocks[0]),
      domainAABB.ySize() / real_c(numBlocks[1]),
      domainAABB.zSize() / real_c(numBlocks[2]));

   real_t maxError = blockSize.min() / real_t(10);

   SetupBlockForest setupBlockforest;
   setupBlockforest.addRootBlockExclusionFunction(F(distanceOctree, maxError));
   setupBlockforest.addWorkloadMemorySUIDAssignmentFunction(blockforest::uniformWorkloadAndMemoryAssignment);


   setupBlockforest.init(domainAABB, numBlocks[0], numBlocks[1], numBlocks[2], false, false, false);
   WALBERLA_LOG_DEVEL(setupBlockforest.toString());


   std::vector< Vector3<real_t> > vertexPositions;
   vertexPositions.reserve(mesh.n_vertices());
   for (auto vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt)
   {
      vertexPositions.push_back(toWalberla(mesh.point(*vIt)));
   }

   std::vector< const blockforest::SetupBlock* > setupBlocks;
   setupBlockforest.getBlocks(setupBlocks);

   // Check wether all vertices are located in allocated blocks
   std::vector< Vector3<real_t> > uncoveredVertices(vertexPositions);

   for (auto bIt = setupBlocks.begin(); bIt != setupBlocks.end(); ++bIt)
   {
      const AABB & aabb = (*bIt)->getAABB();

      uncoveredVertices.erase(std::remove_if(uncoveredVertices.begin(), uncoveredVertices.end(), PointInAABB(aabb)), uncoveredVertices.end());
   }

   WALBERLA_CHECK(uncoveredVertices.empty(), "Not all vertices of the mesh are located in allocated blocks!");

   //setupBlockforest.assignAllBlocksToRootProcess();
   //setupBlockforest.writeVTKOutput( "setupblockforest" );
}


template< typename MeshType >
void run( const std::string & meshFile, const uint_t numTotalBlocks )
{
   auto mesh = make_shared<MeshType>();
   mesh::readAndBroadcast( meshFile, *mesh);

   auto aabb = computeAABB( *mesh );

   auto domainAABB = aabb.getScaled( typename MeshType::Scalar(1.26) ); // AABB containing the test points

   auto triDist = make_shared< mesh::TriangleDistance<MeshType> >( mesh );
   auto  distanceOctree = make_shared< DistanceOctree<MeshType> >( triDist );

   Vector3<uint_t> numBlocks = math::getFactors3D( numTotalBlocks, domainAABB.sizes() );

   test< ExcludeMeshExterior< DistanceOctree< MeshType > > >( distanceOctree, *mesh, domainAABB, numBlocks );
   test< ExcludeMeshInterior< DistanceOctree< MeshType > > >( distanceOctree, *mesh, domainAABB, numBlocks );
}


int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );
   if( args.size() != 3 )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " MESH_FILE NUM_BLOCKS" );

   const std::string & meshFile       = args[1];
   const uint_t        numTotalBlocks = stringToNum< uint_t >( args[2] );

   run< mesh::TriangleMesh >( meshFile, numTotalBlocks );
   run< mesh::FloatTriangleMesh >( meshFile, numTotalBlocks );
   run< mesh::PythonTriangleMesh >( meshFile, numTotalBlocks );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}