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
//! \file MeshDistanceCallbackMollerTrumboreTest.cpp
//! \ingroup mesh
//! \author Brendan Waters <brendan.waters@sydney.edu.au>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "field/AddToStorage.h"

#include "geometry/InitBoundaryHandling.h"

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/DistanceComputations.h"
#include "mesh_common/DistanceFunction.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/MeshBodyWallDistanceCallback.h"
#include "mesh_common/vtk/CommonDataSources.h"

#include "mesh/blockforest/BlockForestInitialization.h"
#include "mesh/boundary/BoundaryInfo.h"
#include "mesh/boundary/BoundaryLocation.h"
#include "mesh/boundary/BoundaryLocationFunction.h"
#include "mesh/boundary/BoundarySetup.h"
#include "mesh/boundary/BoundaryUIDFaceDataSource.h"

#include "stencil/D3Q27.h"

#include <vector>
#include <string>

namespace walberla {
namespace mesh {
   constexpr uint_t FieldGhostLayer{1};
   
   using flag_t      = walberla::uint32_t;
   using FlagField_T = FlagField< flag_t >;

   const FlagUID FluidFlagUID("Fluid");
   const FlagUID ObjectUID("Object");

   template< typename FlagField_T >
   void registerFlagFieldForMeshObject( const std::shared_ptr<blockforest::StructuredBlockForest>& blocks,
                                       const BlockDataID& flagFieldId, const FlagUID & UID ) 
   {
      for (auto& block : *blocks)
      {
            auto * flagField = block.getData<FlagField_T>(flagFieldId);
            if ( !flagField->flagExists(UID))
            flagField->registerFlag(UID);
      }
   }

   class TestMeshDistance
   {
   public:

      TestMeshDistance( std::function<double(const Cell &, const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)>&wallDistanceFct ):
                        elementInitialiser(wallDistanceFct){}

      template<typename FlagField_T>
      void testFromFlagField( const shared_ptr<StructuredBlockForest> & blocks, ConstBlockDataID flagFieldID,
                              FlagUID boundaryFlagUID, FlagUID domainFlagUID)
      {
         for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            testFromFlagField<FlagField_T>(blocks, &*blockIt, flagFieldID, boundaryFlagUID, domainFlagUID );
      }

      template<typename FlagField_T>
      void testFromFlagField( const shared_ptr<StructuredBlockForest> &blocks, IBlock * block, ConstBlockDataID flagFieldID,
                              FlagUID boundaryFlagUID, FlagUID domainFlagUID )
      {
         auto * flagField = block->getData< FlagField_T >(flagFieldID);
         
         if (!(flagField->flagExists(boundaryFlagUID) && flagField->flagExists(domainFlagUID)))
            return;

         auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
         auto domainFlag = flagField->getFlag(domainFlagUID);

         // Loop through the flagField
         for (auto it = flagField->beginWithGhostLayerXYZ(cell_idx_c(flagField->nrOfGhostLayers() - 1)); it != flagField->end(); ++it)
         {
            if (!isFlagSet(it, domainFlag))
                  continue;

            // Check neighbors in all directions
            for( auto dir = stencil::D3Q27::begin(); dir != stencil::D3Q27::end(); ++dir )
            {     
               if (isFlagSet(it.neighbor(dir.cx(), dir.cy(), dir.cz(), 0), boundaryFlag))
               {
                  const real_t q = elementInitialiser(   Cell(it.x(), it.y(), it.z()), 
                                                         Cell(it.x() +  dir.cx(), it.y() + + dir.cy(), it.z() +  dir.cz()), 
                                                         blocks, *block);

                  WALBERLA_CHECK_GREATER_EQUAL(q, real_c(0.0));
                  WALBERLA_CHECK_LESS_EQUAL(q, real_c(1.0));
               }
            }
         }
      }
   private:
      const std::function<double(const Cell &, const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)> elementInitialiser; 
   };


   template<typename MeshType>
   void test( const std::string & meshFile )
   {
      // build grid
      auto mesh = make_shared<MeshType>();
      mesh::readAndBroadcast( meshFile, *mesh);
      
      const Vector3<real_t> mesh_scaling_factor { 20 };
      mesh::scale(*mesh, mesh_scaling_factor );

      auto triDist = make_shared< mesh::TriangleDistance<MeshType> >( mesh );
      auto distanceOctree = make_shared< mesh::DistanceOctree< MeshType > >( triDist );

      std::vector< real_t > test_dx = { 0.33, 1.0, 1.5 };
      std::vector< real_t > mesh_shift = { 0.0, 0.33, 0.5, 1.0 };

      for (const auto dx : test_dx )
      {
         for (const auto shift : mesh_shift )
         {
            const Vector3<real_t> translation_shift { shift * dx };
            AABB aabb = computeAABB( *mesh );
            aabb.scale( 1.2 ); // AABB containing the test points
            aabb.translate(translation_shift); // Shift test points relative to mesh to get arbitrary distances.

            // build blockforest
            mesh::ComplexGeometryStructuredBlockforestCreator bfc( aabb, Vector3< real_t >(dx),
                                                                  mesh::makeExcludeMeshInterior(distanceOctree, dx));

            const shared_ptr<StructuredBlockForest> blocks = bfc.createStructuredBlockForest( Vector3<uint_t>(32,32,32) );        

            // set flags
            const BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", FieldGhostLayer);

            registerFlagFieldForMeshObject<FlagField_T>(blocks, flagFieldId, ObjectUID);

            mesh::BoundarySetup boundarySetup( blocks, makeMeshDistanceFunction( distanceOctree ), FieldGhostLayer );
            boundarySetup.setFlag<FlagField_T>(flagFieldId, ObjectUID, mesh::BoundarySetup::INSIDE);

            // Set remaining cells to fluid
            geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, FluidFlagUID);

            mesh::MeshBodyWallDistance< MeshType > meshWallDistanceCallback( distanceOctree );
            std::function< real_t(const Cell&, const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
               meshWallDistanceFunctor = meshWallDistanceCallback;

            TestMeshDistance testWrapper(meshWallDistanceFunctor);
            testWrapper.testFromFlagField< FlagField_T >(blocks, flagFieldId, ObjectUID, FluidFlagUID);
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
         WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " mesh.obj" );

      const std::string & meshFile = args[1];

      test<mesh::TriangleMesh>( meshFile );
      test<mesh::PythonTriangleMesh>( meshFile );
      

      return EXIT_SUCCESS;
   }


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}