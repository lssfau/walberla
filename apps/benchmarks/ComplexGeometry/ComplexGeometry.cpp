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
//! \file ComplexGeometry.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/loadbalancing/StaticParMetis.h"

#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/logging/Logging.h"
#include "core/math/IntegerFactorization.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"

#include "geometry/InitBoundaryHandling.h"
#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "lbm/BlockForestEvaluation.h"
#include "lbm/PerformanceEvaluation.h"
#include "lbm/PerformanceLogger.h"
#include "lbm/boundary/factories/DefaultBoundaryHandling.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/communication/SparsePdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/refinement/TimeStep.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SplitPureSweep.h"
#include "lbm/vtk/VTKOutput.h"

#include "mesh/blockforest/BlockExclusion.h"
#include "mesh/blockforest/BlockForestInitialization.h"
#include "mesh/blockforest/BlockWorkloadMemory.h"
#include "mesh/blockforest/RefinementSelection.h"
#include "mesh/boundary/BoundaryInfo.h"
#include "mesh/boundary/BoundaryLocation.h"
#include "mesh/boundary/BoundaryLocationFunction.h"
#include "mesh/boundary/BoundarySetup.h"
#include "mesh/boundary/BoundaryUIDFaceDataSource.h"
#include "mesh/boundary/ColorToBoundaryMapper.h"

#include "stencil/D3Q19.h"

#include "timeloop/SweepTimeloop.h"

#include <cmath>
#include <string>
#include <vector>

#include "mesh_common/DistanceComputations.h"
#include "mesh_common/DistanceFunction.h"
#include "mesh_common/MatrixVectorOperations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/vtk/CommonDataSources.h"
#include "mesh_common/vtk/VTKMeshWriter.h"

namespace walberla {

template< typename MeshType >
void vertexToFaceColor( MeshType & mesh, const typename MeshType::Color & defaultColor )
{
   WALBERLA_CHECK( mesh.has_vertex_colors() );
   mesh.request_face_colors();

   for( auto faceIt = mesh.faces_begin(); faceIt != mesh.faces_end(); ++faceIt )
   {
      typename MeshType::Color vertexColor;

      bool useVertexColor = true;

      auto vertexIt = mesh.fv_iter( *faceIt );
      WALBERLA_ASSERT( vertexIt.is_valid() );
      
      vertexColor = mesh.color( *vertexIt );

      ++vertexIt;
      while( vertexIt.is_valid() && useVertexColor )
      {
         if( vertexColor != mesh.color( *vertexIt ) )
            useVertexColor = false;
         ++vertexIt;
      }

      mesh.set_color( *faceIt, useVertexColor ? vertexColor : defaultColor );
   }
}


int main( int argc, char * argv[] )
{
   Environment env( argc, argv );
   if( !env.config() )
   {
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << argv[0] << " INPUT_FILE" );
   }

   mpi::MPIManager::instance()->useWorldComm();

   const auto & config = *( env.config() );

   Config::BlockHandle configBlock = config.getOneBlock( "ComplexGeometry" );

   const std::string     meshFile            = configBlock.getParameter< std::string     >( "meshFile"                                 );
   const real_t          dx                  = configBlock.getParameter< real_t          >( "coarseDx"                                 );
   const real_t          omega               = configBlock.getParameter< real_t          >( "coarseOmega"                              );
   //const uint_t          blockPerProcess     = configBlock.getParameter< uint_t          >( "blocksPerProcess",    uint_t(6)           );
   const uint_t          timeSteps           = configBlock.getParameter< uint_t          >( "coarseTimeSteps"                          );
   const Vector3<real_t> bodyForce           = configBlock.getParameter< Vector3<real_t> >( "bodyForce"                                );
   //const bool            sparseCommunication = configBlock.getParameter< bool            >( "sparseCommunication", true                );
   const Vector3<real_t> domainBlowUp        = configBlock.getParameter< Vector3<real_t> >( "domainBlowUp",        Vector3<real_t>(6)  );
   const Vector3<uint_t> blockSize           = configBlock.getParameter< Vector3<uint_t> >( "blockSize",           Vector3<uint_t>(16) );
         uint_t          numLevels           = configBlock.getParameter< uint_t          >( "numLevels",           uint_t(2)           );

   numLevels = std::max( numLevels, uint_t(1) );

   //uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() );

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( meshFile );

   auto mesh = make_shared< mesh::TriangleMesh >();
   mesh->request_vertex_colors();
   WALBERLA_LOG_DEVEL_ON_ROOT( "Loading mesh" );
   mesh::readAndBroadcast( meshFile, *mesh);
   vertexToFaceColor( *mesh, mesh::TriangleMesh::Color(255,255,255) );

   WALBERLA_LOG_DEVEL_ON_ROOT( "Adding distance info to mesh" );
   auto triDist = make_shared< mesh::TriangleDistance<mesh::TriangleMesh> >( mesh );
   WALBERLA_LOG_DEVEL_ON_ROOT( "Building distance octree" );
   auto distanceOctree = make_shared< mesh::DistanceOctree<mesh::TriangleMesh> >( triDist );
   WALBERLA_LOG_DEVEL_ON_ROOT( "done. Octree has height " << distanceOctree->height() );

   distanceOctree->writeVTKOutput("distanceOctree");

   auto aabb = computeAABB( *mesh );
   aabb.scale(  domainBlowUp );

   mesh::ComplexGeometryStructuredBlockforestCreator bfc( aabb, Vector3<real_t>( dx ), mesh::makeExcludeMeshInterior( distanceOctree, dx ) );
   auto meshWorkloadMemory = mesh::makeMeshWorkloadMemory( distanceOctree, dx );
   meshWorkloadMemory.setInsideCellWorkload(1);
   meshWorkloadMemory.setOutsideCellWorkload(1);
   bfc.setWorkloadMemorySUIDAssignmentFunction( meshWorkloadMemory );
   bfc.setPeriodicity( Vector3<bool>(true) );
   bfc.setRefinementSelectionFunction( makeRefinementSelection( distanceOctree, numLevels - uint_t(1), dx, dx * real_t(1) ) );

   auto structuredBlockforest = bfc.createStructuredBlockForest( blockSize );

   typedef lbm::D3Q19<lbm::collision_model::SRT, false, lbm::force_model::SimpleConstant> LatticeModel_T;

   using flag_t = walberla::uint8_t;
   using FlagField_T = FlagField<flag_t>;
   using PdfField_T = lbm::PdfField<LatticeModel_T>;

   LatticeModel_T latticeModel{ lbm::collision_model::SRT( omega ), lbm::force_model::SimpleConstant( bodyForce ) };

   static const uint_t NUM_GHOSTLAYERS = 4;

   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( structuredBlockforest, "pdf field", latticeModel, Vector3<real_t>(0), real_t(1), NUM_GHOSTLAYERS, field::fzyx );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( structuredBlockforest, "flag field", NUM_GHOSTLAYERS );

   const FlagUID fluidFlagUID( "Fluid" );
   typedef lbm::DefaultBoundaryHandlingFactory< LatticeModel_T, FlagField_T > BHFactory;

   auto boundariesConfig = configBlock.getOneBlock( "Boundaries" );

   BlockDataID boundaryHandlingId = BHFactory::addBoundaryHandlingToStorage( structuredBlockforest, "boundary handling", flagFieldId, pdfFieldId, fluidFlagUID,
                                                                             boundariesConfig.getParameter< Vector3<real_t> >( "velocity0", Vector3<real_t>() ),
                                                                             boundariesConfig.getParameter< Vector3<real_t> >( "velocity1", Vector3<real_t>() ),
                                                                             boundariesConfig.getParameter< real_t > ( "pressure0", real_c( 1.0 ) ),
                                                                             boundariesConfig.getParameter< real_t > ( "pressure1", real_c( 1.001 ) ) );

   mesh::ColorToBoundaryMapper<mesh::TriangleMesh> colorToBoundryMapper(( mesh::BoundaryInfo( BHFactory::getNoSlipBoundaryUID() ) ));

  // colorToBoundryMapper.set( mesh::TriangleMesh::Color(255,0,0), mesh::BoundaryInfo( BHFactory::getPressure0BoundaryUID() ) );
  // colorToBoundryMapper.set( mesh::TriangleMesh::Color(0,0,255), mesh::BoundaryInfo( BHFactory::getPressure1BoundaryUID() ) );
  // colorToBoundryMapper.set( mesh::TriangleMesh::Color(255,255,255), mesh::BoundaryInfo( BHFactory::getNoSlipBoundaryUID() ) );

   auto boundaryLocations = colorToBoundryMapper.addBoundaryInfoToMesh( *mesh );

   mesh::VTKMeshWriter< mesh::TriangleMesh > meshWriter( mesh, "meshBoundaries", 1 );
   meshWriter.addDataSource( make_shared< mesh::BoundaryUIDFaceDataSource< mesh::TriangleMesh > >( boundaryLocations ) );
   meshWriter.addDataSource( make_shared< mesh::ColorFaceDataSource< mesh::TriangleMesh > >() );
   meshWriter.addDataSource( make_shared< mesh::ColorVertexDataSource< mesh::TriangleMesh > >() );
   meshWriter();

   WALBERLA_LOG_DEVEL_ON_ROOT( "Voxelizing mesh" );
   mesh::BoundarySetup boundarySetup( structuredBlockforest, makeMeshDistanceFunction( distanceOctree ), NUM_GHOSTLAYERS );
   //WALBERLA_LOG_DEVEL( "Writing Voxelisation" );
   //boundarySetup.writeVTKVoxelfile();
   WALBERLA_LOG_DEVEL_ON_ROOT( "Setting up fluid cells" );
   boundarySetup.setDomainCells<BHFactory::BoundaryHandling>( boundaryHandlingId, mesh::BoundarySetup::OUTSIDE );
   WALBERLA_LOG_DEVEL_ON_ROOT( "Setting up boundaries" );
   boundarySetup.setBoundaries<BHFactory::BoundaryHandling>( boundaryHandlingId, makeBoundaryLocationFunction( distanceOctree, boundaryLocations ), mesh::BoundarySetup::INSIDE );
   WALBERLA_LOG_DEVEL_ON_ROOT( "done" );

   lbm::BlockForestEvaluation<FlagField_T>( structuredBlockforest, flagFieldId, fluidFlagUID ).logInfoOnRoot();
   lbm::PerformanceLogger<FlagField_T> perfLogger( structuredBlockforest, flagFieldId, fluidFlagUID, 100 );

   SweepTimeloop timeloop( structuredBlockforest->getBlockStorage(), timeSteps );

   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, fluidFlagUID );
   auto refinementTimeStep = lbm::refinement::makeTimeStep<LatticeModel_T, BHFactory::BoundaryHandling > ( structuredBlockforest, sweep, pdfFieldId, boundaryHandlingId );
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( refinementTimeStep ), "Refinement time step" );   

   // log remaining time
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "remaining time logger" );

   // LBM stability check
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< PdfField_T, FlagField_T >( env.config(), structuredBlockforest, pdfFieldId,
                                                                                                             flagFieldId, fluidFlagUID ) ),
                                  "LBM stability check" );

   timeloop.addFuncAfterTimeStep( perfLogger, "Evaluator: performance logging" );


   // add VTK output to time loop
   lbm::VTKOutput< LatticeModel_T, FlagField_T >::addToTimeloop( timeloop, structuredBlockforest, env.config(), pdfFieldId, flagFieldId, fluidFlagUID );

   WcTimingPool timingPool;
   WALBERLA_LOG_INFO_ON_ROOT( "Starting timeloop" );
   timeloop.run( timingPool );
   WALBERLA_LOG_INFO_ON_ROOT( "Timeloop done" );
   timingPool.unifyRegisteredTimersAcrossProcesses();
   timingPool.logResultOnRoot( timing::REDUCE_TOTAL, true );

   return EXIT_SUCCESS;
}


} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::main( argc, argv );
}
