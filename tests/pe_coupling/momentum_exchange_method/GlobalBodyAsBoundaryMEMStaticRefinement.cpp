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
//! \file GlobalBodyAsBoundaryMEMStaticRefinement.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include <blockforest/loadbalancing/StaticCurve.h>
#include <blockforest/SetupBlockForest.h>

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/refinement/all.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "pe/basic.h"
#include "pe/Types.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <functional>

namespace global_body_as_boundary_mem_static_refinement
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

//////////////
// TYPEDEFS //
//////////////

// PDF field, flag field & body field
typedef lbm::D3Q19< lbm::collision_model::TRT >  LatticeModel_T;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

const uint_t FieldGhostLayers = 4;

// boundary handling
typedef pe_coupling::SimpleBB< LatticeModel_T, FlagField_T > MO_SBB_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, MO_SBB_T > BoundaryHandling_T;

using BodyTypeTuple = std::tuple<pe::Plane>;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID MO_SBB_Flag( "moving obstacle SBB" );

/////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void refinementSelection( SetupBlockForest& forest, uint_t levels, AABB refinementBox )
{
   real_t dx = real_t(1); // dx on finest level
   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      uint_t blockLevel = block->getLevel();
      uint_t levelScalingFactor = ( uint_t(1) << (levels - uint_t(1) - blockLevel) );
      real_t dxOnLevel = dx * real_c(levelScalingFactor);
      AABB blockAABB = block->getAABB();

      // extend block AABB by ghostlayers
      AABB extendedBlockAABB = blockAABB.getExtended( dxOnLevel * real_c(FieldGhostLayers) );

      if( extendedBlockAABB.intersects( refinementBox ) )
         if( blockLevel < ( levels - uint_t(1) ) )
            block->setMarker( true );
   }
}

static void workloadAndMemoryAssignment( SetupBlockForest& forest )
{
   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      block->setWorkload( numeric_cast< workload_t >( uint_t(1) << block->getLevel() ) );
      block->setMemory( numeric_cast< memory_t >(1) );
   }
}

static shared_ptr< StructuredBlockForest > createBlockStructure( const AABB & domainAABB, Vector3<uint_t> blockSizeInCells,
                                                                 uint_t numberOfLevels, bool keepGlobalBlockInformation = false )
{
   SetupBlockForest sforest;

   Vector3<uint_t> numberOfFineBlocksPerDirection( uint_c(domainAABB.size(0)) / blockSizeInCells[0],
                                                   uint_c(domainAABB.size(1)) / blockSizeInCells[1],
                                                   uint_c(domainAABB.size(2)) / blockSizeInCells[2] );

   for(uint_t i = 0; i < 3; ++i )
   {
      WALBERLA_CHECK_EQUAL( numberOfFineBlocksPerDirection[i] * blockSizeInCells[i], uint_c(domainAABB.size(i)),
                            "Domain can not be decomposed in direction " << i << " into fine blocks of size " << blockSizeInCells[i] );
   }

   uint_t levelScalingFactor = ( uint_t(1) << ( numberOfLevels - uint_t(1) ) );
   Vector3<uint_t> numberOfCoarseBlocksPerDirection( numberOfFineBlocksPerDirection / levelScalingFactor );

   for(uint_t i = 0; i < 3; ++i )
   {
      WALBERLA_CHECK_EQUAL(numberOfCoarseBlocksPerDirection[i] * levelScalingFactor, numberOfFineBlocksPerDirection[i],
                            "Domain can not be refined in direction " << i << " according to the specified number of levels!" );
   }


   // refinement box is in the left lower corner of the domain
   AABB refinementBox( domainAABB.xMin(), domainAABB.yMin(), domainAABB.zMin(),
                       domainAABB.xMin()+real_t(1), domainAABB.yMin()+real_t(1), domainAABB.zMin()+real_t(1) );

   WALBERLA_LOG_INFO_ON_ROOT(" - refinement box: " << refinementBox);

   sforest.addRefinementSelectionFunction( std::bind( refinementSelection, std::placeholders::_1, numberOfLevels, refinementBox ) );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

   sforest.init( domainAABB, numberOfCoarseBlocksPerDirection[0], numberOfCoarseBlocksPerDirection[1], numberOfCoarseBlocksPerDirection[2], false, false, false );

   // calculate process distribution
   const memory_t memoryLimit = math::Limits< memory_t >::inf();

   sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );

   WALBERLA_LOG_INFO_ON_ROOT( sforest );

   MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)
   shared_ptr< StructuredBlockForest > sbf =
         make_shared< StructuredBlockForest >( make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, keepGlobalBlockInformation ),
                                               blockSizeInCells[0], blockSizeInCells[1], blockSizeInCells[2]);
   sbf->createCellBoundingBoxes();

   return sbf;
}

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyFieldID ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), bodyFieldID_ ( bodyFieldID ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;


}; // class MyBoundaryHandling


BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField       = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField        = block->getData< PdfField_T > ( pdfFieldID_ );
   BodyField_T * bodyField       = block->getData< BodyField_T >( bodyFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "moving obstacle boundary handling", flagField, fluid,
          MO_SBB_T( "MO_SBB", MO_SBB_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   // boundary conditions are set by mapping the (moving) planes into the domain

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}
//*******************************************************************************************************************



//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Short test case that simulates a Lid driven cavity-like setup to check refinement near global bodies.
 *
 * The boundaries are set by mapping PE planes into the domain.
 * The top plane is moving with a given velocity.
 * Static grid refinement is applied in the left part of the domain, resulting in level boarders along the planes.
 *
 * This test case is used to check the applicability and correctness of the boundary conditions used in the
 * momentum exchange method (currently SBB) for usage along global bodies with level boarders.
 *
 * The domain is on purpose very small ( 2x1x1 coarse blocks ), to facilitate debugging.
 *
 * Special care has to be taken when refinement boundaries are on the same PE body:
 *
 *  - the mapping has to be consistent on the coarse and the fine level,
 *    i.e. a boundary cell on the coarse level has to be 8 boundary cells on the fine level.
 *    This implies that only cell-aligned and axis-aligned PE bodies are allowed to have this refinement boundary.
 *    All others (inclined planes, cylinders, spheres) must have and maintain the same refinement level
 *    throughout the entire simulation.
 *
 *  - boundary conditions that access PDF values from cells apart from the near-boundary cell
 *    (usually higher order BCs like CLI) can usually not be used.
 *    Communication patterns in LBM with refinement are complex (see Schornbaum, Ruede - "Massively Parallel Algorithms
 *    for the Lattice Boltzmann Method on NonUniform Grids" (2016)). Thus accesses to ghost layer PDF values that
 *    happen in those BCs are highly dangerous since those are often given NaN values.
 *    Problems are here two-fold:
 *     - boundary handling on fine levels is carried out in two of the four ghost layers. This is done twice on the
 *       finer levels and in this second time, the third and fourth ghost layer feature NaN values only
 *       (due to the swap during the streaming step).
 *     - on coarse blocks, the ghost layer at the coarse-fine boundary does not contain valid PDF values since the
 *       proceeding fine-to-coarse communication sets the PDF values directly inside the coarse block. They are
 *       therefore also not usable for the BC.
 *
 * Due to these restrictions, it is currently only possible to use SimpleBB in those cases.
 * Note, however, that this restriction drops if the a constant block level is maintained along a specific global body
 * throughout the simulation.
 *
 * Run this test with debug functionality to switch on the NaN checks in the refinement functions.
 *
 */
//*******************************************************************************************************************

int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   ///////////////////
   // Customization //
   ///////////////////

   //logging::Logging::instance()->setLogLevel(logging::Logging::DETAIL);

   bool vtkIO = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--vtkIO" ) == 0 ) { vtkIO = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////

   Vector3<uint_t> domainSize( 32, 16, 16 );

   const uint_t numberOfLevels = uint_t(2);
   const real_t relaxationTime = real_t(1);
   const real_t wallVelocity = real_t(0.01);

   std::string baseFolder = "vtk_out";

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t finestLevel = numberOfLevels - uint_t(1);
   const uint_t levelScalingFactor = ( uint_t(1) << finestLevel );

   const uint_t timesteps = uint_t(5);

   Vector3<uint_t> coarseBlocksPerDirection( 2, 1, 1 );
   Vector3<uint_t> blockSizeInCells(domainSize[0] / ( coarseBlocksPerDirection[0] * levelScalingFactor ),
                                    domainSize[1] / ( coarseBlocksPerDirection[1] * levelScalingFactor ),
                                    domainSize[2] / ( coarseBlocksPerDirection[2] * levelScalingFactor ) );

   AABB simulationDomain( real_t(0), real_t(0), real_t(0), real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );
   auto blocks = createBlockStructure( simulationDomain, blockSizeInCells, numberOfLevels );

   //write domain decomposition to file
   if( vtkIO )
   {
      vtk::writeDomainDecomposition( blocks, "domain_decomposition", baseFolder );
   }


   /////////////////
   // PE COUPLING //
   /////////////////

   // set up pe functionality
   shared_ptr<pe::BodyStorage>  globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");

   // create pe bodies

   // bounding planes (global)
   const auto planeMaterial = pe::createMaterial( "myPlaneMat", real_t(8920), real_t(0), real_t(1), real_t(1), real_t(0), real_t(1), real_t(1), real_t(0), real_t(0) );

   // planes in E and W direction
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(1,0,0), Vector3<real_t>(0,0,0), planeMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(-1,0,0), Vector3<real_t>(real_c(domainSize[0]),0,0), planeMaterial );

   // planes in S and N direction
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,1,0), Vector3<real_t>(0,0,0), planeMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,-1,0), Vector3<real_t>(0,real_c(domainSize[1]),0), planeMaterial );

   // planes in B and T direction
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,1), Vector3<real_t>(0,0,0), planeMaterial );
   auto topPlane = pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,real_c(domainSize[2])), planeMaterial );
   topPlane->setLinearVel(wallVelocity, real_t(0), real_t(0));

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( real_t(1) / relaxationTime, lbm::collision_model::TRT::threeSixteenth, finestLevel ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         FieldGhostLayers, field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::zyxf, FieldGhostLayers );

   // add boundary handling
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );

   // map planes into the LBM simulation -> act as no-slip boundaries
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_SBB_Flag );


   ///////////////
   // TIME LOOP //
   ///////////////

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;
   blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
   scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );
   commFunction = scheme;

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );

   auto refinementTimestep = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, sweep, pdfFieldID, boundaryHandlingID );

   if( vtkIO )
   {
      // flag field (written only once in the first time step, ghost layers are also written, written for each block separately)
      auto flagFieldVTK = vtk::createVTKOutput_BlockData(blocks, "flag_field", uint_t(1), FieldGhostLayers, false,
                                                         baseFolder, "simulation_step", false, true, true, false);
      flagFieldVTK->addCellDataWriter(make_shared<field::VTKWriter<FlagField_T> >(flagFieldID, "FlagField"));
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(flagFieldVTK), "VTK (flag field data)");

      // pdf field
      auto pdfFieldVTKfine = vtk::createVTKOutput_BlockData(blocks, "fluid_field_fine_steps", uint_t(1), FieldGhostLayers, false,
                                                            baseFolder, "simulation_step", false, true, true, false);

      pdfFieldVTKfine->addCellDataWriter( make_shared< field::VTKWriter<PdfField_T> >(pdfFieldID, "PdfField"));
      pdfFieldVTKfine->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTKfine->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      refinementTimestep->addPostCollideVoidFunction( lbm::refinement::FunctorWrapper(vtk::writeFiles( pdfFieldVTKfine )), "VTK (fluid field data post collide)", finestLevel);
      refinementTimestep->addPostBoundaryHandlingVoidFunction( lbm::refinement::FunctorWrapper(vtk::writeFiles( pdfFieldVTKfine )), "VTK (fluid field data post bh)", finestLevel);
      refinementTimestep->addPostStreamVoidFunction( lbm::refinement::FunctorWrapper(vtk::writeFiles( pdfFieldVTKfine )), "VTK (fluid field data post stream)", finestLevel);


      // pdf field
      auto pdfFieldVTKcoarse = vtk::createVTKOutput_BlockData(blocks, "fluid_field_coarse_steps", uint_t(1), FieldGhostLayers, false,
                                                              baseFolder, "simulation_step", false, true, true, false);

      pdfFieldVTKcoarse->addCellDataWriter( make_shared< field::VTKWriter<PdfField_T> >(pdfFieldID, "PdfField"));
      pdfFieldVTKcoarse->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTKcoarse->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      uint_t coarseLevel = uint_t(0);
      refinementTimestep->addPostCollideVoidFunction( lbm::refinement::FunctorWrapper(vtk::writeFiles( pdfFieldVTKcoarse )), "VTK (fluid field data collide stream)", coarseLevel);
      refinementTimestep->addPostBoundaryHandlingVoidFunction( lbm::refinement::FunctorWrapper(vtk::writeFiles( pdfFieldVTKcoarse )), "VTK (fluid field data post bh)", coarseLevel);
      refinementTimestep->addPostStreamVoidFunction( lbm::refinement::FunctorWrapper(vtk::writeFiles( pdfFieldVTKcoarse )), "VTK (fluid field data post stream)", coarseLevel);

   }

   // add LBM sweep with refinement
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( refinementTimestep ), "LBM refinement time step" );

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;
   timeloop.run(timeloopTiming);
   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace global_body_as_boundary_mem_static_refinement

int main( int argc, char **argv ){
   global_body_as_boundary_mem_static_refinement::main(argc, argv);
}
