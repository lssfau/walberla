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
//! \file BodyAtBlockBoarderCheck.cpp
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
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/all.h"
#include "core/SharedFunctor.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/refinement/all.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/basic.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include <pe_coupling/utility/all.h>

#include <vector>
#include <iomanip>
#include <iostream>
#include <functional>

namespace body_at_block_boarder_check
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

// PDF field, flag field & body field
using LatticeModel_T = lbm::D3Q19<lbm::collision_model::TRT>;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

const uint_t FieldGhostLayers = 4;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

// boundary handling
typedef pe_coupling::SimpleBB< LatticeModel_T, FlagField_T >  MO_T;
typedef BoundaryHandling< FlagField_T, Stencil_T, MO_T > BoundaryHandling_T;

using BodyTypeTuple = std::tuple<pe::Sphere> ;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag ( "fluid" );
const FlagUID MO_Flag ( "moving obstacle" );
const FlagUID FormerMO_Flag ( "former moving obstacle" );

static void refinementSelection( SetupBlockForest& forest, const uint_t levels, const AABB refinementBox )
{

   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      if( block->getAABB().intersects( refinementBox ) )
         if( block->getLevel() < ( levels - uint_t(1) ) )
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

static shared_ptr< StructuredBlockForest > createBlockStructure( AABB domainAABB, Vector3<uint_t> blockSizeInCells, uint_t numberOfLevels,
                                                                 const bool keepGlobalBlockInformation = false )
{
   SetupBlockForest sforest;

   Vector3<uint_t> numberOfCoarseBlocksPerDirection( uint_c(domainAABB.xMax()) / blockSizeInCells[0],
                                                     uint_c(domainAABB.yMax()) / blockSizeInCells[1],
                                                     uint_c(domainAABB.zMax()) / blockSizeInCells[2] );

   AABB refinementBox( domainAABB.xMin(), domainAABB.yMin(), domainAABB.zMin(), domainAABB.xMax() * real_c(0.5), domainAABB.yMax(), domainAABB.zMax() );
   sforest.addRefinementSelectionFunction( std::bind( refinementSelection, std::placeholders::_1, numberOfLevels, refinementBox ) );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

   sforest.init( domainAABB, numberOfCoarseBlocksPerDirection[0], numberOfCoarseBlocksPerDirection[1], numberOfCoarseBlocksPerDirection[2], true, true, true );

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

   BoundaryHandling_T * handling = new BoundaryHandling_T( "fixed obstacle boundary handling", flagField, fluid,
                                                           MO_T (  "MO_BB",  MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}

class BodyEvaluator
{
   public:
   BodyEvaluator( const shared_ptr< StructuredBlockStorage > & blocks,
                  const BlockDataID & bodyStorageID, const pe::Vec3 & referenceVelocity )
      : blocks_( blocks ), bodyStorageID_( bodyStorageID ), referenceVelocity_( referenceVelocity )
      { }

      void operator()()
      {
         for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
         {
            for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
            {
                std::cout << "Body at position " << bodyIt->getPosition() << " has velocity " << bodyIt->getLinearVel() << ".\n";
                pe::Vec3 velocityDiff = referenceVelocity_ - bodyIt->getLinearVel();
                WALBERLA_CHECK_LESS( std::fabs( velocityDiff[0] ), real_c(1e-5) );
                WALBERLA_CHECK_LESS( std::fabs( velocityDiff[1] ), real_c(1e-5) );
                WALBERLA_CHECK_LESS( std::fabs( velocityDiff[2] ), real_c(1e-5) );
            }
         }
      }
private:
      shared_ptr< StructuredBlockStorage > blocks_;
      const BlockDataID bodyStorageID_;
      const pe::Vec3 referenceVelocity_;
};


/*
 * This test case evaluates the correctness of the pe and the pe-coupling at the boundary between two blocks.
 * To do so, two spheres are created, one completely inside one block and one at the boarder between several blocks.
 * The flow field is initialized with a uniform velocity of <0.1, 0, 0> and so are the spheres.
 * The hydrodynamic force on the spheres should thus be <0,0,0> and the spheres' velocity should remain unchanged.
 * Additionally, the combination with refinement is tested. The total number of used levels can be specified via an input argument, the default being 1.
 */

//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   Vector3<real_t> domainSize( real_c(128), real_c(64), real_c(64) );
   Vector3<uint_t> blockSizeInCells( uint_c(32), uint_c(32), uint_c(32) );
   char *p;
   uint_t numberOfLevels = (argc == 2) ? uint_c(std::strtol(argv[1], &p, 10)) : uint_c(1);

   const real_t omega  = real_c(1);
   const real_t dx     = real_c(1);
   const real_t radius = real_c(5);

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   math::AABB simulationDomain = math::AABB( pe::Vec3( real_c(0), real_c(0), real_c(0) ), domainSize );

   // create fully periodic domain with refined blocks
   auto blocks = createBlockStructure( simulationDomain, blockSizeInCells, numberOfLevels );

   ////////
   // PE //
   ////////

   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "Storage");
   auto ccdID         = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   auto fcdID         = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");
   pe::cr::HCSITS cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID);


   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe
   const real_t overlap = real_c( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );

   auto sphereMaterialID = pe::createMaterial( "sphereMat", real_c(1) , real_c(0.3), real_c(0.2), real_c(0.2), real_c(0.24), real_c(200), real_c(200), real_c(0), real_c(0) );
   // create two spheres: one which overlaps with a block boundary and one inside the block
   pe::Vec3 referenceVelocity( real_c(0.1), real_c(0), real_c(0) );
   auto sphereAtBorder = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_c(31), real_c(31), real_c(31) ), radius, sphereMaterialID );
   if( sphereAtBorder != nullptr )
   {
      sphereAtBorder->setLinearVel( referenceVelocity );
   }
   auto sphereInInner = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_c(7), real_c(7), real_c(7) ), radius, sphereMaterialID );
   if( sphereInInner != nullptr )
   {
      sphereInInner->setLinearVel( referenceVelocity );
   }
   syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );

   // add PDF field ( uInit = <0.1,0,0>, rhoInit = 1 )
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         referenceVelocity, real_c(1),
                                                                         FieldGhostLayers, field::fzyx );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::fzyx, FieldGhostLayers );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
                                    MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );

   // initially map pe bodies into the LBM simulation
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID,  MO_Flag );

   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   uint_t timeSteps( uint_c(10) );
   SweepTimeloop timeloop( blocks->getBlockStorage(), timeSteps );

   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   auto refinementTimestep = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, sweep, pdfFieldID, boundaryHandlingID );

   const uint_t finestLevel = numberOfLevels - uint_c(1);

   BodyEvaluator evalBodies( blocks, bodyStorageID, referenceVelocity );

   // add body evaluator
   refinementTimestep->addPostStreamVoidFunction( lbm::refinement::FunctorWrapper( evalBodies ), "Body Evaluation", finestLevel );

   // add pe timestep
   const real_t dtPE = real_c(1) / real_c( uint_c(1) << finestLevel );

   refinementTimestep->addPostStreamVoidFunction( lbm::refinement::FunctorWrapper( pe_coupling::TimeStep( blocks, bodyStorageID, cr, syncCall, dtPE, 1 ) ),
                                                  "pe Time Step", finestLevel );

   // add sweep for updating the pe body mapping into the LBM simulation
   refinementTimestep->addPostStreamVoidFunction( lbm::refinement::SweepAsFunctorWrapper( pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >( blocks, pdfFieldID, boundaryHandlingID,
                                                  bodyStorageID, globalBodyStorage, bodyFieldID, MO_Flag, FormerMO_Flag ), blocks ),
                                                  "Body Mapping", finestLevel );

   // add sweep for restoring PDFs in cells previously occupied by pe bodies
   typedef pe_coupling::EquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T > Reconstructor_T;
   Reconstructor_T reconstructor( blocks, boundaryHandlingID, bodyFieldID );
   refinementTimestep->addPostStreamVoidFunction( lbm::refinement::SweepAsFunctorWrapper( pe_coupling::PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructor_T > ( blocks, pdfFieldID,
                                                  boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, reconstructor, FormerMO_Flag, Fluid_Flag  ), blocks ),
                                                  "PDF Restore", finestLevel );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( refinementTimestep ), "LBM refinement time step" );

   timeloop.run();

   return 0;

}

} //namespace body_at_block_boarder_check

int main( int argc, char **argv ){
   body_at_block_boarder_check::main(argc, argv);
}
