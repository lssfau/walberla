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
//! \file CommunicationEquivalence.cpp
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/communication/NonUniformBufferedScheme.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "boundary/BoundaryHandling.h"

#include "core/DataTypes.h"
#include "core/SharedFunctor.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Limits.h"
#include "core/mpi/Environment.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "field/AddToStorage.h"
#include "field/iterators/FieldIterator.h"
#include "field/StabilityChecker.h"

#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimpleUBB.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/refinement/TimeStep.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "timeloop/SweepTimeloop.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>

//#define TEST_USES_VTK_OUTPUT
#ifdef TEST_USES_VTK_OUTPUT
#include "vtk/all.h"
#endif

enum TestMode { ENTIRE_TOP, TOP, MIDDLE, ENTIRE_BOTTOM };
static const TestMode testMode = TOP;

namespace walberla{

//////////////
// TYPEDEFS //
//////////////

typedef lbm::D3Q19< lbm::collision_model::SRT,      false > LatticeModel_T;
//typedef lbm::D3Q19< lbm::collision_model::TRT,      false > LatticeModel_T;
//typedef lbm::D3Q19< lbm::collision_model::D3Q19MRT, false > LatticeModel_T;
using Stencil_T = LatticeModel_T::Stencil;

using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers = 4;

typedef lbm::NoSlip< LatticeModel_T, flag_t >     NoSlip_T;
typedef lbm::SimpleUBB< LatticeModel_T, flag_t >  UBB_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, UBB_T > BoundaryHandling_T;

///////////
// FLAGS //
///////////

const FlagUID  Fluid_Flag( "fluid" );
const FlagUID    UBB_Flag( "velocity bounce back" );
const FlagUID NoSlip_Flag( "no slip" );



/////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void refinementSelection( SetupBlockForest& forest, const uint_t levels )
{
   const AABB & domain = forest.getDomain();

   if( testMode == ENTIRE_TOP )
   {
      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         auto aabb = block->getAABB();
         if( realIsEqual( aabb.zMax(), domain.zMax() ) )
            if( block->getLevel() < ( levels - uint_t(1) ) )
               block->setMarker( true );
      }
   }
   else if( testMode == TOP )
   {
      AABB topCorner( domain.xMin(),
                      domain.yMin(),
                      domain.zMax() - domain.zMax() / real_t(14),
                      domain.xMin() + domain.xMax() / real_t(14),
                      domain.yMin() + domain.yMax() / real_t(14),
                      domain.zMax() );

      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         if( block->getAABB().intersects( topCorner ) )
            if( block->getLevel() < ( levels - uint_t(1) ) )
               block->setMarker( true );
      }
   }
   else if( testMode == MIDDLE )
   {
      const real_t xSpan = domain.xSize() / real_t(32);
      const real_t ySpan = domain.ySize() / real_t(32);
      const real_t zSpan = domain.zSize() / real_t(32);

      const real_t xMiddle = ( domain.xMin() + domain.xMax() ) / real_t(2);
      const real_t yMiddle = ( domain.yMin() + domain.yMax() ) / real_t(2);
      const real_t zMiddle = ( domain.zMin() + domain.zMax() ) / real_t(2);

      AABB middleBox( xMiddle - xSpan, yMiddle - ySpan, zMiddle - zSpan,
                      xMiddle + xSpan, yMiddle + ySpan, zMiddle + zSpan );

      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         if( block->getAABB().intersects( middleBox ) )
            if( block->getLevel() < ( levels - uint_t(1) ) )
               block->setMarker( true );
      }
   }
   else if( testMode == ENTIRE_BOTTOM )
   {
      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         auto aabb = block->getAABB();
         if( realIsEqual( aabb.zMin(), domain.zMin() ) )
            if( block->getLevel() < ( levels - uint_t(1) ) )
               block->setMarker( true );
      }
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

static shared_ptr< StructuredBlockForest > createBlockStructure( const uint_t levels,
                                                                 const uint_t numberOfXBlocks,        const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                                                                 const uint_t numberOfXCellsPerBlock, const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                                                                 const bool   xPeriodic = false,      const bool   yPeriodic = false,      const bool   zPeriodic = false,
                                                                 const bool keepGlobalBlockInformation = false )
{
   // initialize SetupBlockForest = determine domain decomposition
   SetupBlockForest sforest;

   sforest.addRefinementSelectionFunction( std::bind( refinementSelection, std::placeholders::_1, levels ) );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

   sforest.init( AABB( real_c(0), real_c(0), real_c(0), real_c( numberOfXBlocks * numberOfXCellsPerBlock ),
                                                        real_c( numberOfYBlocks * numberOfYCellsPerBlock ),
                                                        real_c( numberOfZBlocks * numberOfZCellsPerBlock ) ),
                 numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, xPeriodic, yPeriodic, zPeriodic );

   // calculate process distribution
   const memory_t memoryLimit = math::Limits< memory_t >::inf();

   sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );

   WALBERLA_LOG_INFO_ON_ROOT( sforest );

   MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)
   shared_ptr< StructuredBlockForest > sbf =
         make_shared< StructuredBlockForest >( make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, keepGlobalBlockInformation ),
                                               numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock );
   sbf->createCellBoundingBoxes();
   return sbf;
}



///////////////////////
// BOUNDARY HANDLING //
///////////////////////

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldId, const BlockDataID & pdfFieldId, const real_t topVelocity ) :
      flagFieldId_( flagFieldId ), pdfFieldId_( pdfFieldId ), topVelocity_( topVelocity ) {}

   BoundaryHandling_T * operator()( IBlock* const block ) const;

private:

   const BlockDataID flagFieldId_;
   const BlockDataID  pdfFieldId_;

   const real_t topVelocity_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block ) const
{
   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfFieldId_ );

   const flag_t fluid = flagField->registerFlag( Fluid_Flag );

   return new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                  NoSlip_T( "no slip", NoSlip_Flag, pdfField ),
                                     UBB_T( "velocity bounce back", UBB_Flag, pdfField, topVelocity_, real_c(0), real_c(0) ) );
}



void setFlags( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & boundaryHandlingId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId );

      const uint_t level = blocks->getLevel(*block);

      CellInterval domainBB = blocks->getDomainCellBB( level );
      blocks->transformGlobalToBlockLocalCellInterval( domainBB, *block );

      domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
      domainBB.xMax() += cell_idx_c( FieldGhostLayers );
      domainBB.yMin() -= cell_idx_c( FieldGhostLayers );
      domainBB.yMax() += cell_idx_c( FieldGhostLayers );
      domainBB.zMin() -= cell_idx_c( 1 );
      domainBB.zMax() += cell_idx_c( 1 );

      // TOP
      CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      boundaryHandling->forceBoundary( UBB_Flag, top );

      // BOTTOM
      CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
      boundaryHandling->forceBoundary( NoSlip_Flag, bottom );

      boundaryHandling->fillWithDomain( FieldGhostLayers );
   }
}



/////////
// VTK //
/////////

#ifdef TEST_USES_VTK_OUTPUT
shared_ptr< vtk::VTKOutput> createFluidFieldVTKWriter( shared_ptr< StructuredBlockForest > & blocks,
                                                       const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId )
{
   auto pdfFieldVTKWriter = vtk::createVTKOutput_BlockData( blocks, "fluid_field", uint_t(50), uint_t(1), false );

   blockforest::communication::NonUniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< lbm::refinement::PdfFieldSyncPackInfo< LatticeModel_T > >( pdfFieldId ) );
   pdfFieldVTKWriter->addBeforeFunction( pdfGhostLayerSync );

   field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldId );
   fluidFilter.addFlag( Fluid_Flag );

   const auto & aabb = blocks->getDomain();

   vtk::ChainedFilter combine;
   combine.addFilter( fluidFilter );
   if( testMode == MIDDLE )
   {
      vtk::AABBCellFilter aabbFilter( AABB( aabb.xMin(), real_t(7), aabb.zMin(), aabb.xMax(), real_t(8), aabb.zMax() ) );
      combine.addFilter( aabbFilter );
   }
   else
   {
      vtk::AABBCellFilter aabbFilter( AABB( aabb.xMin(), real_t(1), aabb.zMin(), aabb.xMax(), real_t(2), aabb.zMax() ) );
      combine.addFilter( aabbFilter );
   }
   pdfFieldVTKWriter->addCellInclusionFilter( combine );

   auto velocityWriter = make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldId, "VelocityFromPDF" );
   auto  densityWriter = make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldId, "DensityFromPDF" );
   pdfFieldVTKWriter->addCellDataWriter( velocityWriter );
   pdfFieldVTKWriter->addCellDataWriter( densityWriter );

   return pdfFieldVTKWriter;
}
#endif



class EqualityChecker
{
public:

   EqualityChecker( const shared_ptr< StructuredBlockStorage > & blocks,
                    const BlockDataID & pdfFieldId1, const BlockDataID & pdfFieldId2 ) :
      blocks_( blocks ), pdfFieldId1_( pdfFieldId1 ), pdfFieldId2_( pdfFieldId2 ) {}

   void operator()()
   {
      for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
      {
         PdfField_T * pdfField1 = block->getData< PdfField_T >( pdfFieldId1_ );
         PdfField_T * pdfField2 = block->getData< PdfField_T >( pdfFieldId2_ );

         auto it1 = pdfField1->begin();
         auto it2 = pdfField2->begin();

         while( it1 != pdfField1->end() )
         {
            WALBERLA_CHECK_FLOAT_EQUAL( *it1, *it2 );
            ++it1;
            ++it2;
         }
         WALBERLA_CHECK_EQUAL( it2, pdfField2->end() );
      }
   }

private:

   shared_ptr< StructuredBlockStorage > blocks_;

   BlockDataID pdfFieldId1_;
   BlockDataID pdfFieldId2_;
};



//////////
// MAIN //
//////////

int main( int argc, char ** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   bool shortrun = false;
   for( int i = 1; i < argc; ++i )
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) shortrun = true;

   logging::Logging::printHeaderOnStream();

   uint_t levels = shortrun ? uint_t(3) : uint_t(4);
   if( testMode == ENTIRE_TOP || testMode == ENTIRE_BOTTOM )
     levels = shortrun ? uint_t(3) : uint_t(5);

   uint_t xBlocks = uint_t(1);
   uint_t yBlocks = uint_t(1);
   uint_t zBlocks = uint_t(1);
   if( testMode == ENTIRE_TOP || testMode == ENTIRE_BOTTOM )
   {
      zBlocks = uint_t(2);
   }
   else if( testMode == TOP )
   {
      xBlocks = uint_t(3);
      yBlocks = uint_t(3);
      zBlocks = uint_t(3);
   }
   else if( testMode == MIDDLE )
   {
      xBlocks = uint_t(4);
      yBlocks = uint_t(4);
      zBlocks = uint_t(4);
   }

   const uint_t xCells = shortrun ? uint_t(4) : uint_t(10);
   const uint_t yCells = ( testMode == TOP && !shortrun ) ? uint_t(8) : uint_t(4);
   const uint_t zCells = shortrun ? uint_t(4) : uint_t(10);

   auto blocks = createBlockStructure( levels, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, true, true, false );
   
   const real_t Re = real_t(10);
   const real_t omega = ( testMode == ENTIRE_TOP || testMode == ENTIRE_BOTTOM ) ? real_c(1.9) : real_c(1.3);
   const real_t nu = ( real_t(1) / omega - real_c(0.5) ) / real_t(3);
   const real_t L = real_c( zBlocks * zCells );
   const real_t topVelocity = ( Re * nu ) / L;

   WALBERLA_LOG_INFO_ON_ROOT( "Performing Couette simulation with:"
                              "\n   - " << xBlocks << " x " << yBlocks << " x " << zBlocks << " blocks on the initial grid"
                              "\n   - " << xCells << " x " << yCells << " x " << zCells << " cells per block"
                              "\n   - " << levels << " levels"
                              "\n   - Reynolds number: " << Re <<
                              "\n   - omega (coarsest grid): " << omega <<
                              "\n   - omega (finest grid): " << lbm::collision_model::levelDependentRelaxationParameter( levels - uint_t(1), omega, uint_t(0) ) <<
                              "\n   - top velocity: " << topVelocity );

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( omega ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT( omega, 1.85 ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT( omega ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT::constructTRT( omega, lbm::collision_model::TRT::lambda_d( omega ) ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT( 1.19, 1.4, 1.2, omega, 1.4, 1.98 ) );

   BlockDataID pdfFieldId1 = lbm::addPdfFieldToStorage( blocks, "pdf field (1)", latticeModel,
                                                       Vector3< real_t >( topVelocity / real_c(2), real_c(0), real_c(0) ),
                                                       real_t(1), FieldGhostLayers );

   BlockDataID pdfFieldId2 = lbm::addPdfFieldToStorage( blocks, "pdf field (2)", latticeModel,
                                                       Vector3< real_t >( topVelocity / real_c(2), real_c(0), real_c(0) ),
                                                       real_t(1), FieldGhostLayers );

   BlockDataID flagFieldId1 = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field (1)", FieldGhostLayers );
   BlockDataID flagFieldId2 = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field (2)", FieldGhostLayers );

   BlockDataID boundaryHandlingId1 = blocks->addBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldId1, pdfFieldId1, topVelocity ),
                                                                                 "boundary handling (1)" );
   setFlags( blocks, boundaryHandlingId1 );

   BlockDataID boundaryHandlingId2 = blocks->addBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldId2, pdfFieldId2, topVelocity ),
                                                                                 "boundary handling (2)" );
   setFlags( blocks, boundaryHandlingId2);

   uint_t timeSteps = shortrun ? uint_t(2) : uint_t(101);
      
   SweepTimeloop timeloop( blocks->getBlockStorage(), timeSteps );
   
   auto mySweep1 = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId1, flagFieldId1, Fluid_Flag );
   auto mySweep2 = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId2, flagFieldId2, Fluid_Flag );

   auto tstep1 = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, mySweep1, pdfFieldId1, boundaryHandlingId1 );
   auto tstep2 = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, mySweep2, pdfFieldId2, boundaryHandlingId2 );
   tstep1->optimizeCommunication( true );
   tstep2->optimizeCommunication( false );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( tstep1 ), "LBM refinement time step (1)" );
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( tstep2 ), "LBM refinement time step (2)" );
   
   timeloop.addFuncAfterTimeStep( EqualityChecker( blocks, pdfFieldId1, pdfFieldId2 ), "equivalence checker" );

   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< lbm::PdfField< LatticeModel_T >, FlagField_T >( blocks, pdfFieldId1, flagFieldId1, Fluid_Flag,
                                                                                                                                  uint_t(1), false, true ) ),
                                  "LBM stability check" );

#ifdef TEST_USES_VTK_OUTPUT
   auto pdfFieldVTKWriter = createFluidFieldVTKWriter( blocks, pdfFieldId1, flagFieldId1 );
   timeloop.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTKWriter ), "VTK (fluid field data)" );
#endif

   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), 3.0 ), "Remaining time logger" );

#ifdef TEST_USES_VTK_OUTPUT
   vtk::writeDomainDecomposition( blocks, "domain_decomposition" );
   field::createVTKOutput< FlagField_T >( flagFieldId1, *blocks, "flag_field", uint_t(1), uint_t(1), false )();
#endif

   WcTimingPool timeloopTiming;
   timeloop.run( timeloopTiming );
   timeloopTiming.logResultOnRoot();

   logging::Logging::printFooterOnStream();
   
   return EXIT_SUCCESS;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
