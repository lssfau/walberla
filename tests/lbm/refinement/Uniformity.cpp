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
//! \file Uniformity.cpp
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
#include "core/math/Sample.h"
#include "core/mpi/Environment.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "field/AddToStorage.h"
#include "field/iterators/FieldIterator.h"
#include "field/StabilityChecker.h"

#include "lbm/boundary/NoSlip.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/refinement/TimeStep.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "timeloop/SweepTimeloop.h"

#include <algorithm>
#include <cstdlib>
#include <functional>

//#define TEST_USES_VTK_OUTPUT
#ifdef TEST_USES_VTK_OUTPUT
#include "vtk/all.h"
#endif



namespace uniformity_test {

///////////
// USING //
///////////

using namespace walberla;
using walberla::real_t;
using walberla::uint_t;

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

// dummy boundary handling
typedef lbm::NoSlip< LatticeModel_T, flag_t > NoSlip_T;
typedef BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T > BoundaryHandling_T;

///////////
// FLAGS //
///////////

const FlagUID  Fluid_Flag( "fluid" );
const FlagUID NoSlip_Flag( "no slip" );



/////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void refinementSelection( SetupBlockForest& forest, const uint_t levels )
{
   const AABB & domain = forest.getDomain();

   const real_t xSpan = domain.xSize() / real_t(32);
   const real_t ySpan = domain.ySize() / real_t(32);
   const real_t zSpan = domain.zSize() / real_t(64);

   const real_t xMiddle = ( domain.xMin() + domain.xMax() ) / real_t(2);
   const real_t yMiddle = ( domain.yMin() + domain.yMax() ) / real_t(2);
   const real_t zMiddle = ( domain.zMin() + domain.zMax() ) / real_t(2);

   AABB middleBox( xMiddle - xSpan, yMiddle - ySpan, zMiddle +             zSpan,
                   xMiddle + xSpan, yMiddle + ySpan, zMiddle + real_t(3) * zSpan );

   AABB shiftedBox( xMiddle +             xSpan, yMiddle +             ySpan, zMiddle +             zSpan,
                    xMiddle + real_t(3) * xSpan, yMiddle + real_t(3) * ySpan, zMiddle + real_t(3) * zSpan );

   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      if( block->getAABB().intersects( middleBox ) || block->getAABB().intersects( shiftedBox ) )
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

class MyBoundaryHandling // dummy boundary handling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldId, const BlockDataID & pdfFieldId ) :
      flagFieldId_( flagFieldId ), pdfFieldId_( pdfFieldId ) {}

   BoundaryHandling_T * operator()( IBlock* const block ) const;

private:

   const BlockDataID flagFieldId_;
   const BlockDataID  pdfFieldId_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block ) const
{
   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfFieldId_ );

   const flag_t fluid = flagField->registerFlag( Fluid_Flag );

   flagField->setWithGhostLayer( fluid );

   return new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                  NoSlip_T( "no slip", NoSlip_Flag, pdfField ) );
}



/////////
// VTK //
/////////

#ifdef TEST_USES_VTK_OUTPUT
shared_ptr< vtk::VTKOutput> createFluidFieldVTKWriter( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & pdfFieldId )
{
   auto pdfFieldVTKWriter = vtk::createVTKOutput_BlockData( blocks, "fluid_field", uint_t(500), uint_t(0) );

#ifndef TEST_BASED_ON_UNIFORM_GRID
   blockforest::communication::NonUniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< lbm::refinement::PdfFieldSyncPackInfo< LatticeModel_T > >( pdfFieldId ) );
   pdfFieldVTKWriter->addBeforeFunction( pdfGhostLayerSync );
#endif

   const auto & aabb = blocks->getDomain();
   vtk::AABBCellFilter aabbFilter( AABB( aabb.xMin(), real_t(19), aabb.zMin(), aabb.xMax(), real_t(20), aabb.zMax() ) );
   pdfFieldVTKWriter->addCellInclusionFilter( aabbFilter );

   auto velocityWriter = make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldId, "VelocityFromPDF" );
   auto  densityWriter = make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldId, "DensityFromPDF" );
   pdfFieldVTKWriter->addCellDataWriter( velocityWriter );
   pdfFieldVTKWriter->addCellDataWriter( densityWriter );

   return pdfFieldVTKWriter;
}
#endif



class AccuracyChecker
{
public:

   AccuracyChecker( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & pdfFieldId,
                    const Vector3< real_t > & velocity, const uint_t checkFrequency ) :
      blocks_( blocks ),
      executionCounter_( uint_c(0) ), checkFrequency_( ( checkFrequency > 0 ) ? checkFrequency : uint_c(1) ),
      pdfFieldId_( pdfFieldId ), velocity_( velocity ) {}

   void operator()()
   {
      ++executionCounter_;
      if( ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
         return;

      math::Sample error;

      for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
      {
         PdfField_T * pdfField = block->getData< PdfField_T >( pdfFieldId_ );

         for( auto cell = pdfField->beginXYZ(); cell != pdfField->end(); ++cell )
         {
            Vector3< real_t > velocity = pdfField->getVelocity( cell.x(), cell.y(), cell.z() );
            Vector3< real_t > diff = velocity - velocity_;

            error.insert( diff.length() / velocity_.length() );
         }
      }

      error.mpiGatherRoot();
      WALBERLA_LOG_INFO_ON_ROOT( "error:\n- min:    " << error.min() << "\n- max:    " << error.max() << "\n- avg:    " << error.mean()
                                                      << "\n- median: " << error.median() );
   }

private:

   shared_ptr< StructuredBlockStorage > blocks_;

         uint_t executionCounter_;
   const uint_t checkFrequency_;

   BlockDataID pdfFieldId_;

   Vector3< real_t > velocity_;
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

   const uint_t levels = uint_t(4);

   const uint_t xBlocks = uint_t(4);
   const uint_t yBlocks = uint_t(4);
   const uint_t zBlocks = uint_t(4);

   const uint_t xCells = shortrun ? uint_t(4) : uint_t(10);
   const uint_t yCells = shortrun ? uint_t(4) : uint_t(10);
   const uint_t zCells = shortrun ? uint_t(4) : uint_t(10);

   auto blocks = createBlockStructure( levels, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, true, true, true );

   const real_t omega = real_c(1.3);
   const Vector3< real_t > velocity( real_c(0.01), real_c(0.005), real_c(0.008) );

   WALBERLA_LOG_INFO_ON_ROOT( "Performing simulation with:"
                              "\n   - " << xBlocks << " x " << yBlocks << " x " << zBlocks << " blocks on the initial grid"
                              "\n   - " << xCells << " x " << yCells << " x " << zCells << " cells per block"
                              "\n   - " << levels << " levels"
                              "\n   - omega (coarsest grid): " << omega <<
                              "\n   - omega (finest grid): " << lbm::collision_model::levelDependentRelaxationParameter( levels - uint_t(1), omega, uint_t(0) ) <<
                              "\n   - velocity: " << velocity );

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( omega ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT( omega, 1.85 ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT( omega ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT::constructTRT( omega, lbm::collision_model::TRT::lambda_d( omega ) ) );
   //LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT( 1.19, 1.4, 1.2, omega, 1.4, 1.98 ) );

   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, velocity, real_t(1), FieldGhostLayers );
   BlockDataID tmpFieldId = lbm::addPdfFieldToStorage( blocks, "tmp field", latticeModel, FieldGhostLayers );

   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdfField_T * pdfField = block->getData< PdfField_T >( pdfFieldId );
      for( auto it = pdfField->beginWithGhostLayerXYZ(); it != pdfField->end(); ++it )
      {
         if( !pdfField->isInInnerPart( it.cell() ) )
         {
            for( uint_t i = 0; i < PdfField_T::F_SIZE; ++i )
               it[i] = std::numeric_limits< PdfField_T::value_type >::quiet_NaN();
         }
      }
      PdfField_T * tmpField = block->getData< PdfField_T >( tmpFieldId );
      std::fill( tmpField->beginWithGhostLayer(), tmpField->end(), std::numeric_limits< PdfField_T::value_type >::quiet_NaN() );
   }

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field", FieldGhostLayers );

   BlockDataID boundaryHandlingId = blocks->addBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldId, pdfFieldId ), "dummy boundary handling" );

   const uint_t timeSteps = shortrun ? uint_t(2) : uint_t(1001);
   SweepTimeloop timeloop( blocks->getBlockStorage(), timeSteps );

   auto mySweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, tmpFieldId, flagFieldId, Fluid_Flag );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, mySweep, pdfFieldId, boundaryHandlingId ) ),
                                   "LBM refinement time step" );

   if( !shortrun )
   {
      timeloop.addFuncAfterTimeStep( AccuracyChecker( blocks, pdfFieldId, velocity, uint_t(50) ), "accuracy checker" );

      timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< lbm::PdfField< LatticeModel_T >, FlagField_T >( blocks, pdfFieldId, flagFieldId, Fluid_Flag,
                                                                                                                                     uint_t(1), false, true ) ),
                                     "LBM stability check" );
   }

#ifdef TEST_USES_VTK_OUTPUT
   auto pdfFieldVTKWriter = createFluidFieldVTKWriter( blocks, pdfFieldId );
   timeloop.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTKWriter ), "VTK (fluid field data)" );
#endif

   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), 3.0 ), "Remaining time logger" );

#ifdef TEST_USES_VTK_OUTPUT
   vtk::writeDomainDecomposition( blocks );
   field::createVTKOutput< FlagField_T >( flagFieldId, *blocks, "flag_field", uint_t(1), uint_t(1) )();
#endif

   WcTimingPool timeloopTiming;
   timeloop.run( timeloopTiming );
   timeloopTiming.logResultOnRoot();

   // check constant velocity

   //typedef GhostLayerField<real_t,1> ErrorField;
   //BlockDataID errorFieldId = field::addToStorage< ErrorField >( blocks, "error field", real_t(0), field::fzyx, FieldGhostLayers );

   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdfField_T * pdfField = block->getData< PdfField_T >( pdfFieldId );
      //ErrorField * errorField = block->getData< ErrorField >( errorFieldId );

      for( auto cell = pdfField->beginXYZ(); cell != pdfField->end(); ++cell )
      {
         Vector3< real_t > cellVelocity = pdfField->getVelocity( cell.x(), cell.y(), cell.z() );
         Vector3< real_t > diff = cellVelocity - velocity;

         WALBERLA_CHECK_FLOAT_EQUAL( diff.length() / velocity.length(), real_t(0) );

         //errorField->get( cell.x(), cell.y(), cell.z() ) = diff.length() / velocity.length();
      }
   }
   //field::createVTKOutput< ErrorField >( errorFieldId, *blocks, "error_field" )();

   logging::Logging::printFooterOnStream();

   return EXIT_SUCCESS;
}

} // namespace uniformity_test

int main( int argc, char ** argv )
{
   return uniformity_test::main( argc, argv );
}
