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
//! \file NonUniformGrid.cpp
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/NonUniformBufferedScheme.h"
#include "blockforest/loadbalancing/DynamicCurve.h"
#include "blockforest/loadbalancing/DynamicDiffusive.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/SharedFunctor.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/timing/TimingPool.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/FlagUID.h"
#include "field/iterators/FieldIterator.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"

#include "lbm/BlockForestEvaluation.h"
#include "lbm/PerformanceEvaluation.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimpleUBB.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/refinement/PdfFieldSyncPackInfo.h"
#include "lbm/refinement/TimeStep.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SplitPureSweep.h"
#include "lbm/sweeps/SplitSweep.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"

#include "sqlite/SQLite.h"

#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/BlockCellDataWriter.h"
#include "vtk/Initialization.h"
#include "vtk/VTKOutput.h"

#include <cstdlib>
#include <functional>
#include <iostream>
#include <memory>
#include <type_traits>



namespace non_uniform_grid {

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;
using walberla::real_t;

//////////////
// TYPEDEFS //
//////////////

using D3Q19_SRT_INCOMP = lbm::D3Q19<lbm::collision_model::SRT, false>;
using D3Q19_SRT_COMP = lbm::D3Q19<lbm::collision_model::SRT, true>;
using D3Q19_TRT_INCOMP = lbm::D3Q19<lbm::collision_model::TRT, false>;
using D3Q19_TRT_COMP = lbm::D3Q19<lbm::collision_model::TRT, true>;
using D3Q19_MRT_INCOMP = lbm::D3Q19<lbm::collision_model::D3Q19MRT, false>;

template< typename LatticeModel_T >
struct Types
{
   using Stencil_T = typename LatticeModel_T::Stencil;
   using PdfField_T = lbm::PdfField< LatticeModel_T >;
};

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers  = uint_t(4);
const uint_t BlockForestLevels = uint_t(4);

///////////
// FLAGS //
///////////

const FlagUID  Fluid_Flag( "fluid" );
const FlagUID    UBB_Flag( "velocity bounce back" );
const FlagUID NoSlip_Flag( "no slip" );






/////////////////////
// OUTPUT HELPERS  //
/////////////////////

template< typename LatticeModel_T, class Enable = void >
struct StencilString;

template< typename LatticeModel_T >
struct StencilString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value >::type >
{
   static const char * str() { return "D3Q19"; }
};

template< typename LatticeModel_T, class Enable = void >
struct CollisionModelString;

template< typename LatticeModel_T >
struct CollisionModelString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                          lbm::collision_model::SRT_tag >::value >::type >
{
   static const char * str() { return "SRT"; }
};

template< typename LatticeModel_T >
struct CollisionModelString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                          lbm::collision_model::TRT_tag >::value >::type >
{
   static const char * str() { return "TRT"; }
};

template< typename LatticeModel_T >
struct CollisionModelString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                          lbm::collision_model::MRT_tag >::value >::type >
{
   static const char * str() { return "MRT"; }
};






/////////////
// CONFIG  //
/////////////

static inline void getCells( const Config::BlockHandle & configBlock, uint_t & xCells, uint_t & yCells, uint_t & zCells )
{
   xCells = configBlock.getParameter< uint_t >( "xCells", 10 );
   yCells = configBlock.getParameter< uint_t >( "yCells", 10 );
   zCells = configBlock.getParameter< uint_t >( "zCells", 10 );
}

static void getCellsAndCoarseBlocks( const Config::BlockHandle & configBlock, const uint_t numberOfProcesses, const uint_t fineBlocksPerProcess,
                                     uint_t & xCells,  uint_t & yCells,  uint_t & zCells,
                                     uint_t & xBlocks, uint_t & yBlocks, uint_t & zBlocks )
{
   WALBERLA_CHECK( ( fineBlocksPerProcess == uint_t(1) ) || ( fineBlocksPerProcess == uint_t(2) ) ||
                   ( fineBlocksPerProcess == uint_t(4) ) || ( fineBlocksPerProcess == uint_t(8) ) );
   WALBERLA_CHECK_GREATER_EQUAL( numberOfProcesses, uint_t(64) / fineBlocksPerProcess );
   WALBERLA_CHECK_EQUAL( numberOfProcesses % ( uint_t(64) / fineBlocksPerProcess ), uint_t(0) );

   getCells( configBlock, xCells, yCells, zCells );

   xBlocks = uint_t(3);
   yBlocks = numberOfProcesses / ( uint_t(64) / fineBlocksPerProcess );
   zBlocks = uint_t(3);
}

#ifndef NDEBUG
static inline uint_t numberOfBlocks( const uint_t numberOfProcesses, const uint_t fineBlocksPerProcess )
{
   return ( numberOfProcesses / ( uint_t(64) / fineBlocksPerProcess ) ) * uint_t(107);
}

static uint_t numberOfBlocksOnLevel( const uint_t level, const uint_t numberOfProcesses, const uint_t fineBlocksPerProcess )
{
   const uint_t factor = numberOfProcesses / ( uint_t(64) / fineBlocksPerProcess );

   if( level == uint_t(0) )
      return factor * uint_t(7);
   else if( level == uint_t(1) )
      return factor * uint_t(12);
   else if( level == uint_t(2) )
      return factor * uint_t(24);
   else if( level == uint_t(3) )
      return factor * uint_t(64);

   return uint_t(0);
}
#endif






///////////////////////////
// BLOCK STRUCTURE SETUP //
///////////////////////////

static void refinementSelection( SetupBlockForest& forest )
{
   const AABB & domain = forest.getDomain();

   real_t xSize = ( domain.xSize() / real_t(12) ) * real_c( 0.99 );
   real_t zSize = ( domain.zSize() / real_t(12) ) * real_c( 0.99 );

   AABB leftCorner( domain.xMin(), domain.yMin(), domain.zMax() - zSize,
                    domain.xMin() + xSize, domain.yMax(), domain.zMax() );

   AABB rightCorner( domain.xMax() - xSize, domain.yMin(), domain.zMax() - zSize,
                     domain.xMax(), domain.yMax(), domain.zMax() );

   for(auto & block : forest)
   {
      auto & aabb = block.getAABB();
      if( leftCorner.intersects( aabb ) || rightCorner.intersects( aabb ) )
      {
         if( block.getLevel() < ( BlockForestLevels - uint_t(1) ) )
            block.setMarker( true );
      }
   }
}

static void workloadAndMemoryAssignment( SetupBlockForest & forest, const memory_t memoryPerBlock ) {

   for(auto & block : forest)
   {
      block.setWorkload( numeric_cast< workload_t >( uint_t(1) << block.getLevel() ) );
      block.setMemory( memoryPerBlock );
   }
}

void createSetupBlockForest( blockforest::SetupBlockForest & sforest, const Config::BlockHandle & configBlock, const uint_t numberOfProcesses )
{
   uint_t numberOfXCellsPerBlock;
   uint_t numberOfYCellsPerBlock;
   uint_t numberOfZCellsPerBlock;
   uint_t numberOfXBlocks;
   uint_t numberOfYBlocks;
   uint_t numberOfZBlocks;
   
   const uint_t bufferProcesses = configBlock.getParameter< uint_t >( "bufferProcesses", 0 );
   const uint_t fineBlocksPerProcess = configBlock.getParameter< uint_t >( "fineBlocksPerProcess", 4 );
   
   const uint_t workerProcesses = numberOfProcesses - bufferProcesses;

   getCellsAndCoarseBlocks( configBlock, workerProcesses, fineBlocksPerProcess, numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock,
                            numberOfXBlocks, numberOfYBlocks, numberOfZBlocks );

   const memory_t memoryPerBlock = numeric_cast< memory_t >( ( numberOfXCellsPerBlock + uint_t(2) * FieldGhostLayers ) *
                                                             ( numberOfYCellsPerBlock + uint_t(2) * FieldGhostLayers ) *
                                                             ( numberOfZCellsPerBlock + uint_t(2) * FieldGhostLayers ) *
                                                             uint_c( 19 * sizeof(real_t) ) ) / numeric_cast< memory_t >( 1024 * 1024 );

   sforest.addRefinementSelectionFunction( refinementSelection );
   sforest.addWorkloadMemorySUIDAssignmentFunction( std::bind( workloadAndMemoryAssignment, std::placeholders::_1, memoryPerBlock ) );

   sforest.init( AABB( real_t(0), real_t(0), real_t(0), real_c( numberOfXBlocks * numberOfXCellsPerBlock ),
                                                        real_c( numberOfYBlocks * numberOfYCellsPerBlock ),
                                                        real_c( numberOfZBlocks * numberOfZCellsPerBlock ) ),
                 numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, false, false, false );

#ifndef NDEBUG
   WALBERLA_ASSERT_EQUAL( sforest.getNumberOfBlocks(), numberOfBlocks( workerProcesses, fineBlocksPerProcess ) );
   for( auto i = uint_t(0); i != BlockForestLevels; ++i )
   {
      std::vector< blockforest::SetupBlock * > blocks;
      sforest.getBlocks( blocks, i );
      WALBERLA_ASSERT_EQUAL( blocks.size(), numberOfBlocksOnLevel( i, workerProcesses, fineBlocksPerProcess ) );
   }
#endif

   MPIManager::instance()->useWorldComm();

   const memory_t memoryLimit = configBlock.getParameter< memory_t >( "memoryLimit", numeric_cast< memory_t >(256) );

   sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), numberOfProcesses, bufferProcesses, memoryLimit, true,
                        bufferProcesses != uint_t(0) );

   WALBERLA_LOG_INFO_ON_ROOT( "SetupBlockForest created successfully:\n" << sforest );
}



shared_ptr< blockforest::StructuredBlockForest > createStructuredBlockForest( const Config::BlockHandle & configBlock )
{
   uint_t numberOfXCellsPerBlock;
   uint_t numberOfYCellsPerBlock;
   uint_t numberOfZCellsPerBlock;
   getCells( configBlock, numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock );

   if( configBlock.isDefined( "sbffile" ) )
   {
      std::string sbffile = configBlock.getParameter< std::string >( "sbffile" );

      WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure: loading from file \'" << sbffile << "\' ..." );

      MPIManager::instance()->useWorldComm();

      auto bf = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sbffile.c_str(), true, false );

      auto sbf = std::make_shared< StructuredBlockForest >( bf, numberOfXCellsPerBlock,
                                                                                     numberOfYCellsPerBlock,
                                                                                     numberOfZCellsPerBlock );
      sbf->createCellBoundingBoxes();

      return sbf;
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure ..." );

   blockforest::SetupBlockForest sforest;
   createSetupBlockForest( sforest, configBlock, uint_c( MPIManager::instance()->numProcesses() ) );

   auto bf = std::make_shared< blockforest::BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, false );

   auto sbf = std::make_shared< blockforest::StructuredBlockForest >( bf, numberOfXCellsPerBlock,
                                                                                                            numberOfYCellsPerBlock,
                                                                                                            numberOfZCellsPerBlock );
   sbf->createCellBoundingBoxes();

   return sbf;
}



class ReGrid
{
public:

   ReGrid( const weak_ptr< StructuredBlockForest > & forest, const int type, const uint_t regridAt ) :
      forest_( forest ), type_( type ), executionCounter_( uint_t(0) ), regridAt_( regridAt )
   {}

   void operator()()
   {
      if( ( executionCounter_ + uint_t(1) ) == regridAt_ )
      {
         auto forest = forest_.lock();
         WALBERLA_CHECK_NOT_NULLPTR( forest );
         forest->getBlockForest().refresh();
      }
      ++executionCounter_;
   }
   
   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > &, const BlockForest & forest );

private:

   weak_ptr< StructuredBlockForest > forest_;
   
   int type_;
   
   uint_t executionCounter_;
   uint_t regridAt_;
};

void ReGrid::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                         std::vector< const Block * > &, const BlockForest & forest )
{
   const AABB & domain = forest.getDomain();
   
   const real_t xSize = domain.xSize() / real_t(12);
   const real_t zSize = domain.zSize() / real_t(12);
   
   AABB left;
   AABB right;
   
   if( type_ == 0 )
   {
      left.init( domain.xMin() + real_c(1.01) * xSize, domain.yMin(), domain.zMax() - real_c(0.99) * zSize,
                 domain.xMin() + real_c(1.99) * xSize, domain.yMax(), domain.zMax() );
      right.init( domain.xMax() - real_c(1.99) * xSize, domain.yMin(), domain.zMax() - real_c(0.99) * zSize, 
                  domain.xMax() - real_c(1.01) * xSize, domain.yMax(), domain.zMax() );
   }
   else
   {
      left.init( domain.xMin(), domain.yMin(), domain.zMin(),
                 domain.xMin() + real_c(0.99) * xSize, domain.yMax(), domain.zMin() + real_c(0.99) * zSize );
      right.init( domain.xMax() - real_c(0.99) * xSize, domain.yMin(), domain.zMin(),
                  domain.xMax(), domain.yMax(), domain.zMin() + real_c(0.99) * zSize );
   }
   
   for(auto & minTargetLevel : minTargetLevels)
   {
      auto & aabb = minTargetLevel.first->getAABB();
      if( left.intersects( aabb ) || right.intersects( aabb ) )
         minTargetLevel.second = BlockForestLevels - uint_t(1);
      else
         minTargetLevel.second = uint_t(0);
   }
}






///////////////////////
// BOUNDARY HANDLING //
///////////////////////

template< typename LatticeModel_T >
struct MyBoundaryTypes
{
   using NoSlip_T = lbm::NoSlip<LatticeModel_T, flag_t>;
   using UBB_T = lbm::SimpleUBB<LatticeModel_T, flag_t>;

   using BoundaryHandling_T = BoundaryHandling<FlagField_T, typename Types<LatticeModel_T>::Stencil_T, NoSlip_T, UBB_T>;
};  

template< typename LatticeModel_T >
class MyBoundaryHandling : public blockforest::AlwaysInitializeBlockDataHandling< typename MyBoundaryTypes< LatticeModel_T >::BoundaryHandling_T >
{
public:

   using NoSlip_T = typename MyBoundaryTypes< LatticeModel_T >::NoSlip_T;
   using UBB_T = typename MyBoundaryTypes< LatticeModel_T >::UBB_T;

   using BoundaryHandling_T = typename MyBoundaryTypes< LatticeModel_T >::BoundaryHandling_T;



   MyBoundaryHandling( const weak_ptr< StructuredBlockForest > & forest,
                       const BlockDataID & flagField, const BlockDataID & pdfField, const real_t velocity ) :
      forest_( forest ), flagField_( flagField ), pdfField_( pdfField ), velocity_( velocity ) {}

   BoundaryHandling_T * initialize( IBlock * const block ) override;

private:

   weak_ptr< StructuredBlockForest > forest_;

   const BlockDataID flagField_;
   const BlockDataID  pdfField_;

   const real_t velocity_;

}; // class MyBoundaryHandling

template< typename LatticeModel_T >
typename MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T *
MyBoundaryHandling<LatticeModel_T>::initialize( IBlock * const block )
{
   using PdfField_T = typename Types<LatticeModel_T>::PdfField_T;

   WALBERLA_ASSERT_NOT_NULLPTR( block );

   FlagField_T * flagField = block->getData< FlagField_T >( flagField_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfField_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                                           NoSlip_T( "no slip", NoSlip_Flag, pdfField ),
                                                           UBB_T( "velocity bounce back", UBB_Flag, pdfField, velocity_, real_c(0), real_c(0) ) );

   auto forest = forest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest );
   
   const uint_t level = forest->getLevel( *block );
   CellInterval domainBB = forest->getDomainCellBB( level );
   forest->transformGlobalToBlockLocalCellInterval( domainBB, *block );
   domainBB.expand( cell_idx_t(1) );

   // no slip WEST
   CellInterval west( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( NoSlip_Flag, west );

   // no slip EAST
   CellInterval east( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( NoSlip_Flag, east );

   // no slip SOUTH
   CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
   handling->forceBoundary( NoSlip_Flag, south );

   // no slip NORTH
   CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( NoSlip_Flag, north );

   // no slip BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   handling->forceBoundary( NoSlip_Flag, bottom );

   // velocity bounce back TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( UBB_Flag, top );

   handling->fillWithDomain( domainBB );

   return handling;
}






/////////
// VTK //
/////////

template< typename LatticeModel_T >
class MyVTKOutput {

public:

   MyVTKOutput( const ConstBlockDataID & pdfField, const ConstBlockDataID & flagField,
                const vtk::VTKOutput::BeforeFunction& pdfGhostLayerSync ) :
      pdfField_( pdfField ), flagField_( flagField ), pdfGhostLayerSync_( pdfGhostLayerSync ) {}

   void operator()( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> > & writers,
                    std::map< std::string, vtk::VTKOutput::CellFilter > &          filters,
                    std::map< std::string, vtk::VTKOutput::BeforeFunction > &      beforeFunctions );

private:

   const ConstBlockDataID pdfField_;
   const ConstBlockDataID flagField_;

   vtk::VTKOutput::BeforeFunction pdfGhostLayerSync_;

}; // class MyVTKOutput

template< typename LatticeModel_T >
void MyVTKOutput<LatticeModel_T>::operator()( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> > & writers,
                                              std::map< std::string, vtk::VTKOutput::CellFilter > &          filters,
                                              std::map< std::string, vtk::VTKOutput::BeforeFunction > &      beforeFunctions )
{
   // block data writers

   writers.push_back( make_shared< lbm::VelocityVTKWriter<LatticeModel_T> >( pdfField_, "VelocityFromPDF" ) );
   writers.push_back( make_shared< lbm::DensityVTKWriter<LatticeModel_T> >( pdfField_, "DensityFromPDF" ) );

   writers.push_back( make_shared< field::VTKWriter< FlagField_T > >( flagField_, "FlagField" ) );

   // cell filters

   field::FlagFieldCellFilter<FlagField_T> fluidFilter( flagField_ );
   fluidFilter.addFlag( Fluid_Flag );
   filters[ "FluidFilter" ] = fluidFilter;

   field::FlagFieldCellFilter<FlagField_T> obstacleFilter( flagField_ );
   obstacleFilter.addFlag( NoSlip_Flag );
   obstacleFilter.addFlag(    UBB_Flag );
   filters[ "ObstacleFilter" ] = obstacleFilter;

   // before functions

   beforeFunctions[ "PDFGhostLayerSync" ] = pdfGhostLayerSync_;
}






////////////////////
// THE SIMULATION //
////////////////////

template< typename LatticeModel_T, typename Sweep_T >
void addRefinementTimeStep( SweepTimeloop & timeloop, shared_ptr< blockforest::StructuredBlockForest > & blocks,
                            const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId,
                            const shared_ptr<WcTimingPool> & timingPool, const shared_ptr<WcTimingPool> & levelwiseTimingPool,
                            const bool syncComm, const bool fullComm, const bool linearExplosion,
                            shared_ptr< Sweep_T > & sweep, const std::string & info )
{
   using BH_T = typename MyBoundaryHandling< LatticeModel_T >::BoundaryHandling_T;

   auto ts = lbm::refinement::makeTimeStep< LatticeModel_T, BH_T >( blocks, sweep, pdfFieldId, boundaryHandlingId );
   ts->asynchronousCommunication( !syncComm );
   ts->optimizeCommunication( !fullComm );
   ts->performLinearExplosion( linearExplosion );
   ts->enableTiming( timingPool, levelwiseTimingPool );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor(ts), info );
}

template< typename LatticeModel_T, class Enable = void >
struct AddRefinementTimeStep
{
   static void add( SweepTimeloop & timeloop, shared_ptr< blockforest::StructuredBlockForest > & blocks,
                    const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const BlockDataID & boundaryHandlingId,
                    const shared_ptr<WcTimingPool> & timingPool, const shared_ptr<WcTimingPool> & levelwiseTimingPool,
                    const bool split, const bool pure, const bool syncComm, const bool fullComm, const bool linearExplosion )
   {
      if( split )
      {
         if( pure )
         {
            using Sweep_T = lbm::SplitPureSweep< LatticeModel_T >;
            auto mySweep = make_shared< Sweep_T >( pdfFieldId );

            addRefinementTimeStep< LatticeModel_T, Sweep_T >( timeloop, blocks, pdfFieldId, boundaryHandlingId, timingPool, levelwiseTimingPool,
                                                              syncComm, fullComm, linearExplosion, mySweep, "LBM refinement time step (split pure LB sweep)" );
         }
         else
         {
            using Sweep_T = lbm::SplitSweep<LatticeModel_T, FlagField_T>;
            auto mySweep = make_shared< Sweep_T >( pdfFieldId, flagFieldId, Fluid_Flag );

            addRefinementTimeStep< LatticeModel_T, Sweep_T >( timeloop, blocks, pdfFieldId, boundaryHandlingId, timingPool, levelwiseTimingPool,
                                                              syncComm, fullComm, linearExplosion, mySweep, "LBM refinement time step (split LB sweep)" );
         }
      }
      else
      {
         auto mySweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag );

         addRefinementTimeStep< LatticeModel_T >( timeloop, blocks, pdfFieldId, boundaryHandlingId, timingPool, levelwiseTimingPool,
                                                  syncComm, fullComm, linearExplosion, mySweep, "LBM refinement time step (cell-wise LB sweep)" );
      }
   }
};

template< typename LatticeModel_T  >
struct AddRefinementTimeStep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                           lbm::collision_model::MRT_tag >::value >::type >
{
   static void add( SweepTimeloop & timeloop, shared_ptr< blockforest::StructuredBlockForest > & blocks,
                    const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const BlockDataID & boundaryHandlingId,
                    const shared_ptr<WcTimingPool> & timingPool, const shared_ptr<WcTimingPool> & levelwiseTimingPool,
                    const bool /*split*/, const bool /*pure*/, const bool syncComm, const bool fullComm, const bool linearExplosion )
   {
      auto mySweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag );

      addRefinementTimeStep< LatticeModel_T >( timeloop, blocks, pdfFieldId, boundaryHandlingId, timingPool, levelwiseTimingPool,
                                               syncComm, fullComm, linearExplosion, mySweep, "LBM refinement time step (cell-wise LB sweep)" );
   }
};



template< typename LatticeModel_T >
void run( const shared_ptr< Config > & config, const LatticeModel_T & latticeModel, const bool split, const bool pure,
          const bool fzyx, const bool syncComm, const bool fullComm, const bool linearExplosion )
{
   Config::BlockHandle configBlock = config->getBlock( "NonUniformGrid" );

   // creating the block structure

   auto blocks = createStructuredBlockForest( configBlock );

   // add pdf field to blocks

   BlockDataID pdfFieldId = fzyx ? lbm::addPdfFieldToStorage( blocks, "pdf field (fzyx)", latticeModel,
                                                              Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_t(1),
                                                              FieldGhostLayers, field::fzyx ) :
                                   lbm::addPdfFieldToStorage( blocks, "pdf field (zyxf)", latticeModel,
                                                              Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_t(1),
                                                              FieldGhostLayers, field::zyxf );

   // add flag field to blocks

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field", FieldGhostLayers );

   // add LB boundary handling to blocks

   const real_t velocity = configBlock.getParameter< real_t >( "velocity", real_t(0.05) );

   BlockDataID boundaryHandlingId = blocks->addBlockData( make_shared< MyBoundaryHandling< LatticeModel_T > >( blocks, flagFieldId, pdfFieldId, velocity ),
                                                          "boundary handling" );

   // creating the time loop

   const uint_t outerTimeSteps = configBlock.getParameter< uint_t >( "outerTimeSteps", uint_c(3) );
   const uint_t innerTimeSteps = configBlock.getParameter< uint_t >( "innerTimeSteps", uint_c(4) );
   const uint_t timeSteps = outerTimeSteps * innerTimeSteps;

   SweepTimeloop timeloop( blocks->getBlockStorage(), timeSteps );

   // VTK

   blockforest::communication::NonUniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< lbm::refinement::PdfFieldSyncPackInfo< LatticeModel_T > >( pdfFieldId ) );

   MyVTKOutput< LatticeModel_T > myVTKOutput( pdfFieldId, flagFieldId, pdfGhostLayerSync );

   std::map< std::string, vtk::SelectableOutputFunction > vtkOutputFunctions;
   vtk::initializeVTKOutput( vtkOutputFunctions, myVTKOutput, blocks, config );

   for( auto output = vtkOutputFunctions.begin(); output != vtkOutputFunctions.end(); ++output )
      timeloop.addFuncBeforeTimeStep( output->second.outputFunction, std::string("VTK: ") + output->first,
                                      output->second.requiredGlobalStates, output->second.incompatibleGlobalStates );
                                      
   // add 'refinement' LB time step to time loop

   shared_ptr<WcTimingPool> refinementTimeStepTiming = make_shared<WcTimingPool>();
   shared_ptr<WcTimingPool> refinementTimeStepLevelwiseTiming = make_shared<WcTimingPool>();

   AddRefinementTimeStep< LatticeModel_T >::add( timeloop, blocks, pdfFieldId, flagFieldId, boundaryHandlingId, refinementTimeStepTiming,
                                                 refinementTimeStepLevelwiseTiming, split, pure, syncComm, fullComm, linearExplosion );
                                                 
   // dynamic block structure refresh
   
   auto & blockforest = blocks->getBlockForest();
   
   const uint_t regridAt = configBlock.getParameter< uint_t >( "regridAt", uint_t(0) );
   const int regridType = configBlock.getParameter< int >( "regridType", 0 );
   
   const bool dynamicBlockStructure = ( ( regridAt > uint_t(0) ) && ( regridAt <= timeSteps ) );
   
   const bool allowMultipleRefreshCycles = ( regridType == 1 );
   const bool checkForEarlyOutInRefresh = configBlock.getParameter< bool >( "checkForEarlyOutInRefresh", false );
   const bool checkForLateOutInRefresh = configBlock.getParameter< bool >( "checkForLateOutInRefresh", false );
   
   const int dynamicLoadBalancingType = configBlock.getParameter< int >( "dynamicLoadBalancingType", 0 );
   
   const bool curveHilbert = configBlock.getParameter< bool >( "curveHilbert", false );
   const bool curveAllGather = configBlock.getParameter< bool >( "curveAllGather", true );

   const int diffusionMode = configBlock.getParameter< int >( "diffusionMode", 2 );
   const uint_t diffusionMaxIterations = configBlock.getParameter< uint_t >( "diffusionMaxIterations", uint_t(20) );
   const bool diffusionCheckForEarlyAbort = configBlock.getParameter< bool >( "diffusionCheckForEarlyAbort", true );
   const double diffusionAbortThreshold = configBlock.getParameter< double >( "diffusionAbortThreshold", 1.0 );
   const bool diffusionAdaptOutflow = configBlock.getParameter< bool >( "diffusionAdaptOutflow", true );
   const bool diffusionAdaptInflow = configBlock.getParameter< bool >( "diffusionAdaptInflow", true );
   const uint_t diffusionFlowIterations = configBlock.getParameter< uint_t >( "diffusionFlowIterations", uint_t(10) );
   const uint_t diffusionFlowIterationsIncreaseStart = configBlock.getParameter< uint_t >( "diffusionFlowIterationsIncreaseStart", diffusionMaxIterations );
   const double diffusionFlowIterationsIncreaseFactor = configBlock.getParameter< double >( "diffusionFlowIterationsIncreaseFactor", 0.0 );
   const bool diffusionRegardConnectivity = configBlock.getParameter< bool >( "diffusionRegardConnectivity", true );
   const uint_t diffusionDisregardConnectivityStart = configBlock.getParameter< uint_t >( "diffusionDisregardConnectivityStart", diffusionMaxIterations );
   const double diffusionOutflowExceedFactor = configBlock.getParameter< double >( "diffusionOutflowExceedFactor", 1.0 );
   const double diffusionInflowExceedFactor = configBlock.getParameter< double >( "diffusionInflowExceedFactor", 1.0 );
   
   std::ostringstream loadBalanceLogging;
   std::string refreshFunctorName( "block forest refresh" );
   
   if( dynamicBlockStructure )
   {
      blockforest.recalculateBlockLevelsInRefresh( true );
      blockforest.alwaysRebalanceInRefresh( false );
      blockforest.reevaluateMinTargetLevelsAfterForcedRefinement( false );
      blockforest.allowRefreshChangingDepth( false );
      
      blockforest.allowMultipleRefreshCycles( allowMultipleRefreshCycles );
      blockforest.checkForEarlyOutInRefresh( checkForEarlyOutInRefresh );
      blockforest.checkForLateOutInRefresh( checkForLateOutInRefresh );
         
      ReGrid regrid( blocks, regridType, regridAt );
      
      blockforest.setRefreshMinTargetLevelDeterminationFunction( regrid );
      
      if( dynamicLoadBalancingType == 0 )
      {
         blockforest.setRefreshPhantomBlockMigrationPreparationFunction(
                  blockforest::DynamicCurveBalance< blockforest::NoPhantomData >( curveHilbert, curveAllGather ) );
                  
         loadBalanceLogging << "\n   + type:            " << ( curveHilbert ? "Hilbert order" : "Morton order" ) <<
                               "\n   + parallelization: " << ( curveAllGather? "all gather" : "master-slave" );
      }
      else
      {
         using DLDB = blockforest::DynamicDiffusionBalance<blockforest::NoPhantomData>;
         DLDB balancer( diffusionMaxIterations, diffusionFlowIterations );
         if( diffusionMode == 0 )
            balancer.setMode( DLDB::DIFFUSION_PUSH );
         else if( diffusionMode == 1 )
            balancer.setMode( DLDB::DIFFUSION_PULL );
         else
            balancer.setMode( DLDB::DIFFUSION_PUSHPULL );
         balancer.defineProcessWeightLimitByMultipleOfMaxBlockWeight( true );
         balancer.checkForEarlyAbort( diffusionCheckForEarlyAbort, diffusionAbortThreshold );
         balancer.adaptOutflowWithGlobalInformation( diffusionAdaptOutflow );
         balancer.adaptInflowWithGlobalInformation( diffusionAdaptInflow );
         balancer.setDynamicFlowIterationsIncrease( diffusionFlowIterationsIncreaseStart, diffusionFlowIterationsIncreaseFactor );
         balancer.regardConnectivity( diffusionRegardConnectivity );
         balancer.disregardConnectivity( diffusionDisregardConnectivityStart );
         balancer.setOutflowExceedFactor( diffusionOutflowExceedFactor );
         balancer.setOutflowExceedFactor( diffusionInflowExceedFactor );
         
         blockforest.setRefreshPhantomBlockMigrationPreparationFunction( balancer );
                  
         loadBalanceLogging << "\n   + mode:                                " << ( (diffusionMode == 0) ? "push" : ( (diffusionMode == 1) ? "pull" : "push+pull" ) ) <<
                               "\n   + max iterations:                      " << diffusionMaxIterations <<
                               "\n   + check for early abort:               " << ( diffusionCheckForEarlyAbort ? "yes" : "no" );
         if( diffusionCheckForEarlyAbort )
            loadBalanceLogging << " (abort threshold = " << diffusionAbortThreshold << ")";
         if( diffusionMode != 1 )
            loadBalanceLogging << "\n   + adapt outflow with global knowledge: " << ( diffusionAdaptOutflow ? "yes" : "no" );
         if( diffusionMode != 0 )
            loadBalanceLogging << "\n   + adapt inflow with global knowledge:  " << ( diffusionAdaptInflow ? "yes" : "no" );
         loadBalanceLogging << "\n   + flow iterations:                     " << diffusionFlowIterations <<
                               "\n   + dynamic flow iterations increase:    " << ( ( diffusionFlowIterationsIncreaseStart < diffusionMaxIterations ) ? "yes" : "no" );
         if( diffusionFlowIterationsIncreaseStart < diffusionMaxIterations )
            loadBalanceLogging << " (starting at main iteration " << diffusionFlowIterationsIncreaseStart << " with factor " << diffusionFlowIterationsIncreaseFactor << ")";
         loadBalanceLogging << "\n   + take into account connectivity:      " << ( diffusionRegardConnectivity ? "yes" : "no" );
         if( diffusionRegardConnectivity && diffusionDisregardConnectivityStart < diffusionMaxIterations )
            loadBalanceLogging << " (disregard starting at main iteration " << diffusionDisregardConnectivityStart << ")";
         if( diffusionMode != 1 )
            loadBalanceLogging << "\n   + outflow exceed factor:               " << diffusionOutflowExceedFactor;
         if( diffusionMode != 0 )
            loadBalanceLogging << "\n   + inflow exceed factor:                " << diffusionInflowExceedFactor;
      }
      
      timeloop.addFuncBeforeTimeStep( regrid, refreshFunctorName );
   }

   // logging right before the benchmark starts

   uint_t lbBlockForestEvaluationStamp = blockforest.getModificationStamp();
   lbm::BlockForestEvaluation< FlagField_T > lbBlockForestEvaluation( blocks, flagFieldId, Fluid_Flag );
   lbBlockForestEvaluation.logInfoOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark parameters:"
                              "\n- collision model:  " << CollisionModelString< LatticeModel_T >::str() <<
                              "\n- stencil:          " << StencilString< LatticeModel_T >::str() <<
                              "\n- compressible:     " << ( LatticeModel_T::compressible ? "yes" : "no" ) <<
                              "\n- split kernel:     " << ( split ? "yes" : "no" ) <<
                              "\n- pure kernel:      " << ( pure ? "yes (collision is also performed within obstacle cells)" : "no" ) <<
                              "\n- data layout:      " << ( fzyx ? "fzyx (structure of arrays [SoA])" : "zyxf (array of structures [AoS])" ) <<
                              "\n- communication:    " << ( fullComm ? ( syncComm ? "synchronous, full synchronization" : "asynchronous, full synchronization" ) :
                                                                       ( syncComm ? "synchronous, block neighborhood and direction-aware optimizations" : "asynchronous, block neighborhood and direction-aware optimizations" ) ) <<
                              "\n- linear explosion: " << ( linearExplosion ? "yes" : "no" ) );

   if( dynamicBlockStructure )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Dynamic block structure refresh:"
                                 "\n- refresh at time step:          " << regridAt <<
                                 "\n- refresh type:                  " << ( ( regridType == 0 ) ? "shift finest level inwards" : "move finest level to bottom edges" ) <<
                                 "\n- allow multiple refresh cycles: " << ( allowMultipleRefreshCycles ? "yes" : "no" ) <<
                                 "\n- check for early refresh out:   " << ( checkForEarlyOutInRefresh ? "yes" : "no" ) <<
                                 "\n- check for late refresh out:    " << ( checkForLateOutInRefresh ? "yes" : "no" ) <<
                                 "\n- load balancing algorithm:      " << ( ( dynamicLoadBalancingType == 0 ) ? "space filling curve" : "distributed diffusion" ) <<
                                 loadBalanceLogging.str() );
   }
   
   // run the benchmark

   uint_t performanceEvaluationStamp = blockforest.getModificationStamp();
   lbm::PerformanceEvaluation< FlagField_T > performance( blocks, flagFieldId, Fluid_Flag );

   for( uint_t outerRun = 0; outerRun < outerTimeSteps; ++outerRun )
   {
      const uint_t blockForestStampBeforeInnerLoop = blockforest.getModificationStamp();
      
      WcTimingPool timeloopTiming;

      WALBERLA_MPI_WORLD_BARRIER();
      WcTimer timer;
      timer.start();

      for( uint_t innerRun = 0; innerRun < innerTimeSteps; ++innerRun )
         timeloop.singleStep( timeloopTiming );

      timer.end();

      double time = timer.max();
      mpi::reduceInplace( time, mpi::MAX );
      
      const bool blockStructureRefreshDuringMeasurement = ( blockForestStampBeforeInnerLoop != blockforest.getModificationStamp() );

      const auto reducedTimeloopTiming = timeloopTiming.getReduced();
      const auto reducedRTSTiming  = refinementTimeStepTiming->getReduced();
      const auto reducedRTSLTiming = refinementTimeStepLevelwiseTiming->getReduced();
      refinementTimeStepTiming->clear();
      refinementTimeStepLevelwiseTiming->clear();

      WALBERLA_LOG_RESULT_ON_ROOT( "Time loop timing:\n" << *reducedTimeloopTiming );
      WALBERLA_LOG_RESULT_ON_ROOT( "Refinement time step timing:\n" << *reducedRTSTiming );
      WALBERLA_LOG_RESULT_ON_ROOT( "Refinement time step timing (one timer per level):\n" << *reducedRTSLTiming );

      shared_ptr< WcTimingPool > reducedRefreshTiming;
      
      double refreshPerformance( 0.0 );
      double blockLevelDeterminationPerformance( 0.0 );
      double dataMigrationPerformance( 0.0 );
      double loadBalancePerformance( 0.0 );
      double phantomForestCreationPerformance( 0.0 );
      
      if( blockStructureRefreshDuringMeasurement )
      {
         reducedRefreshTiming = blocks->getBlockForest().getRefreshTiming().getReduced();
         blocks->getBlockForest().clearRefreshTiming();
         WALBERLA_LOG_RESULT_ON_ROOT( "Block structure refresh timing:\n" << *reducedRefreshTiming );
         
         WALBERLA_ROOT_SECTION()
         {
            double restTime( 0.0 );
            double refreshTime( 0.0 );
            for( auto it = reducedTimeloopTiming->begin(); it != reducedTimeloopTiming->end(); ++it )
            {
               if( it->first == refreshFunctorName )
                  refreshTime += it->second.total();
               else
                  restTime += it->second.total();
            }
            const double timeStepTime = restTime / double_c(innerTimeSteps);
            
            double blockLevelDeterminationTime( 0.0 );
            double dataMigrationTime( 0.0 );
            double loadBalanceTime( 0.0 );
            double phantomForestCreationTime( 0.0 );
            
            for( auto it = reducedRefreshTiming->begin(); it != reducedRefreshTiming->end(); ++it )
            {
               if( it->first == "block level determination" )
                  blockLevelDeterminationTime += it->second.total();
               else if( it->first == "block level determination (callback function)" )
                  blockLevelDeterminationTime += it->second.total();
               else if( it->first == "block structure update (includes data migration)" )
                  dataMigrationTime += it->second.total();
               else if( it->first == "phantom block redistribution (= load balancing)" )
                  loadBalanceTime += it->second.total();
               else if( it->first == "phantom forest creation" )
                  phantomForestCreationTime += it->second.total();
            }
            
            refreshPerformance = refreshTime / timeStepTime;
            blockLevelDeterminationPerformance = blockLevelDeterminationTime / timeStepTime;
            dataMigrationPerformance = dataMigrationTime / timeStepTime;
            loadBalancePerformance = loadBalanceTime / timeStepTime;
            phantomForestCreationPerformance = phantomForestCreationTime / timeStepTime;
            
            WALBERLA_LOG_RESULT_ON_ROOT( "Block structure refresh performance (in 'coarse time steps'):" <<
                                         "\n- total refresh:             " << refreshPerformance <<
                                         "\n- block level determination: " << blockLevelDeterminationPerformance <<
                                         "\n- phantom forest setup:      " << phantomForestCreationPerformance <<
                                         "\n- dynamic load balancing:    " << loadBalancePerformance <<
                                         "\n- data migration:            " << dataMigrationPerformance );
         }
      }
      
      if( performanceEvaluationStamp != blockforest.getModificationStamp() )
      {
         performanceEvaluationStamp = blockforest.getModificationStamp();
         performance.refresh();
         WALBERLA_LOG_WARNING_ON_ROOT( "ATTENTION: The following performance statistics may not be entirely correct since the block structure did change during the last time measurement!" );
      }
      performance.logResultOnRoot( innerTimeSteps, time );
      
      if( configBlock.getParameter< bool >( "logToSqlDB", true ) &&
          lbBlockForestEvaluationStamp != blockforest.getModificationStamp() )
      {
         lbBlockForestEvaluationStamp = blockforest.getModificationStamp();
         lbBlockForestEvaluation.refresh();
      }

      WALBERLA_ROOT_SECTION()
      {
         // logging in SQL database

         if( configBlock.getParameter< bool >( "logToSqlDB", true ) )
         {
            const std::string sqlFile = configBlock.getParameter< std::string >( "sqlFile", "performance.sqlite" );

            std::map< std::string, int >         integerProperties;
            std::map< std::string, double >      realProperties;
            std::map< std::string, std::string > stringProperties;
            
            stringProperties[ "blockStructureRefreshDuringMeasurement" ] = ( blockStructureRefreshDuringMeasurement ? "yes" : "no" );

            performance.getResultsForSQLOnRoot( integerProperties, realProperties, stringProperties, innerTimeSteps, time );
            lbBlockForestEvaluation.getResultsForSQLOnRoot( integerProperties, realProperties, stringProperties );

            stringProperties[ "collisionModel" ]    = CollisionModelString< LatticeModel_T >::str();
            stringProperties[ "stencil" ]           = StencilString< LatticeModel_T >::str();
            stringProperties[ "compressible" ]      = ( LatticeModel_T::compressible ? "yes" : "no" );
            stringProperties[ "splitKernel" ]       = ( split ? "yes" : "no" );
            stringProperties[ "pureKernel" ]        = ( pure ? "yes" : "no" );
            stringProperties[ "dataLayout" ]        = ( fzyx ? "fzyx" : "zyxf" );
            stringProperties[ "syncCommunication" ] = ( syncComm ? "yes" : "no" );
            stringProperties[ "fullCommunication" ] = ( fullComm ? "yes" : "no" );
            stringProperties[ "linearExplosion" ]   = ( linearExplosion ? "yes" : "no" );
            
            stringProperties[ "dynamicBlockStructure" ] = ( dynamicBlockStructure ? "yes" : "no" );
            if( dynamicBlockStructure )
            {
               integerProperties[ "regridAt" ] = int_c( regridAt );
               stringProperties[ "regridType" ] = ( ( regridType == 0 ) ? "shift finest level inwards" : "move finest level to bottom edges" );
               stringProperties[ "allowMultipleRefreshCycles" ] = ( allowMultipleRefreshCycles ? "yes" : "no" );
               stringProperties[ "checkForEarlyOutInRefresh" ] = ( checkForEarlyOutInRefresh ? "yes" : "no" );
               stringProperties[ "checkForLateOutInRefresh" ] = ( checkForLateOutInRefresh ? "yes" : "no" );
               stringProperties[ "loadBalancingAlgorithm" ] = ( ( dynamicLoadBalancingType == 0 ) ? "space filling curve" : "distributed diffusion" );
               if( dynamicLoadBalancingType == 0 )
               {
                  stringProperties[ "curveType" ] = ( curveHilbert ? "Hilbert order" : "Morton order" );
                  stringProperties[ "curveParallelization" ] = ( curveAllGather? "all gather" : "master-slave" );
               }
               else
               {
                  integerProperties[ "diffusionMode" ] = diffusionMode;
                  integerProperties[ "diffusionMaxIterations" ] = int_c( diffusionMaxIterations );
                  integerProperties[ "diffusionIterations" ] = int_c( blockforest.phantomBlockMigrationIterations() );
                  stringProperties[ "diffusionCheckForEarlyAbort" ] = ( diffusionCheckForEarlyAbort ? "yes" : "no" );
                  realProperties[ "diffusionAbortThreshold" ] = diffusionAbortThreshold;
                  stringProperties[ "diffusionAdaptOutflow" ] = ( diffusionAdaptOutflow ? "yes" : "no" );
                  stringProperties[ "diffusionAdaptInflow" ] = ( diffusionAdaptInflow ? "yes" : "no" );
                  integerProperties[ "diffusionFlowIterations" ] = int_c( diffusionFlowIterations );
                  integerProperties[ "diffusionFlowIterationsIncreaseStart" ] = int_c( diffusionFlowIterationsIncreaseStart );
                  realProperties[ "diffusionFlowIterationsIncreaseFactor" ] = diffusionFlowIterationsIncreaseFactor;
                  stringProperties[ "diffusionRegardConnectivity" ] = ( diffusionRegardConnectivity ? "yes" : "no" );
                  integerProperties[ "diffusionDisregardConnectivityStart" ] = int_c( diffusionDisregardConnectivityStart );
                  realProperties[ "diffusionOutflowExceedFactor" ] = diffusionOutflowExceedFactor;
                  realProperties[ "diffusionInflowExceedFactor" ] = diffusionInflowExceedFactor;
               }
               if( blockStructureRefreshDuringMeasurement )
               {
                  realProperties[ "refreshPerformance" ] = refreshPerformance;
                  realProperties[ "blockLevelDeterminationPerformance" ] = blockLevelDeterminationPerformance;
                  realProperties[ "dataMigrationPerformance" ] = dataMigrationPerformance;
                  realProperties[ "loadBalancePerformance" ] = loadBalancePerformance;
                  realProperties[ "phantomForestCreationPerformance" ] = phantomForestCreationPerformance;
               }
            }

            auto runId = sqlite::storeRunInSqliteDB( sqlFile, integerProperties, stringProperties, realProperties );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedTimeloopTiming, "Timeloop" );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedRTSTiming, "RefinementTimeStep" );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedRTSLTiming, "RefinementTimeStepLevelwise" );
            if( blockStructureRefreshDuringMeasurement )
               sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedRefreshTiming, "BlockForestRefresh" );
         }
      }
   }

   // logging once again at the end of the simulation, identical to logging at the beginning :-)

   if( lbBlockForestEvaluationStamp != blockforest.getModificationStamp() )
      lbBlockForestEvaluation.refresh();
   lbBlockForestEvaluation.logInfoOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark parameters:"
                              "\n- collision model:  " << CollisionModelString< LatticeModel_T >::str() <<
                              "\n- stencil:          " << StencilString< LatticeModel_T >::str() <<
                              "\n- compressible:     " << ( LatticeModel_T::compressible ? "yes" : "no" ) <<
                              "\n- split kernel:     " << ( split ? "yes" : "no" ) <<
                              "\n- pure kernel:      " << ( pure ? "yes (collision is also performed within obstacle cells)" : "no" ) <<
                              "\n- data layout:      " << ( fzyx ? "fzyx (structure of arrays [SoA])" : "zyxf (array of structures [AoS])" ) <<
                              "\n- communication:    " << ( fullComm ? ( syncComm ? "synchronous, full synchronization" : "asynchronous, full synchronization" ) :
                                                                       ( syncComm ? "synchronous, block neighborhood and direction-aware optimizations" : "asynchronous, block neighborhood and direction-aware optimizations" ) ) <<
                              "\n- linear explosion: " << ( linearExplosion ? "yes" : "no" ) );
                              
   if( dynamicBlockStructure )
   {
      if( dynamicLoadBalancingType != 0 )
      {
         loadBalanceLogging.str("");                               
         loadBalanceLogging << "\n   + mode:                                " << ( (diffusionMode == 0) ? "push" : ( (diffusionMode == 1) ? "pull" : "push+pull" ) ) <<
                               "\n   + max iterations:                      " << diffusionMaxIterations << " (performed iterations: " << blockforest.phantomBlockMigrationIterations() << ")" <<
                               "\n   + check for early abort:               " << ( diffusionCheckForEarlyAbort ? "yes" : "no" );
         if( diffusionCheckForEarlyAbort )
            loadBalanceLogging << " (abort threshold = " << diffusionAbortThreshold << ")";
         if( diffusionMode != 1 )
            loadBalanceLogging << "\n   + adapt outflow with global knowledge: " << ( diffusionAdaptOutflow ? "yes" : "no" );
         if( diffusionMode != 0 )
            loadBalanceLogging << "\n   + adapt inflow with global knowledge:  " << ( diffusionAdaptInflow ? "yes" : "no" );
         loadBalanceLogging << "\n   + flow iterations:                     " << diffusionFlowIterations <<
                               "\n   + dynamic flow iterations increase:    " << ( ( diffusionFlowIterationsIncreaseStart < diffusionMaxIterations ) ? "yes" : "no" );
         if( diffusionFlowIterationsIncreaseStart < diffusionMaxIterations )
            loadBalanceLogging << " (starting at main iteration " << diffusionFlowIterationsIncreaseStart << " with factor " << diffusionFlowIterationsIncreaseFactor << ")";
         loadBalanceLogging << "\n   + take into account connectivity:      " << ( diffusionRegardConnectivity ? "yes" : "no" );
         if( diffusionRegardConnectivity && diffusionDisregardConnectivityStart < diffusionMaxIterations )
            loadBalanceLogging << " (disregard starting at main iteration " << diffusionDisregardConnectivityStart << ")";
         if( diffusionMode != 1 )
            loadBalanceLogging << "\n   + outflow exceed factor:               " << diffusionOutflowExceedFactor;
         if( diffusionMode != 0 )
            loadBalanceLogging << "\n   + inflow exceed factor:                " << diffusionInflowExceedFactor;
      }
      
      WALBERLA_LOG_INFO_ON_ROOT( "Dynamic block structure refresh:"
                                 "\n- refresh at time step:          " << regridAt <<
                                 "\n- refresh type:                  " << ( ( regridType == 0 ) ? "shift finest level inwards" : "move finest level to bottom edges" ) <<
                                 "\n- allow multiple refresh cycles: " << ( allowMultipleRefreshCycles ? "yes" : "no" ) <<
                                 "\n- check for early refresh out:   " << ( checkForEarlyOutInRefresh ? "yes" : "no" ) <<
                                 "\n- check for late refresh out:    " << ( checkForLateOutInRefresh ? "yes" : "no" ) <<
                                 "\n- load balancing algorithm:      " << ( ( dynamicLoadBalancingType == 0 ) ? "space filling curve" : "distributed diffusion" ) <<
                                 loadBalanceLogging.str() );
   }
}






//////////
// MAIN //
//////////

enum CM { CMSRT, CMTRT, CMMRT };

int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );

   if( argc < 2 )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::cout << "Usage: " << argv[0] << " path-to-configuration-file [--trt | --mrt] [--comp] [--split [--pure]] [--fzyx] [--sync-comm] [--full-comm] [--linear-exp]\n"
                      "\n"
                      "By default, SRT is selected as collision model, an asynchronous communication scheme with block neighborhood and\n"
                      "direction-aware optimizations is chosen, and an incompressible, basic LB kernel is executed on a PDF field with\n"
                      "layout 'zyxf' (= array of structures [AoS]).\n"
                      "\n"
                      "Optional arguments:\n"
                      " --trt:        collision model = TRT\n"
                      " --mrt:        collision model = MRT\n"
                      " --comp:       LB kernel is switched from incompressible to compressible\n"
                      " --split:      LB kernel split by PDF direction\n"
                      "               Should always be combined with --fzyx.\n"
                      " --pure:       LB kernel is executed in every cell (including obstacle/boundary cells)\n"
                      "               Only available in combination with --split.\n"
                      " --fzyx:       data layout switched to 'fzyx' (structure of arrays [SoA])\n"
                      " --sync-comm:  A synchronous communication scheme is used instead of an asynchronous scheme\n"
                      "               which is used by default.\n"
                      " --full-comm:  A full synchronization of neighboring blocks is performed instead of using a communication\n"
                      "               that uses block neighborhood and direction-aware optimizations.\n"
                      " --linear-exp: When communicating from coarse to fine grids, a linear interpolation scheme is used\n"
                      "               instead of a uniform distribution of a coarse cell to eight fine cells.\n"
                      "\n"
                      "Please note: Depending on the underlying hardware and the configuration of the benchmark (more precisely: the number of cells\n"
                      "             in each block), the best performance may be achieved with split, pure kernels combined with a structure of arrays\n"
                      "             ('fzyx') data layout!" << std::endl;
      }
      return EXIT_SUCCESS;
   }

   logging::Logging::printHeaderOnStream();

#ifdef _OPENMP
   if( std::getenv( "OMP_NUM_THREADS" ) == nullptr )
      WALBERLA_ABORT( "If you are using a version of the benchmark that was compiled with OpenMP you have to "
                      "specify the environment variable \'OMP_NUM_THREADS\' accordingly!" );
#endif

   // open configuration file

   shared_ptr< Config > config = make_shared< Config >();
   config->readParameterFile( argv[1] );

   Config::BlockHandle configBlock = config->getBlock( "NonUniformGrid" );

   if( !configBlock )
      WALBERLA_ABORT( "You have to specify a \"NonUniformGrid\" block in the configuration file!" );
      
   const uint_t bufferProcesses = configBlock.getParameter< uint_t >( "bufferProcesses", 0 );
   const uint_t fineBlocksPerProcess = configBlock.getParameter< uint_t >( "fineBlocksPerProcess", 4 );
   if( fineBlocksPerProcess != uint_t(1) && fineBlocksPerProcess != uint_t(2) && fineBlocksPerProcess != uint_t(4) && fineBlocksPerProcess != uint_t(8) )
      WALBERLA_ABORT( "The number of fine blocks per process (\"fineBlocksPerProcess\") must be equal to either 1, 2, 4, or 8! (You requested " <<
                      fineBlocksPerProcess << ")" );
   
   const bool progressLoggingOnRoot = configBlock.getParameter< bool >( "progressLoggingOnRoot", false );
   if( progressLoggingOnRoot )
   {
      WALBERLA_ROOT_SECTION() { logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS ); }
   }
   const bool includeLoggingToFile = configBlock.getParameter< bool >( "includeLoggingToFile", false );
   if( includeLoggingToFile )
      logging::Logging::instance()->includeLoggingToFile( "log" );

   // In case 'sbffile' and 'processes' are specified in the configuration file:
   // -> just create the block structure and save it to file for later use

   if( configBlock.isDefined( "sbffile" ) && configBlock.isDefined( "processes" ) )
   {
      std::string  sbffile           = configBlock.getParameter< std::string >( "sbffile" );
      const uint_t numberOfProcesses = configBlock.getParameter< uint_t >( "processes" );

      std::ostringstream infoString;
      infoString << "You have selected the option of just creating the block structure (= domain decomposition) and saving the result to file\n"
                    "by specifying the output file name \'" << sbffile << "\' AND providing a targeted number of processes ("
                 << numberOfProcesses << ").\n";

      if( MPIManager::instance()->numProcesses() > 1 )
         WALBERLA_ABORT( infoString.str() << "In this mode you need to start " << argv[0] << " with just one process!" );

      WALBERLA_CHECK_LESS( bufferProcesses, numberOfProcesses );
      
      const uint_t workerProcesses = numberOfProcesses - bufferProcesses;
      if( workerProcesses < (uint_t(64) / fineBlocksPerProcess) || !( ( workerProcesses % (uint_t(64) / fineBlocksPerProcess) ) == uint_t(0) ) )
         WALBERLA_ABORT( infoString.str() << "You selected " << fineBlocksPerProcess << " fine blocks per process -> The number of worker processes must be divisible by " <<
                         (uint_t(64) / fineBlocksPerProcess) << "!" );

      WALBERLA_LOG_INFO_ON_ROOT( infoString.str() << "Creating the block structure ..." );

      blockforest::SetupBlockForest sforest;
      createSetupBlockForest( sforest, configBlock, numberOfProcesses );

      sforest.saveToFile( sbffile.c_str() );

      logging::Logging::printFooterOnStream();
      return EXIT_SUCCESS;
   }

   WALBERLA_CHECK_LESS( bufferProcesses, MPIManager::instance()->numProcesses() );
      
   const uint_t workerProcesses = uint_c( MPIManager::instance()->numProcesses() ) - bufferProcesses;
   if( workerProcesses < (uint_t(64) / fineBlocksPerProcess ) || !( ( workerProcesses % (uint_t(64) / fineBlocksPerProcess ) ) == uint_t(0) ) )
      WALBERLA_ABORT( "You selected " << fineBlocksPerProcess << " fine blocks per process -> The number of worker processes must be divisible by " <<
                      (64 / int_c(fineBlocksPerProcess) ) << "!\n(You requested " << workerProcesses << " worker processes...)" );

   // reading optional parameters from passed arguments

   CM collisionModel    = CMSRT;
   bool compressible    = false;
   bool split           = false;
   bool pure            = false;
   bool fzyx            = false;
   bool syncComm        = false;
   bool fullComm        = false;
   bool linearExplosion = false;

   for( int i = 2; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--trt" )        == 0 ) collisionModel  = CMTRT;
      if( std::strcmp( argv[i], "--mrt" )        == 0 ) collisionModel  = CMMRT;
      if( std::strcmp( argv[i], "--comp" )       == 0 ) compressible    = true;
      if( std::strcmp( argv[i], "--split" )      == 0 ) split           = true;
      if( std::strcmp( argv[i], "--pure" )       == 0 ) pure            = true;
      if( std::strcmp( argv[i], "--fzyx" )       == 0 ) fzyx            = true;
      if( std::strcmp( argv[i], "--sync-comm" )  == 0 ) syncComm        = true;
      if( std::strcmp( argv[i], "--full-comm" )  == 0 ) fullComm        = true;
      if( std::strcmp( argv[i], "--linear-exp" ) == 0 ) linearExplosion = true;
   }

   if( pure && !split )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "You called the benchmark with \"--pure\" but without \"--split\".\n"
                                    "\"Pure\" kernels are only available for \"split\" kernels! Setting \"pure\" to false ..." );
      pure = pure && split; // pure only works in combination with split
   }

   if( collisionModel == CMMRT && ( compressible || split || pure ) )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Options \"--comp\", \"--split\", and \"--pure\" are not available for MRT!\n"
                                    "Setting \"compressible\", \"split\", and \"pure\" to false ..." );
      compressible = false;
      split        = false;
      pure         = false;
   }

   const real_t omega = configBlock.getParameter< real_t >( "omega", real_t(1.4) ); // on the coarsest grid!

   // executing benchmark

   if( collisionModel == CMSRT ) // SRT
   {
      if( compressible )
      {
         D3Q19_SRT_COMP latticeModel = D3Q19_SRT_COMP( lbm::collision_model::SRT( omega ) );
         run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion );
      }
      else
      {
         D3Q19_SRT_INCOMP latticeModel = D3Q19_SRT_INCOMP( lbm::collision_model::SRT( omega ) );
         run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion );
      }
   }
   else if( collisionModel == CMTRT ) // TRT
   {
      if( compressible )
      {
         D3Q19_TRT_COMP latticeModel = D3Q19_TRT_COMP( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
         run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion );
      }
      else
      {
         D3Q19_TRT_INCOMP latticeModel = D3Q19_TRT_INCOMP( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
         run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion );
      }
   }
   else // MRT
   {
      D3Q19_MRT_INCOMP latticeModel = D3Q19_MRT_INCOMP( lbm::collision_model::D3Q19MRT::constructTRTWithMagicNumber( omega ) );
      run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion );
   }

   logging::Logging::printFooterOnStream();

   return EXIT_SUCCESS;
}

} // namespace non_uniform_grid

int main( int argc, char ** argv )
{
	return non_uniform_grid::main( argc, argv );
}
