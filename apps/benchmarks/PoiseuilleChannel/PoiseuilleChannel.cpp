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
//! \file PoiseuilleChannel.cpp
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/AABBRefinementSelection.h"
#include "blockforest/communication/NonUniformBufferedScheme.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
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
#include "core/math/Constants.h"
#include "core/math/Sample.h"
#include "core/math/Vector3.h"
#include "core/mpi/Environment.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/SendBuffer.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "field/AccuracyEvaluation.h"
#include "field/AccuracyEvaluationLinePlot.h"
#include "field/AddToStorage.h"
#include "field/CellCounter.h"
#include "field/FlagField.h"
#include "field/FlagUID.h"
#include "field/StabilityChecker.h"
#include "field/VolumetricFlowRateEvaluation.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/iterators/FieldIterator.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"

#include "lbm/BlockForestEvaluation.h"
#include "lbm/MassEvaluation.h"
#include "lbm/PerformanceEvaluation.h"
#include "lbm/boundary/Curved.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"
#include "lbm/refinement/BoundarySetup.h"
#include "lbm/refinement/PdfFieldSyncPackInfo.h"
#include "lbm/refinement/TimeStep.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/NonEquilibrium.h"
#include "lbm/vtk/Velocity.h"

#include "sqlite/SQLite.h"

#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/BlockCellDataWriter.h"
#include "vtk/Initialization.h"
#include "vtk/VTKOutput.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>



namespace poiseuille_channel {

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;
using walberla::real_t;

//////////////
// TYPEDEFS //
//////////////

typedef lbm::D3Q19< lbm::collision_model::SRT, false, lbm::force_model::SimpleConstant > D3Q19_SRT_INCOMP;
typedef lbm::D3Q19< lbm::collision_model::SRT, true,  lbm::force_model::SimpleConstant > D3Q19_SRT_COMP;
typedef lbm::D3Q19< lbm::collision_model::TRT, false, lbm::force_model::SimpleConstant > D3Q19_TRT_INCOMP;
//typedef lbm::D3Q19< lbm::collision_model::TRT, true,  lbm::force_model::SimpleConstant > D3Q19_TRT_COMP;

typedef lbm::D3Q27< lbm::collision_model::SRT, false, lbm::force_model::SimpleConstant > D3Q27_SRT_INCOMP;
typedef lbm::D3Q27< lbm::collision_model::SRT, true,  lbm::force_model::SimpleConstant > D3Q27_SRT_COMP;
typedef lbm::D3Q27< lbm::collision_model::TRT, false, lbm::force_model::SimpleConstant > D3Q27_TRT_INCOMP;
typedef lbm::D3Q27< lbm::collision_model::TRT, true,  lbm::force_model::SimpleConstant > D3Q27_TRT_COMP;

template< typename LatticeModel_T >
struct Types
{
   using Stencil_T = typename LatticeModel_T::Stencil;
   using PdfField_T = lbm::PdfField<LatticeModel_T>;
};

using flag_t = walberla::uint16_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers  = uint_t(4);

///////////
// FLAGS //
///////////

const FlagUID  Fluid_Flag( "fluid" );
const FlagUID NoSlip_Flag( "no slip" );
const FlagUID Curved_Flag( "curved" );

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

template< typename LatticeModel_T >
struct StencilString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value >::type >
{
   static const char * str() { return "D3Q27"; }
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

//////////////////////
// Parameter Struct //
//////////////////////

struct Setup
{
   uint_t xBlocks;
   uint_t yBlocks;
   uint_t zBlocks;

   uint_t xCells;
   uint_t yCells;
   uint_t zCells;

   real_t Re;
   
   bool circularProfile;
   real_t radius_L; // on the coarsest grid

   real_t viscosity_L; // on the coarsest grid
   real_t acceleration_L; // on the coarsest grid
   real_t maxVelocity_L;
   real_t meanVelocity_L;
   real_t flowRate_L; // volumetric flow rate - in lattice units of the coarsest grid
};






/////////////////////
// BLOCK STRUCTURE //
/////////////////////

class BorderRefinementSelection
{
public:

   BorderRefinementSelection( const Setup & setup, const uint_t level, const real_t bufferDistance ) :
      setup_( setup ), level_( level ), bufferDistance_( bufferDistance ) {}

   void operator()( SetupBlockForest & forest )
   {
      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         if( block->getLevel() < level_ && !channelContains( forest, *block ) && !outsideChannel( forest, *block ) )
            block->setMarker( true );
      }
   }

private:

   bool channelContains( SetupBlockForest & forest, const SetupBlock & block )
   {
      if( setup_.circularProfile )
      {
         auto corners = block.getAABB().corners();
         for( auto p = corners.begin(); p != corners.end(); ++p )
         {
            const real_t dy = (*p)[1] - forest.getDomain().center()[1];
            const real_t dz = (*p)[2] - forest.getDomain().center()[2];
            const real_t r = setup_.radius_L - bufferDistance_;
            if( (dy * dy + dz * dz) >= (r * r) )
               return false;
         }

         return true;
      }

      AABB domain = forest.getDomain();
      domain.setAxisBounds( uint_t(1), domain.yMin() + bufferDistance_, domain.yMax() - bufferDistance_ );

      return domain.contains( block.getAABB() ) && !(forest.atDomainYMinBorder( block )) && !(forest.atDomainYMaxBorder( block ));
   }

   bool outsideChannel( SetupBlockForest & forest, const SetupBlock & block )
   {
      if( setup_.circularProfile )
      {
         const auto d = block.getAABB().sqSignedDistance( forest.getDomain().center() );
         const real_t r = setup_.radius_L + bufferDistance_;
         if( d > (r * r) )
            return true;
         return false;
      }

      return false;
   }

   const Setup & setup_;
   uint_t level_;
   real_t bufferDistance_;

}; // class BorderRefinementSelection



static void workloadAndMemoryAssignment( SetupBlockForest& forest, const memory_t memoryPerBlock )
{
   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      block->setWorkload( numeric_cast< workload_t >( uint_t(1) << block->getLevel() ) );
      block->setMemory( memoryPerBlock );
   }
}



static shared_ptr< SetupBlockForest > createSetupBlockForest( const blockforest::RefinementSelectionFunctions & refinementSelectionFunctions,
                                                              const Setup & setup,
                                                              uint_t numberOfProcesses, const uint_t blocksPerProcess,
                                                              const memory_t memoryPerCell, const memory_t processMemoryLimit,
                                                              const bool outputSetupForest )
{
   shared_ptr< SetupBlockForest > forest = make_shared< SetupBlockForest >();

   const memory_t memoryPerBlock = numeric_cast< memory_t >( ( setup.xCells + uint_t(2) * FieldGhostLayers ) *
                                                             ( setup.yCells + uint_t(2) * FieldGhostLayers ) *
                                                             ( setup.zCells + uint_t(2) * FieldGhostLayers ) ) * memoryPerCell;

   forest->addRefinementSelectionFunction( refinementSelectionFunctions );
   forest->addWorkloadMemorySUIDAssignmentFunction( std::bind( workloadAndMemoryAssignment, std::placeholders::_1, memoryPerBlock ) );

   forest->init( AABB( real_c(0), real_c(0), real_c(0), real_c( setup.xBlocks * setup.xCells ),
                                                        real_c( setup.yBlocks * setup.yCells ),
                                                        real_c( setup.zBlocks * setup.zCells ) ),
                 setup.xBlocks, setup.yBlocks, setup.zBlocks, true, false, !(setup.circularProfile) );

   MPIManager::instance()->useWorldComm();

   if( blocksPerProcess != 0 )
      numberOfProcesses = uint_c( std::ceil( real_c( forest->getNumberOfBlocks() ) / real_c( blocksPerProcess ) ) );

   forest->balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), numberOfProcesses, real_t(0), processMemoryLimit, true );

   if( outputSetupForest ) 
   {
      forest->writeVTKOutput( "domain_decomposition" );
      forest->writeCSV( "process_distribution" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "SetupBlockForest created successfully:\n" << *forest );

   return forest;
}



shared_ptr< blockforest::StructuredBlockForest > createStructuredBlockForest( const Config::BlockHandle & configBlock,
                                                                              const blockforest::RefinementSelectionFunctions & refinementSelectionFunctions,
                                                                              const Setup & setup,
                                                                              const memory_t memoryPerCell, const memory_t processMemoryLimit )
{
   if( configBlock.isDefined( "sbffile" ) )
   {
      std::string sbffile = configBlock.getParameter< std::string >( "sbffile" );

      WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure: loading from file \'" << sbffile << "\' ..." );

      MPIManager::instance()->useWorldComm();

      auto bf = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sbffile.c_str(), true, false );

      auto sbf = std::make_shared< StructuredBlockForest >( bf, setup.xCells, setup.yCells, setup.zCells );
      sbf->createCellBoundingBoxes();

      return sbf;
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure ..." );

   shared_ptr< SetupBlockForest > sforest = createSetupBlockForest( refinementSelectionFunctions, setup,
                                                                    uint_c( MPIManager::instance()->numProcesses() ), uint_t(0),
                                                                    memoryPerCell, processMemoryLimit,
                                                                    configBlock.getParameter< bool >( "outputSetupForest", false ) );

   auto bf = std::make_shared< blockforest::BlockForest >( uint_c( MPIManager::instance()->rank() ), *sforest, false );

   auto sbf = std::make_shared< blockforest::StructuredBlockForest >( bf, setup.xCells, setup.yCells, setup.zCells );
   sbf->createCellBoundingBoxes();
   
   return sbf;
}






///////////////////////
// BOUNDARY HANDLING //
///////////////////////

class Channel
{
public:

   Channel( const Setup & _setup, const AABB & domain ) : setup_( _setup ), domain_( domain ) {}

   const Setup & setup() const { return setup_; }

   bool wallContains( const Vector3< real_t > & point ) const
   {
      if( setup_.circularProfile )
      {
         const auto domainCenter = domain_.center();
         const auto distance = ( point[1] - domainCenter[1] ) * ( point[1] - domainCenter[1] ) +
                               ( point[2] - domainCenter[2] ) * ( point[2] - domainCenter[2] );
         return distance > ( setup_.radius_L * setup_.radius_L );
      }
      else
      {
         return ( point[1] < domain_.yMin() ) || ( point[1] > domain_.yMax() );
      }
   }

   real_t delta( const Vector3< real_t > & fluid, const Vector3< real_t > & boundary ) const
   {
      WALBERLA_CHECK( !wallContains( fluid ) );
      WALBERLA_CHECK( wallContains( boundary ) );

      if( !setup_.circularProfile )
         return real_c( 0.5 );

      const real_t py = domain_.center()[1];
      const real_t pz = domain_.center()[2];
      const real_t r = setup_.radius_L;

      // http://devmag.org.za/2009/04/17/basic-collision-detection-in-2d-part-2/

      const Vector3< real_t > circle( real_t(0), py, pz );

      const Vector3< real_t > f = fluid - circle;
      const Vector3< real_t > d = ( boundary - circle ) - f;

      const real_t a = d[1] * d[1] + d[2] * d[2];
      const real_t b = real_t(2) * ( d[1] * f[1] + d[2] * f[2] );
      const real_t c = f[1] * f[1] + f[2] * f[2] - r * r;

      const real_t bb4ac = b * b - ( real_t(4) * a * c );
      WALBERLA_CHECK_GREATER_EQUAL( bb4ac, real_t(0) );

      const real_t sqrtbb4ac = std::sqrt( bb4ac );

      real_t alpha = ( -b + sqrtbb4ac ) / ( real_t(2) * a );
      const real_t beta = ( -b - sqrtbb4ac ) / ( real_t(2) * a );
      if( alpha < real_t(0) || ( beta >= real_t(0) && beta < alpha )  )
         alpha = beta;

      WALBERLA_CHECK_GREATER_EQUAL( alpha, real_t(0) );
      WALBERLA_CHECK_LESS_EQUAL( alpha, real_t(1) );

      return alpha;
   }

   bool operator()( const Vector3< real_t > & point ) const { return wallContains( point ); }

private:

   Setup setup_;
   AABB  domain_;

}; // class Channel



template< typename LatticeModel_T >
class MyBoundaryHandling
{
public:

   typedef lbm::NoSlip< LatticeModel_T, flag_t >      NoSlip_T;
   typedef lbm::Curved< LatticeModel_T, FlagField_T > Curved_T;

   typedef BoundaryHandling< FlagField_T, typename Types<LatticeModel_T>::Stencil_T, NoSlip_T, Curved_T > BoundaryHandling_T;



   MyBoundaryHandling( const BlockDataID & flagFieldId, const BlockDataID & pdfFieldId ) :
      flagFieldId_( flagFieldId ), pdfFieldId_( pdfFieldId ) {}

   BoundaryHandling_T * operator()( IBlock* const block ) const;

private:

   const BlockDataID flagFieldId_;
   const BlockDataID  pdfFieldId_;

}; // class MyBoundaryHandling

template< typename LatticeModel_T >
typename MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T *
MyBoundaryHandling<LatticeModel_T>::operator()( IBlock * const block ) const
{
   using PdfField_T = typename Types< LatticeModel_T >::PdfField_T;

   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfFieldId_ );

   const flag_t fluid = flagField->registerFlag( Fluid_Flag );

   return new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                   NoSlip_T( "no slip", NoSlip_Flag, pdfField ),
                                   Curved_T( "curved", Curved_Flag, pdfField, flagField, fluid ) );
}



template< typename LatticeModel_T >
class CurvedDeltaValueCalculation
{
public:

   using Stencil = typename LatticeModel_T::Stencil;

   CurvedDeltaValueCalculation( const shared_ptr< StructuredBlockForest > & blocks, const IBlock & block, const Channel & channel ) :
      blocks_( blocks ), block_( block ), channel_( channel ) {}

   shared_ptr< BoundaryConfiguration > operator()( const Cell & boundaryCell, const Vector3< real_t > & p );

private:

   const shared_ptr< StructuredBlockForest > & blocks_;
   const IBlock & block_;

   const Channel & channel_;

}; // class CurvedDeltaValueCalculation

template< typename LatticeModel_T >
shared_ptr< BoundaryConfiguration > CurvedDeltaValueCalculation< LatticeModel_T >::operator()( const Cell & boundaryCell, const Vector3< real_t > & /*p*/ )
{
   std::vector< real_t > deltas( Stencil::Size, real_c(0.5) );

   if( channel_.setup().circularProfile )
   {
      const Vector3< real_t > boundary = blocks_->getBlockLocalCellCenter( block_, boundaryCell );
      if( channel_.wallContains( boundary ) )
      {
         for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
         {
            const Cell fluidCell = Cell( boundaryCell.x() + dir.cx(), boundaryCell.y() + dir.cy(), boundaryCell.z() + dir.cz() );
            const Vector3< real_t > fluid = blocks_->getBlockLocalCellCenter( block_, fluidCell );
            if( !channel_.wallContains( fluid ) && blocks_->getDomain().contains( fluid ) )
               deltas[ dir.toIdx() ] = channel_.delta( fluid, boundary );
         }
      }
   }

   return shared_ptr< BoundaryConfiguration >( new typename MyBoundaryHandling< LatticeModel_T >::Curved_T::Configuration( deltas ) );
}



template< typename LatticeModel_T, typename BoundaryHandling_T >
void setFlags( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & boundaryHandlingId, const Setup & setup, const bool curvedBC )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      BoundaryHandling_T * boundaryHandling = block->template getData< BoundaryHandling_T >( boundaryHandlingId );

      // channel walls

      Channel channel( setup, blocks->getDomain() );

      if( curvedBC ) // curved boundary condition
      {
         CurvedDeltaValueCalculation< LatticeModel_T > deltaCalculation( blocks, *block, channel );

         lbm::refinement::consistentlyForceBoundary< BoundaryHandling_T >( *blocks, dynamic_cast< blockforest::Block & >(*block),
                                                                           boundaryHandlingId, Curved_Flag, channel, deltaCalculation );
      }
      else // staircase (1st order bounce back no-slip boundary condition)
      {
         lbm::refinement::consistentlyForceBoundary< BoundaryHandling_T >( *blocks, dynamic_cast< blockforest::Block & >(*block),
                                                                           boundaryHandlingId, NoSlip_Flag, channel );
      }

      // fluid
      
      const uint_t level = blocks->getLevel(*block);
      CellInterval domainBB = blocks->getDomainCellBB( level );
      blocks->transformGlobalToBlockLocalCellInterval( domainBB, *block );
      domainBB.expand( cell_idx_c( FieldGhostLayers ) );
      
      boundaryHandling->fillWithDomain( domainBB );
   }
}






/////////
// VTK //
/////////

template< typename LatticeModel_T, typename OutputType = float >
class ErrorVTKWriter : public vtk::BlockCellDataWriter< OutputType, 3 >
{
public:

   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   ErrorVTKWriter( const ConstBlockDataID & pdfFieldId, const std::string & id, const Setup & setup ) :
      vtk::BlockCellDataWriter< OutputType, 3 >( id ), bdid_( pdfFieldId ), pdf_( nullptr ), setup_( setup ) {}

protected:

   void configure() override { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); pdf_ = this->block_->template getData< PdfField_T >( bdid_ ); }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdf_ );

      const auto center = this->blockStorage_->getBlockLocalCellCenter( *(this->block_), Cell(x,y,z) );

      const auto & domain = this->blockStorage_->getDomain();
      const real_t channelRadius = setup_.radius_L;
      const auto   channelMiddle = domain.center();

      const uint_t level = this->blockStorage_->getLevel( *(this->block_) );
      const auto & dy = this->blockStorage_->dy( level );
      WALBERLA_ASSERT_FLOAT_EQUAL( dy, this->blockStorage_->dz( level ) );

      const real_t acceleration_L = pdf_->latticeModel().forceModel().force()[0]; // force = acceleration (in lattice units of the current level!)
      const real_t viscosity_L    = pdf_->latticeModel().collisionModel().viscosity(); // in lattice units on the current level
      const real_t radius_L       = channelRadius / dy; // in lattice units on the current level
      
      real_t refVelocity_x( real_t(0) ); // in lattice units on the current level
      if( setup_.circularProfile )
      {
         const real_t middleDistanceY_L = ( center[1] - channelMiddle[1] ) / dy;
         const real_t middleDistanceZ_L = ( center[2] - channelMiddle[2] ) / dy;
         const real_t middleDistance_L_2 = middleDistanceY_L * middleDistanceY_L + middleDistanceZ_L * middleDistanceZ_L;
         refVelocity_x = ( acceleration_L / ( real_t(4) * viscosity_L ) ) * ( radius_L * radius_L - middleDistance_L_2 );
      }
      else
      {
         const real_t middleDistance_L = ( center[1] - channelMiddle[1] ) / dy;
         refVelocity_x = ( acceleration_L / ( real_t(2) * viscosity_L ) ) * ( radius_L * radius_L - middleDistance_L * middleDistance_L );
      }
      
      const auto velocity = pdf_->getVelocity(x,y,z);

      // error is always a relative error:
      // 1st error component -> base = reference velocity
      // 2nd error component -> base = maximum velocity
      // 3rd error component -> base = mean velocity

      if( f == cell_idx_t(0) )
         return numeric_cast< OutputType >( std::fabs( ( refVelocity_x - velocity[0] ) / refVelocity_x ) );
      else if( f == cell_idx_t(1) )
         return numeric_cast< OutputType >( std::fabs( ( refVelocity_x - velocity[0] ) / setup_.maxVelocity_L ) );
      return numeric_cast< OutputType >( std::fabs( ( refVelocity_x - velocity[0] ) / setup_.meanVelocity_L ) );
   }

   const ConstBlockDataID bdid_;
   const PdfField_T * pdf_;
   
   const Setup & setup_;

}; // class ErrorVTKWriter



template< typename LatticeModel_T >
class MyVTKOutput {

public:

   MyVTKOutput( const ConstBlockDataID & pdfField, const ConstBlockDataID & flagField,
                const vtk::VTKOutput::BeforeFunction& pdfGhostLayerSync, const Setup & setup ) :
      setup_( setup ), pdfField_( pdfField ), flagField_( flagField ), pdfGhostLayerSync_( pdfGhostLayerSync ) {}

   void operator()( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> > & writers,
                    std::map< std::string, vtk::VTKOutput::CellFilter > &          filters,
                    std::map< std::string, vtk::VTKOutput::BeforeFunction > &      beforeFunctions );

private:

   const Setup & setup_;

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

   writers.push_back( make_shared< lbm::VelocityVTKWriter<LatticeModel_T>            >( pdfField_, "VelocityFromPDF" ) );
   writers.push_back( make_shared< lbm::VelocityMagnitudeVTKWriter<LatticeModel_T>   >( pdfField_, "VelocityMagnitudeFromPDF" ) );
   writers.push_back( make_shared< lbm::DensityVTKWriter<LatticeModel_T>             >( pdfField_, "DensityFromPDF" ) );
   writers.push_back( make_shared< lbm::NonEqulibriumVTKWriter<LatticeModel_T>       >( pdfField_, "NonEquPart" ) );
   writers.push_back( make_shared< field::VTKWriter< lbm::PdfField<LatticeModel_T> > >( pdfField_, "PDF" ) );
   writers.push_back( make_shared< ErrorVTKWriter<LatticeModel_T>                    >( pdfField_, "Error", setup_ ) );

   writers.push_back( make_shared< field::VTKWriter< FlagField_T > >( flagField_, "FlagField" ) );

   // cell filters

   field::FlagFieldCellFilter<FlagField_T> fluidFilter( flagField_ );
   fluidFilter.addFlag( Fluid_Flag );
   filters[ "FluidFilter" ] = fluidFilter;

   field::FlagFieldCellFilter<FlagField_T> obstacleFilter( flagField_ );
   obstacleFilter.addFlag( NoSlip_Flag );
   filters[ "ObstacleFilter" ] = obstacleFilter;

   // before functions

   beforeFunctions[ "PDFGhostLayerSync" ] = pdfGhostLayerSync_;
}






////////////////
// Evaluation //
////////////////

real_t exactFlowRate( const real_t flowRate )
{
   return flowRate;
}

Vector3< real_t > exactPlatesVelocity( const Vector3< real_t > & p, const shared_ptr< StructuredBlockStorage > & blocks, const Setup & setup )
{
   const real_t middleDistance_L = ( p[1] - blocks->getDomain().center()[1] ) / blocks->dy();

   return Vector3< real_t >( ( setup.acceleration_L / ( real_t(2) * setup.viscosity_L ) ) * ( setup.radius_L * setup.radius_L - middleDistance_L * middleDistance_L ),
                             real_t(0), real_t(0) );
}

Vector3< real_t > exactPipeVelocity( const Vector3< real_t > & p, const shared_ptr< StructuredBlockStorage > & blocks, const Setup & setup )
{
   const auto & domain = blocks->getDomain();
   const auto   channelMiddle = domain.center();

   const real_t middleDistanceY_L = ( p[1] - channelMiddle[1] ) / blocks->dy();
   const real_t middleDistanceZ_L = ( p[2] - channelMiddle[2] ) / blocks->dz();
   const real_t middleDistance_L_2 = middleDistanceY_L * middleDistanceY_L + middleDistanceZ_L * middleDistanceZ_L;

   return Vector3< real_t >( ( setup.acceleration_L / ( real_t(4) * setup.viscosity_L ) ) * ( setup.radius_L * setup.radius_L - middleDistance_L_2 ),
                             real_t(0), real_t(0) );
}






////////////////////
// THE SIMULATION //
////////////////////

template< typename LatticeModel_T >
void run( const shared_ptr< Config > & config, const LatticeModel_T & latticeModel,
          const bool fzyx, const bool syncComm, const bool fullComm, const bool linearExplosion,
          const blockforest::RefinementSelectionFunctions & refinementSelectionFunctions, Setup & setup,
          const memory_t memoryPerCell, const memory_t processMemoryLimit )
{
   Config::BlockHandle configBlock = config->getBlock( "PoiseuilleChannel" );

   if( setup.circularProfile )
   {
      setup.maxVelocity_L = ( setup.acceleration_L * setup.radius_L * setup.radius_L ) / ( real_t(4) * setup.viscosity_L );
      setup.meanVelocity_L = ( setup.acceleration_L * setup.radius_L * setup.radius_L ) / ( real_t(8) * setup.viscosity_L );
      setup.flowRate_L = setup.meanVelocity_L * math::pi * setup.radius_L * setup.radius_L;
   }
   else
   {
      setup.maxVelocity_L = ( setup.acceleration_L * setup.radius_L * setup.radius_L ) / ( real_t(2) * setup.viscosity_L );
      setup.meanVelocity_L = ( setup.acceleration_L * setup.radius_L * setup.radius_L ) / ( real_t(3) * setup.viscosity_L );
      setup.flowRate_L = setup.meanVelocity_L * real_c( setup.yBlocks * setup.yCells * setup.zBlocks * setup.zCells );
   }

   // creating the block structure

   auto blocks = createStructuredBlockForest( configBlock, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );

   // add pdf field to blocks

   const real_t initVelocity = ( configBlock.getParameter< bool >( "initWithMeanVelocity", false ) ) ? setup.meanVelocity_L : real_t(0);

   BlockDataID pdfFieldId = fzyx ? lbm::addPdfFieldToStorage( blocks, "pdf field (fzyx)", latticeModel,
                                                              Vector3< real_t >( initVelocity, real_c(0), real_c(0) ), real_t(1),
                                                              FieldGhostLayers, field::fzyx ) :
                                   lbm::addPdfFieldToStorage( blocks, "pdf field (zyxf)", latticeModel,
                                                              Vector3< real_t >( initVelocity, real_c(0), real_c(0) ), real_t(1),
                                                              FieldGhostLayers, field::zyxf );

   using VelocityAdaptor_T = typename lbm::Adaptor< LatticeModel_T >::VelocityVector;
   using DensityAdaptor_T = typename lbm::Adaptor< LatticeModel_T >::Density;
   BlockDataID velocityAdaptorId = field::addFieldAdaptor< VelocityAdaptor_T >( blocks, pdfFieldId, "velocity adaptor" );
   BlockDataID  densityAdaptorId = field::addFieldAdaptor<  DensityAdaptor_T >( blocks, pdfFieldId, "density adaptor" );
   
   // add flag field to blocks

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field", FieldGhostLayers );

   // add LB boundary handling to blocks

   BlockDataID boundaryHandlingId = blocks->template addBlockData< typename MyBoundaryHandling< LatticeModel_T >::BoundaryHandling_T >(
            MyBoundaryHandling< LatticeModel_T >( flagFieldId, pdfFieldId ), "boundary handling" );

   const bool curvedBC = configBlock.getParameter< bool >( "curvedBC", false );
   setFlags< LatticeModel_T, typename MyBoundaryHandling< LatticeModel_T >::BoundaryHandling_T >( blocks, boundaryHandlingId, setup, curvedBC );

   // creating the time loop

   const uint_t outerTimeSteps = configBlock.getParameter< uint_t >( "outerTimeSteps", uint_c(10) );
   const uint_t innerTimeSteps = configBlock.getParameter< uint_t >( "innerTimeSteps", uint_c(10) );

   SweepTimeloop timeloop( blocks->getBlockStorage(), ( outerTimeSteps * innerTimeSteps ) + uint_t(1) );

   // VTK

   blockforest::communication::NonUniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< lbm::refinement::PdfFieldSyncPackInfo< LatticeModel_T > >( pdfFieldId ) );

   MyVTKOutput< LatticeModel_T > myVTKOutput( pdfFieldId, flagFieldId, pdfGhostLayerSync, setup );

   std::map< std::string, vtk::SelectableOutputFunction > vtkOutputFunctions;
   vtk::initializeVTKOutput( vtkOutputFunctions, myVTKOutput, blocks, config );   
   
   const bool vtkBeforeTimeStep = configBlock.getParameter< bool >( "vtkBeforeTimeStep", true );

   if( vtkBeforeTimeStep )
   {
      for( auto output = vtkOutputFunctions.begin(); output != vtkOutputFunctions.end(); ++output )
         timeloop.addFuncBeforeTimeStep( output->second.outputFunction, std::string("VTK: ") + output->first,
                                         output->second.requiredGlobalStates, output->second.incompatibleGlobalStates );
   }

   // add 'refinement' LB time step to time loop

   shared_ptr<WcTimingPool> refinementTimeStepTiming = make_shared<WcTimingPool>();
   shared_ptr<WcTimingPool> refinementTimeStepLevelwiseTiming = make_shared<WcTimingPool>();

   using BH_T = typename MyBoundaryHandling< LatticeModel_T >::BoundaryHandling_T;

   auto mySweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag );
   auto ts = lbm::refinement::makeTimeStep< LatticeModel_T, BH_T >( blocks, mySweep, pdfFieldId, boundaryHandlingId );
   ts->asynchronousCommunication( !syncComm );
   ts->optimizeCommunication( !fullComm );
   ts->performLinearExplosion( linearExplosion );
   ts->enableTiming( refinementTimeStepTiming, refinementTimeStepLevelwiseTiming );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor(ts), "LBM refinement time step" );
                                                 
   // evaluation 
   
   const auto exactSolutionFunction = setup.circularProfile ? std::bind( exactPipeVelocity, std::placeholders::_1, blocks, setup ) :
                                                              std::bind( exactPlatesVelocity, std::placeholders::_1, blocks, setup );

   auto volumetricFlowRate = field::makeVolumetricFlowRateEvaluation< VelocityAdaptor_T, FlagField_T >( configBlock, blocks, velocityAdaptorId,
                                                                                                        flagFieldId, Fluid_Flag,
                                                                                                        std::bind( exactFlowRate, setup.flowRate_L ),
                                                                                                        exactSolutionFunction );
   volumetricFlowRate->setNormalizationFactor( real_t(1) / setup.maxVelocity_L );
   volumetricFlowRate->setDomainNormalization( Vector3<real_t>( real_t(1) ) );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( volumetricFlowRate ), "volumetric flow rate evaluation" );

   auto accuracyEvaluation = field::makeAccuracyEvaluation< VelocityAdaptor_T, FlagField_T >( configBlock, blocks, velocityAdaptorId,
                                                                                              flagFieldId, Fluid_Flag, exactSolutionFunction );
   accuracyEvaluation->setNormalizationFactor( real_t(1) / setup.maxVelocity_L );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( accuracyEvaluation ), "accuracy evaluation" );
   
   auto linePlot = field::makeAccuracyEvaluationLinePlot< VelocityAdaptor_T, FlagField_T >( configBlock, blocks, velocityAdaptorId,
                                                                                            flagFieldId, Fluid_Flag, exactSolutionFunction );
   linePlot->setNormalizationFactor( real_t(1) / setup.maxVelocity_L );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( field::makeAccuracyEvaluationLinePlotter( configBlock, linePlot ) ), "accuracy evaluation (line plot)" );
   
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( lbm::makeMassEvaluation< DensityAdaptor_T, FlagField_T >( configBlock, blocks, uint_t(0),
                                                                                                                densityAdaptorId, flagFieldId, Fluid_Flag ) ),
                                   "mass evaluation" );
   
   // stability check (non-finite values in the PDF field?)

   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< lbm::PdfField< LatticeModel_T >, FlagField_T >(
                                                        configBlock, blocks, pdfFieldId, flagFieldId, Fluid_Flag ) ),
                                  "LBM stability check" );

   // VTK

   if( !vtkBeforeTimeStep )
   {
      for( auto output = vtkOutputFunctions.begin(); output != vtkOutputFunctions.end(); ++output )
         timeloop.addFuncAfterTimeStep( output->second.outputFunction, std::string("VTK: ") + output->first,
                                        output->second.requiredGlobalStates, output->second.incompatibleGlobalStates );
   }

   // remaining time logger

   const double remainingTimeLoggerFrequency = configBlock.getParameter< double >( "remainingTimeLoggerFrequency", 3.0 );
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "Remaining time logger" );

   // logging right before the simulation starts

   lbm::BlockForestEvaluation< FlagField_T > blockForest( blocks, flagFieldId, Fluid_Flag );
   blockForest.logInfoOnRoot();

   field::CellCounter< FlagField_T > fluidCells( blocks, flagFieldId, Fluid_Flag );
   fluidCells();

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark run data:"
                              "\n- simulation parameters:"
                              "\n   + collision model:    " << CollisionModelString< LatticeModel_T >::str() <<
                              "\n   + stencil:            " << StencilString< LatticeModel_T >::str() <<
                              "\n   + compressible:       " << ( LatticeModel_T::compressible ? "yes" : "no" ) <<
                              "\n   + data layout:        " << ( fzyx ? "fzyx (structure of arrays [SoA])" : "zyxf (array of structures [AoS])" ) <<
                              "\n   + communication:      " << ( fullComm ? ( syncComm ? "synchronous, full synchronization" : "asynchronous, full synchronization" ) :
                                                                          ( syncComm ? "synchronous, block neighborhood and direction-aware optimizations" : "asynchronous, block neighborhood and direction-aware optimizations" ) ) <<
                              "\n   + linear explosion:   " << ( linearExplosion ? "yes" : "no" ) <<
                              "\n- simulation properties:"
                              "\n   + fluid cells:        " << fluidCells.numberOfCells() << " (in total on all levels)" <<
                              "\n   + circular profile:   " << ( setup.circularProfile ? "yes (= pipe flow)" : "no (= flow between parallel plates)" ) <<
                              "\n   + boundary condition: " << ( curvedBC ? "higher order, curved" : "no slip, bounce back" ) <<
                              "\n   + radius:             " << setup.radius_L <<
                              "\n   + Reynolds number:    " << setup.Re <<
                              "\n   + kin. viscosity:     " << setup.viscosity_L << " (on the coarsest grid)" <<
                              "\n   + acceleration:       " << setup.acceleration_L << " (on the coarsest grid)" <<
                              "\n   + maximum velocity:   " << setup.maxVelocity_L <<
                              "\n   + mean velocity:      " << setup.meanVelocity_L <<
                              "\n   + flow rate:          " << setup.flowRate_L << " (in lattice units of the coarsest grid)" <<
                              "\n   + #time steps:        " << timeloop.getNrOfTimeSteps() << " (on the coarsest grid)" );

   // run the simulation

   lbm::PerformanceEvaluation< FlagField_T > performance( blocks, flagFieldId, Fluid_Flag );

   for( uint_t outerRun = 0; outerRun < outerTimeSteps; ++outerRun )
   {
      WcTimingPool timeloopTiming;

      WALBERLA_MPI_WORLD_BARRIER();
      WcTimer timer;
      timer.start();

      for( uint_t innerRun = 0; innerRun < innerTimeSteps; ++innerRun )
         timeloop.singleStep( timeloopTiming );

      timer.end();

      double time = timer.max();
      mpi::reduceInplace( time, mpi::MAX );

      const auto reducedTimeloopTiming = timeloopTiming.getReduced();
      const auto reducedRTSTiming  = refinementTimeStepTiming->getReduced();
      const auto reducedRTSLTiming = refinementTimeStepLevelwiseTiming->getReduced();
      refinementTimeStepTiming->clear();
      refinementTimeStepLevelwiseTiming->clear();

      WALBERLA_LOG_RESULT_ON_ROOT( "Time loop timing:\n" << *reducedTimeloopTiming );
      WALBERLA_LOG_RESULT_ON_ROOT( "Refinement time step timing:\n" << *reducedRTSTiming );
      WALBERLA_LOG_RESULT_ON_ROOT( "Refinement time step timing (one timer per level):\n" << *reducedRTSLTiming );

      performance.logResultOnRoot( innerTimeSteps, time );

      WALBERLA_ROOT_SECTION()
      {
         // logging in SQL database

         if( configBlock.getParameter< bool >( "logToSqlDB", true ) )
         {
            const std::string sqlFile = configBlock.getParameter< std::string >( "sqlFile", "performance.sqlite" );

            std::map< std::string, int >        integerProperties;
            std::map< std::string, double >        realProperties;
            std::map< std::string, std::string > stringProperties;

            performance.getResultsForSQLOnRoot( integerProperties, realProperties, stringProperties, innerTimeSteps, time );
            blockForest.getResultsForSQLOnRoot( integerProperties, realProperties, stringProperties );

            stringProperties[ "collisionModel" ]    = CollisionModelString< LatticeModel_T >::str();
            stringProperties[ "stencil" ]           = StencilString< LatticeModel_T >::str();
            stringProperties[ "compressible" ]      = ( LatticeModel_T::compressible ? "yes" : "no" );
            stringProperties[ "dataLayout" ]        = ( fzyx ? "fzyx" : "zyxf" );
            stringProperties[ "syncCommunication" ] = ( syncComm ? "yes" : "no" );
            stringProperties[ "fullCommunication" ] = ( fullComm ? "yes" : "no" );
            stringProperties[ "linearExplosion" ]   = ( linearExplosion ? "yes" : "no" );

            realProperties[ "Re" ]             = double_c( setup.Re );
            realProperties[ "radius_L" ]       = double_c( setup.radius_L );
            realProperties[ "viscosity_L" ]    = double_c( setup.viscosity_L );
            realProperties[ "acceleration_L" ] = double_c( setup.acceleration_L );
            realProperties[ "maxVelocity_L" ]  = double_c( setup.maxVelocity_L );
            realProperties[ "meanVelocity_L" ] = double_c( setup.meanVelocity_L );
            realProperties[ "flowRate_L" ]     = double_c( setup.flowRate_L );

            stringProperties[ "circularProfile" ] = ( setup.circularProfile ? "yes" : "no" );
            stringProperties[ "curvedBC" ]        = ( curvedBC ? "yes" : "no" );
            
            integerProperties[ "simulationTimeSteps" ] = int_c( timeloop.getNrOfTimeSteps() );
            
            realProperties[ "simulationProgress" ] = double_c( ( outerRun + uint_t(1) ) * innerTimeSteps ) / double_c( outerTimeSteps * innerTimeSteps );

            auto runId = sqlite::storeRunInSqliteDB( sqlFile, integerProperties, stringProperties, realProperties );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedTimeloopTiming, "Timeloop" );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedRTSTiming, "RefinementTimeStep" );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedRTSLTiming, "RefinementTimeStepLevelwise" );
         }
      }
   }

   // Do one more step so that we see the final output
   timeloop.singleStep();

   // logging once again at the end of the simulation, identical to logging at the beginning :-)

   blockForest.logInfoOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark run data:"
                              "\n- simulation parameters:"
                              "\n   + collision model:    " << CollisionModelString< LatticeModel_T >::str() <<
                              "\n   + stencil:            " << StencilString< LatticeModel_T >::str() <<
                              "\n   + compressible:       " << ( LatticeModel_T::compressible ? "yes" : "no" ) <<
                              "\n   + data layout:        " << ( fzyx ? "fzyx (structure of arrays [SoA])" : "zyxf (array of structures [AoS])" ) <<
                              "\n   + communication:      " << ( fullComm ? ( syncComm ? "synchronous, full synchronization" : "asynchronous, full synchronization" ) :
                                                                          ( syncComm ? "synchronous, block neighborhood and direction-aware optimizations" : "asynchronous, block neighborhood and direction-aware optimizations" ) ) <<
                              "\n   + linear explosion:   " << ( linearExplosion ? "yes" : "no" ) <<
                              "\n- simulation properties:"
                              "\n   + fluid cells:        " << fluidCells.numberOfCells() << " (in total on all levels)" <<
                              "\n   + circular profile:   " << ( setup.circularProfile ? "yes (= pipe flow)" : "no (= flow between parallel plates)" ) <<
                              "\n   + boundary condition: " << ( curvedBC ? "higher order, curved" : "no slip, bounce back" ) <<
                              "\n   + radius:             " << setup.radius_L <<
                              "\n   + Reynolds number:    " << setup.Re <<
                              "\n   + kin. viscosity:     " << setup.viscosity_L << " (on the coarsest grid)" <<
                              "\n   + acceleration:       " << setup.acceleration_L << " (on the coarsest grid)" <<
                              "\n   + maximum velocity:   " << setup.maxVelocity_L <<
                              "\n   + mean velocity:      " << setup.meanVelocity_L <<
                              "\n   + flow rate:          " << setup.flowRate_L << " (in lattice units of the coarsest grid)" <<
                              "\n   + #time steps:        " << timeloop.getNrOfTimeSteps() << " (on the coarsest grid)" );
                              
   auto accuracyEvaluationLinePlotBlock = configBlock.getBlock("AccuracyEvaluationLinePlot");
   if( accuracyEvaluationLinePlotBlock && accuracyEvaluationLinePlotBlock.isDefined("filename") )
      (*linePlot)( accuracyEvaluationLinePlotBlock.template getParameter<std::string>("filename") );

   WALBERLA_ROOT_SECTION()
   {
      if( configBlock.getParameter< bool >( "check", false ) )
      {
         if( configBlock.isDefined( "checkFlowRateError" ) )
         {
            const real_t checkFlowRateError = configBlock.getParameter< real_t >( "checkFlowRateError" );
            WALBERLA_CHECK_LESS( std::abs( ( volumetricFlowRate->flowRate() - volumetricFlowRate->solution() ) / volumetricFlowRate->solution() ),
                                 checkFlowRateError );
         }
         if( configBlock.isDefined( "checkErrorL1" ) )
         {
            const real_t checkErrorL1 = configBlock.getParameter< real_t >( "checkErrorL1" );
            WALBERLA_CHECK_LESS( accuracyEvaluation->L1(), checkErrorL1 );
         }
         if( configBlock.isDefined( "checkErrorL2" ) )
         {
            const real_t checkErrorL2 = configBlock.getParameter< real_t >( "checkErrorL2" );
            WALBERLA_CHECK_LESS( accuracyEvaluation->L2(), checkErrorL2 );
         }
         if( configBlock.isDefined( "checkErrorLmax" ) )
         {
            const real_t checkErrorLmax = configBlock.getParameter< real_t >( "checkErrorLmax" );
            WALBERLA_CHECK_LESS( accuracyEvaluation->Lmax(), checkErrorLmax );
         }
      }
   }
}



//////////
// MAIN //
//////////

enum CM { CMSRT, CMTRT };
enum LM { LMD3Q19, LMD3Q27 };

int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );
   
   if( argc < 2 )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::cout << "Usage: " << argv[0] << " path-to-configuration-file [--trt] [--d3q27] [--comp] [--fzyx] [--sync-comm] [--full-comm] [--linear-exp]\n"
                      "\n"
                      "By default, SRT is selected as collision model, an asynchronous communication scheme with block neighborhood and\n"
                      "direction-aware optimizations is chosen, and an incompressible D3Q19 LB kernel is executed on a PDF field with\n"
                      "layout 'zyxf' (= array of structures [AoS]).\n"
                      "\n"
                      "Optional arguments:\n"
                      " --trt:        collision model = TRT\n"
                      " --d3q27:      A D3Q27 model is used.\n"
                      " --comp:       LB kernel is switched from incompressible to compressible\n"
                      " --fzyx:       data layout switched to 'fzyx' (structure of arrays [SoA])\n"
                      " --sync-comm:  A synchronous communication scheme is used instead of an asynchronous scheme\n"
                      "               which is used by default.\n"
                      " --full-comm:  A full synchronization of neighboring blocks is performed instead of using a communication\n"
                      "               that uses block neighborhood and direction-aware optimizations.\n"
                      " --linear-exp: When communicating from coarse to fine grids, a linear interpolation scheme is used\n"
                      "               instead of a uniform distribution of a coarse cell to eight fine cells." << std::endl;
      }
      return EXIT_SUCCESS;
   }   

   logging::Logging::printHeaderOnStream();
   //WALBERLA_ROOT_SECTION() { logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS ); }

#ifdef _OPENMP
   if( std::getenv( "OMP_NUM_THREADS" ) == nullptr )
      WALBERLA_ABORT( "If you are using a version of the program that was compiled with OpenMP you have to "
                      "specify the environment variable \'OMP_NUM_THREADS\' accordingly!" );
#endif

   // open configuration file

   shared_ptr< Config > config = make_shared< Config >();
   config->readParameterFile( argv[1] );

   Config::BlockHandle configBlock = config->getBlock( "PoiseuilleChannel" );

   if( !configBlock )
      WALBERLA_ABORT( "You have to specify a \"PoiseuilleChannel\" block in the configuration file!" );

   ////////////////
   // PARAMETERS //
   ////////////////

   Setup setup;

   setup.xBlocks = configBlock.getParameter< uint_t >( "xBlocks", uint_t(1) );
   setup.yBlocks = configBlock.getParameter< uint_t >( "yBlocks", uint_t(1) );
   setup.zBlocks = configBlock.getParameter< uint_t >( "zBlocks", uint_t(1) );

   setup.xCells = configBlock.getParameter< uint_t >( "xCells", uint_t(10) );
   setup.yCells = configBlock.getParameter< uint_t >( "yCells", uint_t(50) );
   setup.zCells = configBlock.getParameter< uint_t >( "zCells", uint_t(10) );

   setup.Re = configBlock.getParameter< real_t >( "Re", real_t(10) );
   
   // http://www.ae.metu.edu.tr/~ae244/docs/FluidMechanics-by-JamesFay/2003/Textbook/Nodes/chap06/node9.html
   // http://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node106.html

   setup.circularProfile = configBlock.getParameter< bool >( "circularProfile", true );

   if( setup.circularProfile )
   {
      const real_t ry = real_c(0.5) * real_c(setup.yBlocks * setup.yCells);
      const real_t rz = real_c(0.5) * real_c(setup.zBlocks * setup.zCells);
      setup.radius_L = std::min( ry, rz );
   }
   else
   {
      setup.radius_L = real_c(0.5) * real_c(setup.yBlocks * setup.yCells);
   }
   
   // ... in bytes
   const memory_t memoryPerCell = configBlock.getParameter< memory_t >( "memoryPerCell", memory_t( 19 * 8 + 1 ) );
   // ... in MiB
   const memory_t processMemoryLimit = configBlock.getParameter< memory_t >( "processMemoryLimit", memory_t( 512 ) ) * memory_t( 1024 * 1024  );

   const uint_t blocksPerProcess = configBlock.getParameter< uint_t >( "blocksPerProcess", uint_t(8) );

   ////////////////////////
   // REFINEMENT REGIONS //
   ////////////////////////

   blockforest::RefinementSelectionFunctions refinementSelectionFunctions;

   blockforest::AABBRefinementSelection aabbRefinementSelection( configBlock );
   refinementSelectionFunctions.add( aabbRefinementSelection );

   if( configBlock.getParameter< bool >( "refineOnBorder", false ) )
   {
      if( !configBlock.isDefined("borderRefinementLevel")  )
         WALBERLA_ABORT( "You have to specify \'borderRefinementLevel\' in the \"PoiseuilleChannel\" block of the configuration file (" << argv[1] << ")" );

      const real_t borderRefinementBuffer = configBlock.getParameter< real_t >( "borderRefinementBuffer", real_t(0) );

      BorderRefinementSelection borderRefinementSelection( setup, configBlock.getParameter< uint_t >( "borderRefinementLevel" ),
               borderRefinementBuffer );

      refinementSelectionFunctions.add( borderRefinementSelection );
   }

   // In case 'sbffile' and 'saveToFile' are specified in the configuration file:
   // -> just create the block structure and save it to file for later use

   if( configBlock.isDefined( "sbffile" ) && configBlock.isDefined( "saveToFile" ) )
   {
      std::string  sbffile = configBlock.getParameter< std::string >( "sbffile" );

      std::ostringstream infoString;
      infoString << "You have selected the option of just creating the block structure (= domain decomposition) and saving the result to file\n"
                    "by specifying the output file name \'" << sbffile << "\' AND also specifying \'saveToFile\'.\n";

      if( MPIManager::instance()->numProcesses() > 1 )
         WALBERLA_ABORT( infoString.str() << "In this mode you need to start " << argv[0] << " with just one process!" );

      if( blocksPerProcess == 0 )
         WALBERLA_ABORT( infoString.str() << "In this mode you need to specify \'blocksPerProcess\' in the configuration file "
                         "AND \'blocksPerProcess\' must be greater than 0!" );

      WALBERLA_LOG_INFO_ON_ROOT( infoString.str() << "Creating the block structure ..." );

      shared_ptr< SetupBlockForest > sforest = createSetupBlockForest( refinementSelectionFunctions, setup, uint_t(0), blocksPerProcess,
                                                                       memoryPerCell, processMemoryLimit,
                                                                       configBlock.getParameter< bool >( "outputSetupForest", true ) );
      sforest->saveToFile( sbffile.c_str() );

      logging::Logging::printFooterOnStream();
      return EXIT_SUCCESS;
   }

   // reading optional parameters from passed arguments

   CM   collisionModel  = CMSRT;
   LM   lm              = LMD3Q19;
   bool compressible    = false;
   bool fzyx            = false;
   bool syncComm        = false;
   bool fullComm        = false;
   bool linearExplosion = false;

   for( int i = 2; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--trt" )        == 0 ) collisionModel  = CMTRT;
      if( std::strcmp( argv[i], "--d3q27" )      == 0 ) lm              = LMD3Q27;
      if( std::strcmp( argv[i], "--comp" )       == 0 ) compressible    = true;
      if( std::strcmp( argv[i], "--fzyx" )       == 0 ) fzyx            = true;
      if( std::strcmp( argv[i], "--sync-comm" )  == 0 ) syncComm        = true;
      if( std::strcmp( argv[i], "--full-comm" )  == 0 ) fullComm        = true;
      if( std::strcmp( argv[i], "--linear-exp" ) == 0 ) linearExplosion = true;
   }

   if( collisionModel == CMTRT && lm == LMD3Q19 && compressible )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "For TRT and D3Q19, \"--comp\" cannot be specified!\n"
                                    "Setting \"compressible\" to false ..." );
      compressible = false;
   }

   // executing benchmark

   const real_t omega = configBlock.getParameter< real_t >( "omega", real_t(1.4) );

   const real_t magicNumber = configBlock.getParameter< real_t >( "magicNumber", real_t(3) / real_t(16) );

   const real_t lambda_e = configBlock.getParameter< real_t >( "lambda_e", real_t(1.4) );
   const real_t lambda_d = configBlock.getParameter< real_t >( "lambda_d", real_t(1.4) );

   const uint_t relaxationParametersLevel = configBlock.getParameter< uint_t >( "relaxationParametersLevel", uint_t(0) );

   if( collisionModel == CMSRT ) // SRT
   {
      auto cm = lbm::collision_model::SRT( omega, relaxationParametersLevel );

      setup.viscosity_L = cm.viscosity( uint_t(0) );
      if( setup.circularProfile )
         setup.acceleration_L = ( real_t(4) * setup.viscosity_L * setup.viscosity_L * setup.Re ) / ( setup.radius_L * setup.radius_L * setup.radius_L );
      else
         setup.acceleration_L = ( real_t(3) * setup.viscosity_L * setup.viscosity_L * setup.Re ) / ( real_t(2) * setup.radius_L * setup.radius_L * setup.radius_L );

      if( lm == LMD3Q19 )
      {
         if( compressible )
         {
            D3Q19_SRT_COMP latticeModel = D3Q19_SRT_COMP( cm, lbm::force_model::SimpleConstant( setup.acceleration_L, real_t(0), real_t(0) ) );
            run( config, latticeModel, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q19_SRT_INCOMP latticeModel = D3Q19_SRT_INCOMP( cm, lbm::force_model::SimpleConstant( setup.acceleration_L, real_t(0), real_t(0) ) );
            run( config, latticeModel, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
      else
      {
         if( compressible )
         {
            D3Q27_SRT_COMP latticeModel = D3Q27_SRT_COMP( cm, lbm::force_model::SimpleConstant( setup.acceleration_L, real_t(0), real_t(0) ) );
            run( config, latticeModel, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q27_SRT_INCOMP latticeModel = D3Q27_SRT_INCOMP( cm, lbm::force_model::SimpleConstant( setup.acceleration_L, real_t(0), real_t(0) ) );
            run( config, latticeModel, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
   }
   else // if( collisionModel == CMTRT ) // TRT
   {
      bool useLambdas = false;
      if( configBlock.isDefined("lambda_e") && configBlock.isDefined("lambda_d") )
         useLambdas = true;

      auto cm = useLambdas ? lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) :
                             lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel );

      setup.viscosity_L = cm.viscosity( uint_t(0) );
      if( setup.circularProfile )
         setup.acceleration_L = ( real_t(4) * setup.viscosity_L * setup.viscosity_L * setup.Re ) / ( setup.radius_L * setup.radius_L * setup.radius_L );
      else
         setup.acceleration_L = ( real_t(3) * setup.viscosity_L * setup.viscosity_L * setup.Re ) / ( real_t(2) * setup.radius_L * setup.radius_L * setup.radius_L );

      if( lm == LMD3Q19 )
      {
         WALBERLA_CHECK( !compressible );

         D3Q19_TRT_INCOMP latticeModel = D3Q19_TRT_INCOMP( cm, lbm::force_model::SimpleConstant( setup.acceleration_L, real_t(0), real_t(0) ) );
         run( config, latticeModel, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
      }
      else
      {
         if( compressible )
         {
            D3Q27_TRT_COMP latticeModel = D3Q27_TRT_COMP( cm, lbm::force_model::SimpleConstant( setup.acceleration_L, real_t(0), real_t(0) ) );
            run( config, latticeModel, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q27_TRT_INCOMP latticeModel = D3Q27_TRT_INCOMP( cm, lbm::force_model::SimpleConstant( setup.acceleration_L, real_t(0), real_t(0) ) );
            run( config, latticeModel, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
   }

   logging::Logging::printFooterOnStream();
   return EXIT_SUCCESS;
}

} // namespace poiseuille_channel

int main( int argc, char ** argv )
{
   return poiseuille_channel::main( argc, argv );
}
