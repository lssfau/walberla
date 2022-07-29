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
//! \file SchaeferTurek.cpp
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/AABBRefinementSelection.h"
#include "blockforest/BlockForestEvaluation.h"
#include "blockforest/communication/NonUniformBufferedScheme.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/loadbalancing/DynamicCurve.h"
#include "blockforest/loadbalancing/DynamicDiffusive.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/SharedFunctor.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/CheckFunctions.h"
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

#include "domain_decomposition/BlockSweepWrapper.h"

#include "field/AddToStorage.h"
#include "field/CellCounter.h"
#include "field/FlagField.h"
#include "field/FlagUID.h"
#include "field/StabilityChecker.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/blockforest/GradientRefinement.h"
#include "field/iterators/FieldIterator.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"

#include "lbm/BlockForestEvaluation.h"
#include "lbm/MassEvaluation.h"
#include "lbm/PerformanceEvaluation.h"
#include "lbm/blockforest/PostProcessing.h"
#include "lbm/boundary/Curved.h"
#include "lbm/boundary/DynamicUBB.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/Outlet.h"
#include "lbm/boundary/SimplePressure.h"
#include "lbm/boundary/VelocityBoundary.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/VelocityFieldWriter.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D2Q9.h"
#include "lbm/lattice_model/D3Q15.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"
#include "lbm/refinement/BoundarySetup.h"
#include "lbm/refinement/PdfFieldSyncPackInfo.h"
#include "lbm/refinement/TimeStep.h"
#include "lbm/refinement/VorticityBasedLevelDetermination.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SplitPureSweep.h"
#include "lbm/sweeps/SplitSweep.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/NonEquilibrium.h"
#include "lbm/vtk/Velocity.h"

#include "sqlite/SQLite.h"

#include "stencil/D3Q15.h"
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



namespace schaefer_turek {

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;
using walberla::uint8_t;
using walberla::real_t;

//////////////
// TYPEDEFS //
//////////////

using D2Q9_SRT_INCOMP = lbm::D2Q9<lbm::collision_model::SRT, false>;
using D2Q9_SRT_COMP = lbm::D2Q9<lbm::collision_model::SRT, true>;
using D2Q9_TRT_INCOMP = lbm::D2Q9<lbm::collision_model::TRT, false>;
using D2Q9_TRT_COMP = lbm::D2Q9<lbm::collision_model::TRT, true>;

using D3Q15_SRT_INCOMP = lbm::D3Q15<lbm::collision_model::SRT, false>;
using D3Q15_SRT_COMP = lbm::D3Q15<lbm::collision_model::SRT, true>;
using D3Q15_TRT_INCOMP = lbm::D3Q15<lbm::collision_model::TRT, false>;
using D3Q15_TRT_COMP = lbm::D3Q15<lbm::collision_model::TRT, true>;

using D3Q19_SRT_INCOMP = lbm::D3Q19<lbm::collision_model::SRT, false>;
using D3Q19_SRT_COMP = lbm::D3Q19<lbm::collision_model::SRT, true>;
using D3Q19_TRT_INCOMP = lbm::D3Q19<lbm::collision_model::TRT, false>;
using D3Q19_TRT_COMP = lbm::D3Q19<lbm::collision_model::TRT, true>;
using D3Q19_MRT_INCOMP = lbm::D3Q19<lbm::collision_model::D3Q19MRT, false>;

using D3Q27_SRT_INCOMP = lbm::D3Q27<lbm::collision_model::SRT, false>;
using D3Q27_SRT_COMP = lbm::D3Q27<lbm::collision_model::SRT, true>;
using D3Q27_TRT_INCOMP = lbm::D3Q27<lbm::collision_model::TRT, false>;
using D3Q27_TRT_COMP = lbm::D3Q27<lbm::collision_model::TRT, true>;

template< typename LatticeModel_T >
struct Types
{
   using Stencil_T = typename LatticeModel_T::Stencil;
   using PdfField_T = lbm::PdfField< LatticeModel_T >;
};

using flag_t = walberla::uint16_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers  = uint_t(4);

///////////
// FLAGS //
///////////

const FlagUID          Fluid_Flag( "fluid" );
const FlagUID         NoSlip_Flag( "no slip" );
const FlagUID       Obstacle_Flag( "obstacle (staircase)" );
const FlagUID         Curved_Flag( "obstacle (curved)" );
const FlagUID            UBB_Flag( "velocity bounce back" );
const FlagUID PressureOutlet_Flag( "pressure outlet" );
const FlagUID       Outlet21_Flag( "outlet (2/1)" );
const FlagUID       Outlet43_Flag( "outlet (4/3)" );

///////////
// SUIDS //
///////////

const SUID Empty( "empty" );
const Set<SUID> None( Set<SUID>::emptySet() );

////////
// 2D //
////////

template< typename LatticeModel_T >
struct Is2D
{
   static const bool value = LatticeModel_T::Stencil::D == 2;
};

/////////////////////
// OUTPUT HELPERS  //
/////////////////////

template< typename LatticeModel_T, class Enable = void >
struct StencilString;

template< typename LatticeModel_T >
struct StencilString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D2Q9 >::value >::type >
{
   static const char * str() { return "D2Q9"; }
};

template< typename LatticeModel_T >
struct StencilString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q15 >::value >::type >
{
   static const char * str() { return "D3Q15"; }
};

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

template< typename LatticeModel_T >
struct CollisionModelString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                          lbm::collision_model::MRT_tag >::value >::type >
{
   static const char * str() { return "MRT"; }
};

//////////////////////
// Parameter Struct //
//////////////////////

struct Setup
{
   bool pseudo2D; // false == 3D

   uint_t xBlocks;
   uint_t yzBlocks;

   uint_t xCells;
   uint_t yzCells;

   real_t H;
   real_t L;
   bool   strictlyObeyL;

   real_t cylinderxPosition;
   real_t cylinderyPosition;
   real_t cylinderRadius;
   bool   circularCrossSection;

   bool   evaluateForceComponents;
   uint_t nbrOfEvaluationPointsForCoefficientExtremas;

   bool evaluatePressure;
   Vector3< real_t > pAlpha;
   Vector3< real_t > pOmega;

   bool evaluateStrouhal;
   Vector3< real_t > pStrouhal;

   real_t viscosity;
   real_t rho;
   real_t inflowVelocity;
   real_t dx;
   real_t dt;

   real_t viscosity_L;
   real_t inflowVelocity_L;

   real_t raisingTime;
   real_t raisingTime_L;
   real_t sinPeriod;
   real_t sinPeriod_L;
};






/////////////////////
// BLOCK STRUCTURE //
/////////////////////

class Cylinder
{
public:

   Cylinder( const Setup & setup ) : setup_( setup ) {}
   
   bool operator()( const Vector3< real_t > & point ) const { return contains( point ); }
   
   bool contains( const Vector3< real_t > & point ) const;
   bool contains( const AABB & aabb ) const;
   
   bool intersects( const AABB & aabb, const real_t bufferDistance = real_t(0) ) const;
   
   real_t delta( const Vector3< real_t > & fluid, const Vector3< real_t > & boundary ) const;

private:

   Setup setup_;

}; // class Cylinder

bool Cylinder::contains( const Vector3< real_t > & point ) const
{
   const real_t px = setup_.cylinderxPosition;
   const real_t py = setup_.cylinderyPosition;
   const real_t H = setup_.H;
   const real_t r = setup_.cylinderRadius;

   if( setup_.circularCrossSection )
   {
      const real_t xd = point[0] - px;
      const real_t yd = point[1] - py;
      const real_t d = xd * xd + yd * yd;
      return point[2] > real_t(0) && point[2] < H && d <= ( r * r );
   }
   else
   {
      const AABB cylinder( px - r, py - r, real_t(0), px + r, py + r, H );
      return cylinder.contains( point );
   }
}
   
bool Cylinder::contains( const AABB & aabb ) const
{
   if( setup_.circularCrossSection )
   {
      Vector3< real_t > p[8];
      p[0].set( aabb.xMin(), aabb.yMin(), aabb.zMin() );
      p[1].set( aabb.xMax(), aabb.yMin(), aabb.zMin() );
      p[2].set( aabb.xMin(), aabb.yMax(), aabb.zMin() );
      p[3].set( aabb.xMax(), aabb.yMax(), aabb.zMin() );
      p[4].set( aabb.xMin(), aabb.yMin(), aabb.zMax() );
      p[5].set( aabb.xMax(), aabb.yMin(), aabb.zMax() );
      p[6].set( aabb.xMin(), aabb.yMax(), aabb.zMax() );
      p[7].set( aabb.xMax(), aabb.yMax(), aabb.zMax() );
      return contains( p[0] ) && contains( p[1] ) && contains( p[2] ) && contains( p[3] ) &&
             contains( p[4] ) && contains( p[5] ) && contains( p[6] ) && contains( p[7] );
   }
   else
   {
      return contains( aabb.min() ) && contains( aabb.max() );
   }
}

bool Cylinder::intersects( const AABB & aabb, const real_t bufferDistance ) const
{
   const real_t px = setup_.cylinderxPosition;
   const real_t py = setup_.cylinderyPosition;
   const real_t r = setup_.cylinderRadius;
   
   if( setup_.circularCrossSection )
   {
      Vector3< real_t > p( px, py, real_t(0) );

      if( p[0] < aabb.xMin() ) p[0] = aabb.xMin();
      else if( p[0] > aabb.xMax() ) p[0] = aabb.xMax();
      if( p[1] < aabb.yMin() ) p[1] = aabb.yMin();
      else if( p[1] > aabb.yMax() ) p[1] = aabb.yMax();

      const real_t xd = p[0] - px;
      const real_t yd = p[1] - py;
      const real_t d = xd * xd + yd * yd;
      const real_t rr = ( r + bufferDistance ) * ( r + bufferDistance );
      return d <= rr;
   }
   else
   {
      const AABB cylinder( px - r, py - r, real_t(0), px + r, py + r, setup_.H );
      return cylinder.intersects( aabb, bufferDistance );
   }
}
   
real_t Cylinder::delta( const Vector3< real_t > & fluid, const Vector3< real_t > & boundary ) const
{
   WALBERLA_CHECK( !contains( fluid ) );
   WALBERLA_CHECK( contains( boundary ) );
   
   const real_t px = setup_.cylinderxPosition;
   const real_t py = setup_.cylinderyPosition;
   const real_t r = setup_.cylinderRadius;
   
   if( setup_.circularCrossSection )
   {
      // http://devmag.org.za/2009/04/17/basic-collision-detection-in-2d-part-2/
      
      const Vector3< real_t > circle( px, py, real_t(0) );
      
      const Vector3< real_t > f = fluid - circle;
      const Vector3< real_t > d = ( boundary - circle ) - f;
      
      const real_t a = d[0] * d[0] + d[1] * d[1];
      const real_t b = real_t(2) * ( d[0] * f[0] + d[1] * f[1] );
      const real_t c = f[0] * f[0] + f[1] * f[1] - r * r;
      
      const real_t bb4ac = b * b - ( real_t(4) * a * c );
      WALBERLA_CHECK_GREATER_EQUAL( bb4ac, real_t(0) );
      
      const real_t sqrtbb4ac = std::sqrt( bb4ac );
      
      const real_t alpha = std::min( ( -b + sqrtbb4ac ) / ( real_t(2) * a ), ( -b - sqrtbb4ac ) / ( real_t(2) * a ) );
      
      WALBERLA_CHECK_GREATER_EQUAL( alpha, real_t(0) );
      WALBERLA_CHECK_LESS_EQUAL( alpha, real_t(1) );
      
      return alpha;
   }
   
   const AABB cylinder( px - r, py - r, real_t(0), px + r, py + r, setup_.H );
   
   if( fluid[0] <= cylinder.xMin() )
   {
      const real_t xdiff = cylinder.xMin() - fluid[0];
      
      if( fluid[1] <= cylinder.yMin() )
      {
         const real_t ydiff = cylinder.yMin() - fluid[1];
         if( xdiff >= ydiff )
         {
            WALBERLA_CHECK_LESS_EQUAL( fluid[0], boundary[0] );
            return xdiff / ( boundary[0] - fluid[0] );
         }
         WALBERLA_CHECK_LESS_EQUAL( fluid[1], boundary[1] );
         return ydiff / ( boundary[1] - fluid[1] );
      }
      else if( fluid[1] >= cylinder.yMax() )
      {
         const real_t ydiff = fluid[1] - cylinder.yMax();
         if( xdiff >= ydiff )
         {
            WALBERLA_CHECK_LESS_EQUAL( fluid[0], boundary[0] );
            return xdiff / ( boundary[0] - fluid[0] );
         }
         WALBERLA_CHECK_GREATER_EQUAL( fluid[1], boundary[1] );
         return ydiff / ( fluid[1] - boundary[1] );
      }
   
      WALBERLA_CHECK_LESS_EQUAL( fluid[0], boundary[0] );
      return xdiff / ( boundary[0] - fluid[0] );
   }
   else if( fluid[0] >= cylinder.xMax() )
   {
      const real_t xdiff = fluid[0] - cylinder.xMax();
      
      if( fluid[1] <= cylinder.yMin() )
      {
         const real_t ydiff = cylinder.yMin() - fluid[1];
         if( xdiff >= ydiff )
         {
            WALBERLA_CHECK_GREATER_EQUAL( fluid[0], boundary[0] );
            return xdiff / ( fluid[0] - boundary[0] );
         }
         WALBERLA_CHECK_LESS_EQUAL( fluid[1], boundary[1] );
         return ydiff / ( boundary[1] - fluid[1] );
      }
      else if( fluid[1] >= cylinder.yMax() )
      {
         const real_t ydiff = fluid[1] - cylinder.yMax();
         if( xdiff >= ydiff )
         {
            WALBERLA_CHECK_GREATER_EQUAL( fluid[0], boundary[0] );
            return xdiff / ( fluid[0] - boundary[0] );
         }
         WALBERLA_CHECK_GREATER_EQUAL( fluid[1], boundary[1] );
         return ydiff / ( fluid[1] - boundary[1] );
      }
      
      WALBERLA_CHECK_GREATER_EQUAL( fluid[0], boundary[0] );
      return xdiff / ( fluid[0] - boundary[0] );
   }

   if( fluid[1] <= cylinder.yMin() )
   {
      WALBERLA_CHECK_LESS_EQUAL( fluid[1], boundary[1] );
      return ( cylinder.yMin() - fluid[1] ) / ( boundary[1] - fluid[1] );
   }
   
   WALBERLA_CHECK_GREATER_EQUAL( fluid[1], cylinder.yMax() );
   WALBERLA_CHECK_GREATER_EQUAL( fluid[1], boundary[1] );
   return ( fluid[1] - cylinder.yMax() ) / ( fluid[1] - boundary[1] );
}



class CylinderRefinementSelection
{
public:

   CylinderRefinementSelection( const Cylinder & cylinder, const uint_t level, const real_t bufferDistance ) :
      cylinder_( cylinder ), level_( level ), bufferDistance_( bufferDistance ) {}

   void operator()( SetupBlockForest & forest )
   {
      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         const AABB & aabb = block->getAABB();

         if( block->getLevel() < level_ && cylinder_.intersects( aabb, bufferDistance_ ) && !cylinder_.contains( aabb ) )
            block->setMarker( true );
      }
   }

private:

   Cylinder cylinder_;
   uint_t level_;
   real_t bufferDistance_;

}; // class CylinderRefinementSelection



// check refinement at inflow/outflow and make sure that ALL blocks at these two boundaries are at the same level!
void setInflowOutflowToSameLevel( SetupBlockForest & forest, const Setup & setup )
{
   uint_t maxInflowLevel( uint_t(0) );
   uint_t maxOutflowLevel( uint_t(0) );

   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      if( forest.atDomainXMinBorder(*block) )
      {
         maxInflowLevel = std::max( maxInflowLevel, block->getLevel() + ( block->isMarked() ? uint_t(1) : uint_t(0) ) );
      }
      if( setup.strictlyObeyL )
      {
         AABB aabb = block->getAABB();
         aabb.extend( real_c( FieldGhostLayers ) * ( setup.dx / real_c( math::uintPow2( block->getLevel() ) ) ) );
         auto p = aabb.center();
         p[0] = setup.L;
         if( aabb.contains(p) )
         {
            maxOutflowLevel = std::max( maxOutflowLevel, block->getLevel() + ( block->isMarked() ? uint_t(1) : uint_t(0) ) );
         }
      }
      else
      {
         if( forest.atDomainXMaxBorder(*block) )
            maxOutflowLevel = std::max( maxOutflowLevel, block->getLevel() + ( block->isMarked() ? uint_t(1) : uint_t(0) ) );
      }
   }

   mpi::allReduceInplace( maxInflowLevel, mpi::MAX );
   mpi::allReduceInplace( maxOutflowLevel, mpi::MAX );

   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      if( forest.atDomainXMinBorder(*block) )
      {
         if( block->getLevel() < maxInflowLevel && !(block->isMarked()) )
            block->setMarker( true );
      }
      if( setup.strictlyObeyL )
      {
         AABB aabb = block->getAABB();
         aabb.extend( real_c( FieldGhostLayers ) * ( setup.dx / real_c( math::uintPow2( block->getLevel() ) ) ) );
         auto p = aabb.center();
         p[0] = setup.L;
         if( aabb.contains(p) && block->getLevel() < maxOutflowLevel && !(block->isMarked()) )
         {
            block->setMarker( true );
         }
      }
      else
      {
         if( forest.atDomainXMaxBorder(*block) && block->getLevel() < maxOutflowLevel && !(block->isMarked()) )
            block->setMarker( true );
      }
   }
}



void Pseudo2DRefinementSelectionCorrection( SetupBlockForest & forest )
{
   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      if( ! forest.atDomainZMinBorder(*block) )
         block->setMarker( false );
   }
}



static void workloadMemoryAndSUIDAssignment( SetupBlockForest & forest, const memory_t memoryPerBlock, const Setup & setup )
{
   if( setup.pseudo2D )
   {
      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         if( forest.atDomainZMinBorder(*block) )
         {
            block->setWorkload( numeric_cast< workload_t >( uint_t(1) << block->getLevel() ) );
            block->setMemory( memoryPerBlock );
         }
         else
         {
            block->setWorkload( workload_t(0) );
            block->setMemory( memory_t(0) );
            block->addState( Empty );
         }
      }
   }
   else
   {
      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         block->setWorkload( numeric_cast< workload_t >( uint_t(1) << block->getLevel() ) );
         block->setMemory( memoryPerBlock );
      }
   }
}



static shared_ptr< SetupBlockForest > createSetupBlockForest( const blockforest::RefinementSelectionFunctions & refinementSelectionFunctions,
                                                              const Setup & setup,
                                                              uint_t numberOfProcesses, const uint_t blocksPerProcess,
                                                              const memory_t memoryPerCell, const memory_t processMemoryLimit,
                                                              const bool outputSetupForest )
{
   shared_ptr< SetupBlockForest > forest = make_shared< SetupBlockForest >();

   const memory_t memoryPerBlock = numeric_cast< memory_t >( ( ( setup.pseudo2D ? uint_t(1) : setup.yzCells ) + uint_t(2) * FieldGhostLayers ) *
                                                             ( setup.yzCells + uint_t(2) * FieldGhostLayers ) *
                                                             ( setup.xCells + uint_t(2) * FieldGhostLayers ) ) * memoryPerCell;

   forest->addRefinementSelectionFunction( refinementSelectionFunctions );
   forest->addWorkloadMemorySUIDAssignmentFunction( std::bind( workloadMemoryAndSUIDAssignment, std::placeholders::_1, memoryPerBlock, std::cref( setup ) ) );

   forest->init( AABB( real_c(0), real_c(0), real_c(0),
                       setup.H * ( real_c(setup.xBlocks) * real_c(setup.xCells) ) / ( real_c(setup.yzBlocks) * real_c(setup.yzCells) ), setup.H, setup.H ),
                 setup.xBlocks, setup.yzBlocks, setup.pseudo2D ? uint_t(1) : setup.yzBlocks, false, false, false );

   MPIManager::instance()->useWorldComm();

   if( blocksPerProcess != 0 )
      numberOfProcesses = uint_c( std::ceil( real_c( forest->getNumberOfBlocks() ) / real_c( blocksPerProcess ) ) );

   forest->balanceLoad( blockforest::StaticLevelwiseCurveBalanceWeighted(true), numberOfProcesses, real_t(0), processMemoryLimit, true, false );

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

      auto sbf = std::make_shared< StructuredBlockForest >( bf, setup.xCells, setup.yzCells, setup.pseudo2D ? uint_t(1) : setup.yzCells );
      sbf->createCellBoundingBoxes();

      return sbf;
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure ..." );

   shared_ptr< SetupBlockForest > sforest = createSetupBlockForest( refinementSelectionFunctions, setup,
                                                                    uint_c( MPIManager::instance()->numProcesses() ), uint_t(0),
                                                                    memoryPerCell, processMemoryLimit,
                                                                    configBlock.getParameter< bool >( "outputSetupForest", false ) );

   auto bf = std::make_shared< blockforest::BlockForest >( uint_c( MPIManager::instance()->rank() ), *sforest, false );

   auto sbf = std::make_shared< blockforest::StructuredBlockForest >( bf, setup.xCells, setup.yzCells, setup.pseudo2D ? uint_t(1) : setup.yzCells );
   sbf->createCellBoundingBoxes();
   
   return sbf;
}






///////////////////////
// BOUNDARY HANDLING //
///////////////////////



template< bool Pseudo2D = false >
class SinusInflowVelocity
{
public:

   SinusInflowVelocity( const real_t velocity, const real_t raisingTime, const real_t sinPeriod, const real_t H ) :
      raisingTime_( raisingTime ), sinPeriod_( sinPeriod ), H_( H )
   {
      uTerm_ = Pseudo2D ? ( real_t(4) * velocity ) : ( real_t(16) * velocity );
      HTerm_ = Pseudo2D ? ( real_t(1) / ( H * H ) ) : ( real_t(1) / ( H * H * H * H ) );
      tConstTerm_ = uTerm_ * HTerm_;
   }

   void operator()( const real_t t )
   {
      tConstTerm_ = ( sinPeriod_ > real_t(0) ) ? ( std::abs( std::sin( math::pi * t / sinPeriod_ ) ) ) : real_t(1);
      tConstTerm_ *= uTerm_ * HTerm_;
      tConstTerm_ *= ( raisingTime_ > real_t(0) ) ? std::min( t / raisingTime_, real_t(1) ) : real_t(1);
   }

   Vector3< real_t > operator()( const Vector3< real_t > & x, const real_t )
   {
      return Vector3< real_t >( Pseudo2D ? ( tConstTerm_ * x[1] * ( H_ - x[1] ) ) :
                                           ( tConstTerm_ * x[1] * x[2] * ( H_ - x[1] ) * ( H_ - x[2] ) ), real_t(0), real_t(0) );
   }

private:

   real_t raisingTime_;
   real_t sinPeriod_;
   real_t H_;

   real_t uTerm_;
   real_t HTerm_;

   real_t tConstTerm_;

};






template< typename LatticeModel_T >
struct MyBoundaryTypes
{
   using NoSlip_T = lbm::NoSlip<LatticeModel_T, flag_t>;
   using Obstacle_T = lbm::NoSlip<LatticeModel_T, flag_t>;
   using Curved_T = lbm::Curved<LatticeModel_T, FlagField_T>;
   using DynamicUBB_T = lbm::DynamicUBB<LatticeModel_T, flag_t, SinusInflowVelocity<Is2D<LatticeModel_T>::value>>;
   using Outlet21_T = lbm::Outlet<LatticeModel_T, FlagField_T, 2, 1>;
   using Outlet43_T = lbm::Outlet<LatticeModel_T, FlagField_T, 4, 3>;
   using PressureOutlet_T = lbm::SimplePressure<LatticeModel_T, flag_t>;

   using BoundaryHandling_T = BoundaryHandling<FlagField_T, typename Types<LatticeModel_T>::Stencil_T, NoSlip_T, Obstacle_T, Curved_T, DynamicUBB_T, Outlet21_T, Outlet43_T, PressureOutlet_T>;
};


template< typename LatticeModel_T >
class MyBoundaryHandling : public blockforest::AlwaysInitializeBlockDataHandling< typename MyBoundaryTypes< LatticeModel_T >::BoundaryHandling_T >
{
public:

   using NoSlip_T = typename MyBoundaryTypes<LatticeModel_T>::NoSlip_T;
   using Obstacle_T = typename MyBoundaryTypes<LatticeModel_T>::Obstacle_T;
   using Curved_T = typename MyBoundaryTypes<LatticeModel_T>::Curved_T;
   using DynamicUBB_T = typename MyBoundaryTypes<LatticeModel_T>::DynamicUBB_T;
   using Outlet21_T = typename MyBoundaryTypes<LatticeModel_T>::Outlet21_T;
   using Outlet43_T = typename MyBoundaryTypes<LatticeModel_T>::Outlet43_T;
   using PressureOutlet_T = typename MyBoundaryTypes<LatticeModel_T>::PressureOutlet_T;

   using BoundaryHandling_T = typename MyBoundaryTypes<LatticeModel_T>::BoundaryHandling_T;



   MyBoundaryHandling( const BlockDataID & flagFieldId, const BlockDataID & pdfFieldId,
                       const weak_ptr< StructuredBlockStorage > & blocks,
                       const Setup & setup, const shared_ptr< lbm::TimeTracker > & timeTracker ) :
      flagFieldId_( flagFieldId ), pdfFieldId_( pdfFieldId ), blocks_( blocks ), setup_( setup ), timeTracker_( timeTracker )
   {}

   BoundaryHandling_T * initialize( IBlock * const block ) override;

private:

   const BlockDataID flagFieldId_;
   const BlockDataID  pdfFieldId_;

   weak_ptr< StructuredBlockStorage > blocks_;

   Setup setup_;
   shared_ptr< lbm::TimeTracker > timeTracker_;

}; // class MyBoundaryHandling

template< typename LatticeModel_T >
typename MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T *
MyBoundaryHandling<LatticeModel_T>::initialize( IBlock * const block )
{
   using PdfField_T = typename Types< LatticeModel_T >::PdfField_T;

   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfFieldId_ );

   const flag_t fluid = flagField->registerFlag( Fluid_Flag );

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks );

   SinusInflowVelocity<Is2D< LatticeModel_T >::value> velocity( setup_.inflowVelocity_L, setup_.raisingTime_L, setup_.sinPeriod_L, setup_.H );

   return new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                    NoSlip_T( "no slip", NoSlip_Flag, pdfField ),
                                  Obstacle_T( "obstacle (staircase)", Obstacle_Flag, pdfField ),
                                    Curved_T( "obstacle (curved)", Curved_Flag, pdfField, flagField, fluid ),
                                DynamicUBB_T( "velocity bounce back", UBB_Flag, pdfField, timeTracker_, blocks->getLevel(*block), velocity, block->getAABB() ),
                                  Outlet21_T( "outlet (2/1)", Outlet21_Flag, pdfField, flagField, fluid ),
                                  Outlet43_T( "outlet (4/3)", Outlet43_Flag, pdfField, flagField, fluid ),
                            PressureOutlet_T( "pressure outlet", PressureOutlet_Flag, pdfField, real_t(1) ) );
}



template< typename LatticeModel_T >
class CurvedDeltaValueCalculation
{
public:

   using Stencil = typename LatticeModel_T::Stencil;
   
   CurvedDeltaValueCalculation( const shared_ptr< StructuredBlockForest > & blocks, const IBlock & block, const Cylinder & cylinder ) :
      blocks_( blocks ), block_( block ), cylinder_( cylinder ) {}

   shared_ptr< BoundaryConfiguration > operator()( const Cell & boundaryCell, const Vector3< real_t > & p );

private:

   const shared_ptr< StructuredBlockForest > & blocks_;
   const IBlock & block_;

   const Cylinder & cylinder_;

}; // class CurvedDeltaValueCalculation

template< typename LatticeModel_T >
shared_ptr< BoundaryConfiguration > CurvedDeltaValueCalculation< LatticeModel_T >::operator()( const Cell & boundaryCell, const Vector3< real_t > & /*p*/ )
{
   std::vector< real_t > deltas( Stencil::Size, real_c(0.5) );
   
   const Vector3< real_t > boundary = blocks_->getBlockLocalCellCenter( block_, boundaryCell );
   if( cylinder_.contains( boundary ) )
   {
      for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
      {
         const Cell fluidCell = Cell( boundaryCell.x() + dir.cx(), boundaryCell.y() + dir.cy(), boundaryCell.z() + dir.cz() );
         const Vector3< real_t > fluid = blocks_->getBlockLocalCellCenter( block_, fluidCell );
         if( !cylinder_.contains( fluid ) && blocks_->getDomain().contains( fluid ) )
            deltas[ dir.toIdx() ] = cylinder_.delta( fluid, boundary );
      }
   }
   
   return shared_ptr< BoundaryConfiguration >( new typename MyBoundaryHandling< LatticeModel_T >::Curved_T::Configuration( deltas ) );
}



template< typename LatticeModel_T >
class BoundarySetter
{
public:

   using BoundaryHandling_T = typename MyBoundaryTypes< LatticeModel_T >::BoundaryHandling_T;

   BoundarySetter( const weak_ptr<StructuredBlockForest> & blockForest, const BlockDataID & boundaryHandlingId, const Setup & setup,
                   const int obstacleBoundary, const int outletType,
                   const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors ) :
      blockForest_( blockForest ), boundaryHandlingId_( boundaryHandlingId ), setup_( setup ),
      obstacleBoundary_( obstacleBoundary ), outletType_( outletType ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {}

   void operator()();

private:

   weak_ptr<StructuredBlockForest> blockForest_;
   
   BlockDataID boundaryHandlingId_;
   
   Setup setup_;
   
   int obstacleBoundary_;
   int outletType_;
   
   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;   
};

template< typename LatticeModel_T >
void BoundarySetter< LatticeModel_T >::operator()()
{
   auto blocks = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks );

   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
   {
      BoundaryHandling_T * boundaryHandling = block->template getData< BoundaryHandling_T >( boundaryHandlingId_ );

      const uint_t level = blocks->getLevel(*block);

      CellInterval domainBB = blocks->getDomainCellBB( level );
      domainBB.expand( cell_idx_t(1) );
      if( setup_.strictlyObeyL )
      {
         CellInterval tempBB;
         blocks->getCellBBFromAABB( tempBB, AABB( real_t(0), real_t(0), real_t(0), setup_.L, setup_.H, setup_.H ), uint_t(0) );
         domainBB.xMax() = std::min( ( tempBB.xMax() + cell_idx_t(1) ) << level, domainBB.xMax() );
      }
      blocks->transformGlobalToBlockLocalCellInterval( domainBB, *block );

      boundaryHandling->forceDomain( domainBB );

      // cylinder

      Cylinder cylinder( setup_ );

      if( obstacleBoundary_ == 1 ) // curved boundary
      {
         CurvedDeltaValueCalculation< LatticeModel_T > deltaCalculation( blocks, *block, cylinder );

         lbm::refinement::consistentlyForceBoundary< BoundaryHandling_T >( *blocks, dynamic_cast< blockforest::Block & >(*block),
                                                                           boundaryHandlingId_, Curved_Flag, cylinder, deltaCalculation );
      }
      else // staircase
      {
         lbm::refinement::consistentlyForceBoundary< BoundaryHandling_T >( *blocks, dynamic_cast< blockforest::Block & >(*block),
                                                                           boundaryHandlingId_, Obstacle_Flag, cylinder );
      }

      // inflow WEST

      CellInterval west( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
      boundaryHandling->forceBoundary( UBB_Flag, west );

      // outlet EAST

      CellInterval east( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      if( outletType_ == 0 )
         boundaryHandling->forceBoundary( PressureOutlet_Flag, east );
      else if( outletType_ == 1 )
         boundaryHandling->forceBoundary( Outlet21_Flag, east );
      else
         boundaryHandling->forceBoundary( Outlet43_Flag, east );

      // no slip SOUTH
      CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
      boundaryHandling->forceBoundary( NoSlip_Flag, south );

      // no slip NORTH
      CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      boundaryHandling->forceBoundary( NoSlip_Flag, north );

      // no slip BOTTOM
      CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
      boundaryHandling->forceBoundary( NoSlip_Flag, bottom );

      // no slip TOP
      CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      boundaryHandling->forceBoundary( NoSlip_Flag, top );
   }
}






////////////////////////
// DYNAMIC REFINEMENT //
////////////////////////

class Pseudo2DBlockStateDetermination // used as a 'PhantomBlockForest::BlockStateDeterminationFunction'
{
public:

   Pseudo2DBlockStateDetermination( const BlockForest & forest, const SUID & state, const bool markEmptyBlocks = true ) :
      forest_( forest ), markEmptyBlocks_( markEmptyBlocks ), state_( state )
   {}

   Set<SUID> operator()( const std::vector< std::pair< BlockID, Set<SUID> > > & source, const BlockID & destintation );

private:

   const BlockForest & forest_;
   bool markEmptyBlocks_; // false = non-empty blocks are set to state_
   SUID state_; 
};

Set<SUID> Pseudo2DBlockStateDetermination::operator()( const std::vector< std::pair< BlockID, Set<SUID> > > & source, const BlockID & destintation )
{
   auto & domain = forest_.getDomain();
   auto aabb = forest_.getAABBFromBlockId( destintation );
   
   Set<SUID> state;
   for( auto it = source.begin(); it != source.end(); ++it )
   {
      for( auto suid = it->second.begin(); suid != it->second.end(); ++suid )
         if( *suid != state_ )
            state += *suid;
   }
   
   if( markEmptyBlocks_ )
   {
      if( ! realIsEqual( aabb.zMin(), domain.zMin(), real_c( 1.0E-6 ) * ( aabb.zMax() - aabb.zMin() ) ) )
         state += state_;
   }
   else
   {
      if( realIsEqual( aabb.zMin(), domain.zMin(), real_c( 1.0E-6 ) * ( aabb.zMax() - aabb.zMin() ) ) )
         state += state_;
   }
   
   return state;
}


// used as a 'BlockForest::RefreshMinTargetLevelDeterminationFunction
void keepInflowOutflowAtTheSameLevel( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                      std::vector< const Block * > & blocksAlreadyMarkedForRefinement,
                                      const blockforest::BlockForest & forest, const Setup & setup )
{
   uint_t maxInflowLevel( uint_t(0) );
   uint_t maxOutflowLevel( uint_t(0) );

   // In addition to keeping in- and outflow blocks at the same level, this callback also
   // prevents these blocks from coarsening.

   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      const Block * const block = it->first;
      if( forest.atDomainXMinBorder(*block) )
      {
         it->second = std::max( it->second, block->getLevel() );
         maxInflowLevel = std::max( maxInflowLevel, it->second );
      }
      if( setup.strictlyObeyL )
      {
         AABB aabb = block->getAABB();
         aabb.extend( real_c( FieldGhostLayers ) * ( setup.dx / real_c( math::uintPow2( block->getLevel() ) ) ) );
         auto p = aabb.center();
         p[0] = setup.L;
         if( aabb.contains(p) )
         {
            it->second = std::max( it->second, block->getLevel() );
            maxOutflowLevel = std::max( maxOutflowLevel, it->second );
         }
      }
      else if( forest.atDomainXMaxBorder(*block) )
      {
         it->second = std::max( it->second, block->getLevel() );
         maxOutflowLevel = std::max( maxOutflowLevel, it->second );
      }
   }
   for( auto it = blocksAlreadyMarkedForRefinement.begin(); it != blocksAlreadyMarkedForRefinement.end(); ++it )
   {
      const Block * const block = *it;
      if( forest.atDomainXMinBorder(*block) )
      {
         maxInflowLevel = std::max( maxInflowLevel, block->getTargetLevel() );
      }
      if( setup.strictlyObeyL )
      {
         AABB aabb = block->getAABB();
         aabb.extend( real_c( FieldGhostLayers ) * ( setup.dx / real_c( math::uintPow2( block->getLevel() ) ) ) );
         auto p = aabb.center();
         p[0] = setup.L;
         if( aabb.contains(p) )
            maxOutflowLevel = std::max( maxOutflowLevel, block->getTargetLevel() );
      }
      else if( forest.atDomainXMaxBorder(*block) )
         maxOutflowLevel = std::max( maxOutflowLevel, block->getTargetLevel() );
   }

   mpi::allReduceInplace( maxInflowLevel, mpi::MAX );
   mpi::allReduceInplace( maxOutflowLevel, mpi::MAX );

   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      const Block * const block = it->first;
      if( forest.atDomainXMinBorder(*block) )
      {
         it->second = maxInflowLevel;
      }
      if( setup.strictlyObeyL )
      {
         AABB aabb = block->getAABB();
         aabb.extend( real_c( FieldGhostLayers ) * ( setup.dx / real_c( math::uintPow2( block->getLevel() ) ) ) );
         auto p = aabb.center();
         p[0] = setup.L;
         if( aabb.contains(p) )
            it->second = maxOutflowLevel;
      }
      else if( forest.atDomainXMaxBorder(*block) )
         it->second = maxOutflowLevel;
   }
}



// used as a 'BlockForest::RefreshMinTargetLevelDeterminationFunction
void pseudo2DTargetLevelCorrection( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                    std::vector< const Block * > &, const blockforest::BlockForest & forest )
{
   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      if( ! forest.atDomainZMinBorder( *(it->first ) ) )
         it->second = uint_t(0);
   }
}



class Pseudo2DPhantomWeight // used as a 'PhantomBlockForest::PhantomBlockDataAssignmentFunction'
{
public:

   using weight_t = uint8_t;

   Pseudo2DPhantomWeight( const weight_t _weight ) : weight_( _weight ) {}

   weight_t weight() const { return weight_; }

private:

   weight_t weight_;
};

class Pseudo2DPhantomWeightAssignment
{
public:

   Pseudo2DPhantomWeightAssignment( const SUID & state, const bool markEmptyBlocks = true ) :
      markEmptyBlocks_( markEmptyBlocks ), state_( state )
   {}
   
   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & )
   {
      for( auto it = blockData.begin(); it != blockData.end(); ++it )
      {
         if( it->first->getState().contains( state_ ) )
            it->second = Pseudo2DPhantomWeight( markEmptyBlocks_ ? uint8_t(0) : uint8_t(1) );
         else
            it->second = Pseudo2DPhantomWeight( markEmptyBlocks_ ? uint8_t(1) : uint8_t(0) );
      }
   }

private:

   bool markEmptyBlocks_; // false = state_ marks non-empty blocks
   SUID state_; 
};

struct Pseudo2DPhantomWeightPackUnpack
{
   void operator()( mpi::SendBuffer & buffer, const PhantomBlock & block )
   {
      buffer << block.getData< Pseudo2DPhantomWeight >().weight();
   }

   void operator()( mpi::RecvBuffer & buffer, const PhantomBlock &, walberla::any & data )
   {
      Pseudo2DPhantomWeight::weight_t w;
      buffer >> w;
      data = Pseudo2DPhantomWeight( w );
   }
};






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

   writers.push_back( make_shared< lbm::VelocityVTKWriter<LatticeModel_T>            >( pdfField_, "VelocityFromPDF" ) );
   writers.push_back( make_shared< lbm::VelocityMagnitudeVTKWriter<LatticeModel_T>   >( pdfField_, "VelocityMagnitudeFromPDF" ) );
   writers.push_back( make_shared< lbm::DensityVTKWriter<LatticeModel_T>             >( pdfField_, "DensityFromPDF" ) );
   writers.push_back( make_shared< lbm::NonEqulibriumVTKWriter<LatticeModel_T>       >( pdfField_, "NonEquPart" ) );
   writers.push_back( make_shared< field::VTKWriter< lbm::PdfField<LatticeModel_T> > >( pdfField_, "PDF" ) );

   writers.push_back( make_shared< field::VTKWriter< FlagField_T > >( flagField_, "FlagField" ) );

   // cell filters

   field::FlagFieldCellFilter<FlagField_T> fluidFilter( flagField_ );
   fluidFilter.addFlag( Fluid_Flag );
   filters[ "FluidFilter" ] = fluidFilter;

   field::FlagFieldCellFilter<FlagField_T> obstacleFilter( flagField_ );
   obstacleFilter.addFlag(         NoSlip_Flag );
   obstacleFilter.addFlag(       Obstacle_Flag );
   obstacleFilter.addFlag(         Curved_Flag );
   obstacleFilter.addFlag(            UBB_Flag );
   obstacleFilter.addFlag( PressureOutlet_Flag );
   obstacleFilter.addFlag(       Outlet21_Flag );
   obstacleFilter.addFlag(       Outlet43_Flag );
   filters[ "ObstacleFilter" ] = obstacleFilter;

   // before functions

   beforeFunctions[ "PDFGhostLayerSync" ] = pdfGhostLayerSync_;
}






////////////////
// Evaluation //
////////////////

template< typename LatticeModel_T >
class Evaluation
{
public:

   using PdfField_T = typename Types< LatticeModel_T >::PdfField_T;
   using Stencil_T = typename LatticeModel_T::Stencil;

   Evaluation( const weak_ptr< StructuredBlockStorage > & blocks, const uint_t checkFrequency,
               const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const FlagUID & fluid, const FlagUID & obstacle,
               const Setup & setup,
               const bool logToStream = true, const bool logToFile = true, const std::string& filename = std::string("SchaeferTurek.txt"),
               const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
               const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      initialized_( false ), blocks_( blocks ),
      executionCounter_( uint_t(0) ), checkFrequency_( checkFrequency ), pdfFieldId_( pdfFieldId ), flagFieldId_( flagFieldId ),
      fluid_( fluid ), obstacle_( obstacle ), setup_( setup ), D_( uint_t(0) ), AD_( real_t(0) ), AL_( real_t(0) ), forceEvaluationExecutionCount_( uint_t(0) ),
      strouhalRising_( false ), strouhalNumberRealD_( real_t(0) ), strouhalNumberDiscreteD_( real_t(0) ), strouhalEvaluationExecutionCount_( uint_t(0) ),
      logToStream_( logToStream ), logToFile_( logToFile ), filename_( filename ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {
      forceSample_.resize( uint_t(2) );
      coefficients_.resize( uint_t(4) );
      coefficientExtremas_.resize( uint_t(4) );

      WALBERLA_ROOT_SECTION()
      {
         if( logToFile_ )
         {
            std::ofstream file( filename_.c_str() );
            file << "# time step [1], force (x) [2], force (y) [3], force (z) [4], "
                    "cD (real area) [5], cL (real area) [6], cD (discrete area) [7], cL (discrete area) [8], "
                    "pressure difference (in lattice units) [9], pressure difference (in Pa) [10], vortex velocity (in lattice units) [11], "
                    "Strouhal number (real D) [12], Strouhal number (discrete D) [13]" << std::endl;
            if( !setup_.evaluatePressure )
               file << "# ATTENTION: pressure was not evaluated, pressure difference is set to zero!" << std::endl;
            if( !setup_.evaluateStrouhal )
               file << "# ATTENTION: vortex velocities were not evaluated, Strouhal number is set to zero!" << std::endl;
            file.close();
         }
      }
   }

   void operator()();
   
   void operator()( const uint_t level, const uint_t executionCount ); // for calculating the force
   void operator()( IBlock * block, const uint_t level, const uint_t executionCount ); // for calculating the force

   void prepareResultsForSQL();
   void getResultsForSQLOnRoot( std::map< std::string, double > & realProperties, std::map< std::string, int > & integerProperties ) const;
   
   void check( const shared_ptr< Config > & config );
   
   void refresh();
   
protected:

   void evaluate( real_t & cDRealArea, real_t & cLRealArea, real_t & cDDiscreteArea, real_t & cLDiscreteArea,
                  real_t & pressureDifference_L, real_t & pressureDifference );



   bool initialized_;

   weak_ptr< StructuredBlockStorage > blocks_;
   std::map< IBlock *, std::vector< std::pair< Cell, stencil::Direction > > >  directions_;

   uint_t executionCounter_;
   uint_t checkFrequency_;
   
   BlockDataID pdfFieldId_;
   BlockDataID flagFieldId_;

   FlagUID fluid_;
   FlagUID obstacle_;
   
   Setup setup_;

   uint_t D_;
   real_t AD_;
   real_t AL_;
   
   Vector3< real_t > force_;
   std::vector< math::Sample > forceSample_;
   uint_t forceEvaluationExecutionCount_;
   
   std::vector< std::deque< real_t > > coefficients_;
   std::vector< std::pair< real_t, real_t > > coefficientExtremas_;

   std::vector< real_t > strouhalVelocities_;
   std::vector< uint_t > strouhalTimeStep_;
   bool strouhalRising_;
   real_t strouhalNumberRealD_;
   real_t strouhalNumberDiscreteD_;
   uint_t strouhalEvaluationExecutionCount_;

   bool logToStream_;
   bool logToFile_;
   std::string filename_;
   
   std::map< std::string, double > sqlResults_;
   
   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // class Evaluation



template< typename LatticeModel_T >
void Evaluation< LatticeModel_T >::operator()()
{
   if( checkFrequency_ == uint_t(0) )
      return;

   ++executionCounter_;
   if( ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
      return;

   real_t cDRealArea( real_t(0) );
   real_t cLRealArea( real_t(0) );
   real_t cDDiscreteArea( real_t(0) );
   real_t cLDiscreteArea( real_t(0) );
   
   real_t pressureDifference_L( real_t(0) );
   real_t pressureDifference( real_t(0) );
      
   evaluate( cDRealArea, cLRealArea, cDDiscreteArea, cLDiscreteArea, pressureDifference_L, pressureDifference );

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks );

   // Strouhal number (needs vortex shedding frequency)

   real_t vortexVelocity( real_t(0) );

   if( setup_.evaluateStrouhal )
   {
      auto block = blocks->getBlock( setup_.pStrouhal );
      if( block != nullptr )
      {
         const PdfField_T * const pdfField = block->template getData< PdfField_T >( pdfFieldId_ );
         const auto cell = blocks->getBlockLocalCell( *block, setup_.pStrouhal );
         WALBERLA_ASSERT( pdfField->xyzSize().contains( cell ) );
         vortexVelocity += pdfField->getVelocity( cell )[0];
      }

      mpi::reduceInplace( vortexVelocity, mpi::SUM );
   }

   WALBERLA_ROOT_SECTION()
   {
      coefficients_[0].push_back( cDRealArea );
      coefficients_[1].push_back( cLRealArea );
      coefficients_[2].push_back( cDDiscreteArea );
      coefficients_[3].push_back( cLDiscreteArea );
      
      if( coefficients_[0].size() > setup_.nbrOfEvaluationPointsForCoefficientExtremas )
      {
         for( uint_t i = uint_t(0); i < uint_t(4); ++i )
            coefficients_[i].pop_front();
      }
      
      for( uint_t i = uint_t(0); i < uint_t(4); ++i )
      {
         coefficientExtremas_[i] = std::make_pair( *(coefficients_[i].begin()), *(coefficients_[i].begin()) );
         for( auto v = coefficients_[i].begin(); v != coefficients_[i].end(); ++v )
         {
            coefficientExtremas_[i].first  = std::min( coefficientExtremas_[i].first,  *v );
            coefficientExtremas_[i].second = std::max( coefficientExtremas_[i].second, *v );
         }
      }

      std::ostringstream oss;
      if( setup_.evaluateForceComponents )
      {
         oss << "\nforce components (evaluated in time step " << forceEvaluationExecutionCount_ << "):"
                "\n   x: " << forceSample_[0].min() << " (min), " << forceSample_[0].max() << " (max), "
                           << forceSample_[0].mean() << " (mean), " << forceSample_[0].median() << " (median), "
                           << forceSample_[0].stdDeviation() << " (stdDeviation)" <<
                "\n   y: " << forceSample_[1].min() << " (min), " << forceSample_[1].max() << " (max), "
                           << forceSample_[1].mean() << " (mean), " << forceSample_[1].median() << " (median), "
                           << forceSample_[1].stdDeviation() << " (stdDeviation)";
      }

      if( logToStream_ )
      {
         WALBERLA_LOG_RESULT_ON_ROOT( "force acting on cylinder (in dimensionless lattice units of the coarsest grid - evaluated in time step "
                                      << forceEvaluationExecutionCount_ << "):\n   " << force_ << oss.str() <<
                                      "\ndrag and lift coefficients (including extrema of last " << ( coefficients_[0].size() * checkFrequency_ ) << " time steps):"
                                      "\n   \"real\" area:"
                                      "\n      c_D: " << cDRealArea << " (min = " << coefficientExtremas_[0].first << ", max = " << coefficientExtremas_[0].second << ")" <<
                                      "\n      c_L: " << cLRealArea << " (min = " << coefficientExtremas_[1].first << ", max = " << coefficientExtremas_[1].second << ")" <<
                                      "\n   discrete area:"
                                      "\n      c_D: " << cDDiscreteArea << " (min = " << coefficientExtremas_[2].first << ", max = " << coefficientExtremas_[2].second << ")" <<
                                      "\n      c_L: " << cLDiscreteArea << " (min = " << coefficientExtremas_[3].first << ", max = " << coefficientExtremas_[3].second << ")" );
      }
      
      if( setup_.evaluatePressure && logToStream_ )
      {
         WALBERLA_LOG_RESULT_ON_ROOT( "pressure:\n   difference: " << pressureDifference << " Pa (" << pressureDifference_L  << ")" );
      }
      
      if( setup_.evaluateStrouhal )
      {
         // We evaluate the derivative (-> strouhalRising_) in order to find the local minima and maxima of the velocity over time.
         // If we know the time between a local minimum and a local maximum, we can calculate the frequency.
         // -> "Smooth noise-robust differentiators" (http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/)

         if( strouhalVelocities_.size() < uint_t(11) )
         {
            strouhalVelocities_.push_back( vortexVelocity );
         }
         else
         {
            for( uint_t i = uint_t(0); i < uint_t(10); ++i )
               strouhalVelocities_[i] = strouhalVelocities_[i+1];
            strouhalVelocities_[10] = vortexVelocity;
            
            const real_t f1 = strouhalVelocities_[ 6] - strouhalVelocities_[4];
            const real_t f2 = strouhalVelocities_[ 7] - strouhalVelocities_[3];
            const real_t f3 = strouhalVelocities_[ 8] - strouhalVelocities_[2];
            const real_t f4 = strouhalVelocities_[ 9] - strouhalVelocities_[1];
            const real_t f5 = strouhalVelocities_[10] - strouhalVelocities_[0];
            
            const real_t diff = ( real_c(322) * f1 + real_c(256) * f2 + real_c(39) * f3 - real_c(32) * f4 - real_c(11) * f5 ) / real_c(1536);
            
            if( ( diff > real_t(0) ) != strouhalRising_ )
            {
               strouhalRising_ = ( diff > real_t(0) );
               
               if( strouhalTimeStep_.size() < uint_t(3) )
               {
                  strouhalTimeStep_.push_back( executionCounter_ );
               }
               else
               {
                  strouhalTimeStep_[0] = strouhalTimeStep_[1];
                  strouhalTimeStep_[1] = strouhalTimeStep_[2];
                  strouhalTimeStep_[2] = executionCounter_;
               }
            }
         }

         if( strouhalTimeStep_.size() == uint_t(3) )
         {
            const real_t uMean = Is2D< LatticeModel_T >::value ? ( real_c(2) * setup_.inflowVelocity_L / real_c(3) ) :
                                                                 ( real_c(4) * setup_.inflowVelocity_L / real_c(9) );
            const real_t D = real_t(2) * setup_.cylinderRadius / setup_.dx;
      
            strouhalNumberRealD_     =        D   / ( uMean * real_c( strouhalTimeStep_[2] - strouhalTimeStep_[0] ) );
            strouhalNumberDiscreteD_ = real_c(D_) / ( uMean * real_c( strouhalTimeStep_[2] - strouhalTimeStep_[0] ) );
            
            strouhalEvaluationExecutionCount_ = executionCounter_ - uint_t(1);

            if( logToStream_ )
            {            
               WALBERLA_LOG_RESULT_ON_ROOT( "Strouhal number (evaluated in time step " << strouhalEvaluationExecutionCount_ << "):"
                                            "\n   D/U (in lattice units): " << ( D  / uMean ) << " (\"real\" D), " << ( real_c(D_)  / uMean ) << " (discrete D)"
                                            "\n   T: " << ( real_c( strouhalTimeStep_[2] - strouhalTimeStep_[0] ) * setup_.dt ) << " s ("
                                                       << real_c( strouhalTimeStep_[2] - strouhalTimeStep_[0] ) << ")"
                                            "\n   f: " << ( real_t(1) / ( real_c( strouhalTimeStep_[2] - strouhalTimeStep_[0] ) * setup_.dt ) ) << " Hz ("
                                                       << ( real_t(1) / real_c( strouhalTimeStep_[2] - strouhalTimeStep_[0] ) ) << ")"
                                            "\n   St (\"real\" D):   " << strouhalNumberRealD_ <<
                                            "\n   St (discrete D): " << strouhalNumberDiscreteD_ );
            }
         }
      }
      
      if( logToFile_ )
      {
         std::ofstream file( filename_.c_str(), std::ios_base::app );
         file << executionCounter_ - uint_t(1) << " " << force_[0] << " " << force_[1] << " " << force_[2] <<
                                                  " " << cDRealArea << " " << cLRealArea << " " << cDDiscreteArea << " " << cLDiscreteArea <<
                                                  " " << pressureDifference_L << " " << pressureDifference <<
                                                  " " << vortexVelocity << " " << strouhalNumberRealD_ << " " << strouhalNumberDiscreteD_ << std::endl;
         file.close();
      }
   }

   // WALBERLA_MPI_WORLD_BARRIER();
}



template< typename LatticeModel_T >
void Evaluation< LatticeModel_T >::operator()( const uint_t level, const uint_t executionCount )
{
   // Supposed to be used as a post collide function on every level

   if( checkFrequency_ == uint_t(0) || executionCounter_ % checkFrequency_ != 0 || level > uint_t(0) || executionCount != uint_t(0) )
      return;

   if( !initialized_ )
      refresh();

   force_[0] = real_t(0);
   force_[1] = real_t(0);
   force_[2] = real_t(0);

   forceSample_[0].clear();
   forceSample_[1].clear();
   
   forceEvaluationExecutionCount_ = executionCounter_;
}



template< typename LatticeModel_T >
void Evaluation< LatticeModel_T >::operator()( IBlock * block, const uint_t level, const uint_t executionCount )
{
   // Supposed to be used as a post boundary handling function on every level

   if( checkFrequency_ == uint_t(0) || executionCounter_ % checkFrequency_ != 0 || math::uintPow2(level) != ( executionCount + uint_t(1) ) )
      return;

   if( directions_.find( block ) != directions_.end() )
   {
      const PdfField_T * const pdfField = block->template getData< PdfField_T >( pdfFieldId_ );

      const auto & directions = directions_[ block ];
      for( auto pair = directions.begin(); pair != directions.end(); ++pair )
      {
         const Cell cell( pair->first );
         const stencil::Direction direction( pair->second );

         const real_t scaleFactor = real_t(1) / real_c( uint_t(1) << ( (Is2D< LatticeModel_T >::value ? uint_t(1) : uint_t(2)) * level ) );

         const real_t boundaryValue = pdfField->get( cell.x() + stencil::cx[direction],
                                                     cell.y() + stencil::cy[direction],
                                                     cell.z() + stencil::cz[direction], Stencil_T::idx[ stencil::inverseDir[direction] ] );

         const real_t fluidValue = pdfField->get( cell.x(), cell.y(), cell.z(), Stencil_T::idx[ direction ] );

         const real_t f = scaleFactor * ( boundaryValue + fluidValue );

         force_[0] += real_c( stencil::cx[ direction ] ) * f;
         force_[1] += real_c( stencil::cy[ direction ] ) * f;
         force_[2] += real_c( stencil::cz[ direction ] ) * f;

         if( setup_.evaluateForceComponents )
         {
            forceSample_[0].insert( real_c( stencil::cx[ direction ] ) * f );
            forceSample_[1].insert( real_c( stencil::cy[ direction ] ) * f );
         }
      }
   }
}



template< typename LatticeModel_T >
void Evaluation< LatticeModel_T >::prepareResultsForSQL()
{
   real_t cDRealArea( real_t(0) );
   real_t cLRealArea( real_t(0) );
   real_t cDDiscreteArea( real_t(0) );
   real_t cLDiscreteArea( real_t(0) );
   
   real_t pressureDifference_L( real_t(0) );
   real_t pressureDifference( real_t(0) );
      
   evaluate( cDRealArea, cLRealArea, cDDiscreteArea, cLDiscreteArea, pressureDifference_L, pressureDifference );

   WALBERLA_ROOT_SECTION()
   {   
      sqlResults_[ "forceX_L" ] = double_c( force_[0] );
      sqlResults_[ "forceY_L" ] = double_c( force_[1] );
      sqlResults_[ "forceZ_L" ] = double_c( force_[2] );
      
      if( setup_.evaluateForceComponents )
      {
         sqlResults_[ "forceXMin_L" ]    = double_c( forceSample_[0].min() );
         sqlResults_[ "forceXMax_L" ]    = double_c( forceSample_[0].max() );
         sqlResults_[ "forceXAvg_L" ]    = double_c( forceSample_[0].mean() );
         sqlResults_[ "forceXMedian_L" ] = double_c( forceSample_[0].median() );
         sqlResults_[ "forceXStdDev_L" ] = double_c( forceSample_[0].stdDeviation() );
      
         sqlResults_[ "forceYMin_L" ]    = double_c( forceSample_[1].min() );
         sqlResults_[ "forceYMax_L" ]    = double_c( forceSample_[1].max() );
         sqlResults_[ "forceYAvg_L" ]    = double_c( forceSample_[1].mean() );
         sqlResults_[ "forceYMedian_L" ] = double_c( forceSample_[1].median() );
         sqlResults_[ "forceYStdDev_L" ] = double_c( forceSample_[1].stdDeviation() );
      }
      else
      {
         sqlResults_[ "forceXMin_L" ]    = 0.0;
         sqlResults_[ "forceXMax_L" ]    = 0.0;
         sqlResults_[ "forceXAvg_L" ]    = 0.0;
         sqlResults_[ "forceXMedian_L" ] = 0.0;
         sqlResults_[ "forceXStdDev_L" ] = 0.0;
      
         sqlResults_[ "forceYMin_L" ]    = 0.0;
         sqlResults_[ "forceYMax_L" ]    = 0.0;
         sqlResults_[ "forceYAvg_L" ]    = 0.0;
         sqlResults_[ "forceYMedian_L" ] = 0.0;
         sqlResults_[ "forceYStdDev_L" ] = 0.0;
      }
      
      sqlResults_[ "cDRealArea" ]     = double_c( cDRealArea );
      sqlResults_[ "cLRealArea" ]     = double_c( cLRealArea );
      sqlResults_[ "cDDiscreteArea" ] = double_c( cDDiscreteArea );
      sqlResults_[ "cLDiscreteArea" ] = double_c( cLDiscreteArea );
      
      sqlResults_[ "cDRealAreaMin" ]     = double_c( coefficientExtremas_[0].first  );
      sqlResults_[ "cDRealAreaMax" ]     = double_c( coefficientExtremas_[0].second );
      sqlResults_[ "cLRealAreaMin" ]     = double_c( coefficientExtremas_[1].first  );
      sqlResults_[ "cLRealAreaMax" ]     = double_c( coefficientExtremas_[1].second );
      sqlResults_[ "cDDiscreteAreaMin" ] = double_c( coefficientExtremas_[2].first  );
      sqlResults_[ "cDDiscreteAreaMax" ] = double_c( coefficientExtremas_[2].second );
      sqlResults_[ "cLDiscreteAreaMin" ] = double_c( coefficientExtremas_[3].first  );
      sqlResults_[ "cLDiscreteAreaMax" ] = double_c( coefficientExtremas_[3].second );
      
      sqlResults_[ "pressureDifference_L" ] = double_c( pressureDifference_L );
      sqlResults_[ "pressureDifference" ]   = double_c( pressureDifference );

      sqlResults_[ "strouhalNumberRealD" ]     = double_c( strouhalNumberRealD_ );
      sqlResults_[ "strouhalNumberDiscreteD" ] = double_c( strouhalNumberDiscreteD_ );
   }
}



template< typename LatticeModel_T >
void Evaluation< LatticeModel_T >::getResultsForSQLOnRoot( std::map< std::string, double > & realProperties,
                                                           std::map< std::string, int > & integerProperties ) const
{
   WALBERLA_ROOT_SECTION()
   {
      for( auto result = sqlResults_.begin(); result != sqlResults_.end(); ++result )
         realProperties[ result->first ] = result->second;
         
      integerProperties[ "forceEvaluationTimeStep" ] = int_c( forceEvaluationExecutionCount_ );
      integerProperties[ "strouhalEvaluationTimeStep" ] = int_c( strouhalEvaluationExecutionCount_ );
   }
}



template< typename LatticeModel_T >
void Evaluation< LatticeModel_T >::check( const shared_ptr< Config > & config )
{
   real_t cDRealArea( real_t(0) );
   real_t cLRealArea( real_t(0) );
   real_t cDDiscreteArea( real_t(0) );
   real_t cLDiscreteArea( real_t(0) );
   
   real_t pressureDifference_L( real_t(0) );
   real_t pressureDifference( real_t(0) );
      
   evaluate( cDRealArea, cLRealArea, cDDiscreteArea, cLDiscreteArea, pressureDifference_L, pressureDifference );

   WALBERLA_ROOT_SECTION()
   {   
      Config::BlockHandle configBlock = config->getBlock( "SchaeferTurek" );
      
      const real_t checkCDRealAreaLowerBound = configBlock.getParameter< real_t >( "checkCDRealAreaLowerBound", -real_c(1E6) );
      const real_t checkCDRealAreaUpperBound = configBlock.getParameter< real_t >( "checkCDRealAreaUpperBound",  real_c(1E6) );
      
      const real_t checkCLRealAreaLowerBound = configBlock.getParameter< real_t >( "checkCLRealAreaLowerBound", -real_c(1E6) );
      const real_t checkCLRealAreaUpperBound = configBlock.getParameter< real_t >( "checkCLRealAreaUpperBound",  real_c(1E6) );

      const real_t checkCDDiscreteAreaLowerBound = configBlock.getParameter< real_t >( "checkCDDiscreteAreaLowerBound", -real_c(1E6) );
      const real_t checkCDDiscreteAreaUpperBound = configBlock.getParameter< real_t >( "checkCDDiscreteAreaUpperBound",  real_c(1E6) );
      
      const real_t checkCLDiscreteAreaLowerBound = configBlock.getParameter< real_t >( "checkCLDiscreteAreaLowerBound", -real_c(1E6) );
      const real_t checkCLDiscreteAreaUpperBound = configBlock.getParameter< real_t >( "checkCLDiscreteAreaUpperBound",  real_c(1E6) );
      
      const real_t checkPressureDiffLowerBound = configBlock.getParameter< real_t >( "checkPressureDiffLowerBound", -real_c(1E6) );
      const real_t checkPressureDiffUpperBound = configBlock.getParameter< real_t >( "checkPressureDiffUpperBound",  real_c(1E6) );
      
      const real_t checkStrouhalNbrRealDLowerBound = configBlock.getParameter< real_t >( "checkStrouhalNbrRealDLowerBound", -real_c(1E6) );
      const real_t checkStrouhalNbrRealDUpperBound = configBlock.getParameter< real_t >( "checkStrouhalNbrRealDUpperBound",  real_c(1E6) );
      
      const real_t checkStrouhalNbrDiscreteDLowerBound = configBlock.getParameter< real_t >( "checkStrouhalNbrDiscreteDLowerBound", -real_c(1E6) );
      const real_t checkStrouhalNbrDiscreteDUpperBound = configBlock.getParameter< real_t >( "checkStrouhalNbrDiscreteDUpperBound",  real_c(1E6) );

      WALBERLA_CHECK_GREATER( cDRealArea, checkCDRealAreaLowerBound );
      WALBERLA_CHECK_LESS   ( cDRealArea, checkCDRealAreaUpperBound );
      
      WALBERLA_CHECK_GREATER( cLRealArea, checkCLRealAreaLowerBound );
      WALBERLA_CHECK_LESS   ( cLRealArea, checkCLRealAreaUpperBound );
      
      WALBERLA_CHECK_GREATER( cDDiscreteArea, checkCDDiscreteAreaLowerBound );
      WALBERLA_CHECK_LESS   ( cDDiscreteArea, checkCDDiscreteAreaUpperBound );
      
      WALBERLA_CHECK_GREATER( cLDiscreteArea, checkCLDiscreteAreaLowerBound );
      WALBERLA_CHECK_LESS   ( cLDiscreteArea, checkCLDiscreteAreaUpperBound );
      
      if( setup_.evaluatePressure )
      {
         WALBERLA_CHECK_GREATER( pressureDifference, checkPressureDiffLowerBound );
         WALBERLA_CHECK_LESS   ( pressureDifference, checkPressureDiffUpperBound );
      }
      
      if( setup_.evaluateStrouhal )
      {
         WALBERLA_CHECK_GREATER( strouhalNumberRealD_, checkStrouhalNbrRealDLowerBound );
         WALBERLA_CHECK_LESS   ( strouhalNumberRealD_, checkStrouhalNbrRealDUpperBound );
      
         WALBERLA_CHECK_GREATER( strouhalNumberDiscreteD_, checkStrouhalNbrDiscreteDLowerBound );
         WALBERLA_CHECK_LESS   ( strouhalNumberDiscreteD_, checkStrouhalNbrDiscreteDUpperBound );
      }
   }
}



template< typename LatticeModel_T >
void Evaluation< LatticeModel_T >::refresh()
{
   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks );

   if( Is2D< LatticeModel_T >::value )
   {
      const real_t z = blocks->getDomain().zMin() + real_c(0.5) * blocks->dz( blocks->getNumberOfLevels() - uint_t(1) );
      setup_.pAlpha[2] = z;
      setup_.pOmega[2] = z;
      setup_.pStrouhal[2] = z;
   }

   // Calculate obstacle surface areas required for evaluating drag and lift force

   real_t yMin( setup_.H );
   real_t yMax( real_t(0) );
   
   real_t AD( real_t(0) );
   real_t AL( real_t(0) );
   
   directions_.clear();

   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
   {
      const FlagField_T * const flagField = block->template getData< FlagField_T >( flagFieldId_ );

      auto fluid    = flagField->getFlag( fluid_ );
      auto obstacle = flagField->getFlag( obstacle_ );

      const uint_t level = blocks->getLevel(*block);
      const real_t area = real_t(1) / real_c( uint_t(1) << ( (Is2D< LatticeModel_T >::value ? uint_t(1) : uint_t(2)) * level ) );

      auto xyzSize = flagField->xyzSize();

#ifndef NDEBUG
      if( Is2D< LatticeModel_T >::value )
      {
         WALBERLA_ASSERT( blocks->atDomainZMinBorder( *block ) );
         WALBERLA_ASSERT_EQUAL( xyzSize.zMin(), xyzSize.zMax() );
      }
#endif

      for( auto z = xyzSize.zMin(); z <= xyzSize.zMax(); ++z ) {
         for( auto y = xyzSize.yMin(); y <= xyzSize.yMax(); ++y ) {
            for( auto x = xyzSize.xMin(); x <= xyzSize.xMax(); ++x )
            {
               if( flagField->isFlagSet( x, y, z, fluid ) )
               {
                  for( auto it = Stencil_T::beginNoCenter(); it != Stencil_T::end(); ++it )
                  {
                     auto nx = x + cell_idx_c( it.cx() );
                     auto ny = y + cell_idx_c( it.cy() );
                     auto nz = z + cell_idx_c( it.cz() );

                     if( flagField->isFlagSet( nx, ny, nz, obstacle ) )
                     {
                        directions_[ block.get() ].push_back( std::make_pair( Cell(x,y,z), *it ) );
                        
                        const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, Cell(nx,ny,nz) );
                        
                        if( it.cx() == 1 && it.cy() == 0 && it.cz() == 0 )
                        {
                           AD += area;
                        }
                        else if( it.cx() == 0 && it.cz() == 0 )
                        {
                           if( it.cy() == -1 )
                           {
                              yMax = std::max( yMax, p[1] );
                           }
                           else if( it.cy() == 1 )
                           {
                              yMin = std::min( yMin, p[1] );
                              AL += area;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   mpi::reduceInplace( yMin, mpi::MIN );
   mpi::reduceInplace( yMax, mpi::MAX );
   
   mpi::reduceInplace( AD, mpi::SUM );
   mpi::reduceInplace( AL, mpi::SUM );

   WALBERLA_ROOT_SECTION()
   {
      const Cell yMinCell = blocks->getCell( real_t(0), yMin, real_t(0) );
      const Cell yMaxCell = blocks->getCell( real_t(0), yMax, real_t(0) );
      
      D_ = uint_c( yMaxCell[1] - yMinCell[1] ) + uint_t(1);
      
      AD_ = AD;
      AL_ = AL;
   }
   
   // Check if points alpha and omega (required for evaluating the pressure difference) are located in fluid cells

   if( setup_.evaluatePressure )
   {
      int alpha( 0 );
      int omega( 0 );

      auto block = blocks->getBlock( setup_.pAlpha );
      if( block != nullptr )
      {
         const FlagField_T * const flagField = block->template getData< FlagField_T >( flagFieldId_ );

         const auto cell = blocks->getBlockLocalCell( *block, setup_.pAlpha );
         WALBERLA_ASSERT( flagField->xyzSize().contains( cell ) );

         const auto fluid = flagField->getFlag( fluid_ );

         if( !flagField->isFlagSet( cell, fluid ) )
         {
            const auto aabb = blocks->getBlockLocalCellAABB( *block, cell );
            const Vector3< real_t > pAlpha = setup_.pAlpha;
            setup_.pAlpha[0] = aabb.xMin() - aabb.xSize() / real_t(2);
            WALBERLA_LOG_WARNING( "Cell for evaluating pressure difference at point alpha " << pAlpha << " is not a fluid cell!"
                                  "\nChanging point alpha to " << setup_.pAlpha << " ..." );
         }
         else
         {
            alpha = 1;
         }
      }

      block = blocks->getBlock( setup_.pOmega );
      if( block != nullptr )
      {
         const FlagField_T * const flagField = block->template getData< FlagField_T >( flagFieldId_ );

         const auto cell = blocks->getBlockLocalCell( *block, setup_.pOmega );
         WALBERLA_ASSERT( flagField->xyzSize().contains( cell ) );

         const auto fluid = flagField->getFlag( fluid_ );

         if( !flagField->isFlagSet( cell, fluid ) )
         {
            const auto aabb = blocks->getBlockLocalCellAABB( *block, cell );
            const Vector3< real_t > pOmega = setup_.pOmega;
            setup_.pOmega[0] = aabb.xMax() + aabb.xSize() / real_t(2);
            WALBERLA_LOG_WARNING( "Cell for evaluating pressure difference at point omega " << pOmega << " is not a fluid cell!"
                                  "\nChanging point omega to " << setup_.pOmega << " ..." );
         }
         else
         {
            omega = 1;
         }
      }

      mpi::allReduceInplace( alpha, mpi::SUM );
      mpi::allReduceInplace( omega, mpi::SUM );

      if( alpha == 0 )
      {
         block = blocks->getBlock( setup_.pAlpha );
         if( block != nullptr )
         {
            const FlagField_T * const flagField = block->template getData< FlagField_T >( flagFieldId_ );

            const auto cell = blocks->getBlockLocalCell( *block, setup_.pAlpha );
            WALBERLA_ASSERT( flagField->xyzSize().contains( cell ) );

            const auto fluid = flagField->getFlag( fluid_ );

            if( !flagField->isFlagSet( cell, fluid ) )
               WALBERLA_ABORT( "Cell for evaluating pressure difference at point alpha " << setup_.pAlpha << " is still not a fluid cell!" );

            alpha = 1;
         }
         mpi::reduceInplace( alpha, mpi::SUM );
      }

      if( omega == 0 )
      {
         block = blocks->getBlock( setup_.pOmega );
         if( block != nullptr )
         {
            const FlagField_T * const flagField = block->template getData< FlagField_T >( flagFieldId_ );

            const auto cell = blocks->getBlockLocalCell( *block, setup_.pOmega );
            WALBERLA_ASSERT( flagField->xyzSize().contains( cell ) );

            const auto fluid = flagField->getFlag( fluid_ );

            if( !flagField->isFlagSet( cell, fluid ) )
               WALBERLA_ABORT( "Cell for evaluating pressure difference at point omega " << setup_.pOmega << " is still not a fluid cell!" );

            omega = 1;
         }
         mpi::reduceInplace( omega, mpi::SUM );
      }

      WALBERLA_ROOT_SECTION()
      {
         if( alpha == 0 )
            WALBERLA_ABORT( "Point alpha " << setup_.pAlpha << " (required for evaluating the pressure difference) is not located inside the fluid domain!" );
         WALBERLA_ASSERT_EQUAL( alpha, 1 );

         if( omega == 0 )
            WALBERLA_ABORT( "Point omega " << setup_.pOmega << " (required for evaluating the pressure difference) is not located inside the fluid domain!" );
         WALBERLA_ASSERT_EQUAL( omega, 1 );
      }
   }

   // Check if point for evaluating the Strouhal number is located inside of a fluid cell

   if( setup_.evaluateStrouhal )
   {
      int strouhal( 0 );

      auto block = blocks->getBlock( setup_.pStrouhal );
      if( block != nullptr )
      {
         const FlagField_T * const flagField = block->template getData< FlagField_T >( flagFieldId_ );

         const auto cell = blocks->getBlockLocalCell( *block, setup_.pStrouhal );
         WALBERLA_ASSERT( flagField->xyzSize().contains( cell ) );

         const auto fluid = flagField->getFlag( fluid_ );

         if( !flagField->isFlagSet( cell, fluid ) )
            WALBERLA_ABORT( "Cell for evaluating the Strouhal number at point " << setup_.pStrouhal << " is not a fluid cell!" );

         strouhal = 1;
      }

      mpi::reduceInplace( strouhal, mpi::SUM );

      WALBERLA_ROOT_SECTION()
      {
         if( strouhal == 0 )
            WALBERLA_ABORT( "Point " << setup_.pStrouhal << " (required for evaluating the Strouhal number) is not located inside the fluid domain!" );
         WALBERLA_ASSERT_EQUAL( strouhal, 1 );
      }
   }

   initialized_ = true;
}



template< typename LatticeModel_T >
void Evaluation< LatticeModel_T >::evaluate( real_t & cDRealArea, real_t & cLRealArea, real_t & cDDiscreteArea, real_t & cLDiscreteArea,
                                             real_t & pressureDifference_L, real_t & pressureDifference )
{
   if( !initialized_ )
      refresh();

   // force on obstacle

   mpi::reduceInplace( force_, mpi::SUM );

   if( setup_.evaluateForceComponents )
   {
      forceSample_[0].mpiGatherRoot();
      forceSample_[1].mpiGatherRoot();
   }

   // pressure difference

   real_t pAlpha( real_t(0) );
   real_t pOmega( real_t(0) );

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks );

   if( setup_.evaluatePressure )
   {
      auto block = blocks->getBlock( setup_.pAlpha );
      if( block != nullptr )
      {
         const PdfField_T * const pdfField = block->template getData< PdfField_T >( pdfFieldId_ );
         const auto cell = blocks->getBlockLocalCell( *block, setup_.pAlpha );
         WALBERLA_ASSERT( pdfField->xyzSize().contains( cell ) );
         pAlpha += pdfField->getDensity( cell ) / real_c(3);
      }

      block = blocks->getBlock( setup_.pOmega );
      if( block != nullptr )
      {
         const PdfField_T * const pdfField = block->template getData< PdfField_T >( pdfFieldId_ );
         const auto cell = blocks->getBlockLocalCell( *block, setup_.pOmega );
         WALBERLA_ASSERT( pdfField->xyzSize().contains( cell ) );
         pOmega += pdfField->getDensity( cell ) / real_c(3);
      }

      mpi::reduceInplace( pAlpha, mpi::SUM );
      mpi::reduceInplace( pOmega, mpi::SUM );
   }

   WALBERLA_ROOT_SECTION()
   {
      const real_t uMean = Is2D< LatticeModel_T >::value ? ( real_c(2) * setup_.inflowVelocity_L / real_c(3) ) :
                                                           ( real_c(4) * setup_.inflowVelocity_L / real_c(9) );

      const real_t D = real_t(2) * setup_.cylinderRadius / setup_.dx;
      const real_t H = setup_.H / setup_.dx;
      
      cDRealArea = ( real_t(2) * force_[0] ) / ( uMean * uMean * D * (Is2D< LatticeModel_T >::value ? real_t(1) : H) );
      cLRealArea = ( real_t(2) * force_[1] ) / ( uMean * uMean * D * (Is2D< LatticeModel_T >::value ? real_t(1) : H) );
      
      cDDiscreteArea = ( real_t(2) * force_[0] ) / ( uMean * uMean * AD_ );
      cLDiscreteArea = ( real_t(2) * force_[1] ) / ( uMean * uMean * AL_ );
      
      pressureDifference_L = pAlpha - pOmega;
      pressureDifference   = ( pressureDifference_L * setup_.rho * setup_.dx * setup_.dx ) / ( setup_.dt * setup_.dt );
   }
}






template< typename LatticeModel_T >
class EvaluationRefresh
{
public:

   EvaluationRefresh( const weak_ptr< Evaluation<LatticeModel_T> > & evaluation ) : evaluation_( evaluation ) {}

   void operator()()
   {
      auto evaluation = evaluation_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( evaluation );
      evaluation->refresh();
   }

   void operator()( BlockForest &, const PhantomBlockForest & )
   {
      auto evaluation = evaluation_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( evaluation );
      evaluation->refresh();
   }

private:

   weak_ptr< Evaluation<LatticeModel_T> > evaluation_;
};






/////////////////
// SNAPSHOTING //
/////////////////

class SnapshotSimulator
{
public:

   SnapshotSimulator( const weak_ptr<blockforest::StructuredBlockForest> & blocks,
                      const uint_t failAt, const uint_t failRangeBegin, const uint_t failRangeEnd, const bool failRebalance ) :
      blocks_( blocks ), executionCounter_( 0 ),
      failAt_( failAt ), failRangeBegin_( failRangeBegin ), failRangeEnd_( failRangeEnd ), failRebalance_( failRebalance )
   {}

   void operator()()
   {
      if( executionCounter_ == failAt_ )
      {
         auto blocks = blocks_.lock();
         WALBERLA_CHECK_NOT_NULLPTR( blocks );

         const uint_t process = uint_c( mpi::MPIManager::instance()->rank() );

         if( process >= failRangeBegin_ && process <= failRangeEnd_ )
         {
            for( uint_t i = uint_t(0); i != blocks->numberOfBlockDataItems(); ++i )
               blocks->clearBlockData( BlockDataID(i) );
         }

         blocks->getBlockForest().restoreSnapshot( *this, failRebalance_ );
      }
      ++executionCounter_;
   }

   void operator()( std::vector<uint_t> & sendTo, std::vector<uint_t> & recvFrom )
   {
      const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );
      const uint_t process = uint_c( mpi::MPIManager::instance()->rank() );

      if( processes > uint_t(1) )
      {
         auto blocks = blocks_.lock();
         WALBERLA_CHECK_NOT_NULLPTR( blocks );
         WALBERLA_CHECK_EQUAL( blocks->getProcess(), process );

         const uint_t shift = processes / uint_t(2);

         sendTo.clear();
         recvFrom.clear();

         sendTo.push_back( ( process + shift ) % processes );
         recvFrom.push_back( ( shift > process ) ? ( processes - ( shift - process ) ) : ( process - shift ) );
      }
   }

   uint_t operator()( const uint_t previousProcess )
   {
      if( previousProcess >= failRangeBegin_ && previousProcess <= failRangeEnd_ )
      {
         const uint_t processes = uint_c( mpi::MPIManager::instance()->numProcesses() );
         const uint_t shift = processes / uint_t(2);
         return ( previousProcess + shift ) % processes;
      }

      return previousProcess;
   }

private:

   weak_ptr<blockforest::StructuredBlockForest> blocks_;

   uint_t executionCounter_;
   uint_t failAt_;

   uint_t failRangeBegin_;
   uint_t failRangeEnd_;

   bool failRebalance_;

};






////////////////////
// THE SIMULATION //
////////////////////

template< typename LatticeModel_T, typename Sweep_T >
void addRefinementTimeStep( SweepTimeloop & timeloop, shared_ptr< blockforest::StructuredBlockForest > & blocks,
                            const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId,
                            const shared_ptr<WcTimingPool> & timingPool, const shared_ptr<WcTimingPool> & levelwiseTimingPool,
                            const bool syncComm, const bool fullComm, const bool linearExplosion,
                            shared_ptr< Sweep_T > & sweep, const std::string & info,
                            const shared_ptr< Evaluation< LatticeModel_T > > & evaluation,
                            const shared_ptr< lbm::TimeTracker > & timeTracker )
{
   using BH_T = typename MyBoundaryHandling< LatticeModel_T >::BoundaryHandling_T;

   auto ts = lbm::refinement::makeTimeStep< LatticeModel_T, BH_T >( blocks, sweep, pdfFieldId, boundaryHandlingId, None, Empty );
   ts->asynchronousCommunication( !syncComm );
   ts->optimizeCommunication( !fullComm );
   ts->performLinearExplosion( linearExplosion );
   ts->enableTiming( timingPool, levelwiseTimingPool );

   ts->template addPostCollideVoidFunction< Evaluation< LatticeModel_T > >( evaluation, "obstacle force reset" );
   ts->template addPostBoundaryHandlingBlockFunction< Evaluation< LatticeModel_T > >( evaluation, "obstacle force calculation" );
   ts->template addPostBoundaryHandlingVoidFunction< lbm::TimeTracker >( timeTracker, "time tracker" );

   timeloop.addFuncBeforeTimeStep( makeSharedFunctor(ts), info );
}

template< typename LatticeModel_T, class Enable = void >
struct AddRefinementTimeStep
{
   static void add( SweepTimeloop & timeloop, shared_ptr< blockforest::StructuredBlockForest > & blocks,
                    const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const BlockDataID & boundaryHandlingId,
                    const shared_ptr<WcTimingPool> & timingPool, const shared_ptr<WcTimingPool> & levelwiseTimingPool,
                    const bool split, const bool pure, const bool syncComm, const bool fullComm, const bool linearExplosion,
                    const shared_ptr< Evaluation< LatticeModel_T > > & evaluation,
                    const shared_ptr< lbm::TimeTracker > & timeTracker )
   {
      if( split )
      {
         if( pure )
         {
            using Sweep_T = lbm::SplitPureSweep< LatticeModel_T >;
            auto mySweep = make_shared< Sweep_T >( pdfFieldId );

            addRefinementTimeStep< LatticeModel_T, Sweep_T >( timeloop, blocks, pdfFieldId, boundaryHandlingId, timingPool, levelwiseTimingPool,
                                                              syncComm, fullComm, linearExplosion, mySweep, "LBM refinement time step (split pure LB sweep)",
                                                              evaluation, timeTracker );
         }
         else
         {
            using Sweep_T = lbm::SplitSweep<LatticeModel_T, FlagField_T>;
            auto mySweep = make_shared< Sweep_T >( pdfFieldId, flagFieldId, Fluid_Flag );

            addRefinementTimeStep< LatticeModel_T, Sweep_T >( timeloop, blocks, pdfFieldId, boundaryHandlingId, timingPool, levelwiseTimingPool,
                                                              syncComm, fullComm, linearExplosion, mySweep, "LBM refinement time step (split LB sweep)",
                                                              evaluation, timeTracker );
         }
      }
      else
      {                                         
         auto mySweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag );

         addRefinementTimeStep< LatticeModel_T >( timeloop, blocks, pdfFieldId, boundaryHandlingId, timingPool, levelwiseTimingPool,
                                                  syncComm, fullComm, linearExplosion, mySweep, "LBM refinement time step (cell-wise LB sweep)",
                                                  evaluation, timeTracker );
      }
   }
};

template< typename LatticeModel_T  >
struct AddRefinementTimeStep< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, lbm::collision_model::MRT_tag >::value ||
                                                                       std::is_same< typename LatticeModel_T::Stencil, stencil::D2Q9 >::value ||
                                                                       std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q15 >::value ||
                                                                       std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value
                                                                       >::type >
{
   static void add( SweepTimeloop & timeloop, shared_ptr< blockforest::StructuredBlockForest > & blocks,
                    const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const BlockDataID & boundaryHandlingId,
                    const shared_ptr<WcTimingPool> & timingPool, const shared_ptr<WcTimingPool> & levelwiseTimingPool,
                    const bool /*split*/, const bool /*pure*/, const bool syncComm, const bool fullComm, const bool linearExplosion,
                    const shared_ptr< Evaluation< LatticeModel_T > > & evaluation,
                    const shared_ptr< lbm::TimeTracker > & timeTracker )
   {
      auto mySweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag );

      addRefinementTimeStep< LatticeModel_T >( timeloop, blocks, pdfFieldId, boundaryHandlingId, timingPool, levelwiseTimingPool,
                                               syncComm, fullComm, linearExplosion, mySweep, "LBM refinement time step (cell-wise LB sweep)",
                                               evaluation, timeTracker );
   }
};



template< typename LatticeModel_T >
void run( const shared_ptr< Config > & config, const LatticeModel_T & latticeModel,
          const bool split, const bool pure, const bool fzyx, const bool syncComm, const bool fullComm, const bool linearExplosion,
          const blockforest::RefinementSelectionFunctions & refinementSelectionFunctions, Setup & setup,
          const memory_t memoryPerCell, const memory_t processMemoryLimit )
{
   Config::BlockHandle configBlock = config->getBlock( "SchaeferTurek" );
   
   setup.viscosity_L = latticeModel.collisionModel().viscosity( uint_t(0) );
   setup.dt = ( setup.viscosity_L * setup.dx *setup.dx ) / setup.viscosity; // [s]
   setup.inflowVelocity_L = setup.inflowVelocity * setup.dt / setup.dx;
   setup.raisingTime_L = setup.raisingTime / setup.dt;
   setup.sinPeriod_L = setup.sinPeriod / setup.dt;

   // creating the block structure

   auto blocks = createStructuredBlockForest( configBlock, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );

   // add pdf field to blocks

   const real_t initVelocity = ( configBlock.getParameter< bool >( "initWithVelocity", false ) ) ?
            ( Is2D< LatticeModel_T >::value ? ( real_t(2) * setup.inflowVelocity_L / real_c(3) ) : ( real_t(4) * setup.inflowVelocity_L / real_c(9) ) ) : real_t(0);

   BlockDataID pdfFieldId = fzyx ? lbm::addPdfFieldToStorage( blocks, "pdf field (fzyx)", latticeModel,
                                                              Vector3< real_t >( initVelocity, real_c(0), real_c(0) ), real_t(1),
                                                              FieldGhostLayers, field::fzyx, None, Empty ) :
                                   lbm::addPdfFieldToStorage( blocks, "pdf field (zyxf)", latticeModel,
                                                              Vector3< real_t >( initVelocity, real_c(0), real_c(0) ), real_t(1),
                                                              FieldGhostLayers, field::zyxf, None, Empty );

   // add density adaptor

   using DensityAdaptor_T = typename lbm::Adaptor< LatticeModel_T >::Density;
   BlockDataID densityAdaptorId = field::addFieldAdaptor< DensityAdaptor_T >( blocks, pdfFieldId, "density adaptor", None, Empty );
   
   // add velocity field + initialize velocity field writer (only used for simulations with an adaptive block structure)

   using VelocityField_T = field::GhostLayerField<Vector3<real_t>, 1>;
   BlockDataID velocityFieldId = field::addToStorage< VelocityField_T >( blocks, "velocity", Vector3<real_t>(0), field::zyxf, FieldGhostLayers, true, None, Empty );

   using VelocityFieldWriter_T = lbm::VelocityFieldWriter<typename Types<LatticeModel_T>::PdfField_T, VelocityField_T>;
   BlockSweepWrapper< VelocityFieldWriter_T > velocityFieldWriter( blocks, VelocityFieldWriter_T( pdfFieldId, velocityFieldId ), None, Empty );
   velocityFieldWriter();

   // add flag field to blocks

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field", FieldGhostLayers, true, None, Empty );

   field::FlagFieldEvaluationFilter<FlagField_T> flagFieldFilter( flagFieldId, Fluid_Flag );

   // add LB boundary handling to blocks

   shared_ptr< lbm::TimeTracker > timeTracker = make_shared< lbm::TimeTracker >();

   BlockDataID boundaryHandlingId = blocks->addBlockData( make_shared< MyBoundaryHandling< LatticeModel_T > >( flagFieldId, pdfFieldId, blocks, setup, timeTracker ),
                                                          "boundary handling", None, Empty );

   const int obstacleBoundary = configBlock.getParameter< int >( "obstacleBoundary", 0 );
   const int outletType       = configBlock.getParameter< int >( "outletType", 0 );

   BoundarySetter<LatticeModel_T> boundarySetter( blocks, boundaryHandlingId, setup, obstacleBoundary, outletType, None, Empty );
   boundarySetter();

   // add 'bool' field to every block (required for LB post processing when blocks split/merge in order to keep boundaries consistent)

   BlockDataID markerDataId = blocks->addBlockData( make_shared< lbm::MarkerData< LatticeModel_T, field::FlagFieldEvaluationFilter<FlagField_T> > >( pdfFieldId, flagFieldFilter ),
                                                    "LBM marker data (for dynamic refinement post processing)", None, Empty );

   // creating the time loop

   uint_t outerTimeSteps = configBlock.getParameter< uint_t >( "outerTimeSteps", uint_c(0) );
   uint_t innerTimeSteps = configBlock.getParameter< uint_t >( "innerTimeSteps", uint_c(0) );
   if( configBlock.isDefined( "minSimulationTime" ) )
   {
      const uint_t timeSteps = uint_c( configBlock.getParameter< real_t >( "minSimulationTime" ) / setup.dt + real_c(0.5) );
      if( innerTimeSteps == uint_t(0) )
      {
         outerTimeSteps = uint_t(1);
         innerTimeSteps = timeSteps;
      }
      else
      {
         outerTimeSteps = timeSteps / innerTimeSteps;
         if( timeSteps % innerTimeSteps != uint_t(0) )
            ++outerTimeSteps;
      }
   }

   SweepTimeloop timeloop( blocks->getBlockStorage(), ( outerTimeSteps * innerTimeSteps ) + uint_t(1) );

   // VTK

   blockforest::communication::NonUniformBufferedScheme< typename lbm::NeighborsStencil<LatticeModel_T>::type > pdfGhostLayerSync( blocks, None, Empty );
   pdfGhostLayerSync.addPackInfo( make_shared< lbm::refinement::PdfFieldSyncPackInfo< LatticeModel_T > >( pdfFieldId ) );

   MyVTKOutput< LatticeModel_T > myVTKOutput( pdfFieldId, flagFieldId, pdfGhostLayerSync );

   std::map< std::string, vtk::SelectableOutputFunction > vtkOutputFunctions;
   vtk::initializeVTKOutput( vtkOutputFunctions, myVTKOutput, blocks, config );   
   
   const bool vtkBeforeTimeStep = configBlock.getParameter< bool >( "vtkBeforeTimeStep", true );

   if( vtkBeforeTimeStep )
   {
      for( auto output = vtkOutputFunctions.begin(); output != vtkOutputFunctions.end(); ++output )
         timeloop.addFuncBeforeTimeStep( output->second.outputFunction, std::string("VTK: ") + output->first,
                                         output->second.requiredGlobalStates, output->second.incompatibleGlobalStates );
   }

   // evaluation

   const uint_t      evaluationCheckFrequency = configBlock.getParameter< uint_t >( "evaluationCheckFrequency", uint_t(100) );
   const bool        evaluationLogToStream    = configBlock.getParameter< bool >( "evaluationLogToStream", true );
   const bool        evaluationLogToFile      = configBlock.getParameter< bool >( "evaluationLogToFile", true );
   const std::string evaluationFilename       = configBlock.getParameter< std::string >( "evaluationFilename", std::string("SchaeferTurek.txt") );

   shared_ptr< Evaluation< LatticeModel_T > > evaluation( new Evaluation< LatticeModel_T >( blocks, evaluationCheckFrequency, pdfFieldId, flagFieldId, Fluid_Flag,
                                                                                            ( obstacleBoundary == 1 ) ? Curved_Flag : Obstacle_Flag,
                                                                                            setup, evaluationLogToStream, evaluationLogToFile, evaluationFilename,
                                                                                            None, Empty ) );
   // block structure refresh (rebalance + redistribute blocks -> dynamic load balancing)

   auto & blockforest = blocks->getBlockForest();
   
   const uint_t blockforestRefreshFrequency = configBlock.getParameter< uint_t >( "blockforestRefreshFrequency", uint_t(0) );
   const bool dynamicBlockStructure = ( blockforestRefreshFrequency != uint_t(0) );
   
   std::string adaptiveRefinementLog;

   if( dynamicBlockStructure )
   {
      // block level determination (= adaptive block structure)

      blockforest::MinTargetLevelDeterminationFunctions minTargetLevelDeterminationFunctions;

      const int refinementType = configBlock.getParameter< int >( "refinementType", 1 );
      if( refinementType == 0 )
      {
         const real_t lowerLimit = configBlock.getParameter< real_t >( "curlLowerLimit" );
         const real_t upperLimit = configBlock.getParameter< real_t >( "curlUpperLimit" );

         lbm::refinement::VorticityBasedLevelDetermination< field::FlagFieldEvaluationFilter<FlagField_T>, Is2D<LatticeModel_T>::value > vorticityRefinement(
            velocityFieldId, flagFieldFilter, upperLimit, lowerLimit, configBlock.getParameter< uint_t >( "maxLevel", uint_t(0) ) );

         minTargetLevelDeterminationFunctions.add( vorticityRefinement );
         
         std::ostringstream oss;
         oss << " (vorticity : [" << lowerLimit << "," << upperLimit << "])";
         adaptiveRefinementLog = oss.str();
      }
      else
      {
         const real_t lowerLimit = configBlock.getParameter< real_t >( "gradLowerLimit" );
         const real_t upperLimit = configBlock.getParameter< real_t >( "gradUpperLimit" );
         
         field::GradientRefinement< VelocityField_T, field::FlagFieldEvaluationFilter<FlagField_T>, Is2D< LatticeModel_T >::value > gradientRefinement(
            velocityFieldId, flagFieldFilter, upperLimit, lowerLimit, configBlock.getParameter< uint_t >( "maxLevel", uint_t(0) ) );

         minTargetLevelDeterminationFunctions.add( gradientRefinement );
         
         std::ostringstream oss;
         oss << " (gradient : [" << lowerLimit << "," << upperLimit << "])";
         adaptiveRefinementLog = oss.str();
      }

      minTargetLevelDeterminationFunctions.add( std::bind( keepInflowOutflowAtTheSameLevel, std::placeholders::_1, std::placeholders::_2, 
                                                           std::placeholders::_3, std::cref(setup) ) );

      if( Is2D< LatticeModel_T >::value )
         minTargetLevelDeterminationFunctions.add( pseudo2DTargetLevelCorrection );

      blockforest.setRefreshMinTargetLevelDeterminationFunction( minTargetLevelDeterminationFunctions );

      blockforest.reevaluateMinTargetLevelsAfterForcedRefinement( true ); // required in order to guarantee the same level for all in- and outflow blocks

      blockforest.checkForEarlyOutInRefresh( configBlock.getParameter< bool >( "checkForEarlyOutInRefresh", true ) );
      blockforest.checkForLateOutInRefresh( configBlock.getParameter< bool >( "checkForLateOutInRefresh", true ) );

      // activate adaptive refinement by registering the block structure refresh at the time loop

      timeloop.addFuncBeforeTimeStep( blockforest.getRefreshFunctor( blockforestRefreshFrequency ), "block forest refresh" );
   }

   const uint_t blockforestSnapshotFrequency = configBlock.getParameter< uint_t >( "blockforestSnapshotFrequency", uint_t(0) );
   if( blockforestSnapshotFrequency != uint_t(0) )
   {
      blockforest.addCallbackFunctionAfterBlockDataIsRestored( boundarySetter );
      if( dynamicBlockStructure )
         blockforest.addCallbackFunctionAfterBlockDataIsRestored( velocityFieldWriter );
      if( evaluationCheckFrequency != uint_t(0) )
         blockforest.addCallbackFunctionAfterBlockDataIsRestored( EvaluationRefresh< LatticeModel_T >( evaluation ) );

      const uint_t failAt         = configBlock.getParameter< uint_t >( "failAt" );
      const uint_t failRangeBegin = configBlock.getParameter< uint_t >( "failRangeBegin" );
      const uint_t failRangeEnd   = configBlock.getParameter< uint_t >( "failRangeEnd" );
      const bool   failRebalance  = configBlock.getParameter<  bool  >( "failRebalance" );

      SnapshotSimulator snapshotSimulator( blocks, failAt, failRangeBegin, failRangeEnd, failRebalance );

      timeloop.addFuncBeforeTimeStep( blockforest.getSnapshotCreationFunctor( snapshotSimulator, blockforestSnapshotFrequency ), "block forest snapshot creator" );
      timeloop.addFuncBeforeTimeStep( snapshotSimulator, "block forest snapshot (failure simulator)" );
   }

   if( dynamicBlockStructure || blockforestSnapshotFrequency != uint_t(0) )
   {
      // set load balancing function

      if( configBlock.getParameter< uint_t >( "dynamicLoadBalancing", uint_t(0) ) == uint_t(0) )
      {
         const bool hilbert = configBlock.getParameter< bool >( "curveHilbert" );
         const bool allGather = configBlock.getParameter< bool >( "curveAllGather" );
         
         if( Is2D< LatticeModel_T >::value )
         {
            blockforest.setRefreshBlockStateDeterminationFunction( Pseudo2DBlockStateDetermination( blockforest, Empty ) );
            blockforest.setRefreshPhantomBlockDataAssignmentFunction( Pseudo2DPhantomWeightAssignment( Empty ) );
            blockforest.setRefreshPhantomBlockMigrationPreparationFunction(
                     blockforest::DynamicCurveBalance< Pseudo2DPhantomWeight >( hilbert, allGather ) );
         }
         else
            blockforest.setRefreshPhantomBlockMigrationPreparationFunction(
                     blockforest::DynamicCurveBalance< blockforest::NoPhantomData >( hilbert, allGather ) );
      }
      else
      {
         const uint_t maxIterations = configBlock.getParameter< uint_t >( "diffusionMaxIterations" );
         const uint_t flowIterations = configBlock.getParameter< uint_t >( "diffusionFlowIterations" );

         if( Is2D< LatticeModel_T >::value )
         {
            blockforest.setRefreshBlockStateDeterminationFunction( Pseudo2DBlockStateDetermination( blockforest, Empty ) );
            blockforest.setRefreshPhantomBlockDataAssignmentFunction( Pseudo2DPhantomWeightAssignment( Empty ) );
            blockforest.setRefreshPhantomBlockDataPackFunction( Pseudo2DPhantomWeightPackUnpack() );
            blockforest.setRefreshPhantomBlockDataUnpackFunction( Pseudo2DPhantomWeightPackUnpack() );
            blockforest.setRefreshPhantomBlockMigrationPreparationFunction(
                     blockforest::DynamicDiffusionBalance< Pseudo2DPhantomWeight >( maxIterations, flowIterations ) );
         }
         else
            blockforest.setRefreshPhantomBlockMigrationPreparationFunction(
                     blockforest::DynamicDiffusionBalance< blockforest::NoPhantomData >( maxIterations, flowIterations ) );
      }

      // add callback functions which are executed after all block data was unpacked after the dynamic load balancing

      // for blocks that have *not* migrated: store current flag field state (required for lbm::PostProcessing)
      blockforest.addRefreshCallbackFunctionAfterBlockDataIsUnpacked( lbm::MarkerFieldGenerator< LatticeModel_T, field::FlagFieldEvaluationFilter<FlagField_T> >(
               pdfFieldId, markerDataId, flagFieldFilter ) );
      // (re)set boundaries = (re)initialize flag field for every block with respect to the new block structure (the size of neighbor blocks might have changed)
      blockforest.addRefreshCallbackFunctionAfterBlockDataIsUnpacked( blockforest::BlockForest::RefreshCallbackWrappper( boundarySetter ) );
      // treat boundary-fluid cell conversions
      blockforest.addRefreshCallbackFunctionAfterBlockDataIsUnpacked( lbm::PostProcessing< LatticeModel_T, field::FlagFieldEvaluationFilter<FlagField_T> >(
               pdfFieldId, markerDataId, flagFieldFilter ) );
      // (re)set velocity field (velocity field data is not migrated!)
      blockforest.addRefreshCallbackFunctionAfterBlockDataIsUnpacked( blockforest::BlockForest::RefreshCallbackWrappper( velocityFieldWriter ) );
      // update the state of the evaluation object
      if( evaluationCheckFrequency != uint_t(0) )
         blockforest.addRefreshCallbackFunctionAfterBlockDataIsUnpacked( EvaluationRefresh< LatticeModel_T >( evaluation ) );
      // log the current state of the block structure (involves allreduce!)
      blockforest.addRefreshCallbackFunctionAfterBlockDataIsUnpacked( blockforest::logDuringRefresh );
   }

   // add 'refinement' LB time step to time loop

   shared_ptr<WcTimingPool> refinementTimeStepTiming = make_shared<WcTimingPool>();
   shared_ptr<WcTimingPool> refinementTimeStepLevelwiseTiming = make_shared<WcTimingPool>();

   AddRefinementTimeStep< LatticeModel_T >::add( timeloop, blocks, pdfFieldId, flagFieldId, boundaryHandlingId, refinementTimeStepTiming,
                                                 refinementTimeStepLevelwiseTiming, split, pure, syncComm, fullComm, linearExplosion, evaluation, timeTracker );

   // store velocity in a seperate field

   if( dynamicBlockStructure )
      timeloop.addFuncBeforeTimeStep( velocityFieldWriter, "velocity writer" );

   // evaluation

   timeloop.addFuncBeforeTimeStep( SharedFunctor< Evaluation< LatticeModel_T > >(evaluation), "evaluation" );
                                                                                       
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( lbm::makeMassEvaluation< DensityAdaptor_T, FlagField_T, Is2D< LatticeModel_T >::value >(
            configBlock, blocks, uint_t(0), densityAdaptorId, flagFieldId, Fluid_Flag, "MassEvaluation", None, Empty ) ), "mass evaluation" );

   // VTK

   if( !vtkBeforeTimeStep )
   {
      for( auto output = vtkOutputFunctions.begin(); output != vtkOutputFunctions.end(); ++output )
         timeloop.addFuncAfterTimeStep( output->second.outputFunction, std::string("VTK: ") + output->first,
                                        output->second.requiredGlobalStates, output->second.incompatibleGlobalStates );
   }

   // stability check (non-finite values in the PDF field?)

   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< lbm::PdfField< LatticeModel_T >, FlagField_T >(
                                                        configBlock, blocks, pdfFieldId, flagFieldId, Fluid_Flag, "StabilityChecker", None, Empty ) ),
                                  "LBM stability check" );

   // remaining time logger

   const double remainingTimeLoggerFrequency = configBlock.getParameter< double >( "remainingTimeLoggerFrequency", 3.0 );
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "Remaining time logger" );

   // logging right before the simulation starts

   uint_t lbBlockForestEvaluationStamp = blockforest.getModificationStamp();
   lbm::BlockForestEvaluation< FlagField_T, Is2D< LatticeModel_T >::value > lbBlockForestEvaluation( blocks, flagFieldId, Fluid_Flag, None, Empty );
   lbBlockForestEvaluation.logInfoOnRoot();

   uint_t fluidCellsEvaluationStamp = blockforest.getModificationStamp();
   field::CellCounter< FlagField_T > fluidCells( blocks, flagFieldId, Fluid_Flag, None, Empty );
   fluidCells();

   const real_t Re = Is2D< LatticeModel_T >::value ? ( ( real_c(4) * setup.inflowVelocity * setup.cylinderRadius ) / ( real_c(3) * setup.viscosity ) ) : 
                                                     ( ( real_c(8) * setup.inflowVelocity * setup.cylinderRadius ) / ( real_c(9) * setup.viscosity ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark run data:"
                              "\n- simulation parameters:"
                              "\n   + collision model:  " << CollisionModelString< LatticeModel_T >::str() <<
                              "\n   + stencil:          " << StencilString< LatticeModel_T >::str() <<
                              "\n   + compressible:     " << ( LatticeModel_T::compressible ? "yes" : "no" ) <<
                              "\n   + split kernel:     " << ( split ? "yes" : "no" ) <<
                              "\n   + pure kernel:      " << ( pure ? "yes (collision is also performed within obstacle cells)" : "no" ) <<
                              "\n   + data layout:      " << ( fzyx ? "fzyx (structure of arrays [SoA])" : "zyxf (array of structures [AoS])" ) <<
                              "\n   + communication:    " << ( fullComm ? ( syncComm ? "synchronous, full synchronization" : "asynchronous, full synchronization" ) :
                                                                          ( syncComm ? "synchronous, block neighborhood and direction-aware optimizations" : "asynchronous, block neighborhood and direction-aware optimizations" ) ) <<
                              "\n   + linear explosion: " << ( linearExplosion ? "yes" : "no" ) <<
                              "\n- simulation properties:"
                              "\n   + pseudo 2D:           " << ( Is2D< LatticeModel_T >::value ? "yes" : "no" ) <<
                              "\n   + fluid cells:         " << fluidCells.numberOfCells() << " (in total on all levels)" <<
                              "\n   + H:                   " << setup.H << " [m]" <<
                              "\n   + L:                   " << setup.L << " [m]" <<
                              "\n   + cylinder pos.(x):    " << setup.cylinderxPosition << " [m]" <<
                              "\n   + cylinder pos.(y):    " << setup.cylinderyPosition << " [m]" <<
                              "\n   + cylinder radius:     " << setup.cylinderRadius << " [m]" <<
                              "\n   + circular profile:    " << ( setup.circularCrossSection ? "yes" : "no (= box)" ) <<
                              "\n   + cylinder boundary:   " << ( ( obstacleBoundary == 1 ) ? "curved boundary condition" : "staircase approximation" ) <<
                              "\n   + kin. viscosity:      " << setup.viscosity << " [m^2/s] (" << setup.viscosity_L << " - on the coarsest grid)" <<
                              "\n   + rho:                 " << setup.rho << " [kg/m^3]" <<
                              "\n   + inflow velocity:     " << setup.inflowVelocity << " [m/s] (" << setup.inflowVelocity_L << ")" <<
                              "\n   + outlet type:         " << ( ( outletType == 0 ) ? "pressure" : ( ( outletType == 1 ) ? "2/1" : "4/3" ) ) <<
                              "\n   + Reynolds number:     " << Re <<
                              "\n   + dx (coarsest grid):  " << setup.dx << " [m]" <<
                              "\n   + dt (coarsest grid):  " << setup.dt << " [s]" <<
                              "\n   + #time steps:         " << timeloop.getNrOfTimeSteps() << " (on the coarsest grid, " << ( real_t(1) / setup.dt ) << " for 1s of real time)"
                              "\n   + simulation time:     " << ( real_c( timeloop.getNrOfTimeSteps() ) * setup.dt ) << " [s]"
                              "\n   + adaptive refinement: " << ( dynamicBlockStructure ? "yes" : "no" ) << adaptiveRefinementLog );

   // run the simulation

   uint_t performanceEvaluationStamp = blockforest.getModificationStamp();
   lbm::PerformanceEvaluation< FlagField_T > performance( blocks, flagFieldId, Fluid_Flag, None, Empty );

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
      if( dynamicBlockStructure )
      {
         reducedRefreshTiming = blocks->getBlockForest().getRefreshTiming().getReduced();
         blocks->getBlockForest().clearRefreshTiming();
         WALBERLA_LOG_RESULT_ON_ROOT( "Block structure refresh timing:\n" << *reducedRefreshTiming );
      }

      if( performanceEvaluationStamp != blockforest.getModificationStamp() )
      {
         performanceEvaluationStamp = blockforest.getModificationStamp();
         performance.refresh();
         WALBERLA_LOG_WARNING_ON_ROOT( "ATTENTION: The following performance statistics may not be entirely correct since the block structure did change during the last time measurement!" );
      }
      performance.logResultOnRoot( innerTimeSteps, time );

      if( evaluationCheckFrequency != uint_t(0) )
         evaluation->prepareResultsForSQL();

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

            std::map< std::string, int >        integerProperties;
            std::map< std::string, double >        realProperties;
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

            stringProperties[ "pseudo2D" ]             = ( Is2D< LatticeModel_T >::value ? "yes" : "no" );
            realProperties  [ "H" ]                    = double_c( setup.H );
            realProperties  [ "L" ]                    = double_c( setup.L );
            realProperties  [ "cylinderxPosition" ]    = double_c( setup.cylinderxPosition );
            realProperties  [ "cylinderyPosition" ]    = double_c( setup.cylinderyPosition );
            realProperties  [ "cylinderRadius" ]       = double_c( setup.cylinderRadius );
            stringProperties[ "circularCrossSection" ] = ( setup.circularCrossSection ? "yes" : "no" );
            stringProperties[ "cylinderBoundary" ]     = ( ( obstacleBoundary == 1 ) ? "curved" : "staircase" );
            realProperties  [ "viscosity" ]            = double_c( setup.viscosity );
            realProperties  [ "viscosity_L" ]          = double_c( setup.viscosity_L );
            realProperties  [ "rho" ]                  = double_c( setup.rho );
            realProperties  [ "inflowVelocity" ]       = double_c( setup.inflowVelocity );
            realProperties  [ "inflowVelocity_L" ]     = double_c( setup.inflowVelocity_L );
            stringProperties[ "outletType" ]           = ( outletType == 0 ) ? "pressure" : ( ( outletType == 1 ) ? "2/1" : "4/3" );
            realProperties  [ "Re" ]                   = double_c( Re );
            realProperties  [ "dx" ]                   = double_c( setup.dx );
            realProperties  [ "dt" ]                   = double_c( setup.dt );
            stringProperties[ "adaptiveRefinement" ]   = ( dynamicBlockStructure ? "yes" : "no" );
            
            integerProperties[ "simulationTimeSteps" ] = int_c( timeloop.getNrOfTimeSteps() );
            
            realProperties[ "simulationProgress" ] = double_c( ( outerRun + uint_t(1) ) * innerTimeSteps ) / double_c( outerTimeSteps * innerTimeSteps );
            if( evaluationCheckFrequency != uint_t(0) )
               evaluation->getResultsForSQLOnRoot( realProperties, integerProperties );

            auto runId = sqlite::storeRunInSqliteDB( sqlFile, integerProperties, stringProperties, realProperties );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedTimeloopTiming, "Timeloop" );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedRTSTiming, "RefinementTimeStep" );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedRTSLTiming, "RefinementTimeStepLevelwise" );
            if( dynamicBlockStructure )
               sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedRefreshTiming, "BlockForestRefresh" );
         }
      }
   }

   // Do one more step so that we see the final output
   timeloop.singleStep();

   // logging once again at the end of the simulation, identical to logging at the beginning :-)

   if( lbBlockForestEvaluationStamp != blockforest.getModificationStamp() )
      lbBlockForestEvaluation.refresh();
   lbBlockForestEvaluation.logInfoOnRoot();

   if( fluidCellsEvaluationStamp != blockforest.getModificationStamp() )
      fluidCells();

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark run data:"
                              "\n- simulation parameters:"
                              "\n   + collision model:  " << CollisionModelString< LatticeModel_T >::str() <<
                              "\n   + stencil:          " << StencilString< LatticeModel_T >::str() <<
                              "\n   + compressible:     " << ( LatticeModel_T::compressible ? "yes" : "no" ) <<
                              "\n   + split kernel:     " << ( split ? "yes" : "no" ) <<
                              "\n   + pure kernel:      " << ( pure ? "yes (collision is also performed within obstacle cells)" : "no" ) <<
                              "\n   + data layout:      " << ( fzyx ? "fzyx (structure of arrays [SoA])" : "zyxf (array of structures [AoS])" ) <<
                              "\n   + communication:    " << ( fullComm ? ( syncComm ? "synchronous, full synchronization" : "asynchronous, full synchronization" ) :
                                                                          ( syncComm ? "synchronous, block neighborhood and direction-aware optimizations" : "asynchronous, block neighborhood and direction-aware optimizations" ) ) <<
                              "\n   + linear explosion: " << ( linearExplosion ? "yes" : "no" ) <<
                              "\n- simulation properties:"
                              "\n   + pseudo 2D:           " << ( Is2D< LatticeModel_T >::value ? "yes" : "no" ) <<
                              "\n   + fluid cells:         " << fluidCells.numberOfCells() << " (in total on all levels)" <<
                              "\n   + H:                   " << setup.H << " [m]" <<
                              "\n   + L:                   " << setup.L << " [m]" <<
                              "\n   + cylinder pos.(x):    " << setup.cylinderxPosition << " [m]" <<
                              "\n   + cylinder pos.(y):    " << setup.cylinderyPosition << " [m]" <<
                              "\n   + cylinder radius:     " << setup.cylinderRadius << " [m]" <<
                              "\n   + circular profile:    " << ( setup.circularCrossSection ? "yes" : "no (= box)" ) <<
                              "\n   + cylinder boundary:   " << ( ( obstacleBoundary == 1 ) ? "curved boundary condition" : "staircase approximation" ) <<
                              "\n   + kin. viscosity:      " << setup.viscosity << " [m^2/s] (" << setup.viscosity_L << " - on the coarsest grid)" <<
                              "\n   + rho:                 " << setup.rho << " [kg/m^3]" <<
                              "\n   + inflow velocity:     " << setup.inflowVelocity << " [m/s] (" << setup.inflowVelocity_L << ")" <<
                              "\n   + outlet type:         " << ( ( outletType == 0 ) ? "pressure" : ( ( outletType == 1 ) ? "2/1" : "4/3" ) ) <<
                              "\n   + Reynolds number:     " << Re <<
                              "\n   + dx (coarsest grid):  " << setup.dx << " [m]" <<
                              "\n   + dt (coarsest grid):  " << setup.dt << " [s]" <<
                              "\n   + #time steps:         " << timeloop.getNrOfTimeSteps() << " (on the coarsest grid, " << ( real_t(1) / setup.dt ) << " for 1s of real time)"
                              "\n   + simulation time:     " << ( real_c( timeloop.getNrOfTimeSteps() ) * setup.dt ) << " [s]"
                              "\n   + adaptive refinement: " << ( dynamicBlockStructure ? "yes" : "no" ) << adaptiveRefinementLog );

   if( evaluationCheckFrequency != uint_t(0) && configBlock.getParameter< bool >( "check", false ) )
      evaluation->check( config );
}



//////////
// MAIN //
//////////

enum CM { CMSRT, CMTRT, CMMRT };
enum LM { LMD2Q9, LMD3Q15, LMD3Q19, LMD3Q27 };

int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );
   
   if( argc < 2 )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::cout << "Usage: " << argv[0] << " path-to-configuration-file [--trt | --mrt] [--d2q9 | --d3q15 | --d3q27] [--comp] [--split [--pure]] [--fzyx] [--sync-comm] [--full-comm] [--linear-exp]\n"
                      "\n"
                      "By default, SRT is selected as collision model, an asynchronous communication scheme with block neighborhood and\n"
                      "direction-aware optimizations is chosen, and an incompressible, basic D3Q19 LB kernel is executed on a PDF field with\n"
                      "layout 'zyxf' (= array of structures [AoS]).\n"
                      "\n"
                      "Optional arguments:\n"
                      " --trt:        collision model = TRT\n"
                      " --mrt:        collision model = MRT\n"
                      " --d2q9:       A D2Q9 model is used.\n"
                      " --d3q15:      A D3Q15 model is used.\n"
                      " --d3q27:      A D3Q27 model is used.\n"
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
                      "Please note: Depending on the underlying hardware and the configuration of the application (more precisely: the number of cells\n"
                      "             in each block), the best performance may be achieved with split, pure kernels combined with a structure of arrays\n"
                      "             ('fzyx') data layout!" << std::endl;
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

   Config::BlockHandle configBlock = config->getBlock( "SchaeferTurek" );

   if( !configBlock )
      WALBERLA_ABORT( "You have to specify a \"SchaeferTurek\" block in the configuration file!" );

   //////////
   // INFO //
   //////////

   WALBERLA_LOG_INFO_ON_ROOT( "//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n"
                              "//                                                                                                                      //\n"
                              "//                                               Schaefer Turek Benchmark                                               //\n"
                              "//                                                                                                                      //\n"
                              "// Reference: Schaefer, M. and Turek, S. (1996) Benchmark computations of laminar flow around a cylinder (with support  //\n"
                              "//            by F. Durst, E. Krause and R. Rannacher), in E. Hirschel (Ed.): Flow Simulation with High-Performance     //\n"
                              "//            Computers II. DFG Priority Research Program Results 1993-1995, No. 48 in Notes on Numerical Fluid         //\n"
                              "//            Mechanics, pp.547-566, Vieweg, Weisbaden.                                                                 //\n"
                              "//                                                                                                                      //\n"
                              "//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" );

   ////////////////
   // PARAMETERS //
   ////////////////

   // reading optional parameters from passed arguments

   CM   collisionModel  = CMSRT;
   LM   lm              = LMD3Q19;
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
      if( std::strcmp( argv[i], "--d2q9" )       == 0 ) lm              = LMD2Q9;
      if( std::strcmp( argv[i], "--d3q15" )      == 0 ) lm              = LMD3Q15;
      if( std::strcmp( argv[i], "--d3q27" )      == 0 ) lm              = LMD3Q27;
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
      WALBERLA_LOG_WARNING_ON_ROOT( "You called the simulation with \"--pure\" but without \"--split\".\n"
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

   if( collisionModel == CMMRT && ( lm == LMD3Q27 || lm == LMD3Q15 || lm == LMD2Q9 ) )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "D2Q9, D3Q15, and D3Q27 are not available for MRT!\n"
                                    "Setting lattice model to D3Q19 ..." );
      lm = LMD3Q19;
   }

   if( ( lm == LMD3Q27 || lm == LMD3Q15 || lm == LMD2Q9 ) && split )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "For D2Q9, D3Q15, and D3Q27 there are no \"split\" kernels!\n"
                                    "Setting \"split\" and \"pure\" to false ..." );
      split        = false;
      pure         = false;
   }

   // configuration file

   Setup setup;
   setup.pseudo2D = (lm == LMD2Q9);

   setup.xCells  = configBlock.getParameter< uint_t >( "xCells", uint_t(16) );
   setup.yzCells = configBlock.getParameter< uint_t >( "yzCells", uint_t(16) );

   setup.H = configBlock.getParameter< real_t >( "H", real_c(0.41) ); // [m]
   setup.L = configBlock.getParameter< real_t >( "L", real_c(2.50) ); // [m]
   
   setup.strictlyObeyL = configBlock.getParameter< bool >( "strictlyObeyL", true );

   setup.yzBlocks = configBlock.getParameter< uint_t >( "yzBlocks", uint_t(3) );
   
   const real_t xBlocks = ( setup.L / ( setup.H / ( real_c(setup.yzBlocks) * real_c(setup.yzCells) ) ) ) / real_c(setup.xCells);
   setup.xBlocks = uint_c( xBlocks );
   if( xBlocks > real_c( setup.xBlocks ) ) setup.xBlocks += uint_t(1);
   
   setup.cylinderxPosition = configBlock.getParameter< real_t >( "cylinderxPosition", real_c(0.5) ); // [m]
   setup.cylinderyPosition = configBlock.getParameter< real_t >( "cylinderyPosition", real_c(0.2) ); // [m]
   setup.cylinderRadius = configBlock.getParameter< real_t >( "cylinderRadius", real_c(0.05) ); // [m]
   setup.circularCrossSection = configBlock.getParameter< bool >( "circularCrossSection", true );

   setup.viscosity = configBlock.getParameter< real_t >( "kinViscosity", real_c(0.001) ); // [m^2 / s]
   setup.rho = configBlock.getParameter< real_t >( "rho", real_t(1) ); // [kg / m^3]
   setup.inflowVelocity = configBlock.getParameter< real_t >( "inflowVelocity", real_t(0.45) ); // [m/s]
   setup.raisingTime = configBlock.getParameter< real_t >( "raisingTime", real_t(0) ); // [s]
   setup.sinPeriod = configBlock.getParameter< real_t >( "sinPeriod", real_t(0) ); // [s]
   setup.dx = setup.H / real_c( setup.yzBlocks * setup.yzCells );

   setup.evaluateForceComponents = configBlock.getParameter< bool >( "evaluateForceComponents", false );
   setup.nbrOfEvaluationPointsForCoefficientExtremas = configBlock.getParameter< uint_t >( "nbrOfEvaluationPointsForCoefficientExtremas", uint_t(100) );
   setup.nbrOfEvaluationPointsForCoefficientExtremas = std::max( setup.nbrOfEvaluationPointsForCoefficientExtremas, uint_t(1) );

   setup.evaluatePressure = configBlock.getParameter< bool >( "evaluatePressure", false );
   setup.pAlpha = configBlock.getParameter< Vector3<real_t> >( "pAlpha", Vector3<real_t>( real_c(0.45), real_c(0.2), real_c(0.205) ) );
   setup.pOmega = configBlock.getParameter< Vector3<real_t> >( "pOmega", Vector3<real_t>( real_c(0.55), real_c(0.2), real_c(0.205) ) );

   setup.evaluateStrouhal = configBlock.getParameter< bool >( "evaluateStrouhal", false );
   setup.pStrouhal = configBlock.getParameter< Vector3<real_t> >( "pStrouhal", Vector3<real_t>( real_c(1), real_c(0.325), real_c(0.205) ) );

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

   if( configBlock.getParameter< bool >( "useCylinderForRefinement", false ) )
   {
      if( !configBlock.isDefined("cylinderRefinementLevel")  )
         WALBERLA_ABORT( "You have to specify \'cylinderRefinementLevel\' in the \"SchaeferTurek\" block of the configuration file (" << argv[1] << ")" );

      const real_t cylinderRefinementBuffer = configBlock.getParameter< real_t >( "cylinderRefinementBuffer", real_t(0) );

      Cylinder cylinder( setup );
      CylinderRefinementSelection cylinderRefinementSelection( cylinder, configBlock.getParameter< uint_t >( "cylinderRefinementLevel" ),
               cylinderRefinementBuffer );

      refinementSelectionFunctions.add( cylinderRefinementSelection );
   }

   refinementSelectionFunctions.add( std::bind( setInflowOutflowToSameLevel, std::placeholders::_1, setup ) );

   if( setup.pseudo2D )
      refinementSelectionFunctions.add( Pseudo2DRefinementSelectionCorrection );

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

   // executing benchmark

   real_t omega = configBlock.getParameter< real_t >( "omega", real_t(1.4) );

   const real_t magicNumber = configBlock.getParameter< real_t >( "magicNumber", real_t(3) / real_t(16) );

   const real_t lambda_e = configBlock.getParameter< real_t >( "lambda_e", real_t(1.4) );
   const real_t lambda_d = configBlock.getParameter< real_t >( "lambda_d", real_t(1.4) );

   const real_t s1  = configBlock.getParameter< real_t >( "s1",  real_t(1.4) );
   const real_t s2  = configBlock.getParameter< real_t >( "s2",  real_t(1.4) );
   const real_t s4  = configBlock.getParameter< real_t >( "s4",  real_t(1.4) );
   const real_t s9  = configBlock.getParameter< real_t >( "s9",  real_t(1.4) );
   const real_t s10 = configBlock.getParameter< real_t >( "s10", real_t(1.4) );
   const real_t s16 = configBlock.getParameter< real_t >( "s16", real_t(1.4) );

   uint_t relaxationParametersLevel = configBlock.getParameter< uint_t >( "relaxationParametersLevel", uint_t(0) );

   if( configBlock.isDefined("latticeInflowVelocity") )
   {
      const real_t latticeInflowVelocity = configBlock.getParameter< real_t >( "latticeInflowVelocity" );
      const real_t dt = latticeInflowVelocity * setup.dx / setup.inflowVelocity;
      omega = real_t(1) / ( ( real_c(3) * dt * setup.viscosity ) / ( setup.dx * setup.dx ) + real_c(0.5) );
      relaxationParametersLevel = uint_t(0);
   }

   if( collisionModel == CMSRT ) // SRT
   {
      if( lm == LMD2Q9 )
      {
         if( compressible )
         {
            D2Q9_SRT_COMP latticeModel = D2Q9_SRT_COMP( lbm::collision_model::SRT( omega, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D2Q9_SRT_INCOMP latticeModel = D2Q9_SRT_INCOMP( lbm::collision_model::SRT( omega, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
      else if( lm == LMD3Q15 )
      {
         if( compressible )
         {
            D3Q15_SRT_COMP latticeModel = D3Q15_SRT_COMP( lbm::collision_model::SRT( omega, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q15_SRT_INCOMP latticeModel = D3Q15_SRT_INCOMP( lbm::collision_model::SRT( omega, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
      else if( lm == LMD3Q19 )
      {
         if( compressible )
         {
            D3Q19_SRT_COMP latticeModel = D3Q19_SRT_COMP( lbm::collision_model::SRT( omega, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q19_SRT_INCOMP latticeModel = D3Q19_SRT_INCOMP( lbm::collision_model::SRT( omega, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
      else
      {
         if( compressible )
         {
            D3Q27_SRT_COMP latticeModel = D3Q27_SRT_COMP( lbm::collision_model::SRT( omega, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q27_SRT_INCOMP latticeModel = D3Q27_SRT_INCOMP( lbm::collision_model::SRT( omega, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
   }
   else if( collisionModel == CMTRT ) // TRT
   {
      bool useLambdas = false;
      if( ( ! configBlock.isDefined("latticeInflowVelocity") ) && configBlock.isDefined("lambda_e") && configBlock.isDefined("lambda_d") )
         useLambdas = true;

      if( lm == LMD2Q9 )
      {
         if( compressible )
         {
            D2Q9_TRT_COMP latticeModel = useLambdas ? D2Q9_TRT_COMP( lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) ) :
                                                      D2Q9_TRT_COMP( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D2Q9_TRT_INCOMP latticeModel = useLambdas ? D2Q9_TRT_INCOMP( lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) ) :
                                                        D2Q9_TRT_INCOMP( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
      else if( lm == LMD3Q15 )
      {
         if( compressible )
         {
            D3Q15_TRT_COMP latticeModel = useLambdas ? D3Q15_TRT_COMP( lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) ) :
                                                       D3Q15_TRT_COMP( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q15_TRT_INCOMP latticeModel = useLambdas ? D3Q15_TRT_INCOMP( lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) ) :
                                                         D3Q15_TRT_INCOMP( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
      else if( lm == LMD3Q19 )
      {
         if( compressible )
         {
            D3Q19_TRT_COMP latticeModel = useLambdas ? D3Q19_TRT_COMP( lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) ) :
                                                       D3Q19_TRT_COMP( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q19_TRT_INCOMP latticeModel = useLambdas ? D3Q19_TRT_INCOMP( lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) ) :
                                                         D3Q19_TRT_INCOMP( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
      else
      {
         if( compressible )
         {
            D3Q27_TRT_COMP latticeModel = useLambdas ? D3Q27_TRT_COMP( lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) ) :
                                                       D3Q27_TRT_COMP( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
         else
         {
            D3Q27_TRT_INCOMP latticeModel = useLambdas ? D3Q27_TRT_INCOMP( lbm::collision_model::TRT( lambda_e, lambda_d, relaxationParametersLevel ) ) :
                                                         D3Q27_TRT_INCOMP( lbm::collision_model::TRT::constructWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
            run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
         }
      }
   }
   else // MRT
   {
      WALBERLA_CHECK( !compressible );
      WALBERLA_CHECK( !split );
      WALBERLA_CHECK( !pure );
      WALBERLA_CHECK( lm == LMD3Q19 );

      bool useS = false;
      if( ( ! configBlock.isDefined("latticeInflowVelocity") ) &&
          configBlock.isDefined("s1") && configBlock.isDefined("s2")  && configBlock.isDefined("s4") &&
          configBlock.isDefined("s9") && configBlock.isDefined("s10") && configBlock.isDefined("s16") )
         useS = true;

      D3Q19_MRT_INCOMP latticeModel = useS ? D3Q19_MRT_INCOMP( lbm::collision_model::D3Q19MRT( s1, s2, s4, s9, s10, s16, relaxationParametersLevel ) ) :
                                             D3Q19_MRT_INCOMP( lbm::collision_model::D3Q19MRT::constructTRTWithMagicNumber( omega, magicNumber, relaxationParametersLevel ) );
      run( config, latticeModel, split, pure, fzyx, syncComm, fullComm, linearExplosion, refinementSelectionFunctions, setup, memoryPerCell, processMemoryLimit );
   }

   logging::Logging::printFooterOnStream();
   return EXIT_SUCCESS;
}

} // namespace schaefer_turek

int main( int argc, char ** argv )
{
   return schaefer_turek::main( argc, argv );
}
