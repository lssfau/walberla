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
//! \file AMRSedimentSettling.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/loadbalancing/InfoCollection.h"
#include "blockforest/loadbalancing/DynamicCurve.h"
#include "blockforest/loadbalancing/DynamicDiffusive.h"
#include "blockforest/loadbalancing/DynamicParMetis.h"
#include "blockforest/loadbalancing/StaticCurve.h"
#include "blockforest/loadbalancing/StaticParMetis.h"
#include "blockforest/loadbalancing/weight_assignment/MetisAssignmentFunctor.h"
#include "blockforest/loadbalancing/weight_assignment/WeightAssignmentFunctor.h"
#include "blockforest/AABBRefinementSelection.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/mpi/Broadcast.h"

#include "domain_decomposition/BlockSweepWrapper.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/VelocityFieldWriter.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/refinement/all.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "pe/basic.h"
#include "pe/Types.h"
#include "pe/fcd/GJKEPACollideFunctor.h"
#include "pe/vtk/SphereVtkOutput.h"
#include "pe/vtk/EllipsoidVtkOutput.h"
#include "pe/cr/ICR.h"
#include "pe/synchronization/ClearSynchronization.h"
#include "pe/amr/InfoCollection.h"

#include "pe_coupling/amr/all.h"
#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

namespace amr_sediment_settling
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
using LatticeModel_T = lbm::D3Q19<lbm::collision_model::TRT, false>;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
using BodyField_T = GhostLayerField<pe::BodyID, 1>;
using VelocityField_T = GhostLayerField<Vector3<real_t>, 1>;

const uint_t FieldGhostLayers = 4;

// boundary handling
using NoSlip_T = lbm::NoSlip<LatticeModel_T, flag_t>;

using MO_T = pe_coupling::CurvedLinear<LatticeModel_T, FlagField_T>;

using BoundaryHandling_T = BoundaryHandling<FlagField_T, Stencil_T, NoSlip_T, MO_T>;

using BodyTypeTuple = std::tuple<pe::Sphere, pe::Ellipsoid, pe::Plane>;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID NoSlip_Flag( "no slip" );
const FlagUID MO_Flag( "moving obstacle" );
const FlagUID FormerMO_Flag( "former moving obstacle" );


//////////////////////////////////////
// DYNAMIC REFINEMENT FUNCTIONALITY //
//////////////////////////////////////


/*
 * Refinement check based on gradient magnitude
 * If gradient magnitude is below lowerLimit in all cells of a block, that block could be coarsened.
 * If the gradient value is above the upperLimit for at least one cell, that block gets marked for refinement.
 * Else, the block remains on the current level.
 */
template< typename LatticeModel_T, typename Filter_T >
class VectorGradientRefinement
{
public:
   using VectorField_T = GhostLayerField<Vector3<real_t>, 1>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   VectorGradientRefinement( const ConstBlockDataID & fieldID, const Filter_T & filter,
                             const real_t upperLimit, const real_t lowerLimit, const uint_t maxLevel ) :
         fieldID_( fieldID ), filter_( filter ),
         upperLimit_( upperLimit ), lowerLimit_( lowerLimit ), maxLevel_( maxLevel )
   {}

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > & blocksAlreadyMarkedForRefinement,
                    const BlockForest & forest );

private:

   ConstBlockDataID fieldID_;

   Filter_T filter_;

   real_t upperLimit_;
   real_t lowerLimit_;

   uint_t maxLevel_;

}; // class GradientRefinement

template< typename LatticeModel_T, typename Filter_T >
void VectorGradientRefinement< LatticeModel_T, Filter_T >::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                                                       std::vector< const Block * > &, const BlockForest & )
{
   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      const Block * const block = it->first;

      const uint_t currentLevelOfBlock = block->getLevel();

      const VectorField_T * uField = block->template getData< VectorField_T >( fieldID_ );

      if( uField == nullptr )
      {
         it->second = uint_t(0);
         continue;
      }

      Matrix3<real_t> uGradient( real_t(0) );

      bool refine( false );
      bool coarsen( true );

      filter_( *block );

      WALBERLA_FOR_ALL_CELLS_XYZ( uField,

          std::vector< Vector3<real_t> > uValues( Stencil_T::Size, Vector3<real_t>(real_t(0)) );

          Vector3<real_t> uInCenterCell = uField->get( x,y,z );

          for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
          {
             // check if boundary treatment is necessary
             if( filter_( x+dir.cx(),y+dir.cy(),z+dir.cz() ) )
             {
                // copy from center cell
                uValues[ *dir ] = uInCenterCell;
             } else {
                uValues[ *dir ] = uField->get( x+dir.cx(),y+dir.cy(),z+dir.cz() );
             }
          }

          // obtain the matrix grad(u) with the help of the gradient formula from
          // See: Ramadugu et al - Lattice differential operators for computational physics (2013)
          // with T = c_s**2
          const auto inv_c_s_sqr = real_t(3);
          uGradient = real_t(0);
          for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
          {
             real_t cx = real_c(dir.cx());
             real_t cy = real_c(dir.cy());
             real_t cz = real_c(dir.cz());

             // grad(ux)
             real_t ux = uValues[ *dir ][0];
             uGradient[ 0 ] += LatticeModel_T::w[ dir.toIdx() ] * cx * ux;
             uGradient[ 3 ] += LatticeModel_T::w[ dir.toIdx() ] * cy * ux;
             uGradient[ 6 ] += LatticeModel_T::w[ dir.toIdx() ] * cz * ux;

             // grad(uy)
             real_t uy = uValues[ *dir ][1];
             uGradient[ 1 ] += LatticeModel_T::w[ dir.toIdx() ] * cx * uy;
             uGradient[ 4 ] += LatticeModel_T::w[ dir.toIdx() ] * cy * uy;
             uGradient[ 7 ] += LatticeModel_T::w[ dir.toIdx() ] * cz * uy;

             // grad(uz)
             real_t uz = uValues[ *dir ][2];
             uGradient[ 2 ] += LatticeModel_T::w[ dir.toIdx() ] * cx * uz;
             uGradient[ 5 ] += LatticeModel_T::w[ dir.toIdx() ] * cy * uz;
             uGradient[ 8 ] += LatticeModel_T::w[ dir.toIdx() ] * cz * uz;

          }
          uGradient *= inv_c_s_sqr;

          auto norm = real_t(0);
          //compute maximums norm of 3x3 matrix
          for (auto i = uint_t(0); i < uint_t(3*3); ++i)
             norm = std::max(norm, std::fabs(uGradient[i]));

          if( norm > lowerLimit_ )
          {
             coarsen = false;
             if( norm > upperLimit_ )
                refine = true;
          }

      )

      if( refine && currentLevelOfBlock < maxLevel_ )
      {
         WALBERLA_ASSERT( !coarsen );
         it->second = currentLevelOfBlock + uint_t(1);
      }
      if( coarsen && currentLevelOfBlock > uint_t(0) )
      {
         WALBERLA_ASSERT( !refine );
         it->second = currentLevelOfBlock - uint_t(1);
      }
   }
}


// Load estimators for spheres and ellipsoids, obtained at SuperMUC Phase 2
// See Sec. 3 in the paper for more infos

/////////////
// Spheres //
/////////////
real_t fittedLBMWeightEvaluationFunctionSpheres(const pe_coupling::BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t F  = blockInfo.numberOfFluidCells;

   // from fits with optimized D3Q19 LBM kernel (no forces)
   real_t weight = real_t(9.990957667738165e-06) * real_c(Ce) + real_t(0.00015749920523711047) * real_c(F) + real_t(-0.08232498905584973);
   return weight;
}

real_t fittedBHWeightEvaluationFunctionSpheres(const pe_coupling::BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t NB = blockInfo.numberOfNearBoundaryCells;

   // from fits with optimized D3Q19 LBM kernel (no forces)
   real_t weight = real_t(6.654810986939097e-06) * real_c(Ce) + real_t(0.0007061414693533274) * real_c(NB) + real_t(-0.1094292992294259);
   return weight;
}

real_t fittedCoupling1WeightEvaluationFunctionSpheres(const pe_coupling::BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t F  = blockInfo.numberOfFluidCells;
   uint_t Pl = blockInfo.numberOfLocalBodies;
   uint_t Ps = blockInfo.numberOfShadowBodies;

   // from fits with optimized D3Q19 LBM kernel (no forces)
   real_t weight = real_t(3.07542641675429e-06) * real_c(Ce) + real_t(2.419364600880769e-07) * real_c(F) + real_t(0.01413718259604757) * real_c(Pl) + real_t(0.027761707343462727) * real_c(Ps) + real_t(-0.13991481483939272);
   return weight;
}

real_t fittedCoupling2WeightEvaluationFunctionSpheres(const pe_coupling::BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t F  = blockInfo.numberOfFluidCells;
   uint_t Pl = blockInfo.numberOfLocalBodies;
   uint_t Ps = blockInfo.numberOfShadowBodies;

   // from fits with optimized D3Q19 LBM kernel (no forces)
   real_t weight = real_t(5.988401232749505e-06) * real_c(Ce) + real_t(3.903532223977357e-06) * real_c(F) + real_t(-0.008802674250816316) * real_c(Pl) + real_t(0.02505020738346139) * real_c(Ps) + real_t(-0.12970723676003335);
   return weight;
}

real_t fittedPEWeightEvaluationFunctionSpheres(const pe_coupling::BlockInfo& blockInfo)
{
   uint_t Pl = blockInfo.numberOfLocalBodies;
   uint_t Ps = blockInfo.numberOfShadowBodies;
   uint_t Ct = blockInfo.numberOfContacts;
   uint_t Sc = blockInfo.numberOfPeSubCycles;

   // from fits with optimized D3Q19 LBM kernel (no forces)
   real_t cPlPs2 = real_t(1.1562854544700417e-06);
   real_t cPl    = real_t(0.0009620525068318354);
   real_t cPs    = real_t(0.00027549401081063894);
   real_t cCt    = real_t(0.0014801932788115464);
   real_t c      = real_t(0.01883682418448259);
   real_t weight = real_c(Sc) * ( cPlPs2 * real_c(Pl+Ps) * real_c(Pl+Ps) + cPl * real_c(Pl) + cPs * real_c(Ps) + cCt * real_c(Ct) + c );
   return weight;
}

real_t fittedTotalWeightEvaluationFunctionSpheres(const pe_coupling::BlockInfo& blockInfo)
{
   return fittedLBMWeightEvaluationFunctionSpheres(blockInfo) + fittedBHWeightEvaluationFunctionSpheres(blockInfo) +
          fittedCoupling1WeightEvaluationFunctionSpheres(blockInfo) + fittedCoupling2WeightEvaluationFunctionSpheres(blockInfo) +
          fittedPEWeightEvaluationFunctionSpheres(blockInfo);
}

////////////////
// Ellipsoids //
////////////////
real_t fittedLBMWeightEvaluationFunctionEllipsoids(const pe_coupling::BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t F  = blockInfo.numberOfFluidCells;
   uint_t NB = blockInfo.numberOfNearBoundaryCells;

   // from fits with optimized D3Q19 LBM kernel (no forces)
   real_t cCe = real_t(4.69973868717e-05);
   real_t cF  = real_t(0.000110568537442);
   real_t weight = cCe * real_c(Ce) + cF * real_c(F);
   if( NB > uint_t(0) ) weight += real_t(5.96551488486e-05) * real_c(Ce) + real_t(-5.75351782026e-05) * real_c(F) + real_t(0.000695800745231) * real_c(NB);
   return weight;
}

real_t fittedCouplingWeightEvaluationFunctionEllipsoids(const pe_coupling::BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t F  = blockInfo.numberOfFluidCells;
   uint_t Pl = blockInfo.numberOfLocalBodies;
   uint_t Ps = blockInfo.numberOfShadowBodies;

   // from fits with optimized D3Q19 LBM kernel (no forces)
   real_t cCe = real_t(0.000176674935526);
   real_t cF  = real_t(-0.000170513513027);
   real_t cPl = real_t(0.0252031634776);
   real_t cPs = real_t(0.0356835220918);
   real_t weight = real_t(0);
   if( (Pl + Ps ) > uint_t(0) ) weight = cCe * real_c(Ce) + cF * real_c(F) + cPl * real_c(Pl) + cPs * real_c(Ps);
   return weight;
}

real_t fittedPEWeightEvaluationFunctionEllipsoids(const pe_coupling::BlockInfo& blockInfo)
{
   uint_t Pl = blockInfo.numberOfLocalBodies;
   uint_t Ps = blockInfo.numberOfShadowBodies;
   uint_t Ct = blockInfo.numberOfContacts;
   uint_t Sc = blockInfo.numberOfPeSubCycles;

   // from fits with optimized D3Q19 LBM kernel (no forces)
   real_t cPlPs2 = real_t(8.24153555785e-06);
   real_t cPl    = real_t(0.00135966650494);
   real_t cPs    = real_t(0.00440464092538);
   real_t cCt    = real_t(0.0216278259881);
   real_t weight = real_c(Sc) * ( cPlPs2 * real_c(Pl+Ps) * real_c(Pl+Ps) + cPl * real_c(Pl) + cPs * real_c(Ps) + cCt * real_c(Ct) );
   return weight;
}

real_t fittedTotalWeightEvaluationFunctionEllipsoids(const pe_coupling::BlockInfo& blockInfo)
{
   return fittedLBMWeightEvaluationFunctionEllipsoids(blockInfo) +
          fittedCouplingWeightEvaluationFunctionEllipsoids(blockInfo) +
          fittedPEWeightEvaluationFunctionEllipsoids(blockInfo);
}


struct TimingPoolLogger
{
   TimingPoolLogger( const shared_ptr<WcTimingPool> & timingPool, const shared_ptr<SweepTimeloop> & timeloop, const uint_t interval )
         : timingPool_( timingPool ), timeloop_( timeloop ), interval_( interval )
   {
   }

   void operator()()
   {
      if( interval_ > uint_t(0) && timeloop_->getCurrentTimeStep() % interval_ == uint_t(0) )
      {
         timingPool_->logResultOnRoot();
      }
   }

private:
   shared_ptr<WcTimingPool> timingPool_;
   shared_ptr<SweepTimeloop> timeloop_;
   uint_t interval_;
};

struct TimingTreeLogger
{
   TimingTreeLogger( const shared_ptr<WcTimingTree> & timingTree, const shared_ptr<SweepTimeloop> & timeloop, const uint_t interval )
         : timingTree_( timingTree ), timeloop_( timeloop ), interval_( interval )
   {
   }

   void operator()()
   {
      if( interval_ > uint_t(0) && timeloop_->getCurrentTimeStep() % interval_ == uint_t(0) )
      {
         timingTree_->synchronize();
         auto reducedTimingTree = timingTree_->getReduced();
         WALBERLA_LOG_INFO_ON_ROOT( reducedTimingTree );
      }
   }

private:
   shared_ptr<WcTimingTree> timingTree_;
   shared_ptr<SweepTimeloop> timeloop_;
   uint_t interval_;
};

/////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void refinementSelection( SetupBlockForest& forest, uint_t levels, const AABB & refinementBox )
{
   auto dx = real_t(1); // dx on finest level
   for (auto &block : forest) {
      uint_t blockLevel = block.getLevel();
      uint_t levelScalingFactor = ( uint_t(1) << (levels - uint_t(1) - blockLevel) );
      real_t dxOnLevel = dx * real_c(levelScalingFactor);
      AABB blockAABB = block.getAABB();

      // extend block AABB by ghostlayers
      AABB extendedBlockAABB = blockAABB.getExtended( dxOnLevel * real_c(FieldGhostLayers) );

      if( extendedBlockAABB.intersects( refinementBox ) )
         if( blockLevel < ( levels - uint_t(1) ) )
            block.setMarker( true );
   }
}

static void workloadAndMemoryAssignment( SetupBlockForest& forest )
{
   for (auto &block : forest) {
      block.setWorkload( numeric_cast< workload_t >( uint_t(1) << block.getLevel() ) );
      block.setMemory( numeric_cast< memory_t >(1) );
   }
}

static shared_ptr< StructuredBlockForest > createBlockStructure( const AABB & domainAABB, Vector3<uint_t> blockSizeInCells,
                                                                 uint_t numberOfLevels, const AABB & refinementBox,
                                                                 bool useBox, const std::string & loadDistributionStrategy,
                                                                 bool keepGlobalBlockInformation = false )
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

   WALBERLA_LOG_INFO_ON_ROOT(" - refinement box: " << refinementBox);

   MPIManager::instance()->useWorldComm();

   sforest.addRefinementSelectionFunction( std::bind( refinementSelection, std::placeholders::_1, numberOfLevels, refinementBox ) );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

   Vector3<bool> periodicity( true, true, false);
   if( useBox )
   {
      periodicity[0] = false;
      periodicity[1] = false;
   }
   sforest.init( domainAABB,
                 numberOfCoarseBlocksPerDirection[0], numberOfCoarseBlocksPerDirection[1], numberOfCoarseBlocksPerDirection[2],
                 periodicity[0], periodicity[1], periodicity[2]);

   // calculate process distribution
   const memory_t memoryLimit = math::Limits< memory_t >::inf();

   if( loadDistributionStrategy == "Hilbert" )
   {
      bool useHilbert = true;
      sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(useHilbert), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );
   } else if ( loadDistributionStrategy == "Morton" )
   {
      bool useHilbert = false;
      sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(useHilbert), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );
   } else if ( loadDistributionStrategy == "ParMetis" )
   {
      blockforest::StaticLevelwiseParMetis::Algorithm algorithm = blockforest::StaticLevelwiseParMetis::Algorithm::PARMETIS_PART_GEOM_KWAY;
      blockforest::StaticLevelwiseParMetis staticParMetis(algorithm);
      sforest.balanceLoad( staticParMetis, uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );
   } else if (loadDistributionStrategy == "Diffusive" )
   {
      // also use Hilbert curve here
      bool useHilbert = true;
      sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(useHilbert), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );
   } else
   {
      WALBERLA_ABORT("Load distribution strategy \"" << loadDistributionStrategy << "\t not implemented! - Aborting" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( sforest );


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
class MyBoundaryHandling : public blockforest::AlwaysInitializeBlockDataHandling< BoundaryHandling_T >
{
public:
   MyBoundaryHandling( const weak_ptr< StructuredBlockStorage > & blocks,
                       const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyFieldID ) :
         blocks_( blocks ), flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), bodyFieldID_ ( bodyFieldID )
   {}

   BoundaryHandling_T * initialize( IBlock * const block ) override;

private:

   weak_ptr< StructuredBlockStorage > blocks_;

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;


}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::initialize( IBlock * const block )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   auto * flagField = block->getData< FlagField_T >( flagFieldID_ );
   auto *  pdfField = block->getData< PdfField_T > ( pdfFieldID_ );
   auto * bodyField = block->getData< BodyField_T >( bodyFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   auto blocksPtr = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocksPtr );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "moving obstacle boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                                           MO_T( "MO", MO_Flag, pdfField, flagField, bodyField, fluid, *blocksPtr, *block ),
                                                           BoundaryHandling_T::Mode::ENTIRE_FIELD_TRAVERSAL);

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}


//*******************************************************************************************************************


//*******************************************************************************************************************
/*!\brief Evaluating the position and velocity of the sediments
 *
 */
//*******************************************************************************************************************
class PropertyLogger
{
public:
   PropertyLogger( const shared_ptr<SweepTimeloop> & timeloop, const shared_ptr< StructuredBlockStorage > & blocks,
                   const BlockDataID & bodyStorageID, const std::string & fileName, bool fileIO) :
      timeloop_( timeloop ), blocks_( blocks ), bodyStorageID_( bodyStorageID ), fileName_( fileName ), fileIO_(fileIO),
      meanPos_( real_t(0) ), meanVel_( real_t(0) ), maxVel_( real_t(0) )
   {
      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str() );
            file << "#\t pos\t vel\t maxVel\n";
            file.close();
         }
      }
   }

   void operator()()
   {
      const uint_t timestep (timeloop_->getCurrentTimeStep() );

      auto numSediments = uint_t(0);
      auto meanPos = real_t(0);
      auto meanVel = real_t(0);
      auto maxVel = real_t(0);

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            meanPos += bodyIt->getPosition()[2];
            meanVel += bodyIt->getLinearVel()[2];
            maxVel = std::max(maxVel, std::fabs(bodyIt->getLinearVel()[2]));
            ++numSediments;
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( numSediments, mpi::SUM );
         mpi::allReduceInplace( meanPos, mpi::SUM );
         mpi::allReduceInplace( meanVel, mpi::SUM );
         mpi::allReduceInplace( maxVel, mpi::MAX );
      }

      meanPos /= real_c(numSediments);
      meanVel /= real_c(numSediments);

      meanPos_ = meanPos;
      meanVel_ = meanVel;
      maxVel_ = maxVel;

      if( fileIO_ )
         writeToFile( timestep );
   }

   real_t getMeanPosition() const
   {
      return meanPos_;
   }

   real_t getMaxVelocity() const
   {
      return maxVel_;
   }

   real_t getMeanVelocity() const
   {
      return meanVel_;
   }


private:
   void writeToFile( uint_t timestep )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName_.c_str(), std::ofstream::app );

         file << timestep << "\t" << meanPos_ << "\t" << meanVel_ << "\t" << maxVel_ << "\n";
         file.close();
      }
   }

   shared_ptr<SweepTimeloop> timeloop_;
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;
   std::string fileName_;
   bool fileIO_;

   real_t meanPos_;
   real_t meanVel_;
   real_t maxVel_;
};

void clearBoundaryHandling( BlockForest & forest, const BlockDataID & boundaryHandlingID )
{
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      auto * boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID);
      boundaryHandling->clear( FieldGhostLayers );
   }
}

void clearBodyField( BlockForest & forest, const BlockDataID & bodyFieldID )
{
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      auto * bodyField = blockIt->getData<BodyField_T>(bodyFieldID);
      bodyField->setWithGhostLayer( NULL );
   }
}

void recreateBoundaryHandling( BlockForest & forest, const BlockDataID & boundaryHandlingID )
{
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      auto * boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID);
      boundaryHandling->fillWithDomain( FieldGhostLayers );
   }
}


class TimingEvaluator
{
public:
   TimingEvaluator( const shared_ptr<WcTimingPool> & levelwiseTimingPool, const shared_ptr<WcTimingTree> & peTimingTree, uint_t numberOfLevels)
   : levelwiseTimingPool_( levelwiseTimingPool ), peTimingTree_( peTimingTree ), numberOfLevels_( numberOfLevels )
   {}

   real_t getTimings(const std::vector<std::string> & timerNames, uint_t level )
   {

      auto timing = real_t(0);
      for (const auto &timerName : timerNames)
      {
         std::string timerNameLvlWise = timerName;// +
         // put level between timer string and possible suffix
         auto suffixBegin = timerNameLvlWise.find_first_of('[');
         if( suffixBegin != std::string::npos)
         {
            // suffix detected
            auto suffixEnd = timerNameLvlWise.find_last_of(']');
            if( suffixEnd != std::string::npos)
            {
               auto timerString = timerNameLvlWise.substr(0,suffixBegin);
               auto suffixString = timerNameLvlWise.substr(suffixBegin,suffixEnd-suffixBegin+1);

               timerNameLvlWise = timerString + "(" + std::to_string(level) + ") " + suffixString; // NOLINT

            }
            else
            {
               WALBERLA_ABORT("Invalid timer string");
            }
         }
         else
         {
            timerNameLvlWise += " (" + std::to_string(level) + ")";;
         }

         if( levelwiseTimingPool_->timerExists(timerNameLvlWise))
            timing += real_c((*levelwiseTimingPool_)[timerNameLvlWise].total());

         if( level == numberOfLevels_- 1)
         {
            if( peTimingTree_->timerExists(timerName))
               timing += real_c((*peTimingTree_)[timerName].total());
         }
      }

      return timing;
   }


private:

   shared_ptr<WcTimingPool> levelwiseTimingPool_;
   shared_ptr<WcTimingTree> peTimingTree_;
   uint_t numberOfLevels_;
};


real_t weightEvaluation(BlockForest & forest,
                        const shared_ptr<pe_coupling::InfoCollection>& couplingInfoCollection,
                        const shared_ptr<blockforest::InfoCollection> & peInfoCollection,
                        real_t peBlockBaseWeight,
                        const std::string & loadEvaluationStrategy,
                        uint_t level,
                        bool useEllipsoids )
{
   auto weight = real_t(0);
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      if( forest.getLevel(*blockIt) != level) continue;


      auto * block = static_cast<blockforest::Block*> (&(*blockIt));
      const auto &blockID = block->getId();

      if(loadEvaluationStrategy == "LBM")
      {
         auto infoIt = couplingInfoCollection->find( blockID );
         weight += pe_coupling::amr::defaultWeightEvaluationFunction(infoIt->second);

      }else if(loadEvaluationStrategy == "PE")
      {
         auto infoIt = peInfoCollection->find( blockID );
         weight += real_c(infoIt->second.computationalWeight) + peBlockBaseWeight;
      }else if(loadEvaluationStrategy == "Fit" || loadEvaluationStrategy == "FitMulti")
      {
         auto infoIt = couplingInfoCollection->find( blockID );
         if( useEllipsoids )
         {
            weight += fittedTotalWeightEvaluationFunctionEllipsoids(infoIt->second);
         }
         else
         {
            weight += fittedTotalWeightEvaluationFunctionSpheres(infoIt->second);
         }
      }else
      {
         WALBERLA_ABORT("Load balancing strategy not defined");
      }
   }
   return weight;
}


uint_t evaluateEdgeCut(BlockForest & forest)
{

   //note: only works for edges in uniform grids

   auto edgecut = uint_t(0); // = edge weights between processes

   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      auto * block = static_cast<blockforest::Block*> (&(*blockIt));

      real_t blockVolume = block->getAABB().volume();
      real_t approximateEdgeLength = std::cbrt( blockVolume );

      uint_t faceNeighborWeight = uint_c(approximateEdgeLength * approximateEdgeLength ); //common face
      uint_t edgeNeighborWeight = uint_c(approximateEdgeLength); //common edge
      uint_t cornerNeighborWeight = uint_c( 1 ); //common corner


      for( const uint_t idx : blockforest::getFaceNeighborhoodSectionIndices() )
      {
         for (auto nb = uint_t(0); nb < block->getNeighborhoodSectionSize(idx); ++nb)
         {
            if( block->neighborExistsRemotely(idx,nb) ) edgecut += faceNeighborWeight;
         }
      }

      for( const uint_t idx : blockforest::getEdgeNeighborhoodSectionIndices() )
      {
         for (auto nb = uint_t(0); nb < block->getNeighborhoodSectionSize(idx); ++nb)
         {
            if( block->neighborExistsRemotely(idx,nb) ) edgecut += edgeNeighborWeight;
         }
      }

      for( const uint_t idx : blockforest::getCornerNeighborhoodSectionIndices() )
      {
         for (auto nb = uint_t(0); nb < block->getNeighborhoodSectionSize(idx); ++nb)
         {
            if( block->neighborExistsRemotely(idx,nb) ) edgecut += cornerNeighborWeight;
         }
      }
   }
   return edgecut;
}


void evaluateTotalSimulationTimePassed(WcTimingPool & timeloopTimingPool, real_t & totalSimTime, real_t & totalLBTime)
{
   shared_ptr< WcTimingPool> reduced = timeloopTimingPool.getReduced(timing::REDUCE_TOTAL, 0);

   std::string simulationString("LBM refinement time step");
   auto totalTime = real_t(0);
   WALBERLA_ROOT_SECTION(){
      totalTime = real_c((*reduced)[simulationString].total());
   }
   totalSimTime = totalTime;

   std::string lbString("refinement checking");
   auto lbTime = real_t(0);
   WALBERLA_ROOT_SECTION(){
      lbTime = real_c((*reduced)[lbString].total());
   }
   totalLBTime = lbTime;

}

void createSedimentLayer(uint_t numberOfSediments, const AABB & generationDomain, real_t diameter, real_t heightBorder,
                         pe::MaterialID peMaterial,
                         pe::cr::HCSITS & cr, const std::function<void(void)> & syncCall,
                         const shared_ptr< StructuredBlockForest > & blocks,
                         const shared_ptr<pe::BodyStorage> & globalBodyStorage, BlockDataID bodyStorageID,
                         real_t gravitationalAcceleration, bool useEllipsoids, bool shortRun)
{
   WALBERLA_LOG_INFO_ON_ROOT("Starting creation of sediments");

   auto xParticle = real_t(0);
   auto yParticle = real_t(0);
   auto zParticle = real_t(0);

   for( uint_t nSed = 0; nSed < numberOfSediments; ++nSed )
   {

      WALBERLA_ROOT_SECTION()
      {
         xParticle = math::realRandom<real_t>(generationDomain.xMin(), generationDomain.xMax());
         yParticle = math::realRandom<real_t>(generationDomain.yMin(), generationDomain.yMax());
         zParticle = math::realRandom<real_t>(generationDomain.zMin(), generationDomain.zMax());
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::broadcastObject( xParticle );
         mpi::broadcastObject( yParticle );
         mpi::broadcastObject( zParticle );
      }

      if( useEllipsoids )
      {
         // prolate ellipsoids
         auto axisFactor = real_t(1.5);
         real_t axisFactor2 = std::sqrt(real_t(1)/axisFactor);
         real_t radius = diameter * real_t(0.5);
         pe::createEllipsoid( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( xParticle, yParticle, zParticle ), Vector3<real_t>(axisFactor*radius, axisFactor2*radius, axisFactor2*radius), peMaterial );
      }
      else
      {
         pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( xParticle, yParticle, zParticle ), diameter * real_t(0.5), peMaterial );
      }

   }

   syncCall();

   // carry out 100 simulations to resolve all overlaps
   for (auto pet = uint_t(1); pet <= uint_t(100); ++pet)
   {
      cr.timestep( real_t(1) );
      syncCall();

      // reset all velocities to zero
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->setLinearVel(Vector3<real_t>(real_t(0)));
            bodyIt->setAngularVel(Vector3<real_t>(real_t(0)));
         }
      }
   }


   const auto maxInitialPeSteps = (shortRun) ? uint_t(10) : uint_t(200000);
   const auto dt_PE_init = real_t(1);

   real_t gravityGeneration = real_t(0.1) * gravitationalAcceleration;
   cr.setGlobalLinearAcceleration(Vector3<real_t>(real_t(0), real_t(0), gravityGeneration));

   auto oldMinBodyPosition = real_t(0);
   real_t convergenceLimit = std::fabs(gravityGeneration);
   for (auto pet = uint_t(1); pet <= maxInitialPeSteps; ++pet)
   {
      cr.timestep( dt_PE_init );
      syncCall();

      real_t minBodyPosition = generationDomain.zMax();
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            minBodyPosition = std::min(bodyIt->getPosition()[2], minBodyPosition);
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace(minBodyPosition, mpi::MIN);
      }

      if( minBodyPosition > heightBorder ) break;

      if( pet % 500 == 0)
      {
         if( std::fabs(minBodyPosition - oldMinBodyPosition) / minBodyPosition  < convergenceLimit ) break;
         oldMinBodyPosition = minBodyPosition;
      }

      WALBERLA_ROOT_SECTION()
      {
         if( pet % 100 == 0)
         {
            WALBERLA_LOG_INFO("[" << pet << "] Min position of all bodies = " << minBodyPosition << " with goal height " << heightBorder);
         }
      }

   }

   // revert gravitational acceleration to 'real' direction
   cr.setGlobalLinearAcceleration(Vector3<real_t>(real_t(0), real_t(0), -gravityGeneration));

   // carry out a few time steps to relax the system towards the real condition
   const auto relaxationTimeSteps = uint_t(std::sqrt(real_t(2)/std::fabs(gravitationalAcceleration)));
   WALBERLA_LOG_INFO_ON_ROOT("Carrying out " << relaxationTimeSteps << " more time steps with correct gravity");
   for (auto pet = uint_t(1); pet <= relaxationTimeSteps; ++pet)
   {
      cr.timestep(dt_PE_init);
      syncCall();
   }

   WALBERLA_LOG_INFO_ON_ROOT("Sediment layer creation done!");

   // reset all velocities to zero
   Vector3<real_t> initialBodyVelocity(real_t(0));
   WALBERLA_LOG_INFO_ON_ROOT("Setting initial velocity " << initialBodyVelocity << " of all bodies");
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         bodyIt->setLinearVel(initialBodyVelocity);
         bodyIt->setAngularVel(Vector3<real_t>(real_t(0)));
      }
   }

   cr.setGlobalLinearAcceleration(Vector3<real_t>(real_t(0)));
}


//*******************************************************************************************************************
/*!\brief Simulation of settling particles inside a rectangular column filled with viscous fluid
 *
 * This application is used in the paper
 *  Rettinger, Ruede - "Dynamic Load Balancing Techniques for Particulate Flow Simulations", submitted to Computation
 * in Section 4 to apply the load estimator and to evaluate different load distribution strategies.
 *
 * It, however, features several different command line arguments that can be used to tweak the simulation.
 * The setup can be horizontally period, a box or a hopper geometry (configurable, as in the paper).
 * The size, resolution and used blocks for the domain partitioning can be changed.
 * It even features adaptive mesh refinement, with different refinement criteria:
 *  - particle based (always on, also for global bodies like bounding planes)
 *  - optionally: vorticity- or gradient-based (with lower and upper limits)
 * Since the paper, however, uses a uniform grid, many evaluation functionalities might not work properly for this case.
 * Initially, all particles are pushed upwards to obtain a dense packing at the top plane.
 *
 * Most importantly, the load balancing can be modified:
 *  - load estimation strategies:
 *    - pure LBM = number of cells per block = constant workload per block
 *    - pure PE = number of local particles + baseweight
 *    - coupling based load estimator = use fitted function from Sec. 3 of paper
 *  - load distribution strategies:
 *    - space-filling curves: Hilbert and Morton
 *    - ParMETIS (and several algorithms and parameters, also multiple constraints possible)
 *    - diffusive (and options)
 *  - load balancing (/refinement check ) frequency
 */
//*******************************************************************************************************************
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   ///////////////////
   // Customization //
   ///////////////////

   // simulation control
   bool shortRun = false;
   bool funcTest = false;
   bool fileIO = true;
   bool logging = false; // logging of physical components
   uint_t vtkWriteFreqDD = 0; //domain decomposition
   uint_t vtkWriteFreqBo = 0; //bodies
   uint_t vtkWriteFreqFl = 0; //fluid
   uint_t vtkWriteFreq = 0; //general
   std::string baseFolder = "vtk_out_AMRSedimentSettling"; // folder for vtk and file output

   // physical setup
   auto GalileoNumber = real_t(50);
   auto densityRatio = real_t(1.5);
   auto diameter = real_t(15);
   auto solidVolumeFraction = real_t(0.1);
   auto blockSize = uint_t(32);
   auto XBlocks = uint_t(12);
   auto YBlocks = uint_t(12);
   auto ZBlocks = uint_t(16);
   bool useBox = false;
   bool useHopper = false;
   bool useEllipsoids = false;
   auto hopperRelHeight = real_t(0.5); // for hopper setup
   auto hopperRelOpening = real_t(0.3); // for hopper setup

   auto timestepsOnFinestLevel = uint_t(80000);

   //numerical parameters
   bool averageForceTorqueOverTwoTimSteps = true;
   auto numberOfLevels = uint_t(1);
   auto refinementCheckFrequency = uint_t(100);
   auto numPeSubCycles = uint_t(10);

   // refinement criteria
   auto lowerFluidRefinementLimit = real_t(0);
   auto upperFluidRefinementLimit = std::numeric_limits<real_t>::infinity();
   bool useVorticityCriterion = false;
   bool useGradientCriterion = false;

   // load balancing
   std::string loadEvaluationStrategy = "LBM"; //LBM, PE, Fit
   std::string loadDistributionStrategy = "Hilbert"; //Morton, Hilbert, ParMetis, Diffusive

   auto parMetis_ipc2redist = real_t(1000);
   auto parMetisTolerance = real_t(-1);
   std::string parMetisAlgorithmString = "ADAPTIVE_REPART";

   auto diffusionFlowIterations = uint_t(15);
   auto diffusionMaxIterations = uint_t(20);


   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortRun" )                 == 0 ) { shortRun = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" )                 == 0 ) { funcTest = true; continue; }
      if( std::strcmp( argv[i], "--fileIO" )                   == 0 ) { fileIO = true; continue; }
      if( std::strcmp( argv[i], "--logging" )                  == 0 ) { logging = true; continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqDD" )           == 0 ) { vtkWriteFreqDD = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqBo" )           == 0 ) { vtkWriteFreqBo = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqFl" )           == 0 ) { vtkWriteFreqFl = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreq" )             == 0 ) { vtkWriteFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--baseFolder" )               == 0 ) { baseFolder = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--densityRatio" )             == 0 ) { densityRatio = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--Ga" )                       == 0 ) { GalileoNumber = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--diameter" )                 == 0 ) { diameter = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--blockSize" )                == 0 ) { blockSize = uint_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--XBlocks" )                  == 0 ) { XBlocks = uint_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--YBlocks" )                  == 0 ) { YBlocks = uint_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--ZBlocks" )                  == 0 ) { ZBlocks = uint_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--useBox" )                   == 0 ) { useBox = true; continue; }
      if( std::strcmp( argv[i], "--useHopper" )                == 0 ) { useHopper = true; continue; }
      if( std::strcmp( argv[i], "--hopperHeight" )             == 0 ) { hopperRelHeight = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--hopperOpening" )            == 0 ) { hopperRelOpening = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--timesteps" )                == 0 ) { timestepsOnFinestLevel = uint_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--noForceAveraging" )         == 0 ) { averageForceTorqueOverTwoTimSteps = false; continue; }
      if( std::strcmp( argv[i], "--numPeSubCycles" )           == 0 ) { numPeSubCycles = uint_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--numLevels" )                == 0 ) { numberOfLevels = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--refinementCheckFrequency" ) == 0 ) { refinementCheckFrequency = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--lowerLimit" )               == 0 ) { lowerFluidRefinementLimit = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--upperLimit" )               == 0 ) { upperFluidRefinementLimit = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--useVorticityCriterion" )    == 0 ) { useVorticityCriterion = true; continue; }
      if( std::strcmp( argv[i], "--useGradientCriterion" )     == 0 ) { useGradientCriterion = true; continue; }
      if( std::strcmp( argv[i], "--loadEvaluationStrategy" )   == 0 ) { loadEvaluationStrategy = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--loadDistributionStrategy" ) == 0 ) { loadDistributionStrategy = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--ipc2redist" )               == 0 ) { parMetis_ipc2redist = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--parMetisTolerance" )        == 0 ) { parMetisTolerance = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--parMetisAlgorithm" )        == 0 ) { parMetisAlgorithmString = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--diffusionFlowIterations" )  == 0 ) { diffusionFlowIterations = uint_c(std::atof(argv[++i])); continue; }
      if( std::strcmp( argv[i], "--diffusionMaxIterations" )   == 0 ) { diffusionMaxIterations = uint_c(std::atof(argv[++i])); continue; }
      if( std::strcmp( argv[i], "--useEllipsoids" )            == 0 ) { useEllipsoids = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   if( funcTest )
   {
      walberla::logging::Logging::instance()->setLogLevel(logging::Logging::LogLevel::WARNING);
   }

   if( fileIO || logging )
   {
      WALBERLA_ROOT_SECTION(){
         // create base directory if it does not yet exist
         filesystem::path tpath( baseFolder );
         if( !filesystem::exists( tpath ) )
            filesystem::create_directory( tpath );
      }
   }

   if( useVorticityCriterion && useGradientCriterion )
   {
      WALBERLA_ABORT("Use either vorticity or gradient criterion for refinement!");
   }

   if( loadEvaluationStrategy != "LBM" && loadEvaluationStrategy != "PE" && loadEvaluationStrategy != "Fit" && loadEvaluationStrategy != "FitMulti")
   {
      WALBERLA_ABORT("Invalid load evaluation strategy: " << loadEvaluationStrategy);
   }

   if( vtkWriteFreq != 0 )
   {
      vtkWriteFreqDD = vtkWriteFreq;
      vtkWriteFreqBo = vtkWriteFreq;
      vtkWriteFreqFl = vtkWriteFreq;
   }

   if( diameter > real_c(blockSize) )
   {
      WALBERLA_LOG_WARNING("PE Body Synchronization might not work since bodies are large compared to block size!");
   }

   if( useHopper )
   {
      WALBERLA_CHECK(hopperRelHeight >= real_t(0) && hopperRelHeight <= real_t(1), "Invalid relative hopper height of " << hopperRelHeight);
      WALBERLA_CHECK(hopperRelOpening >= real_t(0) && hopperRelOpening <= real_t(1), "Invalid relative hopper opening of " << hopperRelOpening);
   }


   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////

   const Vector3<uint_t> domainSize( XBlocks * blockSize, YBlocks * blockSize, ZBlocks * blockSize );
   const auto domainVolume = real_t(domainSize[0] * domainSize[1] * domainSize[2]);
   const real_t sphereVolume = math::pi / real_t(6) * diameter * diameter * diameter;
   const uint_t numberOfSediments = uint_c(std::ceil(solidVolumeFraction * domainVolume / sphereVolume));

   real_t expectedSedimentVolumeFraction = (useBox||useHopper) ? real_t(0.45) : real_t(0.52);
   const real_t expectedSedimentedVolume = real_t(1)/expectedSedimentVolumeFraction * real_c(numberOfSediments) * sphereVolume;
   const real_t expectedSedimentedHeight = std::max(diameter, expectedSedimentedVolume / real_c(domainSize[0] * domainSize[1]));

   const auto uRef = real_t(0.02);
   const real_t xRef = diameter;
   const real_t tRef = xRef / uRef;

   const real_t gravitationalAcceleration = uRef * uRef / ( (densityRatio-real_t(1)) * diameter );
   const real_t viscosity = uRef * diameter / GalileoNumber;
   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);
   const real_t tau = real_t(1) / omega;

   const auto loggingDisplayFrequency = uint_t(100);

   const auto dx = real_t(1);
   const real_t overlap = real_t( 1.5 ) * dx;

   if( useVorticityCriterion && floatIsEqual(lowerFluidRefinementLimit, real_t(0)) && std::isinf(upperFluidRefinementLimit) )
   {
      // use computed criterion instead of user input
      lowerFluidRefinementLimit = real_t(0.05) * uRef;
      upperFluidRefinementLimit = real_t(0.1) * uRef;
   }

   const uint_t finestLevel = numberOfLevels - uint_t(1);
   std::stringstream omega_msg;
   for( uint_t i = 0; i < numberOfLevels; ++i )
   {
      real_t omegaLvl = lbm::collision_model::levelDependentRelaxationParameter( i, omega, finestLevel );
      omega_msg << omegaLvl << " ( on level " << i << ", tau = " << real_t(1)/omega << " ), ";
   }

   const uint_t levelScalingFactor = ( uint_t(1) << finestLevel );
   const uint_t lbmTimeStepsPerTimeLoopIteration = levelScalingFactor;

   const uint_t timesteps = funcTest ? 1 : ( shortRun ? uint_t(100) : uint_t( timestepsOnFinestLevel / lbmTimeStepsPerTimeLoopIteration ) );

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - sediment diameter = " << diameter );
   WALBERLA_LOG_INFO_ON_ROOT(" - Galileo number = " << GalileoNumber );
   WALBERLA_LOG_INFO_ON_ROOT(" - number of sediments: " << numberOfSediments);
   WALBERLA_LOG_INFO_ON_ROOT(" - densityRatio = " << densityRatio );
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: relaxation time (tau) = " << tau << ", kin. visc = " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration );
   WALBERLA_LOG_INFO_ON_ROOT(" - reference values: x = " << xRef << ", t = " << tRef << ", vel = " << uRef);
   WALBERLA_LOG_INFO_ON_ROOT(" - omega: " << omega_msg.str());
   WALBERLA_LOG_INFO_ON_ROOT(" - number of levels: " << numberOfLevels);
   WALBERLA_LOG_INFO_ON_ROOT(" - number of pe sub cycles: " << numPeSubCycles);
   if( useVorticityCriterion )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using vorticity criterion with lower limit = " << lowerFluidRefinementLimit << " and upper limit = " << upperFluidRefinementLimit );
   }
   if( useGradientCriterion )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using gradient criterion with lower limit = " << lowerFluidRefinementLimit << " and upper limit = " << upperFluidRefinementLimit );
   }
   if( vtkWriteFreqDD > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files of domain decomposition to folder \"" << baseFolder << "\" with frequency " << vtkWriteFreqDD);
   }
   if( vtkWriteFreqBo > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files of bodies data to folder \"" << baseFolder << "\" with frequency " << vtkWriteFreqBo);
   }
   if( vtkWriteFreqFl > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files of fluid data to folder \"" << baseFolder << "\" with frequency " << vtkWriteFreqFl);
   }
   if( useEllipsoids )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using (prolate) ellipsoids as sediments");
   }
   if( useBox )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using box setup");
   }
   else if ( useHopper )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using hopper setup");
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using horizontally periodic domain");
   }

   if( refinementCheckFrequency == 0 && numberOfLevels != 1 )
   {
      // determine check frequency automatically based on maximum admissible velocity and block sizes
      auto uMax = real_t(0.1);
      refinementCheckFrequency = uint_c(( overlap + real_c(blockSize) - real_t(2) * real_t(FieldGhostLayers) * dx) / uMax) / lbmTimeStepsPerTimeLoopIteration;
   }
   WALBERLA_LOG_INFO_ON_ROOT(" - refinement / load balancing check frequency (coarse time steps): " << refinementCheckFrequency);
   WALBERLA_LOG_INFO_ON_ROOT(" - load evaluation strategy: " << loadEvaluationStrategy);
   WALBERLA_LOG_INFO_ON_ROOT(" - load distribution strategy: " << loadDistributionStrategy);

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3<uint_t> blockSizeInCells( blockSize );

   AABB simulationDomain( real_t(0), real_t(0), real_t(0), real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );
   AABB sedimentDomain( real_t(0), real_t(0), real_c(domainSize[2]) - expectedSedimentedHeight, real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );

   AABB initialRefinementDomain = sedimentDomain;
   if( useBox || useHopper )
   {
      // require finest levels also along bounding planes -> initially refine everywhere
      initialRefinementDomain = simulationDomain;
   }

   auto blocks = createBlockStructure( simulationDomain, blockSizeInCells, numberOfLevels, initialRefinementDomain, (useBox||useHopper), loadDistributionStrategy );

   //write initial domain decomposition to file
   if( vtkWriteFreqDD > 0 )
   {
      vtk::writeDomainDecomposition( blocks, "initial_domain_decomposition", baseFolder );
   }


   /////////////////
   // PE COUPLING //
   /////////////////

   // set up pe functionality
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();

   auto bodyStorageID  = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID          = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   BlockDataID fcdID   = (useEllipsoids) ? blocks->addBlockData( pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::GJKEPACollideFunctor>(), "FCD" )
                                         : blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   shared_ptr<WcTimingTree> timingTreePE = make_shared<WcTimingTree>();

   // set up collision response
   pe::cr::HCSITS cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, &(*timingTreePE) );
   cr.setMaxIterations(10);
   cr.setRelaxationModel( pe::cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling );

   // set up synchronization procedure
   std::function<void(void)> syncCall = std::bind( pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, &(*timingTreePE), overlap, false );

   // create pe bodies

   // add the sediments
   auto peMaterial = pe::createMaterial( "mat", densityRatio, real_t(1), real_t(0.25), real_t(0.25), real_t(0), real_t(200), real_t(100), real_t(100), real_t(100) );

   // create two planes at bottom and top of domain for a horizontally periodic box
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,1), Vector3<real_t>(0,0,0), peMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,simulationDomain.zMax()), peMaterial );
   if( useBox )
   {
      // add four more planes to obtain a closed box
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(1,0,0), Vector3<real_t>(0,0,0), peMaterial );
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(-1,0,0), Vector3<real_t>(simulationDomain.xMax(),0,0), peMaterial );
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,1,0), Vector3<real_t>(0,0,0), peMaterial );
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,-1,0), Vector3<real_t>(0,simulationDomain.yMax(),0), peMaterial );
   }
   else if ( useHopper )
   {
      // box bounding planes
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(1,0,0), Vector3<real_t>(0,0,0), peMaterial );
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(-1,0,0), Vector3<real_t>(simulationDomain.xMax(),0,0), peMaterial );
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,1,0), Vector3<real_t>(0,0,0), peMaterial );
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,-1,0), Vector3<real_t>(0,simulationDomain.yMax(),0), peMaterial );

      //hopper planes
      real_t xMax = simulationDomain.xMax();
      real_t yMax = simulationDomain.yMax();
      real_t zMax = simulationDomain.zMax();
      Vector3<real_t> p1(0,0,hopperRelHeight*zMax);
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(p1[2],0,hopperRelOpening*xMax-p1[0]), p1, peMaterial );
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,p1[2],hopperRelOpening*yMax-p1[0]), p1, peMaterial );

      Vector3<real_t> p2(xMax,yMax,hopperRelHeight*zMax);
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(-p2[2],0,-((real_t(1)-hopperRelOpening)*xMax-p2[0])), p2, peMaterial );
      pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,-p2[2],-((real_t(1)-hopperRelOpening)*yMax-p2[1])), p2, peMaterial );
   }

   AABB sedimentGenerationDomain( real_t(0), real_t(0), real_t(0.5)*real_c(domainSize[2]), real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );
   createSedimentLayer(numberOfSediments, sedimentGenerationDomain, diameter, sedimentDomain.zMin(), peMaterial, cr, syncCall, blocks, globalBodyStorage, bodyStorageID, gravitationalAcceleration, useEllipsoids, shortRun );

   // reset timer to not cover init stats
   timingTreePE.reset();

   // now we can use the information about the body positions to adapt the refinement

   ///////////////////////////
   // DYNAMIC REFINEMENT, 1 //
   ///////////////////////////

   auto & blockforest = blocks->getBlockForest();
   blockforest.recalculateBlockLevelsInRefresh( true );
   blockforest.alwaysRebalanceInRefresh( true ); //load balancing every time refresh is triggered
   blockforest.reevaluateMinTargetLevelsAfterForcedRefinement( false );
   blockforest.allowRefreshChangingDepth( false );
   blockforest.allowMultipleRefreshCycles( false ); // otherwise info collections are invalid

   {
      blockforest::CombinedMinTargetLevelDeterminationFunctions initialMinTargetLevelDeterminationFunctions;

      blockforest::AABBRefinementSelection aabbRefinementSelection;
      aabbRefinementSelection.addAABB(sedimentDomain,finestLevel );
      initialMinTargetLevelDeterminationFunctions.add( aabbRefinementSelection );

      // refinement along global bodies (bounding planes) to have consistent mapping (required for CLI always, or SimpleBB with non-AABB planes)
      real_t blockExtension = real_c(FieldGhostLayers);
      pe_coupling::amr::GlobalBodyPresenceLevelDetermination globalBodyPresenceRefinement( globalBodyStorage, finestLevel, blockExtension, pe_coupling::selectGlobalBodies );
      initialMinTargetLevelDeterminationFunctions.add(globalBodyPresenceRefinement);

      blockforest.setRefreshMinTargetLevelDeterminationFunction( initialMinTargetLevelDeterminationFunctions );

      for ( auto refreshCycle = uint_t(0); refreshCycle < finestLevel; ++refreshCycle)
      {

         WALBERLA_LOG_INFO_ON_ROOT("Refreshing blockforest...")

         // check refinement criteria and refine/coarsen if necessary
         uint_t stampBefore = blocks->getBlockForest().getModificationStamp();
         blocks->refresh();
         uint_t stampAfter = blocks->getBlockForest().getModificationStamp();

         if( stampBefore == stampAfter )
         {
            break;
         }

         WALBERLA_LOG_INFO_ON_ROOT("Recreating data structures..");

         // rebuild PE data structures
         pe::clearSynchronization( blockforest, bodyStorageID);

         syncCall();

         for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
         {
            auto * ccd = blockIt->getData< pe::ccd::ICCD >( ccdID );
            ccd->reloadBodies();
         }
      }
   }

   uint_t numberOfInitialFineBlocks = blockforest.getNumberOfBlocks(finestLevel);
   mpi::allReduceInplace(numberOfInitialFineBlocks, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("Total number of initial fine blocks in simulation: " << numberOfInitialFineBlocks);

   uint_t numberOfProcesses = uint_c(MPIManager::instance()->numProcesses());


   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega, lbm::collision_model::TRT::threeSixteenth, finestLevel ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         FieldGhostLayers, field::fzyx );
   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::fzyx, FieldGhostLayers );

   // add velocity field and utility
   BlockDataID velocityFieldID = field::addToStorage<VelocityField_T>( blocks, "velocity field", Vector3<real_t>(real_t(0)), field::fzyx, uint_t(2) );

   using VelocityFieldWriter_T = lbm::VelocityFieldWriter<PdfField_T, VelocityField_T>;
   BlockSweepWrapper< VelocityFieldWriter_T > velocityFieldWriter( blocks, VelocityFieldWriter_T( pdfFieldID, velocityFieldID ) );


   shared_ptr<blockforest::communication::NonUniformBufferedScheme<stencil::D3Q27> > velocityCommunicationScheme = make_shared<blockforest::communication::NonUniformBufferedScheme<stencil::D3Q27> >( blocks );
   velocityCommunicationScheme->addPackInfo( make_shared< field::refinement::PackInfo<VelocityField_T, stencil::D3Q27> >( velocityFieldID ) );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addBlockData( make_shared< MyBoundaryHandling >( blocks, flagFieldID, pdfFieldID, bodyFieldID ),
                                                          "boundary handling" );

   // map planes into the LBM simulation -> act as no-slip boundaries
   //pe_coupling::mapBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectGlobalBodies );
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectGlobalBodies );

   // map pe bodies into the LBM simulation
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectRegularBodies );


   // force averaging functionality
   shared_ptr<pe_coupling::BodiesForceTorqueContainer> bodiesFTContainer1 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
   std::function<void(void)> storeForceTorqueInCont1 = std::bind(&pe_coupling::BodiesForceTorqueContainer::store, bodiesFTContainer1);

   shared_ptr<pe_coupling::BodiesForceTorqueContainer> bodiesFTContainer2 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
   std::function<void(void)> setForceTorqueOnBodiesFromCont2 = std::bind(&pe_coupling::BodiesForceTorqueContainer::setOnBodies, bodiesFTContainer2);

   shared_ptr<pe_coupling::ForceTorqueOnBodiesScaler> forceScaler = make_shared<pe_coupling::ForceTorqueOnBodiesScaler>(blocks, bodyStorageID, real_t(0.5));
   std::function<void(void)> setForceScalingFactorToOne = std::bind(&pe_coupling::ForceTorqueOnBodiesScaler::resetScalingFactor,forceScaler,real_t(1));
   std::function<void(void)> setForceScalingFactorToHalf = std::bind(&pe_coupling::ForceTorqueOnBodiesScaler::resetScalingFactor,forceScaler,real_t(0.5));

   if( averageForceTorqueOverTwoTimSteps ) {
      bodiesFTContainer2->store();

      setForceScalingFactorToOne();
   }

   ///////////////////////////
   // DYNAMIC REFINEMENT, 2 //
   ///////////////////////////

   blockforest::CombinedMinTargetLevelDeterminationFunctions minTargetLevelDeterminationFunctions;

   // add refinement criterion based on particle presence
   shared_ptr<pe_coupling::InfoCollection> couplingInfoCollection = walberla::make_shared<pe_coupling::InfoCollection>();
   pe_coupling::amr::BodyPresenceLevelDetermination particlePresenceRefinement( couplingInfoCollection, finestLevel );

   minTargetLevelDeterminationFunctions.add( particlePresenceRefinement );

   // also add (possible) refinement criteria based on fluid quantities

   if( useVorticityCriterion )
   {
      // add refinement criterion based on vorticity magnitude
      field::FlagFieldEvaluationFilter<FlagField_T> flagFieldFilter( flagFieldID, Fluid_Flag );
      lbm::refinement::VorticityBasedLevelDetermination< field::FlagFieldEvaluationFilter<FlagField_T> > vorticityRefinement(
            velocityFieldID, flagFieldFilter, upperFluidRefinementLimit, lowerFluidRefinementLimit, finestLevel );

      minTargetLevelDeterminationFunctions.add( vorticityRefinement );
   }

   if( useGradientCriterion )
   {
      // add refinement criterion based on velocity gradient magnitude
      field::FlagFieldEvaluationFilter<FlagField_T> flagFieldFilter( flagFieldID, Fluid_Flag );
      VectorGradientRefinement< LatticeModel_T, field::FlagFieldEvaluationFilter<FlagField_T> > gradientRefinement(
            velocityFieldID, flagFieldFilter, upperFluidRefinementLimit, lowerFluidRefinementLimit, finestLevel );

      minTargetLevelDeterminationFunctions.add( gradientRefinement );
   }

   // refinement along global bodies (bounding planes) to have consistent mapping (required for CLI always, or SimpleBB with non-AABB planes)
   real_t blockExtension = real_c(FieldGhostLayers);
   pe_coupling::amr::GlobalBodyPresenceLevelDetermination globalBodyPresenceRefinement( globalBodyStorage, finestLevel, blockExtension, pe_coupling::selectGlobalBodies );
   minTargetLevelDeterminationFunctions.add(globalBodyPresenceRefinement);

   blockforest.setRefreshMinTargetLevelDeterminationFunction( minTargetLevelDeterminationFunctions );

   bool curveAllGather = true;
   bool balanceLevelwise = true;

   auto peBlockBaseWeight = real_t(1); //default value, might not be the best
   shared_ptr<blockforest::InfoCollection> peInfoCollection = walberla::make_shared<blockforest::InfoCollection>();

   if( loadDistributionStrategy == "Hilbert" || loadDistributionStrategy == "Morton")
   {
      if( loadDistributionStrategy == "Hilbert")
      {
         bool useHilbert = true;
         blockforest.setRefreshPhantomBlockMigrationPreparationFunction( blockforest::DynamicCurveBalance< blockforest::PODPhantomWeight<real_t> >( useHilbert, curveAllGather, balanceLevelwise ) );
      }
      else if (loadDistributionStrategy == "Morton" )
      {
         bool useHilbert = false;
         blockforest.setRefreshPhantomBlockMigrationPreparationFunction( blockforest::DynamicCurveBalance< blockforest::PODPhantomWeight<real_t> >( useHilbert, curveAllGather, balanceLevelwise ) );
      }

      blockforest.setRefreshPhantomBlockDataPackFunction(blockforest::PODPhantomWeightPackUnpack<real_t>());
      blockforest.setRefreshPhantomBlockDataUnpackFunction(blockforest::PODPhantomWeightPackUnpack<real_t>());

      if( loadEvaluationStrategy == "Fit" )
      {
         if( useEllipsoids )
         {
            pe_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, fittedTotalWeightEvaluationFunctionEllipsoids);
            blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
         } else{
            pe_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, fittedTotalWeightEvaluationFunctionSpheres);
            blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
         }
      }
      else if( loadEvaluationStrategy == "PE" )
      {
         blockforest::WeightAssignmentFunctor weightAssignmentFunctor(peInfoCollection, peBlockBaseWeight );
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else if( loadEvaluationStrategy == "LBM" )
      {
         pe_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, pe_coupling::amr::defaultWeightEvaluationFunction);
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else
      {
         WALBERLA_ABORT("Invalid load evaluation strategy: " << loadEvaluationStrategy);
      }

   }
   else if( loadDistributionStrategy == "ParMetis")
   {

#ifndef WALBERLA_BUILD_WITH_PARMETIS
      WALBERLA_ABORT( "You are trying to use ParMetis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_PARMETIS' to 'ON' in your CMake cache to build against an installed version of ParMetis!" );
#endif

      uint_t ncon = 1;
      if( loadEvaluationStrategy == "FitMulti")
      {
         ncon = 2;
      }

      blockforest::DynamicParMetis::Algorithm parMetisAlgorithm = blockforest::DynamicParMetis::stringToAlgorithm(parMetisAlgorithmString);
      blockforest::DynamicParMetis::WeightsToUse parMetisWeightsToUse = blockforest::DynamicParMetis::WeightsToUse::PARMETIS_BOTH_WEIGHTS;
      blockforest::DynamicParMetis::EdgeSource parMetisEdgeSource = blockforest::DynamicParMetis::EdgeSource::PARMETIS_EDGES_FROM_EDGE_WEIGHTS;

      blockforest::DynamicParMetis dynamicParMetis(parMetisAlgorithm, parMetisWeightsToUse, parMetisEdgeSource, ncon);
      dynamicParMetis.setipc2redist(parMetis_ipc2redist);

      real_t loadImbalanceTolerance = (parMetisTolerance < real_t(1)) ? std::max(real_t(1.05), real_t(1) + real_t(1) / ( real_c(numberOfInitialFineBlocks) / real_c(numberOfProcesses) ) ) : parMetisTolerance;
      std::vector<double> parMetisLoadImbalanceTolerance(ncon, double(loadImbalanceTolerance));
      dynamicParMetis.setImbalanceTolerance(parMetisLoadImbalanceTolerance[0], 0);

      WALBERLA_LOG_INFO_ON_ROOT(" - ParMetis configuration: ");
      WALBERLA_LOG_INFO_ON_ROOT("   - algorithm = " << dynamicParMetis.algorithmToString() );
      WALBERLA_LOG_INFO_ON_ROOT("   - weights to use = " << dynamicParMetis.weightsToUseToString() );
      WALBERLA_LOG_INFO_ON_ROOT("   - edge source = " << dynamicParMetis.edgeSourceToString() );
      WALBERLA_LOG_INFO_ON_ROOT("   - ncon = " << ncon );
      WALBERLA_LOG_INFO_ON_ROOT("   - ipc2redist parameter = " << dynamicParMetis.getipc2redist() );

      blockforest.setRefreshPhantomBlockDataPackFunction(blockforest::DynamicParMetisBlockInfoPackUnpack());
      blockforest.setRefreshPhantomBlockDataUnpackFunction(blockforest::DynamicParMetisBlockInfoPackUnpack());

      if( loadEvaluationStrategy == "Fit" )
      {
         WALBERLA_LOG_INFO_ON_ROOT("   - load imbalance tolerance = <" << parMetisLoadImbalanceTolerance[0] << ">" );
         if( useEllipsoids )
         {
            pe_coupling::amr::MetisAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, fittedTotalWeightEvaluationFunctionEllipsoids);
            blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
         } else{
            pe_coupling::amr::MetisAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, fittedTotalWeightEvaluationFunctionSpheres);
            blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
         }
      }
      else if( loadEvaluationStrategy == "FitMulti" )
      {
         double imbalanceTolerancePE = 10.;
         parMetisLoadImbalanceTolerance[1] = std::min(imbalanceTolerancePE, static_cast<double>(MPIManager::instance()->numProcesses()));
         WALBERLA_LOG_INFO_ON_ROOT("   - load imbalance tolerances = <" << parMetisLoadImbalanceTolerance[0] << ", " << parMetisLoadImbalanceTolerance[1] << ">" );
         dynamicParMetis.setImbalanceTolerance(parMetisLoadImbalanceTolerance[1], 1);

         if( useEllipsoids )
         {
            std::vector< std::function<real_t(const pe_coupling::BlockInfo&)> > weightEvaluationFunctions(ncon);
            weightEvaluationFunctions[0] = fittedLBMWeightEvaluationFunctionEllipsoids;
            weightEvaluationFunctions[1] = fittedPEWeightEvaluationFunctionEllipsoids;
            pe_coupling::amr::MetisAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, weightEvaluationFunctions);
            blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
         } else{
            std::vector< std::function<real_t(const pe_coupling::BlockInfo&)> > weightEvaluationFunctions(ncon);
            weightEvaluationFunctions[0] = fittedLBMWeightEvaluationFunctionSpheres;
            weightEvaluationFunctions[1] = fittedPEWeightEvaluationFunctionSpheres;
            pe_coupling::amr::MetisAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, weightEvaluationFunctions);
            blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
         }
      }
      else if( loadEvaluationStrategy == "PE" )
      {
         blockforest::MetisAssignmentFunctor weightAssignmentFunctor(peInfoCollection, peBlockBaseWeight );
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else if( loadEvaluationStrategy == "LBM" )
      {
         pe_coupling::amr::MetisAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, pe_coupling::amr::defaultWeightEvaluationFunction);
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else
      {
         WALBERLA_ABORT("Invalid load evaluation strategy: " << loadEvaluationStrategy);
      }

      blockforest.setRefreshPhantomBlockMigrationPreparationFunction( dynamicParMetis );

   }
   else if( loadDistributionStrategy == "Diffusive")
   {
      using DB_T = blockforest::DynamicDiffusionBalance< blockforest::PODPhantomWeight<real_t> >;
      DB_T dynamicDiffusion(diffusionMaxIterations, diffusionFlowIterations );
      dynamicDiffusion.setMode(DB_T::Mode::DIFFUSION_PUSH);

      WALBERLA_LOG_INFO_ON_ROOT(" - Dynamic diffusion configuration: ");
      WALBERLA_LOG_INFO_ON_ROOT("   - max iterations = " << dynamicDiffusion.getMaxIterations() );
      WALBERLA_LOG_INFO_ON_ROOT("   - flow iterations = " << dynamicDiffusion.getFlowIterations());

      blockforest.setRefreshPhantomBlockDataPackFunction(blockforest::PODPhantomWeightPackUnpack<real_t>());
      blockforest.setRefreshPhantomBlockDataUnpackFunction(blockforest::PODPhantomWeightPackUnpack<real_t>());
      blockforest.setRefreshPhantomBlockMigrationPreparationFunction( dynamicDiffusion );

      if( loadEvaluationStrategy == "Fit" )
      {
         if( useEllipsoids )
         {
            pe_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, fittedTotalWeightEvaluationFunctionEllipsoids);
            blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
         } else{
            pe_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, fittedTotalWeightEvaluationFunctionSpheres);
            blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
         }
      }
      else if( loadEvaluationStrategy == "PE" )
      {
         blockforest::WeightAssignmentFunctor weightAssignmentFunctor(peInfoCollection, peBlockBaseWeight );
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else if( loadEvaluationStrategy == "LBM" )
      {
         pe_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, pe_coupling::amr::defaultWeightEvaluationFunction);
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else
      {
         WALBERLA_ABORT("Invalid load evaluation strategy: " << loadEvaluationStrategy);
      }

   } else
   {
      WALBERLA_ABORT("Load distribution strategy \"" << loadDistributionStrategy << "\t not implemented! - Aborting" );
   }


   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   auto timeloop = make_shared<SweepTimeloop>( blocks->getBlockStorage(), timesteps );

   if( vtkWriteFreqBo != uint_t(0) ) {

      // pe bodies
      if (useEllipsoids) {
         auto bodyVtkOutput = make_shared<pe::EllipsoidVtkOutput>(bodyStorageID, blocks->getBlockStorage());
         auto bodyVTK = vtk::createVTKOutput_PointData(bodyVtkOutput, "bodies", vtkWriteFreqBo, baseFolder);
         timeloop->addFuncBeforeTimeStep(vtk::writeFiles(bodyVTK), "VTK (sediment data)");

      } else {
         auto bodyVtkOutput = make_shared<pe::SphereVtkOutput>(bodyStorageID, blocks->getBlockStorage());
         auto bodyVTK = vtk::createVTKOutput_PointData(bodyVtkOutput, "bodies", vtkWriteFreqBo, baseFolder);
         timeloop->addFuncBeforeTimeStep(vtk::writeFiles(bodyVTK), "VTK (sediment data)");
      }
   }

   if( vtkWriteFreqFl != uint_t(0) ) {

      // pdf field
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData(blocks, "fluid_field", vtkWriteFreqFl, 0, false, baseFolder);

      field::FlagFieldCellFilter<FlagField_T> fluidFilter(flagFieldID);
      fluidFilter.addFlag(Fluid_Flag);
      pdfFieldVTK->addCellInclusionFilter(fluidFilter);

      pdfFieldVTK->addCellDataWriter(
            make_shared<lbm::VelocityVTKWriter<LatticeModel_T, float> >(pdfFieldID, "VelocityFromPDF"));
      pdfFieldVTK->addCellDataWriter(
            make_shared<lbm::DensityVTKWriter<LatticeModel_T, float> >(pdfFieldID, "DensityFromPDF"));

      timeloop->addFuncBeforeTimeStep(vtk::writeFiles(pdfFieldVTK), "VTK (fluid field data)");
   }

   if( vtkWriteFreqDD != uint_t(0) ) {
      auto domainDecompVTK = vtk::createVTKOutput_DomainDecomposition(blocks, "domain_decomposition", vtkWriteFreqDD, baseFolder );
      timeloop->addFuncBeforeTimeStep( vtk::writeFiles(domainDecompVTK), "VTK (domain decomposition)");
   }

   shared_ptr<WcTimingPool> timeloopTiming = make_shared<WcTimingPool>();
   shared_ptr<WcTimingPool> timeloopRefinementTiming = make_shared<WcTimingPool>();
   shared_ptr<WcTimingPool> timeloopRefinementTimingLevelwise = make_shared<WcTimingPool>();


   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   auto refinementTimestep = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, sweep, pdfFieldID, boundaryHandlingID );

   refinementTimestep->enableTiming( timeloopRefinementTiming, timeloopRefinementTimingLevelwise );

   // Averaging the force/torque over two time steps is said to damp oscillations of the interaction force/torque.
   // See Ladd - " Numerical simulations of particulate suspensions via a discretized Boltzmann equation. Part 1. Theoretical foundation", 1994, p. 302
   if( averageForceTorqueOverTwoTimSteps ) {

      // store force/torque from hydrodynamic interactions in container1
      refinementTimestep->addPostStreamVoidFunction(lbm::refinement::FunctorWrapper(storeForceTorqueInCont1), "Force Storing", finestLevel);

      // set force/torque from previous time step (in container2)
      refinementTimestep->addPostStreamVoidFunction(lbm::refinement::FunctorWrapper(setForceTorqueOnBodiesFromCont2), "Force setting", finestLevel);

      // average the force/torque by scaling it with factor 1/2 (except in first timestep and directly after refinement, there it is 1)
      refinementTimestep->addPostStreamVoidFunction(lbm::refinement::FunctorWrapper(SharedFunctor<pe_coupling::ForceTorqueOnBodiesScaler>(forceScaler)), "Force averaging", finestLevel);
      refinementTimestep->addPostStreamVoidFunction(lbm::refinement::FunctorWrapper(setForceScalingFactorToHalf), "Force scaling adjustment", finestLevel);

      // swap containers
      refinementTimestep->addPostStreamVoidFunction(lbm::refinement::FunctorWrapper(pe_coupling::BodyContainerSwapper(bodiesFTContainer1, bodiesFTContainer2)), "Swap FT container", finestLevel);

   }

   Vector3<real_t> gravitationalForce( real_t(0), real_t(0), -(densityRatio - real_t(1)) * gravitationalAcceleration * sphereVolume );
   refinementTimestep->addPostStreamVoidFunction(lbm::refinement::FunctorWrapper(pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, gravitationalForce )), "Gravitational force", finestLevel );

   // add pe timesteps
   refinementTimestep->addPostStreamVoidFunction(lbm::refinement::FunctorWrapper(pe_coupling::TimeStep( blocks, bodyStorageID, cr, syncCall, real_t(1), numPeSubCycles)),
                                                  "pe Time Step", finestLevel );

   // add sweep for updating the pe body mapping into the LBM simulation
   refinementTimestep->addPostStreamVoidFunction(lbm::refinement::SweepAsFunctorWrapper( pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID,  MO_Flag, FormerMO_Flag, pe_coupling::selectRegularBodies ), blocks ),
                                                 "Body Mapping", finestLevel );

   // add sweep for restoring PDFs in cells previously occupied by pe bodies
   using Reconstructor_T = pe_coupling::EquilibriumReconstructor<LatticeModel_T, BoundaryHandling_T>;
   Reconstructor_T reconstructor( blocks, boundaryHandlingID, bodyFieldID );
   refinementTimestep->addPostStreamVoidFunction(lbm::refinement::SweepAsFunctorWrapper( pe_coupling::PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructor_T > ( blocks, pdfFieldID,
                                                 boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, reconstructor, FormerMO_Flag, Fluid_Flag ), blocks ),
                                                 "PDF Restore", finestLevel );


   // add LBM sweep with refinement
   timeloop->addFuncBeforeTimeStep( makeSharedFunctor( refinementTimestep ), "LBM refinement time step" );

   std::string loggingFileName( baseFolder + "/Logging_Ga");
   loggingFileName += std::to_string(uint_c(GalileoNumber));
   loggingFileName += "_lvl";
   loggingFileName += std::to_string(numberOfLevels);
   loggingFileName += ".txt";
   if( logging  )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing logging output to file \"" << loggingFileName << "\"");
   }
   shared_ptr< PropertyLogger > logger = walberla::make_shared< PropertyLogger >( timeloop, blocks, bodyStorageID,
                                                                                  loggingFileName, fileIO );
   if(logging)
   {
      timeloop->addFuncAfterTimeStep( SharedFunctor< PropertyLogger >( logger ), "Property logger" );
   }


   timeloop->addFuncAfterTimeStep( RemainingTimeLogger( timeloop->getNrOfTimeSteps() ), "Remaining Time Logger" );


   // add top level timing pool output
   timeloop->addFuncAfterTimeStep( TimingPoolLogger( timeloopTiming, timeloop, loggingDisplayFrequency ), "Regular Timing Logger" );

   // add regular refinement timing pool output
   timeloop->addFuncAfterTimeStep( TimingPoolLogger( timeloopRefinementTiming, timeloop, loggingDisplayFrequency ), "Refinement Timing Logger" );

   // add level wise timing pool output
   //if( numberOfLevels != uint_t(1))
   timeloop->addFuncAfterTimeStep( TimingPoolLogger( timeloopRefinementTimingLevelwise, timeloop, loggingDisplayFrequency ), "Refinement Levelwise Timing Logger" );

   // add PE timing tree output
   timeloop->addFuncAfterTimeStep( TimingTreeLogger( timingTreePE, timeloop, loggingDisplayFrequency ), "PE Timing Tree Timing Logger" );


   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   uint_t loadEvaluationFrequency = refinementCheckFrequency;
   TimingEvaluator timingEvaluator( timeloopRefinementTimingLevelwise, timingTreePE, numberOfLevels );

   // file for simulation infos
   std::string infoFileName( baseFolder + "/simulation_info.txt");
   WALBERLA_ROOT_SECTION()
   {
      std::ofstream file;
      file.open( infoFileName.c_str(), std::fstream::out | std::fstream::trunc );
      file << "#i\t t\t tSim\t tLB\t numProcs\t levelwise blocks (min/max/sum)\n";
      file.close();
   }

   // process local timing measurements and predicted loads
   std::string processLocalFiles(baseFolder + "/processLocalFiles");
   WALBERLA_ROOT_SECTION()
   {
      filesystem::path tpath( processLocalFiles );
      if( !filesystem::exists( tpath ) )
         filesystem::create_directory( tpath );
   }
   std::string measurementFileProcessName(processLocalFiles + "/measurements_" + std::to_string(MPIManager::instance()->rank()) + ".txt");
   {
      std::ofstream file;
      file.open( measurementFileProcessName.c_str(), std::fstream::out | std::fstream::trunc );
      file << "#i\t t\t mTotSim\t mLB\t mLBM\t mBH\t mCoup1\t mCoup2\t mRB\t cLBM\t cRB\t numBlocks\n";
      file.close();
   }

   std::string predictionFileProcessName(processLocalFiles + "/predictions_" + std::to_string(MPIManager::instance()->rank()) + ".txt");
   {
      std::ofstream file;
      file.open( predictionFileProcessName.c_str(), std::fstream::out | std::fstream::trunc );
      file << "#i\t t\t wlLBM\t wlBH\t wlCoup1\t wlCoup2\t wlRB\t edgecut\t numBlocks\n";
      file.close();
   }

   std::vector<std::string> LBMTimer;
   LBMTimer.emplace_back("collide");
   LBMTimer.emplace_back("stream");
   LBMTimer.emplace_back("stream & collide");

   std::vector<std::string> bhTimer;
   bhTimer.emplace_back("boundary handling");

   std::vector<std::string> couplingTimer1;
   couplingTimer1.emplace_back("Body Mapping");
   std::vector<std::string> couplingTimer2;
   couplingTimer2.emplace_back("PDF Restore");

   std::vector<std::string> peTimer;
   peTimer.emplace_back("Simulation Step.Collision Detection");
   peTimer.emplace_back("Simulation Step.Collision Response Integration");
   peTimer.emplace_back("Simulation Step.Collision Response Resolution.Collision Response Solving");

   std::vector<std::string> LBMCommTimer;
   LBMCommTimer.emplace_back("communication equal level [pack & send]");
   LBMCommTimer.emplace_back("communication equal level [wait & unpack]");

   std::vector<std::string> peCommTimer;
   //Adapt if using different collision response (like DEM!)
   peCommTimer.emplace_back("Simulation Step.Collision Response Resolution.Velocity Sync");
   peCommTimer.emplace_back("Sync");


   real_t terminationPosition = expectedSedimentedHeight;
   real_t terminationVelocity = real_t(0.05) * uRef;

   auto oldmTotSim = real_t(0);
   auto oldmLB = real_t(0);

   auto measurementFileCounter = uint_t(0);
   auto predictionFileCounter = uint_t(0);

   std::string loadEvaluationStep("load evaluation");

   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {

      // evaluate measurements (note: reflect simulation behavior BEFORE the evaluation)
      if( loadEvaluationFrequency > 0 && i % loadEvaluationFrequency == 0 && i > 0 && fileIO)
      {

         (*timeloopTiming)[loadEvaluationStep].start();

         // write process local timing measurements to files (per process, per load balancing step)
         {

            // evaluate all required timers
            uint_t evalLevel = finestLevel;

            real_t mTotSim = ( (*timeloopTiming).timerExists("LBM refinement time step") ) ? real_c((*timeloopTiming)["LBM refinement time step"].total()) : real_t(0);

            real_t mLB = ( (*timeloopTiming).timerExists("refinement checking") ) ? real_c((*timeloopTiming)["refinement checking"].total()) : real_t(0);

            real_t mLBM = timingEvaluator.getTimings(LBMTimer, evalLevel);
            real_t mBH  = timingEvaluator.getTimings(bhTimer, evalLevel);
            real_t mCoup1 = timingEvaluator.getTimings(couplingTimer1, evalLevel);
            real_t mCoup2 = timingEvaluator.getTimings(couplingTimer2, evalLevel);
            real_t mPE = timingEvaluator.getTimings(peTimer, evalLevel);

            real_t cLBM = timingEvaluator.getTimings(LBMCommTimer, evalLevel);
            real_t cRB = timingEvaluator.getTimings(peCommTimer, evalLevel);

            auto & forest = blocks->getBlockForest();
            uint_t numBlocks = forest.getNumberOfBlocks(finestLevel);

            // write to process local file
            std::ofstream file;
            file.open( measurementFileProcessName.c_str(), std::ofstream::app  );
            file << measurementFileCounter << "\t " << real_c(i) / tRef << "\t"
                 << mTotSim - oldmTotSim << "\t" << mLB - oldmLB << "\t" << mLBM << "\t" << mBH << "\t" << mCoup1 << "\t"
                 << mCoup2 << "\t" << mPE << "\t" << cLBM << "\t" << cRB << "\t" << numBlocks << "\n";
            file.close();

            oldmTotSim = mTotSim;
            oldmLB = mLB;
            measurementFileCounter++;

            // reset timers to have measurement from evaluation to evaluation point
            timeloopRefinementTimingLevelwise->clear();
            timingTreePE.reset();

         }

         // evaluate general simulation infos (on root)
         {
            real_t totalTimeToCurrentTimestep;
            real_t totalLBTimeToCurrentTimestep;
            evaluateTotalSimulationTimePassed(*timeloopTiming, totalTimeToCurrentTimestep, totalLBTimeToCurrentTimestep);
            std::vector<math::DistributedSample> numberOfBlocksPerLevel(numberOfLevels);

            auto & forest = blocks->getBlockForest();
            for( uint_t lvl = 0; lvl < numberOfLevels; ++lvl)
            {
               uint_t numBlocks = forest.getNumberOfBlocks(lvl);
               numberOfBlocksPerLevel[lvl].castToRealAndInsert(numBlocks);
            }

            for( uint_t lvl = 0; lvl < numberOfLevels; ++lvl)
            {
               numberOfBlocksPerLevel[lvl].mpiGatherRoot();
            }

            WALBERLA_ROOT_SECTION()
            {
               std::ofstream file;
               file.open( infoFileName.c_str(), std::ofstream::app  );
               file << i << "\t " << real_c(i) / tRef << "\t"
                    << totalTimeToCurrentTimestep << "\t " << totalLBTimeToCurrentTimestep << "\t " << numberOfProcesses << "\t ";

               for( uint_t lvl = 0; lvl < numberOfLevels; ++lvl)
               {
                  file << uint_c(numberOfBlocksPerLevel[numberOfLevels-1-lvl].min()) << "\t ";
                  file << uint_c(numberOfBlocksPerLevel[numberOfLevels-1-lvl].max()) << "\t ";
                  file << uint_c(numberOfBlocksPerLevel[numberOfLevels-1-lvl].sum()) << "\t ";
               }
               file << "\n";

               file.close();
            }
         }

         (*timeloopTiming)[loadEvaluationStep].end();

      }


      if( refinementCheckFrequency != 0 && i % refinementCheckFrequency == 0)
      {

         WALBERLA_LOG_INFO_ON_ROOT("Checking for refinement and load balancing...")

         std::string refinementCheckStep("refinement checking");
         (*timeloopTiming)[refinementCheckStep].start();

         if( loadEvaluationStrategy != "LBM" ) {

            // first evaluate all data that is required for the refinement checks

            // update info collections for the particle presence based check and the load balancing:
            auto &forest = blocks->getBlockForest();
            pe_coupling::createWithNeighborhood<BoundaryHandling_T>(forest, boundaryHandlingID, bodyStorageID, ccdID,
                                                                    fcdID, numPeSubCycles, *couplingInfoCollection);
            pe::createWithNeighborhoodLocalShadow(forest, bodyStorageID, *peInfoCollection);

            // for the fluid property based check:
            if (useVorticityCriterion || useGradientCriterion) {
               velocityFieldWriter();
               (*velocityCommunicationScheme)();
            }

            WALBERLA_LOG_INFO_ON_ROOT("Refreshing blockforest...")

            // check refinement criteria and refine/coarsen if necessary
            uint_t stampBefore = blocks->getBlockForest().getModificationStamp();
            blocks->refresh();
            uint_t stampAfter = blocks->getBlockForest().getModificationStamp();

            bool recreatingNecessary = false;

            if (stampBefore != stampAfter) {
               recreatingNecessary = true;
            }

            bool reducedRecreationFlag = mpi::allReduce(recreatingNecessary, mpi::LOGICAL_OR);

            if (reducedRecreationFlag != recreatingNecessary) {
               WALBERLA_LOG_INFO("Reduced recreation flag different from individual one");
            }

            recreatingNecessary = reducedRecreationFlag;

            if (recreatingNecessary) {

               WALBERLA_LOG_INFO_ON_ROOT("Recreating data structures..");

               // rebuild PE data structures
               pe::clearSynchronization(blockforest, bodyStorageID);

               syncCall();

               for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
                  auto * ccd = blockIt->getData<pe::ccd::ICCD>(ccdID);
                  ccd->reloadBodies();
               }

               clearBoundaryHandling(forest, boundaryHandlingID);
               clearBodyField(forest, bodyFieldID);

               if (averageForceTorqueOverTwoTimSteps) {

                  // clear containers from old values
                  bodiesFTContainer1->clear();
                  bodiesFTContainer2->clear();

                  // initialize FT container on all blocks anew, i.e. with the currently acting force/torque, which is zero after the refinement step
                  bodiesFTContainer2->store();

                  // set force scaling factor to one after refinement since force history is not present on blocks after refinement
                  // thus the usual averaging of 1/2 (over two time steps) can not be carried out, i.e. it would lead to 1/2 of the acting force
                  // the scaling factor is thus adapted for the next timestep to 1, and then changed back to 1/2 (in the timeloop)
                  setForceScalingFactorToOne();
               }

               recreateBoundaryHandling(forest, boundaryHandlingID);

               // re-set the no-slip flags along the walls
               pe_coupling::mapMovingBodies<BoundaryHandling_T>(*blocks, boundaryHandlingID, bodyStorageID,
                                                                *globalBodyStorage, bodyFieldID, MO_Flag,
                                                                pe_coupling::selectGlobalBodies);

               // re-map the body into the domain (initializing the bodyField as well)
               pe_coupling::mapMovingBodies<BoundaryHandling_T>(*blocks, boundaryHandlingID, bodyStorageID,
                                                                *globalBodyStorage, bodyFieldID, MO_Flag,
                                                                pe_coupling::selectRegularBodies);
            }

         }

         (*timeloopTiming)[refinementCheckStep].end();
      }

      // evaluate predictions (note: reflect the predictions for all upcoming simulations, thus the corresponding measurements have to be taken afterwards)
      if( loadEvaluationFrequency > 0 && i % loadEvaluationFrequency == 0 && fileIO)
      {

         (*timeloopTiming)[loadEvaluationStep].start();

         // write process local load predictions to files (per process, per load balancing step)
         {

            auto wlLBM = real_t(0);
            auto wlBH = real_t(0);
            auto wlCoup1 = real_t(0);
            auto wlCoup2 = real_t(0);
            auto wlRB = real_t(0);

            auto & forest = blocks->getBlockForest();
            pe_coupling::createWithNeighborhood<BoundaryHandling_T>(forest, boundaryHandlingID, bodyStorageID, ccdID, fcdID, numPeSubCycles, *couplingInfoCollection);

            for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt ) {
               auto * block = static_cast<blockforest::Block *> (&(*blockIt));
               const auto &blockID = block->getId();
               auto infoIt = couplingInfoCollection->find(blockID);
               auto blockInfo = infoIt->second;

               if( useEllipsoids )
               {
                  WALBERLA_ABORT("Not yet implemented!");
               }
               else
               {
                  wlLBM   += fittedLBMWeightEvaluationFunctionSpheres(blockInfo);
                  wlBH    += fittedBHWeightEvaluationFunctionSpheres(blockInfo);
                  wlCoup1 += fittedCoupling1WeightEvaluationFunctionSpheres(blockInfo);
                  wlCoup2 += fittedCoupling2WeightEvaluationFunctionSpheres(blockInfo);
                  wlRB    += fittedPEWeightEvaluationFunctionSpheres(blockInfo);
               }

            }

            // note: we count the edge weight doubled here in total (to and from the other process). ParMetis only counts one direction.
            uint_t edgecut = evaluateEdgeCut(forest);

            uint_t numBlocks = forest.getNumberOfBlocks(finestLevel);

            std::ofstream file;
            file.open( predictionFileProcessName.c_str(), std::ofstream::app  );
            file << predictionFileCounter << "\t " << real_c(i) / tRef << "\t"
                 << wlLBM << "\t" << wlBH << "\t" << wlCoup1 << "\t" << wlCoup2 << "\t" << wlRB << "\t"
                 << edgecut << "\t" << numBlocks << "\n";
            file.close();

            predictionFileCounter++;;
         }

         (*timeloopTiming)[loadEvaluationStep].end();

      }

      // perform a single simulation step
      timeloop->singleStep( *timeloopTiming );


      if( logging )
      {
         real_t curMeanPos = logger->getMeanPosition();
         real_t curMeanVel = logger->getMeanVelocity();
         if( curMeanPos <= terminationPosition && std::fabs(curMeanVel) < terminationVelocity )
         {
            WALBERLA_LOG_INFO_ON_ROOT("Sediments passed terminal mean position " << terminationPosition << " and reached velocity " << curMeanVel << " - terminating simulation!");
            break;
         }
      }

   }

   (*timeloopTiming).logResultOnRoot();


   return EXIT_SUCCESS;
}

} // namespace amr_sediment_settling

int main( int argc, char **argv ){
   amr_sediment_settling::main(argc, argv);
}
