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
//! \file AMRSettlingSphere.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/loadbalancing/all.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"
#include "domain_decomposition/BlockSweepWrapper.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/VelocityFieldWriter.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/refinement/all.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "pe/amr/InfoCollection.h"
#include "pe/basic.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"
#include "pe/cr/ICR.h"
#include "pe/Types.h"
#include "pe/synchronization/ClearSynchronization.h"

#include "pe_coupling/amr/all.h"
#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

namespace amr_settling_sphere
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
                                        const real_t inv_c_s_sqr = real_t(3);
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

                                        real_t norm( real_t(0) );
                                        //compute maximums norm of 3x3 matrix
                                        for( uint_t i = uint_t(0); i < uint_t(3*3); ++i )
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
                                                                 uint_t numberOfLevels, real_t diameter, Vector3<real_t> spherePosition,
                                                                 bool useStaticRefinement,
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

   AABB refinementBox;
   if(useStaticRefinement)
   {
      refinementBox = AABB( std::floor(spherePosition[0] - real_t(0.5) * diameter),
                            std::floor(spherePosition[1] - real_t(0.5) * diameter),
                            domainAABB.zMin(),
                            std::ceil( spherePosition[0] + real_t(0.5) * diameter),
                            std::ceil( spherePosition[1] + real_t(0.5) * diameter),
                            domainAABB.zMax() );
   }else{
      refinementBox = AABB( std::floor(spherePosition[0] - real_t(0.5) * diameter),
                            std::floor(spherePosition[1] - real_t(0.5) * diameter),
                            std::floor(spherePosition[2] - real_t(0.5) * diameter),
                            std::ceil( spherePosition[0] + real_t(0.5) * diameter),
                            std::ceil( spherePosition[1] + real_t(0.5) * diameter),
                            std::ceil( spherePosition[2] + real_t(0.5) * diameter) );
   }

   WALBERLA_LOG_INFO_ON_ROOT(" - refinement box: " << refinementBox);

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

   FlagField_T * flagField       = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField        = block->getData< PdfField_T > ( pdfFieldID_ );
   BodyField_T * bodyField       = block->getData< BodyField_T >( bodyFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   auto blocksPtr = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocksPtr );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "moving obstacle boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                                           MO_T( "MO", MO_Flag, pdfField, flagField, bodyField, fluid, *blocksPtr, *block ) );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}


//*******************************************************************************************************************


//*******************************************************************************************************************
/*!\brief Evaluating the position and velocity of the sphere
 *
 */
//*******************************************************************************************************************
class SpherePropertyLogger
{
public:
   SpherePropertyLogger( const shared_ptr<SweepTimeloop> & timeloop, const shared_ptr< StructuredBlockStorage > & blocks,
                         const BlockDataID & bodyStorageID, const std::string & fileName, bool fileIO,
                         real_t xRef, real_t tRef, uint_t lbmTimeStepsPerTimeLoopIteration,
                         real_t diameter, real_t viscosity) :
         timeloop_( timeloop ), blocks_( blocks ), bodyStorageID_( bodyStorageID ), fileName_( fileName ), fileIO_(fileIO),
         xRef_( xRef ), tRef_( tRef ), lbmTimeStepsPerTimeLoopIteration_( lbmTimeStepsPerTimeLoopIteration ),
         diameter_( diameter ), viscosity_( viscosity ),
         position_( real_t(0) ), maxVelocity_( real_t(0) )
   {
      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str() );
            file << "#\t t\t posX\t posY\t posZ\t velX\t velY\t velZ\t posX*\t posY*\t posZ*\t velX*\t velY*\t velZ*\t Re\n";
            file.close();
         }
      }
   }

   void operator()()
   {
      const uint_t timestep (timeloop_->getCurrentTimeStep() * lbmTimeStepsPerTimeLoopIteration_ );

      Vector3<real_t> pos(real_t(0));
      Vector3<real_t> transVel(real_t(0));

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            pos = bodyIt->getPosition();
            transVel = bodyIt->getLinearVel();
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( pos[0], mpi::SUM );
         mpi::allReduceInplace( pos[1], mpi::SUM );
         mpi::allReduceInplace( pos[2], mpi::SUM );

         mpi::allReduceInplace( transVel[0], mpi::SUM );
         mpi::allReduceInplace( transVel[1], mpi::SUM );
         mpi::allReduceInplace( transVel[2], mpi::SUM );
      }

      position_ = pos[2];
      maxVelocity_ = std::max(maxVelocity_, -transVel[2]);

      if( fileIO_ )
         writeToFile( timestep, pos, transVel);
   }

   real_t getPosition() const
   {
      return position_;
   }

   real_t getMaxVelocity() const
   {
      return maxVelocity_;
   }

private:
   void writeToFile( uint_t timestep, const Vector3<real_t> & position, const Vector3<real_t> & velocity )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName_.c_str(), std::ofstream::app );

         auto scaledPosition = position / xRef_;
         auto scaledVelocity = velocity / (xRef_ / tRef_);
         real_t Re = std::fabs(velocity[2]) * diameter_ / viscosity_;

         file << timestep << "\t" << real_c(timestep) / tRef_ << "\t"
              << "\t" << position[0] << "\t" << position[1] << "\t" << position[2]
              << "\t" << velocity[0] << "\t" << velocity[1] << "\t" << velocity[2]
              << "\t" << scaledPosition[0] << "\t" << scaledPosition[1] << "\t" << scaledPosition[2]
              << "\t" << scaledVelocity[0] << "\t" << scaledVelocity[1] << "\t" << scaledVelocity[2]
              << "\t" << Re << "\n";
         file.close();
      }
   }

   shared_ptr<SweepTimeloop> timeloop_;
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;
   std::string fileName_;
   bool fileIO_;
   real_t xRef_, tRef_;
   uint_t lbmTimeStepsPerTimeLoopIteration_;
   real_t diameter_, viscosity_;

   real_t position_;
   real_t maxVelocity_;
};

void clearBoundaryHandling( BlockForest & forest, const BlockDataID & boundaryHandlingID )
{
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      BoundaryHandling_T * boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID);
      boundaryHandling->clear( FieldGhostLayers );
   }
}

void clearBodyField( BlockForest & forest, const BlockDataID & bodyFieldID )
{
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      BodyField_T * bodyField = blockIt->getData<BodyField_T>(bodyFieldID);
      bodyField->setWithGhostLayer( NULL );
   }
}

void recreateBoundaryHandling( BlockForest & forest, const BlockDataID & boundaryHandlingID )
{
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      BoundaryHandling_T * boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID);
      boundaryHandling->fillWithDomain( FieldGhostLayers );
   }
}


class LoadImbalanceEvaluator
{
public:
   LoadImbalanceEvaluator( const shared_ptr<WcTimingPool> & levelwiseTimingPool, const shared_ptr<WcTimingTree> & peTimingTree, uint_t numberOfLevels)
         : levelwiseTimingPool_( levelwiseTimingPool ), peTimingTree_( peTimingTree ), numberOfLevels_( numberOfLevels )
   {
      // workload timer that are present on each level
      timerOnEachLevel_.emplace_back("collide");
      timerOnEachLevel_.emplace_back("boundary handling");
      timerOnEachLevel_.emplace_back("linear explosion");
      timerOnEachLevel_.emplace_back("stream");

      // workload timer only present on finest level
      timerOnFinestLevel_.emplace_back("Body Mapping");
      timerOnFinestLevel_.emplace_back("PDF Restore");
      timerOnFinestLevel_.emplace_back("stream & collide");
      //timerOnFinestLevel_.push_back("pe Time Step");

      // workload timer used in PE
      timerInPE_.emplace_back("CCD");
      timerInPE_.emplace_back("FCD");
      timerInPE_.emplace_back("Integration");

   }

   void operator()()
   {
      std::vector<real_t> minTimingsPerLevel(numberOfLevels_);
      std::vector<real_t> maxTimingsPerLevel(numberOfLevels_);
      std::vector<real_t> avgTimingsPerLevel(numberOfLevels_);
      uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());

      levelwiseTimingPool_->unifyRegisteredTimersAcrossProcesses();
      peTimingTree_->synchronize();

      for( uint_t level = 0; level < numberOfLevels_; ++level)
      {
         real_t timeOnLevelProcessLocal = real_t(0);
         for( auto timerIt = timerOnEachLevel_.begin(); timerIt != timerOnEachLevel_.end(); ++timerIt )
         {
            std::string timerName = *timerIt + " (" + std::to_string(level) + ")";
            timeOnLevelProcessLocal += real_c((*levelwiseTimingPool_)[timerName].total());
         }

         if( level == numberOfLevels_- 1)
         {
            // evaluate more timers on finest level

            for( auto timerIt = timerOnFinestLevel_.begin(); timerIt != timerOnFinestLevel_.end(); ++timerIt )
            {
               std::string timerName = *timerIt + " (" + std::to_string(level) + ")";
               timeOnLevelProcessLocal += real_c((*levelwiseTimingPool_)[timerName].total());
            }
            for( auto timerIt = timerInPE_.begin(); timerIt != timerInPE_.end(); ++timerIt )
            {
               std::string timerName = *timerIt;
               timeOnLevelProcessLocal += real_c((*peTimingTree_)[timerName].total());
            }
         }

         minTimingsPerLevel[level] = mpi::reduce( timeOnLevelProcessLocal, mpi::MIN );
         maxTimingsPerLevel[level] = mpi::reduce( timeOnLevelProcessLocal, mpi::MAX );
         avgTimingsPerLevel[level] = mpi::reduce( timeOnLevelProcessLocal, mpi::SUM ) / real_c(numProcesses);

      }

      for( uint_t level = 0; level < numberOfLevels_; ++level)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Level " << level << ": min = " << minTimingsPerLevel[level] << ", max = " << maxTimingsPerLevel[level] << ", avg = " << avgTimingsPerLevel[level]);
         WALBERLA_LOG_INFO_ON_ROOT("==> Imbalance(max/min) = " <<  maxTimingsPerLevel[level] / minTimingsPerLevel[level] << ", imbalance(max/avg) = " << maxTimingsPerLevel[level]/avgTimingsPerLevel[level]);
      }
   }

private:

   shared_ptr<WcTimingPool> levelwiseTimingPool_;
   shared_ptr<WcTimingTree> peTimingTree_;
   uint_t numberOfLevels_;
   std::vector<std::string> timerOnEachLevel_;
   std::vector<std::string> timerOnFinestLevel_;
   std::vector<std::string> timerInPE_;
};

class TimingResetter
{
public:
   TimingResetter(const shared_ptr<WcTimingPool> & levelwiseTimingPool, const shared_ptr<WcTimingTree> & peTimingTree )
         : levelwiseTimingPool_( levelwiseTimingPool ), peTimingTree_( peTimingTree )
   {}

   void operator()()
   {
      levelwiseTimingPool_->clear();
      peTimingTree_->reset();

   }

private:

   shared_ptr<WcTimingPool> levelwiseTimingPool_;
   shared_ptr<WcTimingTree> peTimingTree_;
};

//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Setup that simulates the settling of a sphere inside a rectangular column filled with viscous fluid
 *
 * The rectangular column is fully periodic and of size [xSizeNonDim x ySizeNonDim x zSizeNonDim] * diameter.
 *
 * The settling behavior can be modified via the Galileo number and the density ratio.
 * Numerical parameters are the diameter (i.e. resolution) and the characteristic settling velocity ug.
 *
 * If numLevels = 0, then a uniform grid is used.
 * Else, adaptive grid refinement according to the specified criteria is applied.
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

   // simulation control
   bool shortrun = false;
   bool funcTest = false;
   bool fileIO = true;
   uint_t vtkWriteFreqDD = 0; //domain decomposition
   uint_t vtkWriteFreqBo = 0; //bodies
   uint_t vtkWriteFreqFl = 0; //fluid
   uint_t vtkWriteFreq = 0; //general
   bool vtkWriteFluidSlice = false;
   std::string baseFolder = "vtk_out_AMRSettlingSphere"; // folder for vtk and file output

   // physical setup
   real_t GalileoNumber = real_t(200);
   real_t densityRatio = real_t(1.1);
   real_t diameter = real_t(20);
   real_t ug = real_t(0.03); // characteristic settling velocity
   uint_t xSizeNonDim = uint_t(16);
   uint_t ySizeNonDim = uint_t(16);
   uint_t zSizeNonDim = uint_t(32);

   //numerical parameters
   bool averageForceTorqueOverTwoTimSteps = true;
   uint_t numberOfLevels = uint_t(3);
   uint_t refinementCheckFrequency = uint_t(0);
   bool useStaticRefinement = false;

   real_t lowerFluidRefinementLimit = real_t(0);
   real_t upperFluidRefinementLimit = std::numeric_limits<real_t>::infinity();

   bool useVorticityCriterion = false;
   bool useGradientCriterion = false;

   // initialize the horizontal sphere velocity to avoid ambiguity due to physical instability that determines the horizontal movement direction
   bool initializeSphereVelocity = false;
   // add small offset to initial sphere position to break numerical symmetry of setup
   bool offsetSphere = false;

   // evaluate and print current imbalances in the workload
   bool evaluateLoadImbalance = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" )         == 0 ) { shortrun = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" )         == 0 ) { funcTest = true; continue; }
      if( std::strcmp( argv[i], "--fileIO" )           == 0 ) { fileIO = true; continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqDD" )   == 0 ) { vtkWriteFreqDD = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqBo" )   == 0 ) { vtkWriteFreqBo = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqFl" )   == 0 ) { vtkWriteFreqFl = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFluidSlice" ) == 0 ) { vtkWriteFluidSlice = true; continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreq" )     == 0 ) { vtkWriteFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--densityRatio" )     == 0 ) { densityRatio = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--Ga" )               == 0 ) { GalileoNumber = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--ug" )               == 0 ) { ug = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--diameter" )         == 0 ) { diameter = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--xSizeNonDim" )      == 0 ) { xSizeNonDim = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--ySizeNonDim" )      == 0 ) { ySizeNonDim = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--zSizeNonDim" )      == 0 ) { zSizeNonDim = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--noForceAveraging" ) == 0 ) { averageForceTorqueOverTwoTimSteps = false; continue; }
      if( std::strcmp( argv[i], "--numLevels" )        == 0 ) { numberOfLevels = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--refinementCheckFrequency" ) == 0 ) { refinementCheckFrequency = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--useStaticRefinement" ) == 0 ) { useStaticRefinement = true; continue; }
      if( std::strcmp( argv[i], "--lowerLimit" )       == 0 ) { lowerFluidRefinementLimit = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--upperLimit" )       == 0 ) { upperFluidRefinementLimit = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--baseFolder" )       == 0 ) { baseFolder = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--useVorticityCriterion" ) == 0 ) { useVorticityCriterion = true; continue; }
      if( std::strcmp( argv[i], "--useGradientCriterion" )  == 0 ) { useGradientCriterion = true; continue; }
      if( std::strcmp( argv[i], "--initializeSphereVelocity" )  == 0 ) { initializeSphereVelocity = true; continue; }
      if( std::strcmp( argv[i], "--offsetSphere" )  == 0 ) { offsetSphere = true; continue; }
      if( std::strcmp( argv[i], "--evaluateLoadImbalance" )  == 0 ) { evaluateLoadImbalance = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   if( funcTest )
   {
      walberla::logging::Logging::instance()->setLogLevel(logging::Logging::LogLevel::WARNING);
   }

   if( fileIO )
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

   if( vtkWriteFreq != 0 )
   {
      vtkWriteFreqDD = vtkWriteFreq;
      vtkWriteFreqBo = vtkWriteFreq;
      vtkWriteFreqFl = vtkWriteFreq;
   }


   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////

   const Vector3<uint_t> domainSize( uint_c(real_t(xSizeNonDim) * diameter ),
                                     uint_c(real_t(ySizeNonDim) * diameter ),
                                     uint_c(real_t(zSizeNonDim) * diameter ) );
   const real_t sphereVolume = math::pi / real_t(6) * diameter * diameter * diameter;

   const real_t xRef = diameter;
   const real_t tRef = xRef / ug;

   const real_t gravitationalAcceleration = ug * ug / ( (densityRatio-real_t(1)) * diameter );
   const real_t viscosity = ug * diameter / GalileoNumber;
   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);
   const real_t tau = real_t(1) / omega;

   const uint_t timesteps = funcTest ? 1 : ( shortrun ? uint_t(2000) : uint_t( 250000 ) );
   const uint_t numPeSubCycles = uint_t(1);

   const uint_t loggingDisplayFrequency = uint_t(100);

   const real_t dx = real_t(1);
   const real_t overlap = real_t( 1.5 ) * dx;

   Vector3<real_t> initialSpherePosition( real_t(0.5) * real_c(domainSize[0]),
                                          real_t(0.5) * real_c(domainSize[1]),
                                          real_t(0.5) * real_c(domainSize[2]));
   if( offsetSphere )
   {
      Vector3<real_t> offset( real_t(0.3), real_t(0.2), real_t(0));
      initialSpherePosition += offset;
   }

   if( useVorticityCriterion && floatIsEqual(lowerFluidRefinementLimit, real_t(0)) && std::isinf(upperFluidRefinementLimit) )
   {
      // use computed criterion instead of user input
      lowerFluidRefinementLimit = real_t(0.05) * ug;
      upperFluidRefinementLimit = real_t(0.1) * ug;
   }

   const uint_t finestLevel = numberOfLevels - uint_t(1);
   std::stringstream omega_msg;
   for( uint_t i = 0; i < numberOfLevels; ++i )
   {
      omega_msg << lbm::collision_model::levelDependentRelaxationParameter( i, omega, finestLevel ) << " ( on level " << i << " ), ";
   }


   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere diameter = " << diameter );
   WALBERLA_LOG_INFO_ON_ROOT(" - Galileo number = " << GalileoNumber );
   WALBERLA_LOG_INFO_ON_ROOT(" - densityRatio = " << densityRatio );
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: relaxation time (tau) = " << tau << ", kin. visc = " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration );
   WALBERLA_LOG_INFO_ON_ROOT(" - initial sphere position = " << initialSpherePosition );
   WALBERLA_LOG_INFO_ON_ROOT(" - reference values: x = " << xRef << ", t = " << tRef << ", vel = " << ug);
   WALBERLA_LOG_INFO_ON_ROOT(" - omega: " << omega_msg.str());
   WALBERLA_LOG_INFO_ON_ROOT(" - number of levels: " << numberOfLevels);
   if(useStaticRefinement)
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - using static refinement");
   }
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
      if( vtkWriteFluidSlice ){
         WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files of sliced fluid data to folder \"" << baseFolder << "\" with frequency " << vtkWriteFreqFl);
      }
      else{
         WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files of full fluid data to folder \"" << baseFolder << "\" with frequency " << vtkWriteFreqFl);
      }
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t levelScalingFactor = ( uint_t(1) << finestLevel );
   const uint_t lbmTimeStepsPerTimeLoopIteration = levelScalingFactor;

   uint_t blockSize = std::max(uint_t(16), uint_c(diameter) );
   Vector3<uint_t> blockSizeInCells( blockSize );

   AABB simulationDomain( real_t(0), real_t(0), real_t(0), real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );
   auto blocks = createBlockStructure( simulationDomain, blockSizeInCells, numberOfLevels, diameter, initialSpherePosition, useStaticRefinement );

   //write domain decomposition to file
   if( vtkWriteFreqDD > 0 )
   {
      vtk::writeDomainDecomposition( blocks, "initial_domain_decomposition", baseFolder );
   }

   if( !useStaticRefinement && refinementCheckFrequency == 0 && numberOfLevels != 1 )
   {
      // determine check frequency automatically based on maximum admissible velocity and block sizes
      real_t uMax = real_t(0.1);
      real_t refinementCheckFrequencyFinestLevel = ( overlap + real_c(blockSize) - real_t(2) * real_t(FieldGhostLayers) * dx) / uMax;
      refinementCheckFrequency = uint_c( refinementCheckFrequencyFinestLevel / real_t(lbmTimeStepsPerTimeLoopIteration));
   }
   WALBERLA_LOG_INFO_ON_ROOT(" - refinement check frequency: " << refinementCheckFrequency);


   /////////////////
   // PE COUPLING //
   /////////////////

   // set up pe functionality
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();

   auto bodyStorageID  = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID          = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   auto fcdID          = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   shared_ptr<WcTimingTree> timingTreePE = make_shared<WcTimingTree>();

   // set up collision response, here DEM solver
   pe::cr::DEM cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, &(*timingTreePE));

   // set up synchronization procedure
   std::function<void(void)> syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, &(*timingTreePE), overlap, false );

   // create pe bodies

   // add the sphere
   const auto sphereMaterial = pe::createMaterial( "mySphereMat", densityRatio , real_t(0.5), real_t(0.1), real_t(0.1), real_t(0.24), real_t(200), real_t(200), real_t(0), real_t(0) );
   auto sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, initialSpherePosition, real_t(0.5) * diameter, sphereMaterial );

   if( initializeSphereVelocity && sphere != nullptr )
   {
      Vector3<real_t> initialSphereVelocity( real_t(0.01) * ug, real_t(0.01) * ug, real_t(0));
      sphere->setLinearVel(initialSphereVelocity);
      WALBERLA_LOG_INFO(" - setting initial sphere velocity of " << initialSphereVelocity);
   }

   uint_t minBlockSizeInCells = blockSizeInCells.min();
   for( uint_t i = 0; i < uint_c(diameter / real_c(minBlockSizeInCells)) + 1; ++i)
      syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega, lbm::collision_model::TRT::threeSixteenth, finestLevel ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         FieldGhostLayers, field::zyxf );
   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::zyxf, FieldGhostLayers );

   // add velocity field and utility
   BlockDataID velocityFieldID = field::addToStorage<VelocityField_T>( blocks, "velocity field", Vector3<real_t>(real_t(0)), field::zyxf, uint_t(2) );

   using VelocityFieldWriter_T = lbm::VelocityFieldWriter<PdfField_T, VelocityField_T>;
   BlockSweepWrapper< VelocityFieldWriter_T > velocityFieldWriter( blocks, VelocityFieldWriter_T( pdfFieldID, velocityFieldID ) );


   shared_ptr<blockforest::communication::NonUniformBufferedScheme<stencil::D3Q27> > velocityCommunicationScheme = make_shared<blockforest::communication::NonUniformBufferedScheme<stencil::D3Q27> >( blocks );
   velocityCommunicationScheme->addPackInfo( make_shared< field::refinement::PackInfo<VelocityField_T, stencil::D3Q27> >( velocityFieldID ) );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addBlockData( make_shared< MyBoundaryHandling >( blocks, flagFieldID, pdfFieldID, bodyFieldID ),
                                                          "boundary handling" );

   // map planes into the LBM simulation -> act as no-slip boundaries
   pe_coupling::mapBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectGlobalBodies );

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

   ////////////////////////
   // DYNAMIC REFINEMENT //
   ////////////////////////

   auto & blockforest = blocks->getBlockForest();
   blockforest.recalculateBlockLevelsInRefresh( true );
   blockforest.alwaysRebalanceInRefresh( false );
   blockforest.reevaluateMinTargetLevelsAfterForcedRefinement( false );
   blockforest.allowRefreshChangingDepth( false );
   blockforest.allowMultipleRefreshCycles( false ); // otherwise info collections are invalid

   blockforest::CombinedMinTargetLevelDeterminationFunctions minTargetLevelDeterminationFunctions;

   // add refinement criterion based on particle presence
   shared_ptr<pe_coupling::InfoCollection> couplingInfoCollection = walberla::make_shared<pe_coupling::InfoCollection>();
   pe_coupling::amr::BodyPresenceLevelDetermination particlePresenceRefinement( couplingInfoCollection, finestLevel );

   minTargetLevelDeterminationFunctions.add( particlePresenceRefinement );

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

   blockforest.setRefreshMinTargetLevelDeterminationFunction( minTargetLevelDeterminationFunctions );

   bool curveHilbert = true; //false = use Morton
   bool curveAllGather = true;
   bool balanceLevelwise = true;
   blockforest.setRefreshPhantomBlockMigrationPreparationFunction( blockforest::DynamicCurveBalance< blockforest::PODPhantomWeight<real_t> >( curveHilbert, curveAllGather, balanceLevelwise ) );

   blockforest.setRefreshPhantomBlockDataPackFunction(blockforest::PODPhantomWeightPackUnpack<real_t>());
   blockforest.setRefreshPhantomBlockDataUnpackFunction(blockforest::PODPhantomWeightPackUnpack<real_t>());

   pe_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(couplingInfoCollection, pe_coupling::amr::defaultWeightEvaluationFunction);
   blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);


   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   shared_ptr<SweepTimeloop> timeloop = make_shared<SweepTimeloop>( blocks->getBlockStorage(), timesteps );

   shared_ptr<WcTimingPool> timeloopTiming = make_shared<WcTimingPool>();
   shared_ptr<WcTimingPool> timeloopRefinementTiming = make_shared<WcTimingPool>();
   shared_ptr<WcTimingPool> timeloopRefinementTimingLevelwise = make_shared<WcTimingPool>();

   if( vtkWriteFreqDD != uint_t(0) ) {
      auto domainDecompVTK = vtk::createVTKOutput_DomainDecomposition(blocks, "domain_decomposition", vtkWriteFreqDD, baseFolder );
      timeloop->addFuncBeforeTimeStep( vtk::writeFiles(domainDecompVTK), "VTK (domain decomposition)");
   }

   if( vtkWriteFreqBo != uint_t(0) ) {
      // pe bodies
      auto bodyVtkOutput = make_shared<pe::SphereVtkOutput>(bodyStorageID, blocks->getBlockStorage());
      auto bodyVTK = vtk::createVTKOutput_PointData(bodyVtkOutput, "bodies", vtkWriteFreqBo, baseFolder);
      timeloop->addFuncBeforeTimeStep(vtk::writeFiles(bodyVTK), "VTK (sphere data)");
   }

   if( vtkWriteFreqFl != uint_t(0) ) {
      // flag field
      //auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", vtkIOFreq, 0, false, baseFolder );
      //flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      //timeloop.addFuncAfterTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      // pdf field
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData(blocks, "fluid_field", vtkWriteFreqFl, 0, false, baseFolder);

      field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
      fluidFilter.addFlag( Fluid_Flag );

      if(vtkWriteFluidSlice)
      {
         AABB sliceAABB( real_t(0), real_c(domainSize[1])*real_t(0.5)-real_t(1), real_t(0),
                         real_c(domainSize[0]), real_c(domainSize[1])*real_t(0.5)+real_t(1), real_c(domainSize[2]) );
         vtk::AABBCellFilter aabbSliceFilter( sliceAABB );

         vtk::ChainedFilter combinedSliceFilter;
         combinedSliceFilter.addFilter( fluidFilter );
         combinedSliceFilter.addFilter( aabbSliceFilter );

         pdfFieldVTK->addCellInclusionFilter( combinedSliceFilter );
      }
      else {
         pdfFieldVTK->addCellInclusionFilter( fluidFilter );
      }

      pdfFieldVTK->addCellDataWriter(
            make_shared<lbm::VelocityVTKWriter<LatticeModel_T, float> >(pdfFieldID, "VelocityFromPDF"));
      pdfFieldVTK->addCellDataWriter(
            make_shared<lbm::DensityVTKWriter<LatticeModel_T, float> >(pdfFieldID, "DensityFromPDF"));

      timeloop->addFuncBeforeTimeStep(vtk::writeFiles(pdfFieldVTK), "VTK (fluid field data)");
   }



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

   // check for convergence of the particle position
   std::string loggingFileName( baseFolder + "/LoggingAMRSettlingSphere_Ga");
   loggingFileName += std::to_string(uint_c(GalileoNumber));
   loggingFileName += "_lvl";
   loggingFileName += std::to_string(numberOfLevels);
   loggingFileName += ".txt";
   if( fileIO  )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing logging output to file \"" << loggingFileName << "\"");
   }
   shared_ptr< SpherePropertyLogger > logger = walberla::make_shared< SpherePropertyLogger >( timeloop, blocks, bodyStorageID,
                                                                                              loggingFileName, fileIO, xRef, tRef,
                                                                                              lbmTimeStepsPerTimeLoopIteration,
                                                                                              diameter, viscosity);
   timeloop->addFuncAfterTimeStep( SharedFunctor< SpherePropertyLogger >( logger ), "Sphere property logger" );

   timeloop->addFuncAfterTimeStep( RemainingTimeLogger( timeloop->getNrOfTimeSteps() ), "Remaining Time Logger" );

   if( evaluateLoadImbalance ) timeloop->addFuncAfterTimeStep( LoadImbalanceEvaluator( timeloopRefinementTimingLevelwise, timingTreePE, numberOfLevels ), "Load Imbalance Evaluator" );

   // add level wise timing pool output
   timeloop->addFuncAfterTimeStep( TimingPoolLogger( timeloopTiming, timeloop, loggingDisplayFrequency ), "Regular Timing Logger" );

   // add regular refinement timing pool output
   timeloop->addFuncAfterTimeStep( TimingPoolLogger( timeloopRefinementTiming, timeloop, loggingDisplayFrequency ), "Refinement Timing Logger" );

   // add level wise timing pool output
   timeloop->addFuncAfterTimeStep( TimingPoolLogger( timeloopRefinementTimingLevelwise, timeloop, loggingDisplayFrequency ), "Refinement Levelwise Timing Logger" );

   // add PE timing tree output
   timeloop->addFuncAfterTimeStep( TimingTreeLogger( timingTreePE, timeloop, loggingDisplayFrequency ), "PE Timing Tree Timing Logger" );


   timeloop->addFuncAfterTimeStep( TimingResetter( timeloopRefinementTimingLevelwise, timingTreePE ), "Timing Resetter" );

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   real_t terminationPosition = diameter;
   real_t curPos = initialSpherePosition[2];
   real_t oldPos = initialSpherePosition[2];

   uint_t numberOfPassesThroughTerminationPosition = 3;
   uint_t passCounter = uint_t(0);

   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {

      if( refinementCheckFrequency != 0 && i % refinementCheckFrequency == 0)
      {
         // first evaluate all data that is required for the refinement checks

         // for the particle presence based check:
         auto & forest = blocks->getBlockForest();
         pe_coupling::createWithNeighborhood<BoundaryHandling_T>(forest, boundaryHandlingID, bodyStorageID, ccdID, fcdID, numPeSubCycles, *couplingInfoCollection);

         // for the fluid property based check:
         if( useVorticityCriterion || useGradientCriterion )
         {
            velocityFieldWriter();
            (*velocityCommunicationScheme)();
         }

         // check refinement criteria and refine/coarsen if necessary
         uint_t stampBefore = blocks->getBlockForest().getModificationStamp();
         blocks->refresh();
         uint_t stampAfter = blocks->getBlockForest().getModificationStamp();

         if(stampBefore == stampAfter)
         {
            // nothing has changed
            continue;
         }

         WALBERLA_LOG_INFO_ON_ROOT("Adapting grid and reinitializing data structures");

         // rebuild PE data structures
         pe::clearSynchronization( blockforest, bodyStorageID);

         for( uint_t syncStep = 0; syncStep < uint_c(diameter / real_c(minBlockSizeInCells)) + 1; ++syncStep)
            syncCall();

         for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
         {
            pe::ccd::ICCD* ccd = blockIt->getData< pe::ccd::ICCD >( ccdID );
            ccd->reloadBodies();
         }

         clearBoundaryHandling(forest, boundaryHandlingID);
         clearBodyField(forest, bodyFieldID);

         if( averageForceTorqueOverTwoTimSteps ) {

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
         pe_coupling::mapBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectGlobalBodies );

         // re-map the body into the domain (initializing the bodyField as well)
         pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectRegularBodies );

         // some evaluation
         uint_t numBlocksFinestLevel = forest.getNumberOfBlocks(numberOfLevels-1);
         mpi::allReduceInplace(numBlocksFinestLevel, mpi::SUM);
         WALBERLA_LOG_INFO_ON_ROOT("Total number of blocks on finest level = " << numBlocksFinestLevel);

      }

      // perform a single simulation step
      timeloop->singleStep( *timeloopTiming );

      oldPos = curPos;
      curPos = logger->getPosition();

      if( curPos <= terminationPosition && oldPos > terminationPosition )
      {
         ++passCounter;
         if( passCounter == numberOfPassesThroughTerminationPosition )
         {
            WALBERLA_LOG_INFO_ON_ROOT("Sphere passed terminal position " << terminationPosition << " for the " << passCounter << ". time - terminating simulation!");
            break;
         }
      }
   }

   timeloopTiming->logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace amr_settling_sphere

int main( int argc, char **argv ){
   amr_settling_sphere::main(argc, argv);
}
