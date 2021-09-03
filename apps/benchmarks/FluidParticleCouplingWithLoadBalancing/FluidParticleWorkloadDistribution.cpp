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
//! \file FluidParticleWorkloadDistribution.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/loadbalancing/all.h"
#include "blockforest/AABBRefinementSelection.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/mpi/Broadcast.h"

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
#include "lbm/sweeps/CellwiseSweep.h"

#include "lbm_mesapd_coupling/amr/BlockInfo.h"
#include "lbm_mesapd_coupling/amr/InfoCollection.h"
#include "lbm_mesapd_coupling/amr/weight_assignment/WeightEvaluationFunctions.h"
#include "lbm_mesapd_coupling/amr/weight_assignment/WeightAssignmentFunctor.h"
#include "lbm_mesapd_coupling/amr/weight_assignment/MetisAssignmentFunctor.h"
#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/InspectionProbe.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/kernel/AssocToBlock.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ClearNextNeighborSync.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/SyncNextNeighborsBlockForest.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include "Utility.h"

namespace fluid_particle_workload_distribution
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

using LatticeModel_T = lbm::D3Q19< lbm::collision_model::TRT, false >;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;
using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

using ParticleField_T = lbm_mesapd_coupling::ParticleField_T;
using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;

const uint_t FieldGhostLayers = 1;

using NoSlip_T = lbm::NoSlip<LatticeModel_T, flag_t> ;
using MO_T = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
using BoundaryHandling_T = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, MO_T >;


///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID NoSlip_Flag( "no slip" );
const FlagUID MO_Flag( "moving obstacle" );
const FlagUID FormerMO_Flag( "former moving obstacle" );


///////////////////////////////////
// LOAD BALANCING FUNCTIONALITY //
//////////////////////////////////



/////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void workloadAndMemoryAssignment( SetupBlockForest& forest )
{
   for (auto &block : forest) {
      block.setWorkload( numeric_cast< workload_t >( uint_t(1) << block.getLevel() ) );
      block.setMemory( numeric_cast< memory_t >(1) );
   }
}

static shared_ptr< StructuredBlockForest > createBlockStructure( const AABB & domainAABB, Vector3<uint_t> blockSizeInCells,
                                                                 bool useBox, const std::string & loadDistributionStrategy,
                                                                 bool keepGlobalBlockInformation = false )
{
   SetupBlockForest sforest;

   Vector3<uint_t> numberOfBlocksPerDirection( uint_c(domainAABB.size(0)) / blockSizeInCells[0],
                                               uint_c(domainAABB.size(1)) / blockSizeInCells[1],
                                               uint_c(domainAABB.size(2)) / blockSizeInCells[2] );

   for(uint_t i = 0; i < 3; ++i )
   {
      WALBERLA_CHECK_EQUAL( numberOfBlocksPerDirection[i] * blockSizeInCells[i], uint_c(domainAABB.size(i)),
                            "Domain can not be decomposed in direction " << i << " into fine blocks of size " << blockSizeInCells[i] );
   }

   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

   Vector3<bool> periodicity( true, true, false);
   if( useBox )
   {
      periodicity[0] = false;
      periodicity[1] = false;
   }
   sforest.init( domainAABB,
                 numberOfBlocksPerDirection[0], numberOfBlocksPerDirection[1], numberOfBlocksPerDirection[2],
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
                       const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID,
                       const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor_T>& ac ) :
         blocks_( blocks ), flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), particleFieldID_( particleFieldID ), ac_( ac )
   {}

   BoundaryHandling_T * initialize( IBlock * const block ) override;

private:

   weak_ptr< StructuredBlockStorage > blocks_;

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;
   shared_ptr<ParticleAccessor_T> ac_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::initialize( IBlock * const block )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   auto * flagField = block->getData< FlagField_T >( flagFieldID_ );
   auto *  pdfField = block->getData< PdfField_T > ( pdfFieldID_ );
   auto* particleField = block->getData<lbm_mesapd_coupling::ParticleField_T>(particleFieldID_);

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   auto blocksPtr = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocksPtr );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "moving obstacle boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                                           MO_T( "MO", MO_Flag, pdfField, flagField, particleField, ac_, fluid, *blocksPtr, *block ),
                                                           BoundaryHandling_T::Mode::ENTIRE_FIELD_TRAVERSAL);

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}

void clearBoundaryHandling( BlockForest & forest, const BlockDataID & boundaryHandlingID ) {
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      BoundaryHandling_T * boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID);
      boundaryHandling->clear( FieldGhostLayers );
   }
}

void recreateBoundaryHandling( BlockForest & forest, const BlockDataID & boundaryHandlingID ) {
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      BoundaryHandling_T * boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID);
      boundaryHandling->fillWithDomain( FieldGhostLayers );
   }
}

void clearParticleField( BlockForest & forest, const BlockDataID & particleFieldID, const ParticleAccessor_T & accessor ) {
   for (auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt) {
      ParticleField_T *particleField = blockIt->getData<ParticleField_T>(particleFieldID);
      particleField->setWithGhostLayer(accessor.getInvalidUid());
   }
}
//*******************************************************************************************************************

void evaluateTimers(WcTimingPool & timingPool,
                    const std::vector<std::vector<std::string> > & timerKeys,
                    std::vector<real_t> & timings )
{

   for (auto & timingsIt : timings)
   {
      timingsIt = real_t(0);
   }

   timingPool.unifyRegisteredTimersAcrossProcesses();

   for (auto i = uint_t(0); i < timerKeys.size(); ++i )
   {
      auto keys = timerKeys[i];
      for (const auto &timerName : keys)
      {
         if(timingPool.timerExists(timerName))
         {
            timings[i] += real_c(timingPool[timerName].total());
         }
      }

   }
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


void evaluateTotalSimulationTimePassed(WcTimingPool & timeloopTimingPool, const std::string & simulationString, const std::string & loadbalancingString,
                                       double & totalSimTime, double & totalLBTime)
{
   shared_ptr< WcTimingPool> reduced = timeloopTimingPool.getReduced(timing::REDUCE_TOTAL, 0);

   double totalTime = 0.0;
   WALBERLA_ROOT_SECTION(){
      totalTime = (*reduced)[simulationString].total();
   }
   totalSimTime = totalTime;

   double lbTime = 0.0;
   WALBERLA_ROOT_SECTION(){
      lbTime = (*reduced)[loadbalancingString].total();
   }
   totalLBTime = lbTime;
}

void createPlane( const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss,
                  const Vector3<real_t> position, const Vector3<real_t> normal)
{
   mesa_pd::data::Particle&& p0 = *ps->create(true);
   p0.setPosition(position);
   p0.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p0.setShapeID(ss->create<mesa_pd::data::HalfSpace>( normal.getNormalized() ));
   p0.setOwner(mpi::MPIManager::instance()->rank());
   p0.setType(0);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   WALBERLA_LOG_INFO_ON_ROOT("Created plane at position " << position << " with normal " << normal.getNormalized  ());
}

void createBasicPlaneSetup(const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss, const math::AABB & simulationDomain)
{
   createPlane(ps, ss, simulationDomain.minCorner(), Vector3<real_t>(0,0, 1));
   //createPlane(ps, ss, Vector3<real_t>(simulationDomain.xMax() * 0.5, simulationDomain.yMax() * 0.5, simulationDomain.zMin() + std::max(simulationDomain.xMax(), simulationDomain.yMax()) * 0.5 ), Vector3<real_t>(0,1, 1));
   createPlane(ps, ss, simulationDomain.maxCorner(), Vector3<real_t>(0,0,-1));
}

void addBoxPlaneSetup(const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss, const math::AABB & simulationDomain)
{
   // add bounding planes (four horizontal directions)
   createPlane(ps, ss, simulationDomain.minCorner(), Vector3<real_t>( 1, 0,0));
   createPlane(ps, ss, simulationDomain.maxCorner(), Vector3<real_t>(-1, 0,0));
   createPlane(ps, ss, simulationDomain.minCorner(), Vector3<real_t>( 0, 1,0));
   createPlane(ps, ss, simulationDomain.maxCorner(), Vector3<real_t>( 0,-1,0));
}

void addHopperPlaneSetup(const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss, const math::AABB & simulationDomain,
                         real_t hopperRelHeight, real_t hopperRelOpening)
{
   //hopper planes
   real_t xMax = simulationDomain.xMax();
   real_t yMax = simulationDomain.yMax();
   real_t zMax = simulationDomain.zMax();
   Vector3<real_t> p1(0,0,hopperRelHeight*zMax);
   Vector3<real_t> n1(p1[2],0,hopperRelOpening*xMax-p1[0]);
   Vector3<real_t> n2(0,p1[2],hopperRelOpening*yMax-p1[0]);
   createPlane(ps, ss, p1, n1);
   createPlane(ps, ss, p1, n2);

   Vector3<real_t> p2(xMax,yMax,hopperRelHeight*zMax);
   Vector3<real_t> n3(-p2[2],0,-((real_t(1)-hopperRelOpening)*xMax-p2[0]));
   Vector3<real_t> n4(0,-p2[2],-((real_t(1)-hopperRelOpening)*yMax-p2[1]));
   createPlane(ps, ss, p2, n3);
   createPlane(ps, ss, p2, n4);
}

void evaluateParticleSimulation(const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss,
                                uint_t numParticles, uint_t numGlobalParticles)
{
   auto numShapes = ss->shapes.size();
   uint_t numLocalParticles = 0;
   uint_t numGhostParticles = 0;
   uint_t numGlobalParticlesOfRank = 0;
   for(auto p = ps->begin(); p != ps->end(); ++p)
   {
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(p->getFlags(), GHOST))
      {
         ++numGhostParticles;
      } else if (isSet(p->getFlags(), GLOBAL))
      {
         ++numGlobalParticlesOfRank;
      } else
      {
         ++numLocalParticles;
      }
      if(p->getShapeID() >= numShapes)
      {
         WALBERLA_LOG_INFO("Found invalid shape id " << *p);
      }
   }
   //WALBERLA_LOG_INFO(numShapes << " " << numLocalParticles << " " << numGhostParticles);

   if(numGlobalParticlesOfRank != numGlobalParticles)
   {
      WALBERLA_LOG_INFO("Number of global particles has changed to " << numGlobalParticlesOfRank);
   }

   mpi::reduceInplace(numLocalParticles, mpi::SUM);
   if(numLocalParticles != numParticles)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Number of particles has changed to " << numLocalParticles);
   }
}

bool isnan(Vector3<real_t> vec)
{
   return std::isnan(vec[0]) || std::isnan(vec[1]) || std::isnan(vec[2]);
}

void checkParticleProperties(const shared_ptr<mesa_pd::data::ParticleStorage> & ps)
{
   for(auto p = ps->begin(); p != ps->end(); ++p)
   {
      if(isnan(p->getHydrodynamicForce())) WALBERLA_LOG_INFO("Found nan in hyd Force " << *p);
      if(isnan(p->getOldHydrodynamicForce())) WALBERLA_LOG_INFO("Found nan in old hyd Force " << *p);
      if(isnan(p->getForce())) WALBERLA_LOG_INFO("Found nan in Force " << *p);
      if(isnan(p->getOldForce())) WALBERLA_LOG_INFO("Found nan in Force " << *p);
   }
}

template< typename Probe_T>
void checkMapping(const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID pdfFieldID,
                  const BlockDataID boundaryHandlingID, Probe_T & probe )
{

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto * pdfField = blockIt->getData< PdfField_T >( pdfFieldID );
      auto * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID );

      WALBERLA_FOR_ALL_CELLS_XYZ(pdfField,
          if (boundaryHandling->isDomain(x, y, z))
          {
             uint_t f = uint_t(0);
             //for( uint_t f = uint_t(0); f < PdfField_T::F_SIZE; ++f )
             //{
                if( !walberla::field::internal::stabilityCheckerIsFinite( pdfField->get( x, y, z, cell_idx_c(f) ) ) )
                {

                   Vector3< real_t > center;
                   blocks->getBlockLocalCellCenter( *blockIt, Cell(x,y,z), center );

                   Cell gCell(x,y,z);
                   blocks->transformBlockLocalToGlobalCell( gCell, *blockIt);

                   WALBERLA_LOG_INFO("Instability found in block local cell( " << x << ", " << y << ", " << z << " ) at index " << f
                                         << " = global cell ( " << gCell.x() << ", " << gCell.y() << ", " << gCell.z() << " ) with cell center ( " << center[0] << ", " << center[1] << ", " << center[2] << " )");

                   probe.setPosition(center);
                   real_t rho;
                   Vector3<real_t> velocity;
                   probe(rho, velocity);

                }
             //}
          }
          );

   }
}

//*******************************************************************************************************************
/*!\brief Simulation of settling particles inside a rectangular column filled with viscous fluid
 *
 * This application is used in the paper
 *  Rettinger, Ruede - "Dynamic Load Balancing Techniques for Particulate Flow Simulations", 2019, Computation
 * in Section 4 to apply the load estimator and to evaluate different load distribution strategies.
 * It has here been adapted to the new mesapd and its coupling.
 * More infos can be found in the PhD thesis by Christoph Rettinger.
 *
 * It, however, features several different command line arguments that can be used to tweak the simulation.
 * The setup can be horizontally period, a box or a hopper geometry (configurable, as in the paper).
 * The size, resolution and used blocks for the domain partitioning can be changed.
 * Initially, all particles are pushed upwards to obtain a dense packing at the top plane.
 *
 * Most importantly, the load balancing can be modified:
 *  - load estimation strategies:
 *    - pure LBM = number of cells per block = constant workload per block
 *    - coupling based load estimator = use fitted function from Sec. 3 of paper
 *  - load distribution strategies:
 *    - space-filling curves: Hilbert and Morton
 *    - ParMETIS (and several algorithms and parameters)
 *    - diffusive (and options)
 *  - load balancing frequency
 */
//*******************************************************************************************************************
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   MPIManager::instance()->useWorldComm();

   ///////////////////
   // Customization //
   ///////////////////

   // simulation control
   bool shortRun = false;
   bool funcTest = false;
   bool fileIO = true;
   bool checkSimulation = false;
   uint_t vtkWriteFreqDD = 0; //domain decomposition
   uint_t vtkWriteFreqPa = 0; //particles
   uint_t vtkWriteFreqFl = 0; //fluid
   uint_t vtkWriteFreq = 0; //general
   uint_t vtkWriteFreqInit = 0; //initial (particle-only) simulation
   std::string baseFolder = "vtk_out_WorkloadDistribution"; // folder for vtk and file output
   bool useProgressLogging = false;

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
   auto hopperRelHeight = real_t(0.5); // for hopper setup
   auto hopperRelOpening = real_t(0.3); // for hopper setup

   auto timesteps = uint_t(80000);

   //numerical parameters
   auto loadBalancingCheckFrequency = uint_t(100);
   auto numRPDSubCycles = uint_t(10);
   bool useBlockForestSync = false;

   // load balancing
   std::string loadEvaluationStrategy = "LBM"; //LBM, Fit
   std::string loadDistributionStrategy = "Hilbert"; //Morton, Hilbert, ParMetis, Diffusive
   real_t blockBaseWeight = real_t(1);

   auto parMetis_ipc2redist = real_t(1000);
   auto parMetisTolerance = real_t(-1);
   std::string parMetisAlgorithmString = "ADAPTIVE_REPART";

   auto diffusionFlowIterations = uint_t(15);
   auto diffusionMaxIterations = uint_t(20);

   bool useNoSlipForPlanes = false;


   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortRun" )                 == 0 ) { shortRun = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" )                 == 0 ) { funcTest = true; continue; }
      if( std::strcmp( argv[i], "--fileIO" )                   == 0 ) { fileIO = true; continue; }
      if( std::strcmp( argv[i], "--useProgressLogging" )       == 0 ) { useProgressLogging = true; continue; }
      if( std::strcmp( argv[i], "--checkSimulation" )          == 0 ) { checkSimulation = true; continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqDD" )           == 0 ) { vtkWriteFreqDD = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqPa" )           == 0 ) { vtkWriteFreqPa = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqFl" )           == 0 ) { vtkWriteFreqFl = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreq" )             == 0 ) { vtkWriteFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--vtkWriteFreqInit" )         == 0 ) { vtkWriteFreqInit = uint_c( std::atof( argv[++i] ) ); continue; }
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
      if( std::strcmp( argv[i], "--timesteps" )                == 0 ) { timesteps = uint_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--numRPDSubCycles" )          == 0 ) { numRPDSubCycles = uint_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--useBlockForestSync" )       == 0 ) { useBlockForestSync = true; continue; }
      if( std::strcmp( argv[i], "--useNoSlipForPlanes" )       == 0 ) { useNoSlipForPlanes = true; continue; }
      if( std::strcmp( argv[i], "--blockBaseWeight" )          == 0 ) { blockBaseWeight = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--loadBalancingCheckFrequency" ) == 0 ) { loadBalancingCheckFrequency = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--loadEvaluationStrategy" )   == 0 ) { loadEvaluationStrategy = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--loadDistributionStrategy" ) == 0 ) { loadDistributionStrategy = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--ipc2redist" )               == 0 ) { parMetis_ipc2redist = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--parMetisTolerance" )        == 0 ) { parMetisTolerance = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--parMetisAlgorithm" )        == 0 ) { parMetisAlgorithmString = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--diffusionFlowIterations" )  == 0 ) { diffusionFlowIterations = uint_c(std::atof(argv[++i])); continue; }
      if( std::strcmp( argv[i], "--diffusionMaxIterations" )   == 0 ) { diffusionMaxIterations = uint_c(std::atof(argv[++i])); continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   if( fileIO )
   {
      WALBERLA_ROOT_SECTION(){
         // create base directory if it does not yet exist
         filesystem::path tpath( baseFolder );
         if( !filesystem::exists( tpath ) )
            filesystem::create_directory( tpath );
      }
      WALBERLA_MPI_BARRIER();
   }

   if(useProgressLogging) logging::Logging::instance()->includeLoggingToFile(baseFolder + "/progress_logging.txt");

   if( loadEvaluationStrategy != "LBM" && loadEvaluationStrategy != "Fit" && loadEvaluationStrategy != "FitMulti")
   {
      WALBERLA_ABORT("Invalid load evaluation strategy: " << loadEvaluationStrategy);
   }

   if( vtkWriteFreq != 0 )
   {
      vtkWriteFreqDD = vtkWriteFreq;
      vtkWriteFreqPa = vtkWriteFreq;
      vtkWriteFreqFl = vtkWriteFreq;
   }

   if( diameter > real_c(blockSize) )
   {
      WALBERLA_LOG_WARNING("Particle Synchronization might not work since bodies are large compared to block size!");
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

   const auto dx = real_t(1);
   const real_t overlap = real_t( 1.5 ) * dx;

   timesteps = funcTest ? 1 : ( shortRun ? uint_t(100) : timesteps );

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - sediment diameter = " << diameter );
   WALBERLA_LOG_INFO_ON_ROOT(" - Galileo number = " << GalileoNumber );
   WALBERLA_LOG_INFO_ON_ROOT(" - number of sediments: " << numberOfSediments);
   WALBERLA_LOG_INFO_ON_ROOT(" - densityRatio = " << densityRatio );
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: relaxation time (tau) = " << tau << ", kin. visc = " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration );
   WALBERLA_LOG_INFO_ON_ROOT(" - reference values: x = " << xRef << ", t = " << tRef << ", vel = " << uRef);
   WALBERLA_LOG_INFO_ON_ROOT(" - omega: " << omega);
   if( vtkWriteFreqDD > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files of domain decomposition to folder \"" << baseFolder << "\" with frequency " << vtkWriteFreqDD);
   }
   if( vtkWriteFreqPa > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files of bodies data to folder \"" << baseFolder << "\" with frequency " << vtkWriteFreqPa);
   }
   if( vtkWriteFreqFl > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files of fluid data to folder \"" << baseFolder << "\" with frequency " << vtkWriteFreqFl);
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

   WALBERLA_LOG_INFO_ON_ROOT(" - refinement / load balancing check frequency: " << loadBalancingCheckFrequency);
   WALBERLA_LOG_INFO_ON_ROOT(" - load evaluation strategy: " << loadEvaluationStrategy);
   WALBERLA_LOG_INFO_ON_ROOT(" - load distribution strategy: " << loadDistributionStrategy);

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3<uint_t> blockSizeInCells( blockSize );

   AABB simulationDomain( real_t(0), real_t(0), real_t(0), real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );
   AABB sedimentDomain( real_t(0), real_t(0), real_c(domainSize[2]) - expectedSedimentedHeight, real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );

   auto blocks = createBlockStructure( simulationDomain, blockSizeInCells, (useBox||useHopper), loadDistributionStrategy );

   //write initial domain decomposition to file
   if( vtkWriteFreqDD > 0 )
   {
      vtk::writeDomainDecomposition( blocks, "initial_domain_decomposition", baseFolder );
   }

   /////////
   // RPD //
   /////////

   const real_t restitutionCoeff = real_t(0.97);
   const real_t frictionCoeffStatic = real_t(0.8);
   const real_t frictionCoeffDynamic = real_t(0.15);
   const real_t collisionTime = real_t(4) * diameter; // from my paper
   const real_t poissonsRatio = real_t(0.22);
   const real_t kappa = real_t(2) * ( real_t(1) - poissonsRatio ) / ( real_t(2) - poissonsRatio ) ;

   auto rpdDomain = std::make_shared<mesa_pd::domain::BlockForestDomain>(blocks->getBlockForestPointer());

   //init data structures
   auto ps = walberla::make_shared<mesa_pd::data::ParticleStorage>(1);
   blocks->addBlockData(mesa_pd::domain::createBlockForestDataHandling(ps), "Particle Storage"); // returned ID is not used, but ps has to be known to blockforest
   auto ss = walberla::make_shared<mesa_pd::data::ShapeStorage>();
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);

   real_t timeStepSizeRPD = real_t(1)/real_t(numRPDSubCycles);
   mesa_pd::kernel::VelocityVerletPreForceUpdate  vvIntegratorPreForce(timeStepSizeRPD);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(timeStepSizeRPD);

   // types: 0 = wall, 1: sphere
   mesa_pd::kernel::LinearSpringDashpot collisionResponse(2);
   collisionResponse.setFrictionCoefficientDynamic(0,1,frictionCoeffDynamic);
   collisionResponse.setFrictionCoefficientDynamic(1,1,frictionCoeffDynamic);
   collisionResponse.setFrictionCoefficientStatic(0,1,frictionCoeffStatic);
   collisionResponse.setFrictionCoefficientStatic(1,1,frictionCoeffStatic);

   const real_t particleMass = densityRatio * sphereVolume;
   const real_t effMass_sphereWall = particleMass;
   const real_t effMass_sphereSphere = particleMass * particleMass / ( real_t(2) * particleMass );
   collisionResponse.setStiffnessAndDamping(0,1,restitutionCoeff,collisionTime,kappa,effMass_sphereWall);
   collisionResponse.setStiffnessAndDamping(1,1,restitutionCoeff,collisionTime,kappa,effMass_sphereSphere);

   mesa_pd::mpi::ReduceProperty reduceProperty;
   mesa_pd::mpi::ReduceContactHistory reduceAndSwapContactHistory;
   mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
   mesa_pd::kernel::AssocToBlock associateToBlock(blocks->getBlockForestPointer());
   mesa_pd::mpi::SyncNextNeighborsBlockForest syncNextNeighborBlockForestFunc;


   // create bounding planes

   createBasicPlaneSetup(ps, ss, simulationDomain);
   if(useBox || useHopper) addBoxPlaneSetup(ps, ss, simulationDomain);
   if(useHopper) addHopperPlaneSetup(ps, ss, simulationDomain, hopperRelHeight, hopperRelOpening);

   auto numGlobalParticles = ps->size();
   WALBERLA_LOG_INFO_ON_ROOT("Created " << numGlobalParticles << " global particles");

   // add the sediments

   WALBERLA_LOG_INFO_ON_ROOT("Starting creation of sediments");

   AABB sedimentGenerationDomain( real_t(0), real_t(0), real_t(0.5)*real_c(domainSize[2]), real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );

   auto xParticle = real_t(0);
   auto yParticle = real_t(0);
   auto zParticle = real_t(0);

   auto rank = mpi::MPIManager::instance()->rank();

   auto sphereShape = ss->create<mesa_pd::data::Sphere>( diameter * real_t(0.5) );
   ss->shapes[sphereShape]->updateMassAndInertia(densityRatio);

   std::mt19937 randomGenerator (static_cast<unsigned int>(2610)); // fixed seed: quasi-random and reproducible

   for( uint_t nSed = 0; nSed < numberOfSediments; ++nSed )
   {

      WALBERLA_ROOT_SECTION()
      {
         xParticle = math::realRandom<real_t>(sedimentGenerationDomain.xMin(), sedimentGenerationDomain.xMax(),randomGenerator);
         yParticle = math::realRandom<real_t>(sedimentGenerationDomain.yMin(), sedimentGenerationDomain.yMax(),randomGenerator);
         zParticle = math::realRandom<real_t>(sedimentGenerationDomain.zMin(), sedimentGenerationDomain.zMax(),randomGenerator);
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::broadcastObject( xParticle );
         mpi::broadcastObject( yParticle );
         mpi::broadcastObject( zParticle );
      }

      auto position = Vector3<real_t>( xParticle, yParticle, zParticle );

      if (!rpdDomain->isContainedInProcessSubdomain(uint_c(rank), position)) continue;
      auto p = ps->create();
      p->setPosition(position);
      p->setInteractionRadius(diameter * real_t(0.5));
      p->setShapeID(sphereShape);
      p->setType(1);
      p->setOwner(rank);

   }

   if(useBlockForestSync)
   {
      ps->forEachParticle(false, mesa_pd::kernel::SelectLocal(), *accessor, associateToBlock, *accessor);
      syncNextNeighborBlockForestFunc(*ps, blocks->getBlockForestPointer(), rpdDomain, overlap);

   } else
   {
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   }


   // Carry out particle-only simulation to obtain dense packing at top plane
   // consists of three phases:
   // 1: carefully resolve initial overlaps due to random generation, no gravity
   // 2: apply low gravity in positive z-direction to create packing until convergence or targeted packing height reached
   // 3: carry out a few time steps with gravity in negative direction to relay the system towards the real setup

   const bool useOpenMP = false;
   const real_t dt_RPD_Init = real_t(1);
   const auto particleSimStepsPhase1 = uint_t(1000);
   const auto maxParticleSimStepsPhase2 = (shortRun) ? uint_t(10) : uint_t(200000);
   const auto particleSimStepsPhase3 = uint_t(std::sqrt(real_t(2)/std::fabs(gravitationalAcceleration)));

   uint_t maxInitialParticleSimSteps = particleSimStepsPhase1 + maxParticleSimStepsPhase2 + particleSimStepsPhase3;

   auto particleVtkOutput = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
   particleVtkOutput->addOutput<mesa_pd::data::SelectParticleLinearVelocity>("velocity");
   particleVtkOutput->addOutput<mesa_pd::data::SelectParticleOwner>("owner");
   particleVtkOutput->setParticleSelector( [sphereShape](const mesa_pd::data::ParticleStorage::iterator& pIt) {return !mesa_pd::data::particle_flags::isSet(pIt->getFlags(), mesa_pd::data::particle_flags::GHOST) && pIt->getShapeID() == sphereShape;} ); //limit output to local sphere
   auto particleVtkWriterInit = vtk::createVTKOutput_PointData(particleVtkOutput, "Particles_init", 1, baseFolder, "simulation_step");

   real_t gravitationalAccelerationGeneration = gravitationalAcceleration;
   auto oldMinParticlePosition = real_t(0);
   real_t phase2ConvergenceLimit = std::fabs(gravitationalAccelerationGeneration);
   real_t heightConvergenceThreshold = sedimentDomain.zMin();

   uint_t beginOfPhase3SimStep = uint_t(0);

   uint_t currentPhase = 1;

   for(auto pet = uint_t(0); pet <= maxInitialParticleSimSteps; ++pet )
   {

      real_t maxPenetrationDepth = 0;
      ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                  [&collisionResponse, &rpdDomain, &maxPenetrationDepth, dt_RPD_Init]
                                        (const size_t idx1, const size_t idx2, auto& ac)
                                  {
                                     mesa_pd::collision_detection::AnalyticContactDetection acd;
                                     mesa_pd::kernel::DoubleCast double_cast;
                                     mesa_pd::mpi::ContactFilter contact_filter;
                                     if (double_cast(idx1, idx2, ac, acd, ac ))
                                     {
                                        if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                        {
                                           maxPenetrationDepth = std::max(maxPenetrationDepth, std::abs(acd.getPenetrationDepth()));
                                           collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(),
                                                             acd.getContactNormal(), acd.getPenetrationDepth(), dt_RPD_Init);
                                        }
                                     }
                                  },
                                  *accessor );

      reduceAndSwapContactHistory(*ps);

      mpi::allReduceInplace(maxPenetrationDepth, mpi::MAX);

      reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);

      ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, mesa_pd::kernel::ExplicitEuler(dt_RPD_Init), *accessor);
      if(useBlockForestSync)
      {
         syncNextNeighborBlockForestFunc(*ps, blocks->getBlockForestPointer(), rpdDomain, overlap);
         ps->forEachParticle(false, mesa_pd::kernel::SelectLocal(), *accessor, associateToBlock, *accessor);

      } else
      {
         syncNextNeighborFunc(*ps, *rpdDomain, overlap);
      }

      if( vtkWriteFreqInit > uint_t(0) && pet % vtkWriteFreqInit == uint_t(0) )
      {
         particleVtkWriterInit->write();
      }


      if(currentPhase == 1)
      {
         // damp velocities to avoid too large ones
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
                              [](const size_t idx, ParticleAccessor_T& ac){
                                 ac.setLinearVelocity(idx, ac.getLinearVelocity(idx) * real_t(0.5));
                                 ac.setAngularVelocity(idx, ac.getAngularVelocity(idx) * real_t(0.5));
                              }, *accessor);

         if(pet > particleSimStepsPhase1)
         {
            WALBERLA_LOG_INFO_ON_ROOT("Starting phase 2 of initial particle simulation, with height threshold = " << heightConvergenceThreshold);
            currentPhase = 2;

            ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
                                 [](const size_t idx, ParticleAccessor_T& ac){
                                    ac.setLinearVelocity(idx, Vector3<real_t>(0.0));
                                    ac.setAngularVelocity(idx, Vector3<real_t>(0.0));
                                 }, *accessor);
         }
      } else if(currentPhase == 2)
      {

         Vector3<real_t> gravitationalForce( real_t(0), real_t(0), (densityRatio - real_t(1)) * gravitationalAccelerationGeneration * sphereVolume );
         lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor );


         real_t minParticlePosition = sedimentGenerationDomain.zMax();
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
                              [&minParticlePosition](const size_t idx, ParticleAccessor_T& ac){
                                 minParticlePosition = std::min(ac.getPosition(idx)[2], minParticlePosition);
                              }, *accessor);

         WALBERLA_MPI_SECTION()
         {
            mpi::allReduceInplace(minParticlePosition, mpi::MIN);
         }

         WALBERLA_ROOT_SECTION()
         {
            if( pet % 100 == 0)
            {
               WALBERLA_LOG_INFO("[" << pet << "] Min position of all particles = " << minParticlePosition << " with goal height " << heightConvergenceThreshold);
            }
         }

         if( minParticlePosition > heightConvergenceThreshold ) currentPhase = 3;

         if( pet % 500 == 0)
         {
            if( std::fabs(minParticlePosition - oldMinParticlePosition) / minParticlePosition  < phase2ConvergenceLimit ) currentPhase = 3;
            oldMinParticlePosition = minParticlePosition;
         }

         if( currentPhase == 3)
         {
            WALBERLA_LOG_INFO_ON_ROOT("Starting phase 3 of initial particle simulation");
            beginOfPhase3SimStep = pet;
         }

      } else if(currentPhase == 3)
      {
         Vector3<real_t> gravitationalForce( real_t(0), real_t(0), -(densityRatio - real_t(1)) * gravitationalAccelerationGeneration * sphereVolume );
         lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor );

         if(pet - beginOfPhase3SimStep > particleSimStepsPhase3)
         {
            Vector3<real_t> initialParticleVelocity(real_t(0));
            WALBERLA_LOG_INFO_ON_ROOT("Setting initial velocity " << initialParticleVelocity << " of all particles");
            // reset velocities to avoid too large ones
            ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
                                 [initialParticleVelocity](const size_t idx, ParticleAccessor_T& ac){
                                    ac.setLinearVelocity(idx, initialParticleVelocity);
                                    ac.setAngularVelocity(idx, Vector3<real_t>(0));
                                 }, *accessor);
            break;
         }
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT("Sediment layer creation done!");

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega, lbm::collision_model::TRT::threeSixteenth ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         FieldGhostLayers, field::zyxf );
   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::zyxf, FieldGhostLayers );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addBlockData( make_shared< MyBoundaryHandling >( blocks, flagFieldID, pdfFieldID, particleFieldID, accessor ), "boundary handling" );

   Vector3<real_t> gravitationalForce( real_t(0), real_t(0), -(densityRatio - real_t(1)) * gravitationalAcceleration * sphereVolume );
   lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel(viscosity, [](real_t r){return (real_t(0.001 + real_t(0.00007)*r))*r;});
   lbm_mesapd_coupling::ParticleMappingKernel<BoundaryHandling_T> particleMappingKernel(blocks, boundaryHandlingID);
   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);
   lbm_mesapd_coupling::ParticleMappingKernel<BoundaryHandling_T> fixedParticleMappingKernel(blocks, boundaryHandlingID);
   lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel initializeHydrodynamicForceTorqueForAveragingKernel;

   WALBERLA_LOG_INFO_ON_ROOT(" - Lubrication correction:");
   WALBERLA_LOG_INFO_ON_ROOT("   - normal cut off distance = " << lubricationCorrectionKernel.getNormalCutOffDistance());
   WALBERLA_LOG_INFO_ON_ROOT("   - tangential translational cut off distance = " << lubricationCorrectionKernel.getTangentialTranslationalCutOffDistance());
   WALBERLA_LOG_INFO_ON_ROOT("   - tangential rotational cut off distance = " << lubricationCorrectionKernel.getTangentialRotationalCutOffDistance());
   const real_t maximumLubricationCutOffDistance = std::max(lubricationCorrectionKernel.getNormalCutOffDistance(), std::max(lubricationCorrectionKernel.getTangentialRotationalCutOffDistance(), lubricationCorrectionKernel.getTangentialTranslationalCutOffDistance()));

   std::function<void(void)> particleMappingCall = [ps, accessor, &movingParticleMappingKernel, &fixedParticleMappingKernel, useNoSlipForPlanes]()
   {
      // map planes into the LBM simulation -> act as no-slip boundaries
      if(useNoSlipForPlanes)
      {
         ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), *accessor, fixedParticleMappingKernel, *accessor, NoSlip_Flag);
      }
      else {
         ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);
      }
      // map particles into the LBM simulation
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);
   };

   particleMappingCall();


   //////////////////////////
   // LOAD BALANCING UTILs //
   //////////////////////////

   auto & blockforest = blocks->getBlockForest();
   blockforest.recalculateBlockLevelsInRefresh( false ); // = only load balancing, no refinement checking
   blockforest.alwaysRebalanceInRefresh( true ); //load balancing every time refresh is triggered
   blockforest.reevaluateMinTargetLevelsAfterForcedRefinement( false );
   blockforest.allowRefreshChangingDepth( false );
   blockforest.allowMultipleRefreshCycles( false ); // otherwise info collections are invalid

   shared_ptr<lbm_mesapd_coupling::amr::InfoCollection> couplingInfoCollection = walberla::make_shared<lbm_mesapd_coupling::amr::InfoCollection>();
   uint_t numberOfProcesses = uint_c(MPIManager::instance()->numProcesses());

   if( loadDistributionStrategy == "Hilbert" || loadDistributionStrategy == "Morton")
   {

      bool curveAllGather = true;
      bool balanceLevelwise = true;

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

      if( loadEvaluationStrategy == "Fit" )
      {
         lbm_mesapd_coupling::amr::WeightEvaluationFunctor weightEvaluationFunctor(couplingInfoCollection, lbm_mesapd_coupling::amr::fittedTotalWeightEvaluationFunction);
         lbm_mesapd_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(weightEvaluationFunctor);
         weightAssignmentFunctor.setBlockBaseWeight(blockBaseWeight);
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else if( loadEvaluationStrategy == "LBM" )
      {
         lbm_mesapd_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(lbm_mesapd_coupling::amr::defaultWeightEvaluationFunction);
         weightAssignmentFunctor.setBlockBaseWeight(blockBaseWeight);
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else
      {
         WALBERLA_ABORT("Invalid load evaluation strategy: " << loadEvaluationStrategy);
      }

      blockforest.setRefreshPhantomBlockDataPackFunction(blockforest::PODPhantomWeightPackUnpack<real_t>());
      blockforest.setRefreshPhantomBlockDataUnpackFunction(blockforest::PODPhantomWeightPackUnpack<real_t>());

   }
   else if( loadDistributionStrategy == "ParMetis")
   {

#ifndef WALBERLA_BUILD_WITH_PARMETIS
      WALBERLA_ABORT( "You are trying to use ParMetis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_PARMETIS' to 'ON' in your CMake cache to build against an installed version of ParMetis!" );
#endif

      blockforest::DynamicParMetis::Algorithm parMetisAlgorithm = blockforest::DynamicParMetis::stringToAlgorithm(parMetisAlgorithmString);
      blockforest::DynamicParMetis::WeightsToUse parMetisWeightsToUse = blockforest::DynamicParMetis::WeightsToUse::PARMETIS_BOTH_WEIGHTS;
      blockforest::DynamicParMetis::EdgeSource parMetisEdgeSource = blockforest::DynamicParMetis::EdgeSource::PARMETIS_EDGES_FROM_EDGE_WEIGHTS;

      blockforest::DynamicParMetis dynamicParMetis(parMetisAlgorithm, parMetisWeightsToUse, parMetisEdgeSource);
      dynamicParMetis.setipc2redist(parMetis_ipc2redist);

      auto numberOfBlocks = XBlocks * YBlocks * ZBlocks;

      real_t loadImbalanceTolerance = (parMetisTolerance < real_t(1)) ? std::max(real_t(1.05), real_t(1) + real_t(1) / ( real_c(numberOfBlocks) / real_c(numberOfProcesses) ) ) : parMetisTolerance;
      dynamicParMetis.setImbalanceTolerance(double(loadImbalanceTolerance), 0);

      WALBERLA_LOG_INFO_ON_ROOT(" - ParMetis configuration: ");
      WALBERLA_LOG_INFO_ON_ROOT("   - algorithm = " << dynamicParMetis.algorithmToString() );
      WALBERLA_LOG_INFO_ON_ROOT("   - weights to use = " << dynamicParMetis.weightsToUseToString() );
      WALBERLA_LOG_INFO_ON_ROOT("   - edge source = " << dynamicParMetis.edgeSourceToString() );
      WALBERLA_LOG_INFO_ON_ROOT("   - ipc2redist parameter = " << dynamicParMetis.getipc2redist() );
      WALBERLA_LOG_INFO_ON_ROOT("   - imbalance tolerance = " << dynamicParMetis.getImbalanceTolerance() );

      blockforest.setRefreshPhantomBlockDataPackFunction(blockforest::DynamicParMetisBlockInfoPackUnpack());
      blockforest.setRefreshPhantomBlockDataUnpackFunction(blockforest::DynamicParMetisBlockInfoPackUnpack());
      blockforest.setRefreshPhantomBlockMigrationPreparationFunction( dynamicParMetis );

      if( loadEvaluationStrategy == "Fit" )
      {
         lbm_mesapd_coupling::amr::WeightEvaluationFunctor weightEvaluationFunctor(couplingInfoCollection, lbm_mesapd_coupling::amr::fittedTotalWeightEvaluationFunction);
         lbm_mesapd_coupling::amr::MetisAssignmentFunctor weightAssignmentFunctor(weightEvaluationFunctor); //attention: special METIS assignment functor!
         weightAssignmentFunctor.setBlockBaseWeight(blockBaseWeight);
         real_t weightMultiplicator = real_t(1000); // values from predictor are in range [0-5] which is too coarse when cast to int as done in parmetis
         weightAssignmentFunctor.setWeightMultiplicator(weightMultiplicator);
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else if( loadEvaluationStrategy == "LBM" )
      {
         lbm_mesapd_coupling::amr::MetisAssignmentFunctor weightAssignmentFunctor(lbm_mesapd_coupling::amr::defaultWeightEvaluationFunction);
         weightAssignmentFunctor.setBlockBaseWeight(blockBaseWeight);
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else
      {
         WALBERLA_ABORT("Invalid load evaluation strategy: " << loadEvaluationStrategy);
      }

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
         lbm_mesapd_coupling::amr::WeightEvaluationFunctor weightEvaluationFunctor(couplingInfoCollection, lbm_mesapd_coupling::amr::fittedTotalWeightEvaluationFunction);
         lbm_mesapd_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(weightEvaluationFunctor);
         weightAssignmentFunctor.setBlockBaseWeight(blockBaseWeight);
         blockforest.setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);
      }
      else if( loadEvaluationStrategy == "LBM" )
      {
         lbm_mesapd_coupling::amr::WeightAssignmentFunctor weightAssignmentFunctor(lbm_mesapd_coupling::amr::defaultWeightEvaluationFunction);
         weightAssignmentFunctor.setBlockBaseWeight(blockBaseWeight);
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

   lbm_mesapd_coupling::amr::BlockInfo emptyExampleBlock(blockSize*blockSize*blockSize, 0, 0, 0, 0, numRPDSubCycles);
   WALBERLA_LOG_INFO_ON_ROOT("An example empty block has the weight: " << lbm_mesapd_coupling::amr::fittedTotalWeightEvaluationFunction(emptyExampleBlock));


   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   timeloop.addFuncBeforeTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

   if( vtkWriteFreqPa != uint_t(0) ) {
      auto particleVtkWriter = vtk::createVTKOutput_PointData(particleVtkOutput, "Particles", vtkWriteFreqPa, baseFolder, "simulation_step");
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( particleVtkWriter ), "VTK (sphere data)" );
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

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(pdfFieldVTK), "VTK (fluid field data)");
   }

   if( vtkWriteFreqDD != uint_t(0) ) {
      auto domainDecompVTK = vtk::createVTKOutput_DomainDecomposition(blocks, "domain_decomposition", vtkWriteFreqDD, baseFolder );
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles(domainDecompVTK), "VTK (domain decomposition)");
   }


   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   blockforest::communication::UniformBufferedScheme< Stencil_T > optimizedPDFCommunicationScheme( blocks );
   optimizedPDFCommunicationScheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) ); // optimized sync

   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   timeloop.add() << BeforeFunction( optimizedPDFCommunicationScheme, "LBM Communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   // streaming & collide
   timeloop.add() << Sweep( makeSharedSweep(sweep), "Stream&Collide" );


   SweepTimeloop timeloopAfterParticles( blocks->getBlockStorage(), timesteps );

   // sweep for updating the particle mapping into the LBM simulation
   bool conserveMomentum = false;
   timeloopAfterParticles.add() << Sweep( lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag,
                                                                                                                        lbm_mesapd_coupling::RegularParticlesSelector(), conserveMomentum), "Particle Mapping" );

   bool recomputeTargetDensity = false;
   auto gradReconstructor = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, omega, recomputeTargetDensity,true);
   blockforest::communication::UniformBufferedScheme< Stencil_T > fullPDFCommunicationScheme( blocks );
   fullPDFCommunicationScheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) ); // full sync
   timeloopAfterParticles.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                                << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag,
                                                                                                                                           gradReconstructor, conserveMomentum) ), "PDF Restore" );

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////


   uint_t loadEvaluationFrequency = loadBalancingCheckFrequency;

   // file for simulation infos
   std::string infoFileName( baseFolder + "/simulation_info.txt");
   WALBERLA_ROOT_SECTION()
   {
      std::ofstream file;
      file.open( infoFileName.c_str(), std::fstream::out | std::fstream::trunc );
      file << "#i\t t\t tSim\t tLoadBal\t numProcs\t blocks (min/max/sum)\n";
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
      file << "#i\t t\t mTotSim\t mLB\t mLBM\t mBH\t mCoup1\t mCoup2\t mRPD\t cLBM\t cRB\t numBlocks\n";
      file.close();
   }

   std::string predictionFileProcessName(processLocalFiles + "/predictions_" + std::to_string(MPIManager::instance()->rank()) + ".txt");
   {
      std::ofstream file;
      file.open( predictionFileProcessName.c_str(), std::fstream::out | std::fstream::trunc );
      file << "#i\t t\t wlLBM\t wlBH\t wlCoup1\t wlCoup2\t wlRPD\t edgecut\t numBlocks\n";
      file.close();
   }

   std::vector< std::vector<std::string> > timerKeys;
   std::vector<std::string> LBMTimer;
   LBMTimer.emplace_back("Stream&Collide");
   timerKeys.push_back(LBMTimer);

   std::vector<std::string> bhTimer;
   bhTimer.emplace_back("Boundary Handling");
   timerKeys.push_back(bhTimer);

   std::vector<std::string> couplingTimer1;
   couplingTimer1.emplace_back("Particle Mapping");
   std::vector<std::string> couplingTimer2;
   couplingTimer2.emplace_back("PDF Restore");
   timerKeys.push_back(couplingTimer1);
   timerKeys.push_back(couplingTimer2);

   std::vector<std::string> rpdTimer;
   rpdTimer.emplace_back("RPD Force");
   rpdTimer.emplace_back("RPD VV1");
   rpdTimer.emplace_back("RPD VV2");
   rpdTimer.emplace_back("RPD Lub");
   rpdTimer.emplace_back("RPD Collision");
   timerKeys.push_back(rpdTimer);

   std::vector<std::string> LBMCommTimer;
   LBMCommTimer.emplace_back("LBM Communication");
   LBMCommTimer.emplace_back("PDF Communication");
   timerKeys.push_back(LBMCommTimer);


   std::vector<std::string> rpdCommTimer;
   rpdCommTimer.emplace_back("Reduce Force Torque");
   rpdCommTimer.emplace_back("Reduce Hyd Force Torque");
   rpdCommTimer.emplace_back("Reduce Contact History");
   rpdCommTimer.emplace_back("Sync");
   timerKeys.push_back(rpdCommTimer);

   std::vector<real_t> timings(timerKeys.size());

   double oldmTotSim = 0.0;
   double oldmLB = 0.0;

   auto measurementFileCounter = uint_t(0);
   auto predictionFileCounter = uint_t(0);

   std::string loadEvaluationStep("load evaluation");
   std::string loadBalancingStep("load balancing");
   std::string simulationStep("simulation");

   WcTimingPool timeloopTiming;
   WcTimingPool simulationTiming;

   lbm_mesapd_coupling::InspectionProbe<PdfField_T,BoundaryHandling_T,ParticleAccessor_T> probeForNaNs( Vector3<real_t>(0, 0, 0),
                                                                                                  blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor,
                                                                                                  true, true, "" );

   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {

      // evaluate measurements (note: reflect simulation behavior BEFORE the evaluation)
      if( loadEvaluationFrequency > 0 && i % loadEvaluationFrequency == 0 && i > 0 && fileIO)
      {

         simulationTiming[loadEvaluationStep].start();

         // write process local timing measurements to files (per process, per load balancing step)
         {
            auto mTotSim = simulationTiming[simulationStep].total();
            auto mLB = simulationTiming[loadBalancingStep].total();

            evaluateTimers(timeloopTiming, timerKeys, timings);

            auto & forest = blocks->getBlockForest();
            uint_t numBlocks = forest.getNumberOfBlocks();

            // write to process local file
            std::ofstream file;
            file.open( measurementFileProcessName.c_str(), std::ofstream::app  );
            file << measurementFileCounter << "\t " << real_c(i) / tRef << "\t"
                 << mTotSim - oldmTotSim << "\t" << mLB - oldmLB << "\t";
            for (real_t timing : timings) {
               file << "\t" << timing;
            }
            file << "\t" << numBlocks << "\n";
            file.close();

            oldmTotSim = mTotSim;
            oldmLB = mLB;
            measurementFileCounter++;

            // reset timer to have measurement from evaluation to evaluation point
            timeloopTiming.clear();

         }

         // evaluate general simulation infos (on root)
         {
            double totalTimeToCurrentTimestep(0.0);
            double totalLBTimeToCurrentTimestep(0.0);
            evaluateTotalSimulationTimePassed(simulationTiming, simulationStep, loadBalancingStep,
                                              totalTimeToCurrentTimestep, totalLBTimeToCurrentTimestep);
            math::DistributedSample numberOfBlocks;

            auto & forest = blocks->getBlockForest();
            uint_t numBlocks = forest.getNumberOfBlocks();
            numberOfBlocks.castToRealAndInsert(numBlocks);
            numberOfBlocks.mpiGatherRoot();

            WALBERLA_ROOT_SECTION()
            {
               std::ofstream file;
               file.open( infoFileName.c_str(), std::ofstream::app  );
               file << i << "\t " << real_c(i) / tRef << "\t"
                    << totalTimeToCurrentTimestep << "\t " << totalLBTimeToCurrentTimestep << "\t "
                    << numberOfProcesses << "\t ";
               file << uint_c(numberOfBlocks.min()) << "\t ";
               file << uint_c(numberOfBlocks.max()) << "\t ";
               file << uint_c(numberOfBlocks.sum()) << "\n ";
               file.close();
            }
         }

         simulationTiming[loadEvaluationStep].end();

      }

      if( loadBalancingCheckFrequency != 0 && i % loadBalancingCheckFrequency == 0)
      {

         if(useProgressLogging) walberla::logging::Logging::instance()->setFileLogLevel(logging::Logging::LogLevel::PROGRESS);

         simulationTiming[loadBalancingStep].start();

         if( loadEvaluationStrategy != "LBM" ) {

            WALBERLA_LOG_INFO_ON_ROOT("Checking for load balancing...");

            // update info collections for the particle presence based check and the load balancing:
            auto &forest = blocks->getBlockForest();
            lbm_mesapd_coupling::amr::updateAndSyncInfoCollection<BoundaryHandling_T,ParticleAccessor_T >(forest, boundaryHandlingID, *accessor, numRPDSubCycles, *couplingInfoCollection);

            if(useProgressLogging)  WALBERLA_LOG_INFO("Created info collection with size " << couplingInfoCollection->size());

            auto numBlocksBefore = forest.getNumberOfBlocks();
            if(useProgressLogging)  WALBERLA_LOG_INFO("Number of blocks before refresh " << numBlocksBefore);

            mesa_pd::mpi::ClearNextNeighborSync CNNS;
            CNNS(*accessor);

            if(useProgressLogging) WALBERLA_LOG_INFO_ON_ROOT("Cleared particle sync and starting refresh");

            blocks->refresh();

            auto numBlocksAfter = forest.getNumberOfBlocks();
            if(useProgressLogging)  WALBERLA_LOG_INFO("Number of blocks after refresh " << numBlocksAfter);

            if(useProgressLogging) WALBERLA_LOG_INFO_ON_ROOT("Refresh finished, recreating all datastructures");

            //WALBERLA_LOG_INFO(rank << ", " << numBlocksBefore << " -> " << numBlocksAfter);
            rpdDomain->refresh();
            if(useBlockForestSync)
            {
               ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, associateToBlock, *accessor);
               syncNextNeighborBlockForestFunc(*ps, blocks->getBlockForestPointer(), rpdDomain, overlap);

            } else
            {
               syncNextNeighborFunc(*ps, *rpdDomain, overlap);
            }

            clearBoundaryHandling(forest, boundaryHandlingID);
            clearParticleField(forest, particleFieldID, *accessor);

            recreateBoundaryHandling(forest, boundaryHandlingID);

            //NOTE: order in mapping is important: first all global/fixed particles,
            // only then moving ones, which do not overwrite the mapping of the global/fixed ones
            particleMappingCall();

         }
         simulationTiming[loadBalancingStep].end();

         if(useProgressLogging) walberla::logging::Logging::instance()->setFileLogLevel(logging::Logging::LogLevel::INFO);

      }

      if(checkSimulation)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Checking time step " << i );
         evaluateParticleSimulation(ps, ss, numberOfSediments, numGlobalParticles);
         checkMapping(blocks, pdfFieldID, boundaryHandlingID, probeForNaNs);
      }

      // evaluate predictions (note: reflect the predictions for all upcoming simulations, thus the corresponding measurements have to be taken afterwards)
      if( loadEvaluationFrequency > 0 && i % loadEvaluationFrequency == 0 && fileIO)
      {

         simulationTiming[loadEvaluationStep].start();

         // write process local load predictions to files (per process, per load balancing step)
         {

            auto wlLBM = real_t(0);
            auto wlBH = real_t(0);
            auto wlCoup1 = real_t(0);
            auto wlCoup2 = real_t(0);
            auto wlRPD = real_t(0);

            auto & forest = blocks->getBlockForest();
            lbm_mesapd_coupling::amr::updateAndSyncInfoCollection<BoundaryHandling_T,ParticleAccessor_T >(forest, boundaryHandlingID, *accessor, numRPDSubCycles, *couplingInfoCollection);

            for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt ) {
               auto * block = static_cast<blockforest::Block *> (&(*blockIt));
               const auto &blockID = block->getId();
               auto infoIt = couplingInfoCollection->find(blockID);
               auto blockInfo = infoIt->second;

               wlLBM   += lbm_mesapd_coupling::amr::fittedLBMWeightEvaluationFunction(blockInfo);
               wlBH    += lbm_mesapd_coupling::amr::fittedBHWeightEvaluationFunction(blockInfo);
               wlCoup1 += lbm_mesapd_coupling::amr::fittedCoup1WeightEvaluationFunction(blockInfo);
               wlCoup2 += lbm_mesapd_coupling::amr::fittedCoup2WeightEvaluationFunction(blockInfo);
               wlRPD   += lbm_mesapd_coupling::amr::fittedRPDWeightEvaluationFunction(blockInfo);

            }

            // note: we count the edge weight doubled here in total (to and from the other process). ParMetis only counts one direction.
            uint_t edgecut = evaluateEdgeCut(forest);
            uint_t numBlocks = forest.getNumberOfBlocks();

            std::ofstream file;
            file.open( predictionFileProcessName.c_str(), std::ofstream::app  );
            file << predictionFileCounter << "\t " << real_c(i) / tRef << "\t"
                 << wlLBM << "\t" << wlBH << "\t" << wlCoup1 << "\t" << wlCoup2 << "\t" << wlRPD << "\t"
                 << edgecut << "\t" << numBlocks << "\n";
            file.close();

            predictionFileCounter++;;
         }

         simulationTiming[loadEvaluationStep].end();

      }

      simulationTiming[simulationStep].start();

      // perform a single simulation step

      timeloop.singleStep( timeloopTiming );

      timeloopTiming["Reduce Hyd Force Torque"].start();
      reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);
      timeloopTiming["Reduce Hyd Force Torque"].end();

      timeloopTiming["RPD Force"].start();
      if( i == 0 )
      {
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, initializeHydrodynamicForceTorqueForAveragingKernel, *accessor );
      }
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque, *accessor );
      timeloopTiming["RPD Force"].end();

      if(checkSimulation)
      {
         checkParticleProperties(ps);
      }

      for(auto subCycle = uint_t(0); subCycle < numRPDSubCycles; ++subCycle )
      {

         timeloopTiming["RPD VV1"].start();
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
         timeloopTiming["RPD VV1"].end();

         timeloopTiming["Sync"].start();
         if(useBlockForestSync)
         {
            syncNextNeighborBlockForestFunc(*ps, blocks->getBlockForestPointer(), rpdDomain, overlap);
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, associateToBlock, *accessor);
         } else
         {
            syncNextNeighborFunc(*ps, *rpdDomain, overlap);
         }
         timeloopTiming["Sync"].end();

         // lubrication correction
         timeloopTiming["RPD Lub"].start();
         ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                     [&lubricationCorrectionKernel,maximumLubricationCutOffDistance, &rpdDomain]
                                           (const size_t idx1, const size_t idx2, auto& ac)
                                     {
                                        mesa_pd::collision_detection::AnalyticContactDetection acd;
                                        acd.getContactThreshold() = maximumLubricationCutOffDistance;
                                        mesa_pd::kernel::DoubleCast double_cast;
                                        mesa_pd::mpi::ContactFilter contact_filter;
                                        if (double_cast(idx1, idx2, ac, acd, ac ))
                                        {
                                           if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                           {
                                              double_cast(acd.getIdx1(), acd.getIdx2(), ac, lubricationCorrectionKernel, ac, acd.getContactNormal(), acd.getPenetrationDepth());
                                           }
                                        }
                                     },
                                     *accessor );
         timeloopTiming["RPD Lub"].end();

         // collision response
         timeloopTiming["RPD Collision"].start();
         ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                     [&collisionResponse, &rpdDomain, timeStepSizeRPD]
                                           (const size_t idx1, const size_t idx2, auto& ac)
                                     {
                                        mesa_pd::collision_detection::AnalyticContactDetection acd;
                                        mesa_pd::kernel::DoubleCast double_cast;
                                        mesa_pd::mpi::ContactFilter contact_filter;
                                        if (double_cast(idx1, idx2, ac, acd, ac ))
                                        {
                                           if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                           {
                                              collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth(), timeStepSizeRPD);
                                           }
                                        }
                                     },
                                     *accessor );

         timeloopTiming["RPD Collision"].end();

         timeloopTiming["Reduce Contact History"].start();
         reduceAndSwapContactHistory(*ps);
         timeloopTiming["Reduce Contact History"].end();

         timeloopTiming["RPD Force"].start();
         lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction, *accessor );
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor );
         timeloopTiming["RPD Force"].end();

         timeloopTiming["Reduce Force Torque"].start();
         reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);
         timeloopTiming["Reduce Force Torque"].end();

         // integration
         timeloopTiming["RPD VV2"].start();
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
         timeloopTiming["RPD VV2"].end();

         timeloopTiming["Sync"].start();
         if(useBlockForestSync)
         {
            syncNextNeighborBlockForestFunc(*ps, blocks->getBlockForestPointer(), rpdDomain, overlap);
         } else
         {
            syncNextNeighborFunc(*ps, *rpdDomain, overlap);
         }
         timeloopTiming["Sync"].end();


      }

      timeloopTiming["RPD Force"].start();
      ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );
      timeloopTiming["RPD Force"].end();

      // update particle mapping
      timeloopAfterParticles.singleStep(timeloopTiming);

      simulationTiming[simulationStep].end();

   }

   simulationTiming.logResultOnRoot();


   return EXIT_SUCCESS;
}

} // namespace fluid_particle_workload_distribution

int main( int argc, char **argv ){
   fluid_particle_workload_distribution::main(argc, argv);
}
