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
//! \file WorkLoadEvaluation.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

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
#include "core/mpi/Broadcast.h"

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
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"
#include "lbm/BlockForestEvaluation.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/basic.h"
#include "pe/cr/ICR.h"
#include "pe/fcd/GJKEPACollideFunctor.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/EllipsoidVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <vector>
#include <iomanip>
#include <iostream>

namespace workload_evaluation
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

// PDF field, flag field & body field
using LatticeModel_T = lbm::D3Q19<lbm::collision_model::TRT, false>;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
using BodyField_T = GhostLayerField<pe::BodyID, 1>;

const uint_t FieldGhostLayers = 1;

// boundary handling
using MO_CLI_T = pe_coupling::CurvedLinear<LatticeModel_T, FlagField_T>;

using BoundaryHandling_T = BoundaryHandling<FlagField_T, Stencil_T, MO_CLI_T>;

using BodyTypeTuple = std::tuple<pe::Sphere, pe::Ellipsoid, pe::Plane>;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );
const FlagUID MO_CLI_Flag  ( "moving obstacle CLI" );
const FlagUID FormerMO_Flag( "former moving obstacle" );


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyFieldID, bool useEntireFieldTraversal ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), bodyFieldID_ ( bodyFieldID ), useEntireFieldTraversal_( useEntireFieldTraversal ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;
   bool useEntireFieldTraversal_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   auto * flagField = block->getData< FlagField_T >( flagFieldID_ );
   auto *  pdfField = block->getData< PdfField_T > ( pdfFieldID_ );
   auto * bodyField = block->getData< BodyField_T >( bodyFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T::Mode mode = (useEntireFieldTraversal_) ? BoundaryHandling_T::Mode::ENTIRE_FIELD_TRAVERSAL : BoundaryHandling_T::Mode::OPTIMIZED_SPARSE_TRAVERSAL ;

   BoundaryHandling_T * handling = new BoundaryHandling_T( "fixed obstacle boundary handling", flagField, fluid,
                                                           MO_CLI_T ( "MO_CLI", MO_CLI_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ),
                                                           mode);

   // boundaries are set by mapping the planes into the domain

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}


class CollisionPropertiesEvaluator
{
public:
   explicit CollisionPropertiesEvaluator( pe::cr::ICR & collisionResponse ) : collisionResponse_( collisionResponse )
   {}

   real_t get()
   {
      real_t maxPen = std::fabs(collisionResponse_.getMaximumPenetration());
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( maxPen, mpi::MAX );
      }
      return maxPen;
   }

private:
   pe::cr::ICR & collisionResponse_;
};

class ContactDistanceEvaluator
{
public:
   ContactDistanceEvaluator( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID ccdID, const BlockDataID fcdID ) :
   blocks_( blocks ), ccdID_(ccdID), fcdID_(fcdID)
   {}

   real_t get()
   {
      auto maximumPenetration = real_t(0);
      for (auto it = blocks_->begin(); it != blocks_->end(); ++it) {
         IBlock &currentBlock = *it;

         auto *ccd = currentBlock.getData<pe::ccd::ICCD>(ccdID_);
         auto *fcd = currentBlock.getData<pe::fcd::IFCD>(fcdID_);
         ccd->generatePossibleContacts();
         pe::Contacts& contacts = fcd->generateContacts( ccd->getPossibleContacts() );
         size_t numContacts( contacts.size() );

         for( size_t i = 0; i < numContacts; ++i )
         {
            const pe::ContactID c( &contacts[i] );
            maximumPenetration = std::max( maximumPenetration, std::fabs(c->getDistance()));
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( maximumPenetration, mpi::MAX );
      }
      return maximumPenetration;
   }
private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID ccdID_;
   const BlockDataID fcdID_;
};

class MaxVelocityEvaluator
{
public:
   MaxVelocityEvaluator( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID bodyStorageID ) :
         blocks_( blocks ), bodyStorageID_(bodyStorageID)
   {}

   Vector3<real_t> get()
   {
      auto maxVelX = real_t(0);
      auto maxVelY = real_t(0);
      auto maxVelZ = real_t(0);

      for (auto it = blocks_->begin(); it != blocks_->end(); ++it) {

         for( auto bodyIt = pe::LocalBodyIterator::begin(*it, bodyStorageID_ ); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt ) {
            auto vel = bodyIt->getLinearVel();
            maxVelX = std::max(maxVelX, std::fabs(vel[0]));
            maxVelY = std::max(maxVelY, std::fabs(vel[1]));
            maxVelZ = std::max(maxVelZ, std::fabs(vel[2]));
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( maxVelX, mpi::MAX );
         mpi::allReduceInplace( maxVelY, mpi::MAX );
         mpi::allReduceInplace( maxVelZ, mpi::MAX );
      }
      return Vector3<real_t>(maxVelX, maxVelY, maxVelZ);
   }

   real_t getMagnitude()
   {
      auto magnitude = real_t(0);

      for (auto it = blocks_->begin(); it != blocks_->end(); ++it) {

         for( auto bodyIt = pe::LocalBodyIterator::begin(*it, bodyStorageID_ ); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt ) {
            magnitude = std::max(magnitude, bodyIt->getLinearVel().length());
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( magnitude, mpi::MAX );
      }
      return magnitude;

   }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;
};

void evaluateFluidQuantities(const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID boundaryHandlingID,
                             uint_t & numCells, uint_t & numFluidCells, uint_t & numNBCells )
{

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID );
      auto xyzSize = boundaryHandling->getFlagField()->xyzSize();
      numCells += xyzSize.numCells();

      for(auto z = cell_idx_t(xyzSize.zMin()); z <= cell_idx_t(xyzSize.zMax()); ++z ){
         for(auto y = cell_idx_t(xyzSize.yMin()); y <= cell_idx_t(xyzSize.yMax()); ++y ){
            for(auto x = cell_idx_t(xyzSize.xMin()); x <= cell_idx_t(xyzSize.xMax()); ++x ) {
               if (boundaryHandling->isDomain(x, y, z)) {
                  ++numFluidCells;
               }
               if (boundaryHandling->isNearBoundary(x, y, z)) {
                  ++numNBCells;
               }
            }
         }
      }
   }
}

void evaluatePEQuantities( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID bodyStorageID,
                           const pe::cr::ICR & cr,
                           uint_t & numLocalParticles, uint_t & numShadowParticles, uint_t & numContacts)
{

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
      auto * bodyStorage = blockIt->getData<pe::Storage>(bodyStorageID);
      pe::BodyStorage const & localStorage  = (*bodyStorage)[pe::StorageType::LOCAL];
      pe::BodyStorage const & shadowStorage = (*bodyStorage)[pe::StorageType::SHADOW];
      numLocalParticles += localStorage.size();
      numShadowParticles += shadowStorage.size();

      numContacts += cr.getNumberOfContactsTreated();
   }
}

void evaluateTimers(WcTimingPool & timingPool, WcTimingTree & peTimingTree,
                    const std::vector<std::vector<std::string> > & timerKeys,
                    std::vector<real_t> & timings )
{

   for (auto & timingsIt : timings)
   {
      timingsIt = real_t(0);
   }

   timingPool.unifyRegisteredTimersAcrossProcesses();
   peTimingTree.synchronize();

   auto scalingFactor = real_t(1000); // milliseconds

   for (auto i = uint_t(0); i < timerKeys.size(); ++i )
   {
      auto keys = timerKeys[i];
      for (const auto &timerName : keys)
      {
         if(timingPool.timerExists(timerName))
         {
            timings[i] += real_c(timingPool[timerName].total()) * scalingFactor;
         }
         if(peTimingTree.timerExists(timerName))
         {
            timings[i] += real_c(peTimingTree[timerName].total()) * scalingFactor;
         }
      }

   }
}

void resetTimers(WcTimingPool & timingPool, WcTimingTree & peTimingTree)
{
   timingPool.clear();
   peTimingTree.reset();
}

//*******************************************************************************************************************
/*! Application to evaluate the workload (time measurements) for a fluid-particle simulation
 *
 * This application is used in the paper
 *  Rettinger, Ruede - "Dynamic Load Balancing Techniques for Particulate Flow Simulations", submitted to Computation
 * in Section 3 to develop and calibrate the workload estimator.
 * The setup features settling particle inside a horizontally periodic box.
 * A comprehensive description is given in Sec. 3.3 of the paper.
 * It uses 4 x 4 x 5 blocks for domain partitioning.
 * For each block ( = each process), the block local quantities are evaluated as well as the timing infos of
 * the fluid-particle interaction algorithm. Those infos are then written to files that can be used later on
 * for function fitting to obtain a workload estimator.
 *
 * NOTE: Since this estimator relies on timing measurements, this evaluation procedure should be carried out everytime
 * a different implementation, hardware or algorithm is used.
 *
 */
//*******************************************************************************************************************
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );


   auto solidVolumeFraction = real_t(0.2);

   // LBM / numerical parameters
   auto blockSize  = uint_t(32);
   auto uSettling = real_t(0.1); // characteristic settling velocity
   auto diameter = real_t(10);

   auto Ga = real_t(30); //Galileo number
   auto numPeSubCycles = uint_t(10);

   auto vtkIOFreq = uint_t(0);
   auto timestepsNonDim = real_t(2.5);
   auto numSamples = uint_t(2000);
   std::string baseFolder = "workload_files"; // folder for vtk and file output

   bool useEllipsoids = false;
   bool optimizeForSmallObstacleFraction = false;
   bool noFileOutput = false;
   bool fixBodies = false;
   bool useEntireFieldTraversal = true;
   bool averageForceTorqueOverTwoTimSteps = true;
   bool useFusedStreamCollide = false;


   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--vtkIOFreq"               ) == 0 ) { vtkIOFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--noFileOutput"            ) == 0 ) { noFileOutput = true; continue; }
      if( std::strcmp( argv[i], "--basefolder"              ) == 0 ) { baseFolder = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--solidVolumeFraction"     ) == 0 ) { solidVolumeFraction = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--diameter"                ) == 0 ) { diameter = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--blockSize"               ) == 0 ) { blockSize = uint_c(std::atof( argv[++i]) ); continue; }
      if( std::strcmp( argv[i], "--uSettling"               ) == 0 ) { uSettling = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--Ga"                      ) == 0 ) { Ga = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--timestepsNonDim"         ) == 0 ) { timestepsNonDim = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--numPeSubCycles"          ) == 0 ) { numPeSubCycles = uint_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--useEllipsoids"           ) == 0 ) { useEllipsoids = true; continue; }
      if( std::strcmp( argv[i], "--optSmallSVF"             ) == 0 ) { optimizeForSmallObstacleFraction = true; continue; }
      if( std::strcmp( argv[i], "--fixBodies"               ) == 0 ) { fixBodies = true; continue; }
      if( std::strcmp( argv[i], "--useEntireFieldTraversal" ) == 0 ) { useEntireFieldTraversal = true; continue; }
      if( std::strcmp( argv[i], "--numSamples"              ) == 0 ) { numSamples = uint_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--noForceAveraging"        ) == 0 ) { averageForceTorqueOverTwoTimSteps = false; continue; }
      if( std::strcmp( argv[i], "--useFusedStreamCollide"   ) == 0 ) { useFusedStreamCollide = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   WALBERLA_CHECK(diameter > real_t(1));
   WALBERLA_CHECK(uSettling > real_t(0));
   WALBERLA_CHECK(Ga > real_t(0));
   WALBERLA_CHECK(solidVolumeFraction > real_t(0));
   WALBERLA_CHECK(solidVolumeFraction < real_t(0.65));

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   const auto XBlocks = uint_t(4);
   const auto YBlocks = uint_t(4);
   const auto ZBlocks = uint_t(5);

   if( MPIManager::instance()->numProcesses() != XBlocks * YBlocks * ZBlocks )
   {
      WALBERLA_LOG_WARNING_ON_ROOT("WARNING! You have specified less or more processes than number of blocks -> the time measurements are no longer blockwise!")
   }

   if( diameter > real_c(blockSize) )
   {
      WALBERLA_LOG_WARNING_ON_ROOT("The bodies might be too large to work with the currently used synchronization!");
   }


   WALBERLA_LOG_INFO_ON_ROOT("Using setup with sedimenting particles -> creating two planes and applying gravitational force")
   if( useEllipsoids ){ WALBERLA_LOG_INFO_ON_ROOT("using ELLIPSOIDS"); }
   else{ WALBERLA_LOG_INFO_ON_ROOT("using SPHERES"); }


   const uint_t XCells = blockSize * XBlocks;
   const uint_t YCells = blockSize * YBlocks;
   const uint_t ZCells = blockSize * ZBlocks;

   const real_t topWallOffset = real_t(1.05) * real_t(blockSize); // move the top wall downwards to take away a certain portion of the overall domain


   // determine number of spheres to generate, if necessary scale diameter a bit to reach desired solid volume fraction
   real_t domainHeight = real_c(ZCells) - topWallOffset;
   real_t fluidVolume =  real_c( XCells * YCells ) * domainHeight;
   real_t solidVolume = solidVolumeFraction * fluidVolume;
   uint_t numberOfParticles = uint_c(std::ceil(solidVolume / ( math::pi / real_t(6) * diameter * diameter * diameter )));
   diameter = std::cbrt( solidVolume / ( real_c(numberOfParticles) * math::pi / real_t(6) ) );

   auto densityRatio = real_t(2.5);

   real_t viscosity = uSettling * diameter / Ga;
   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);

   const real_t gravitationalAcceleration = uSettling * uSettling / ( (densityRatio-real_t(1)) * diameter );

   real_t tref = diameter / uSettling;
   real_t Tref = domainHeight / uSettling;

   uint_t timesteps = uint_c(timestepsNonDim * Tref);

   const real_t dx = real_c(1);
   WALBERLA_LOG_INFO_ON_ROOT("viscosity = " << viscosity);
   WALBERLA_LOG_INFO_ON_ROOT("tau = " << real_t(1)/omega);
   WALBERLA_LOG_INFO_ON_ROOT("diameter = " << diameter);
   WALBERLA_LOG_INFO_ON_ROOT("solid volume fraction = " << solidVolumeFraction);
   WALBERLA_LOG_INFO_ON_ROOT("domain size (in cells) = " << XCells << " x " << ZCells << " x " << ZCells);
   WALBERLA_LOG_INFO_ON_ROOT("number of bodies = " << numberOfParticles);
   WALBERLA_LOG_INFO_ON_ROOT("gravitational acceleration = " << gravitationalAcceleration);
   WALBERLA_LOG_INFO_ON_ROOT("Ga = " << Ga);
   WALBERLA_LOG_INFO_ON_ROOT("uSettling = " << uSettling);
   WALBERLA_LOG_INFO_ON_ROOT("tref = " << tref);
   WALBERLA_LOG_INFO_ON_ROOT("Tref = " << Tref);
   WALBERLA_LOG_INFO_ON_ROOT("timesteps = " << timesteps);
   WALBERLA_LOG_INFO_ON_ROOT("number of workload samples = " << numSamples);

   // create folder to store logging files
   WALBERLA_ROOT_SECTION()
   {
      walberla::filesystem::path path1( baseFolder );
      if( !walberla::filesystem::exists( path1 ) )
         walberla::filesystem::create_directory( path1 );
   }


   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3<bool> periodicity( true );
   periodicity[2] = false;

   // create domain
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid( XBlocks, YBlocks, ZBlocks, blockSize, blockSize, blockSize, dx,
                                                    0, false, false, //one block per process!
                                                    periodicity[0], periodicity[1], periodicity[2], //periodicity
                                                    false );

   ////////
   // PE //
   ////////

   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   BlockDataID fcdID   = (useEllipsoids) ? blocks->addBlockData( pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::GJKEPACollideFunctor>(), "FCD" )
                                         : blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   WcTimingTree timingTreePE;

   pe::cr::HCSITS cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, &timingTreePE );
   cr.setMaxIterations(10);
   cr.setRelaxationModel( pe::cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling );
   cr.setErrorReductionParameter(real_t(0.8));

   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe
   const real_t overlap = real_c( 1.5 ) * dx;

   std::function<void(void)> syncCall = std::bind( pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, &timingTreePE, overlap, false );

   auto generationDomain = AABB( real_t(0), real_t(0), real_t(0), real_c(XCells), real_c(YCells), real_c(ZCells) - topWallOffset);
   auto peMaterial = pe::createMaterial( "mat", densityRatio, real_t(1), real_t(0.25), real_t(0.25), real_t(0), real_t(200), real_t(100), real_t(100), real_t(100) );

   // create two planes at bottom and top of domain
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,1), Vector3<real_t>(0,0,0), peMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,real_c(ZCells)-topWallOffset), peMaterial );

   auto xParticle = real_t(0);
   auto yParticle = real_t(0);
   auto zParticle = real_t(0);

   for( uint_t nPart = 0; nPart < numberOfParticles; ++nPart )
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

      } else{
         pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( xParticle, yParticle, zParticle ), diameter * real_t(0.5), peMaterial );
      }

   }

   syncCall();

   // resolve possible overlaps of the particles due to the random initialization

   // 100 iterations of solver to resolve all major overlaps
   {
      for (auto pet = uint_t(1); pet <= uint_t(100); ++pet )
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
   }


   // resolve remaining overlaps via particle simulation
   {
      const auto initialPeSteps = uint_t(2000);
      const auto dt_PE_init = real_t(1);
      const real_t overlapLimit = real_t(0.001) * diameter;

      WALBERLA_LOG_INFO_ON_ROOT("Particle creation done --- resolving overlaps with goal all < " << overlapLimit / diameter * real_t(100) << "%");

      CollisionPropertiesEvaluator collisionPropertiesEvaluator( cr );

      ContactDistanceEvaluator contactDistanceEvaluator(blocks, ccdID, fcdID);
      MaxVelocityEvaluator maxVelEvaluator(blocks, bodyStorageID);

      for(auto pet = uint_t(1); pet <= initialPeSteps; ++pet )
      {
         cr.timestep( dt_PE_init );
         syncCall();
         real_t maxPen = collisionPropertiesEvaluator.get();

         if( maxPen < overlapLimit )
         {
            WALBERLA_LOG_INFO_ON_ROOT("Carried out " << pet << " PE-only time steps to resolve initial overlaps");
            WALBERLA_LOG_INFO_ON_ROOT("Final max penetration from cr is " << maxPen << " = " << maxPen / diameter * real_t(100) << "%");

            break;
         }

         real_t maxMagnitude = maxVelEvaluator.getMagnitude();

         if( maxMagnitude * dt_PE_init > overlapLimit)
         {
            // avoid too large response velocities by setting them to zero
            for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            {
               for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
               {
                  bodyIt->setLinearVel(Vector3<real_t>(real_t(0)));
                  bodyIt->setAngularVel(Vector3<real_t>(real_t(0)));
               }
            }
         }
         else
         {
            cr.setErrorReductionParameter(real_t(0.8));
         }

         if( pet % uint_t(20) == uint_t(0) )
         {
            WALBERLA_LOG_INFO_ON_ROOT(pet << " - current max overlap = " << maxPen / diameter * real_t(100) << "%, max vel magnitude = " << maxMagnitude );
         }
      }
   }

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

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         uint_t(1), field::fzyx );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::fzyx );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
                                    MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID, useEntireFieldTraversal ), "boundary handling" );


   // initially map pe bodies into the LBM simulation
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_CLI_Flag );

   lbm::BlockForestEvaluation<FlagField_T> bfEval(blocks, flagFieldID, Fluid_Flag);

   WALBERLA_LOG_INFO_ON_ROOT(bfEval.loggingString());

   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   if( vtkIOFreq != uint_t(0) )
   {
      // pe bodies
      if(useEllipsoids)
      {
         auto bodyVtkOutput = make_shared<pe::EllipsoidVtkOutput>( bodyStorageID, blocks->getBlockStorage() );
         auto bodyVTK = vtk::createVTKOutput_PointData( bodyVtkOutput, "bodies", vtkIOFreq, baseFolder );
         timeloop.addFuncBeforeTimeStep( vtk::writeFiles( bodyVTK ), "VTK (body data)" );

      }else
      {
         auto bodyVtkOutput = make_shared<pe::SphereVtkOutput>( bodyStorageID, blocks->getBlockStorage() );
         auto bodyVTK = vtk::createVTKOutput_PointData( bodyVtkOutput, "bodies", vtkIOFreq, baseFolder );
         timeloop.addFuncBeforeTimeStep( vtk::writeFiles( bodyVTK ), "VTK (body data)" );
      }


      // flag field
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", vtkIOFreq, 1, false, baseFolder );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      // pdf field
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", vtkIOFreq, 0, false, baseFolder );

      field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
      fluidFilter.addFlag( Fluid_Flag );
      pdfFieldVTK->addCellInclusionFilter( fluidFilter );

      pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( pdfFieldVTK ), "VTK (fluid field data)" );

      auto domainDecompVTK = vtk::createVTKOutput_DomainDecomposition(blocks, "domain_decomposition", vtkIOFreq, baseFolder );
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles(domainDecompVTK), "VTK (domain decomposition)");
   }

   // sweep for updating the pe body mapping into the LBM simulation
   timeloop.add()
         << Sweep( pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MO_CLI_Flag, FormerMO_Flag, pe_coupling::selectRegularBodies ), "Body Mapping" );

   // sweep for restoring PDFs in cells previously occupied by pe bodies
   using Reconstructor_T = pe_coupling::EquilibriumReconstructor<LatticeModel_T, BoundaryHandling_T>;
   Reconstructor_T reconstructor( blocks, boundaryHandlingID, bodyFieldID);
   timeloop.add()
         << Sweep( pe_coupling::PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructor_T >
                         ( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, reconstructor, FormerMO_Flag, Fluid_Flag, pe_coupling::selectRegularBodies, optimizeForSmallObstacleFraction ), "PDF Restore" );

   shared_ptr<pe_coupling::BodiesForceTorqueContainer> bodiesFTContainer1 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
   std::function<void(void)> storeForceTorqueInCont1 = std::bind(&pe_coupling::BodiesForceTorqueContainer::store, bodiesFTContainer1);
   shared_ptr<pe_coupling::BodiesForceTorqueContainer> bodiesFTContainer2 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
   std::function<void(void)> setForceTorqueOnBodiesFromCont2 = std::bind(&pe_coupling::BodiesForceTorqueContainer::setOnBodies, bodiesFTContainer2);
   shared_ptr<pe_coupling::ForceTorqueOnBodiesScaler> forceScaler = make_shared<pe_coupling::ForceTorqueOnBodiesScaler>(blocks, bodyStorageID, real_t(1));
   std::function<void(void)> setForceScalingFactorToHalf = std::bind(&pe_coupling::ForceTorqueOnBodiesScaler::resetScalingFactor,forceScaler,real_t(0.5));

   if( averageForceTorqueOverTwoTimSteps ) {
      bodiesFTContainer2->store();
   }


   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > scheme( blocks );
   scheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) );
   commFunction = scheme;

   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );

   if( !useFusedStreamCollide )
   {
      // streaming & collide
      timeloop.add() << Sweep( makeCollideSweep(sweep), "Collide" );
   }

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   if( useFusedStreamCollide )
   {
      // streaming & collide
      timeloop.add() << Sweep( makeSharedSweep(sweep), "Stream&Collide" );
   } else
   {
      // streaming & collide
      timeloop.add() << Sweep( makeStreamSweep(sweep), "Stream" );

   }

   // Averaging the force/torque over two time steps is said to damp oscillations of the interaction force/torque.
   // See Ladd - " Numerical simulations of particulate suspensions via a discretized Boltzmann equation. Part 1. Theoretical foundation", 1994, p. 302
   if( averageForceTorqueOverTwoTimSteps ) {

      // store force/torque from hydrodynamic interactions in container1
      timeloop.addFuncAfterTimeStep(storeForceTorqueInCont1, "Force Storing");

      // set force/torque from previous time step (in container2)
      timeloop.addFuncAfterTimeStep(setForceTorqueOnBodiesFromCont2, "Force setting");

      // average the force/torque by scaling it with factor 1/2 (except in first timestep, there it is 1, which it is initially)
      timeloop.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesScaler(blocks, bodyStorageID, real_t(0.5)),  "Force averaging");
      timeloop.addFuncAfterTimeStep( setForceScalingFactorToHalf, "Force scaling adjustment" );

      // swap containers
      timeloop.addFuncAfterTimeStep( pe_coupling::BodyContainerSwapper( bodiesFTContainer1, bodiesFTContainer2 ), "Swap FT container" );

   }

   real_t sphereVolume = diameter * diameter * diameter * math::pi / real_t(6);
   Vector3<real_t> gravitationalForce( real_t(0), real_t(0), -(densityRatio - real_t(1)) * gravitationalAcceleration * sphereVolume );
   timeloop.addFuncAfterTimeStep(pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, gravitationalForce ), "Gravitational force" );

   if( fixBodies ) {
      // reset all forces
      timeloop.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesResetter(blocks, bodyStorageID), "Force Resetting");
   } else{
      // add pe timesteps
      timeloop.addFuncAfterTimeStep( pe_coupling::TimeStep( blocks, bodyStorageID, cr, syncCall, real_t(1), numPeSubCycles ), "pe Time Step" );
   }

   timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );



   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   std::vector< std::vector<std::string> > timerKeys;
   std::vector<std::string> LBMTimer;
   LBMTimer.emplace_back("Stream&Collide");
   LBMTimer.emplace_back("Stream");
   LBMTimer.emplace_back("Collide");
   timerKeys.push_back(LBMTimer);

   std::vector<std::string> bhTimer;
   bhTimer.emplace_back("Boundary Handling");
   timerKeys.push_back(bhTimer);

   std::vector<std::string> couplingTimer1;
   couplingTimer1.emplace_back("Body Mapping");
   std::vector<std::string> couplingTimer2;
   couplingTimer2.emplace_back("PDF Restore");
   timerKeys.push_back(couplingTimer1);
   timerKeys.push_back(couplingTimer2);

   std::vector<std::string> peTimer;
   peTimer.emplace_back("Simulation Step.Collision Detection");
   peTimer.emplace_back("Simulation Step.Collision Response Integration");
   peTimer.emplace_back("Simulation Step.Collision Response Resolution.Collision Response Solving");
   timerKeys.push_back(peTimer);

   uint_t numCells;
   uint_t numFluidCells;
   uint_t numNBCells;
   uint_t numLocalParticles;
   uint_t numShadowParticles;
   uint_t numContacts;
   numCells = uint_t(0);
   numFluidCells = uint_t(0);
   numNBCells = uint_t(0);
   numLocalParticles = uint_t(0);
   numShadowParticles = uint_t(0);
   numContacts = uint_t(0);

   std::vector<real_t> timings(timerKeys.size());

   resetTimers(timeloopTiming, timingTreePE);

   // every rank writes its own file -> numProcesses number of samples!
   int myRank = MPIManager::instance()->rank();

   std::string logFileName = baseFolder + "/load";
   logFileName += "_settling";
   if( useEllipsoids)
   {
      logFileName += "_ellipsoids";
   }
   else
   {
      logFileName += "_spheres";
   }
   logFileName += "_d" + std::to_string(int_c(std::ceil(diameter)));
   logFileName += "_bs" + std::to_string(blockSize);
   logFileName += "_" + std::to_string(myRank) + ".txt";


   std::ofstream file;

   if(!noFileOutput)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Writing load info to file " << logFileName);
      file.open( logFileName.c_str(), std::ofstream::app );
      file << "# svf = " << solidVolumeFraction << ", d = " << diameter << ", domain = " << XCells << "x" << YCells << "x" << ZCells << "\n";
   }


   auto timeStepOfFirstTiming = uint_t(50);

   if( timesteps - timeStepOfFirstTiming < numSamples )
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Less actual time steps than number of required samples!");
   }

   uint_t nSample( 0 ); // number of current sample
   real_t samplingFrequency = real_c(timesteps - timeStepOfFirstTiming) / real_c(numSamples);

   // time loop
   for (uint_t i = 1; i <= timesteps; ++i )
   {
      // perform a single simulation step
      timeloop.singleStep( timeloopTiming );

      // check if current time step should be included in sample
      if( i >= uint_c( samplingFrequency * real_c(nSample) ) + timeStepOfFirstTiming )
      {
         // include -> evaluate all timers and quantities

         evaluateFluidQuantities(blocks, boundaryHandlingID, numCells, numFluidCells, numNBCells);
         evaluatePEQuantities(blocks, bodyStorageID, cr, numLocalParticles, numShadowParticles, numContacts);

         evaluateTimers(timeloopTiming, timingTreePE, timerKeys, timings);

         if(!noFileOutput)
         {
            real_t totalTime = std::accumulate(timings.begin(), timings.end(), real_t(0) );

            file << timeloop.getCurrentTimeStep() << " " << real_c(timeloop.getCurrentTimeStep()) / Tref << " "
                 << numCells << " " << numFluidCells << " " << numNBCells << " "
                 << numLocalParticles << " " << numShadowParticles << " " << numContacts << " " << numPeSubCycles;
            for (real_t timing : timings) {
               file << " " << timing;
            }
            file << " " << totalTime << "\n";
         }

         numCells = uint_t(0);
         numFluidCells = uint_t(0);
         numNBCells = uint_t(0);
         numLocalParticles = uint_t(0);
         numShadowParticles = uint_t(0);
         numContacts = uint_t(0);

         ++nSample;
      }

      // reset timers to always include only a single time step in them
      resetTimers(timeloopTiming, timingTreePE);
   }

   if(!noFileOutput) {
      file.close();
   }

   //timeloopTiming.logResultOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT("Simulation finished!");

   return 0;

}

} //namespace workload_evaluation

int main( int argc, char **argv ){
   workload_evaluation::main(argc, argv);
}
