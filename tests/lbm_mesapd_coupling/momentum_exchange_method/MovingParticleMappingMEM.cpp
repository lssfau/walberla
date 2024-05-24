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
//! \file MovingParticleMapping.cpp
//! \ingroup lbm_mesapd_coupling
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
#include "core/math/all.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"

#include "field/AddToStorage.h"

#include "lbm/boundary/all.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"

#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/DataTypes.h"

#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/kernel/SingleCast.h"

#include "vtk/all.h"
#include "field/vtk/all.h"

namespace moving_particle_mapping
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRT>;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

const uint_t FieldGhostLayers = 1;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;


///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag ( "fluid" );
const FlagUID MO_Flag ( "moving obstacle" );
const FlagUID FormerMO_Flag ( "former moving obstacle" );


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

template <typename ParticleAccessor>
class MyBoundaryHandling
{
public:

   using MO_T = lbm_mesapd_coupling::SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor >;
   using Type = BoundaryHandling< FlagField_T, Stencil_T, MO_T >;

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID,
                       const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor>& ac) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), particleFieldID_( particleFieldID ), ac_( ac ) {}

   Type * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( storage );

      auto * flagField     = block->getData< FlagField_T >( flagFieldID_ );
      auto *  pdfField     = block->getData< PdfField_T > ( pdfFieldID_ );
      auto * particleField = block->getData< lbm_mesapd_coupling::ParticleField_T > ( particleFieldID_ );

      const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

      Type * handling = new Type( "moving obstacle boundary handling", flagField, fluid,
                                  MO_T("MO_BB",  MO_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ) );

      handling->fillWithDomain( FieldGhostLayers );

      return handling;
   }

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor> ac_;

};

template <typename BoundaryHandling_T>
class MappingChecker
{
public:
   MappingChecker(const shared_ptr< StructuredBlockStorage > & blocks,
                  const BlockDataID & boundaryHandlingID, real_t sphereRadius) :
         blocks_( blocks ), boundaryHandlingID_( boundaryHandlingID ),
         sphereRadius_( sphereRadius ), sphereVolume_( math::pi * real_t(4) / real_t(3) * sphereRadius * sphereRadius * sphereRadius )
   {
      WALBERLA_ASSERT(blocks->isXPeriodic());
   }

   // check the mapping in the inner domain of the block and check mapped volume against real sphere volume
   void operator()(std::string & testIdentifier, const Vector3<real_t> & pos, bool periodic )
   {
      uint_t cellCounter( uint_t(0) );
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         auto * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         auto *        flagField = boundaryHandling->getFlagField();

         auto xyzSize = flagField->xyzSize();

         for (auto cellIt : xyzSize) {
            Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter(*blockIt, cellIt);
            real_t distance = (cellCenter - pos).length();
            if( periodic )
            {
               Vector3<real_t> periodicOffset(blocks_->getDomain().xSize(),0,0);
               // check if other (periodic) copies are closer
               distance = std::min(distance, (cellCenter - (pos-periodicOffset)).length());
               distance = std::min(distance, (cellCenter - (pos+periodicOffset)).length());
            }

            if (distance < sphereRadius_) {
               WALBERLA_CHECK(boundaryHandling->isBoundary(cellIt),
                              testIdentifier << "Invalid mapping in cell " << cellIt
                                             << " with center at " << cellCenter
                                             << " - expected boundary cell. Distance cell center - particle center = "
                                             << distance << ".");
               ++cellCounter;
            } else {
               WALBERLA_CHECK(boundaryHandling->isDomain(cellIt),
                              testIdentifier << "Invalid mapping in cell " << cellIt
                                             << " with center at " << cellCenter
                                             << " - expected domain cell. Distance cell center - particle center = "
                                             << distance << ".");
            }
         }
      }

      mpi::allReduceInplace(cellCounter, mpi::SUM);

      // mapped volume has to be - approximately - the same as the real volume
      real_t mappedVolume = real_c(cellCounter); // dx=1
      WALBERLA_CHECK(std::fabs( mappedVolume - sphereVolume_ ) / sphereVolume_ <= real_t(0.1),
                     "Mapped volume " << mappedVolume << " does not fit to real sphere volume " << sphereVolume_ << ".");
   }

   // checks only the mapping in the ghost layers
   void checkGhostLayer(std::string & testIdentifier, const Vector3<real_t> & pos, bool periodic )
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         auto * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         auto *        flagField = boundaryHandling->getFlagField();

         auto xyzSizeWGl = flagField->xyzSizeWithGhostLayer();

         for (auto cellIt : xyzSizeWGl) {
            if( flagField->isInInnerPart(cellIt) ) continue;

            Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter( *blockIt, cellIt);
            real_t distance = (cellCenter - pos).length();
            if( periodic )
            {
               Vector3<real_t> periodicOffset(blocks_->getDomain().xSize(),0,0);
               // check if other (periodic) copies are closer
               distance = std::min(distance, (cellCenter - (pos-periodicOffset)).length());
               distance = std::min(distance, (cellCenter - (pos+periodicOffset)).length());
            }

            if( distance < sphereRadius_ )
            {
               WALBERLA_CHECK( boundaryHandling->isBoundary(cellIt),
                               testIdentifier << ": Invalid mapping in ghost layer cell " << cellIt
                                              << " with center at " << cellCenter
                                              << " - expected boundary cell. Distance cell center - body center = "
                                              << distance << ".");
            }
            else
            {
               WALBERLA_CHECK( boundaryHandling->isDomain(cellIt),
                               testIdentifier << ": Invalid mapping in ghost layer cell " << cellIt
                                              << " with center at " << cellCenter
                                              << " - expected domain cell. Distance cell center - body center = "
                                              << distance << "." );
            }
         }
      }
   }


private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID boundaryHandlingID_;
   real_t sphereRadius_, sphereVolume_;

};

template< typename BoundaryHandling_T>
class MappingResetter
{
public:
   MappingResetter(const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & boundaryHandlingID, const BlockDataID & particleFieldID, walberla::id_t invalidUID) :
         blocks_( blocks ), boundaryHandlingID_( boundaryHandlingID ), particleFieldID_( particleFieldID ), invalidUID_( invalidUID )
   { }

   void operator()()
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         auto * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         auto *        flagField = boundaryHandling->getFlagField();
         auto *    particleField = blockIt->getData< lbm_mesapd_coupling::ParticleField_T >( particleFieldID_ );

         auto xyzSizeWGl = flagField->xyzSizeWithGhostLayer();
         // reset to domain (fluid)
         boundaryHandling->forceDomain(xyzSizeWGl);

         for( auto cellIt = xyzSizeWGl.begin(); cellIt != xyzSizeWGl.end(); ++cellIt )
         {
            // reset body field
            particleField->get(*cellIt) = invalidUID_;
         }
      }
   }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID particleFieldID_;
   walberla::id_t invalidUID_;
};



/*!\brief Test case for the particle mapping functionality
 *
 * The following mapping functions are tested:
 *  - RegularParticlesSelector
 *  - GlobalParticlesSelector
 *  - FixedParticlesSelector
 * The mapping inside a block, at a regular block boarder and at a periodic block boarder is tested.
 * Regular, global and infinite-mass spheres are tested.
 * After each test, the sphere is destroyed and the mapping is erased from the datastructures.
 *
 * The setup is as follows (x marks the different sphere positions, # is a periodic border)
 * ----------------------------------------------------------------------
 * #                      |                      |                      #
 * #                      |                      |                      #
 * #                      |                      |                      #
 * #                      |                      |                      #
 * #X         X          X|                      |                      #
 * #                      |                      |                      #
 * #                      |                      |                      #
 * #                      |                      |                      #
 * #                      |                      |                      #
 * ----------------------------------------------------------------------
 */


//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   WALBERLA_CHECK(mpi::MPIManager::instance()->numProcesses()==3, "Due to periodicity, three processes have to be used!");

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   bool writeVTK = false;
   const real_t omega  = real_t(1);
   const real_t dx     = real_t(1);
   const real_t radius = real_t(5);

   ///////////////////////////
   // DATA STRUCTURES SETUP //
   ///////////////////////////

   Vector3<uint_t> blocksPerDirection(uint_t(3), uint_t(1), uint_t(1));
   Vector3<uint_t> cellsPerBlock(uint_t(20), uint_t(20), uint_t(20));
   Vector3<bool> periodicity(true, false, false);

   auto blocks = blockforest::createUniformBlockGrid( blocksPerDirection[0], blocksPerDirection[1], blocksPerDirection[2],
                                                      cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
                                                      dx,
                                                      1, false, false,
                                                      periodicity[0], periodicity[1], periodicity[2],
                                                      false );

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( omega );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         Vector3<real_t>(real_t(0)), real_t(1),
                                                                         FieldGhostLayers, field::fzyx );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // rpd setup
   mesa_pd::domain::BlockForestDomain domain(blocks->getBlockForestPointer());

   //init data structures
   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   shared_ptr<ParticleAccessor_T> accessor = make_shared<ParticleAccessor_T>(ps, ss);
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );
   mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;

   // coupling
   const real_t overlap = real_t( 1.5 ) * dx;

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::fzyx, FieldGhostLayers );

   // add boundary handling
   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor), "boundary handling" );

   // vtk
   auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field" );
   flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );

   // testing utilities
   MappingChecker<BoundaryHandling_T> mappingChecker(blocks, boundaryHandlingID, radius);
   MappingResetter<BoundaryHandling_T> mappingResetter(blocks, boundaryHandlingID, particleFieldID, accessor->getInvalidUid());

   // mapping functors
   auto regularParticleMapper = lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag, lbm_mesapd_coupling::RegularParticlesSelector(), true );
   auto globalParticleMapper = lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag, lbm_mesapd_coupling::GlobalParticlesSelector(), true );

   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);

   // sphere positions for test scenarios
   Vector3<real_t> positionInsideBlock(real_t(10), real_t(10), real_t(10));
   Vector3<real_t> positionAtBlockBoarder(real_t(19), real_t(10), real_t(10));
   Vector3<real_t> positionAtPeriodicBoarder(real_t(1), real_t(10), real_t(10));

   //////////////////////////
   // MOVING BODY MAPPING //
   /////////////////////////

   //////////////////////////////////
   // TEST SPHERE INSIDE ONE BLOCK //
   //////////////////////////////////

   //////////////////////
   // REGULAR SPHERE 1 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere inside block with moving particle mapping, case 1 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionInsideBlock ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionInsideBlock);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionInsideBlock,false);
      mappingChecker.checkGhostLayer(testIdentifier,positionInsideBlock,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   //////////////////////
   // REGULAR SPHERE 2 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere inside block with moving particle mapping, case 2 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionInsideBlock ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionInsideBlock);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularParticleMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionInsideBlock,false);
      mappingChecker.checkGhostLayer(testIdentifier,positionInsideBlock,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere inside block with moving particle mapping, case 1 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(positionInsideBlock);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);

      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionInsideBlock,false);
      mappingChecker.checkGhostLayer(testIdentifier,positionInsideBlock,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   /////////////////////
   // GLOBAL SPHERE 2 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere inside block with moving particle mapping, case 2 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(positionInsideBlock);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);

      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) globalParticleMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionInsideBlock,false);
      mappingChecker.checkGhostLayer(testIdentifier,positionInsideBlock,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   //////////////////////////////////
   // TEST SPHERE AT BLOCK BOARDER //
   //////////////////////////////////

   //////////////////////
   // REGULAR SPHERE 1 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere at block boarder with moving particle mapping, case 1 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionAtBlockBoarder ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionAtBlockBoarder);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionAtBlockBoarder,false);
      mappingChecker.checkGhostLayer(testIdentifier,positionAtBlockBoarder,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   //////////////////////
   // REGULAR SPHERE 2 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere at block boarder with moving particle mapping, case 2 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionAtBlockBoarder ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionAtBlockBoarder);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularParticleMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionAtBlockBoarder,false);
      mappingChecker.checkGhostLayer(testIdentifier,positionAtBlockBoarder,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere at block boarder with moving particle mapping, case 1 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(positionAtBlockBoarder);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);

      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionAtBlockBoarder,false);
      mappingChecker.checkGhostLayer(testIdentifier,positionAtBlockBoarder,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   /////////////////////
   // GLOBAL SPHERE 2 //
   /////////////////////
   {
      std::string testIdentifier("Test: global sphere at block boarder with moving particle mapping, case 2 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(positionAtBlockBoarder);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);

      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) globalParticleMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionAtBlockBoarder,false);
      mappingChecker.checkGhostLayer(testIdentifier,positionAtBlockBoarder,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   /////////////////////////////////////
   // TEST SPHERE AT PERIODIC BOARDER //
   /////////////////////////////////////

   //////////////////////
   // REGULAR SPHERE 1 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere at periodic boarder with moving particle mapping, case 1 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionAtPeriodicBoarder ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionAtPeriodicBoarder);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionAtPeriodicBoarder,true);
      mappingChecker.checkGhostLayer(testIdentifier,positionAtPeriodicBoarder,true);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   //////////////////////
   // REGULAR SPHERE 2 //
   //////////////////////
   {
      std::string testIdentifier("Test: regular sphere at periodic boarder with moving particle mapping, case 2 ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionAtPeriodicBoarder ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionAtPeriodicBoarder);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularParticleMapper(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      mappingChecker(testIdentifier,positionAtPeriodicBoarder,true);
      mappingChecker.checkGhostLayer(testIdentifier,positionAtPeriodicBoarder,true);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   /////////////////////
   // GLOBAL SPHERE 1 //
   /////////////////////

   // left out since periodicity does not affect global particles

   return 0;

}

} //namespace moving_particle_mapping

int main( int argc, char **argv ){
   moving_particle_mapping::main(argc, argv);
}
