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
//! \file ParticleMapping.cpp
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

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"

#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/SingleCast.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"

#include "vtk/all.h"
#include "field/vtk/all.h"

namespace particle_mapping
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

// boundary handling
using NoSlip_T = lbm::NoSlip< LatticeModel_T, flag_t >;
using BoundaryHandling_T = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T >;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag ( "fluid" );
const FlagUID NoSlip_Flag  ( "no slip" );


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const /*storage*/ ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   //WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField       = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField        = block->getData< PdfField_T > ( pdfFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "fixed obstacle boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ) );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}



class SphereMappingChecker
{
public:
   SphereMappingChecker(const shared_ptr< StructuredBlockStorage > & blocks,
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
         BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         FlagField_T *        flagField        = boundaryHandling->getFlagField();

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
         BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         FlagField_T *        flagField        = boundaryHandling->getFlagField();

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
                                              << " - expected boundary cell. Distance cell center - particle center = "
                                              << distance << ".");
            }
            else
            {
               WALBERLA_CHECK( boundaryHandling->isDomain(cellIt),
                               testIdentifier << ": Invalid mapping in ghost layer cell " << cellIt
                                              << " with center at " << cellCenter
                                              << " - expected domain cell. Distance cell center - particle center = "
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


class HalfSpaceMappingChecker {
public:
   HalfSpaceMappingChecker(const shared_ptr<StructuredBlockStorage> &blocks,
                           const BlockDataID &boundaryHandlingID) :
         blocks_(blocks), boundaryHandlingID_(boundaryHandlingID) {
      WALBERLA_ASSERT(blocks->isXPeriodic());
   }

   // check the mapping in the inner domain of the block
   void operator()(std::string & testIdentifier, const Vector3<real_t> & pos, const Vector3<real_t> & normal )
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         FlagField_T *        flagField        = boundaryHandling->getFlagField();

         auto xyzSize = flagField->xyzSize();

         for (auto cellIt : xyzSize) {
            Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter(*blockIt, cellIt);
            real_t distance = (cellCenter - pos) * normal;

            if (distance <= real_t(0)) {
               WALBERLA_CHECK(boundaryHandling->isBoundary(cellIt),
                              testIdentifier << "Invalid mapping in cell " << cellIt
                                             << " with center at " << cellCenter
                                             << " - expected boundary cell. Distance cell center - particle center = "
                                             << distance << ".");
            } else {
               WALBERLA_CHECK(boundaryHandling->isDomain(cellIt),
                              testIdentifier << "Invalid mapping in cell " << cellIt
                                             << " with center at " << cellCenter
                                             << " - expected domain cell. Distance cell center - particle center = "
                                             << distance << ".");
            }
         }
      }
   }

   // checks only the mapping in the ghost layers
   void checkGhostLayer(std::string & testIdentifier, const Vector3<real_t> & pos, const Vector3<real_t> & normal  )
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         FlagField_T *        flagField        = boundaryHandling->getFlagField();

         auto xyzSizeWGl = flagField->xyzSizeWithGhostLayer();

         for (auto cellIt : xyzSizeWGl) {
            if( flagField->isInInnerPart(cellIt) ) continue;

            Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter(*blockIt, cellIt);
            real_t distance = (cellCenter - pos) * normal;

            if (distance <= real_t(0))
            {
               WALBERLA_CHECK( boundaryHandling->isBoundary(cellIt),
                               testIdentifier << ": Invalid mapping in ghost layer cell " << cellIt
                                              << " with center at " << cellCenter
                                              << " - expected boundary cell. Distance cell center - particle center = "
                                              << distance << ".");
            }
            else
            {
               WALBERLA_CHECK( boundaryHandling->isDomain(cellIt),
                               testIdentifier << ": Invalid mapping in ghost layer cell " << cellIt
                                              << " with center at " << cellCenter
                                              << " - expected domain cell. Distance cell center - particle center = "
                                              << distance << "." );
            }
         }
      }
   }


private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID boundaryHandlingID_;
};

class MappingResetter
{
public:
   MappingResetter(const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & boundaryHandlingID) :
         blocks_( blocks ), boundaryHandlingID_( boundaryHandlingID )
   { }

   void operator()()
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         FlagField_T *        flagField        = boundaryHandling->getFlagField();

         auto xyzSizeWGl = flagField->xyzSizeWithGhostLayer();
         // reset to domain (fluid)
         boundaryHandling->forceDomain(xyzSizeWGl);
      }
   }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID boundaryHandlingID_;
};



/*!\brief Test case for the particle mapping functionality
 *
 * The following mapping functions are tested for spheres:
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
 *
 * Furthermore, the mapping of a plane (halfspace) is tested which resides at the bottom of the domain.
 * It is by definition global and thus uses the GlobalParticlesSelector.
 *
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

   // add boundary handling
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling( flagFieldID, pdfFieldID), "boundary handling" );

   // rpd setup
   mesa_pd::domain::BlockForestDomain domain(blocks->getBlockForestPointer());

   //init data structures
   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor = mesa_pd::data::ParticleAccessorWithShape;
   ParticleAccessor accessor(ps, ss);
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );
   mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;

   // coupling
   const real_t overlap = real_t( 1.5 ) * dx;
   lbm_mesapd_coupling::ParticleMappingKernel<BoundaryHandling_T> particleMappingKernel(blocks, boundaryHandlingID);

   // vtk
   auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field" );
   flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );

   // testing utilities
   SphereMappingChecker sphereMappingChecker(blocks, boundaryHandlingID, radius);
   MappingResetter mappingResetter(blocks, boundaryHandlingID);

   // sphere positions for test scenarios
   Vector3<real_t> positionInsideBlock(real_t(10), real_t(10), real_t(10));
   Vector3<real_t> positionAtBlockBoarder(real_t(19), real_t(10), real_t(10));
   Vector3<real_t> positionAtPeriodicBoarder(real_t(1), real_t(10), real_t(10));

   /////////////////////
   // NO SLIP MAPPING //
   /////////////////////

   //////////////////////////////////
   // TEST SPHERE INSIDE ONE BLOCK //
   //////////////////////////////////

   ////////////////////
   // REGULAR SPHERE //
   ////////////////////
   {
      std::string testIdentifier("Test: regular sphere inside block with no slip mapping ");
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
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      sphereMappingChecker(testIdentifier,positionInsideBlock,false);
      sphereMappingChecker.checkGhostLayer(testIdentifier,positionInsideBlock,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   ///////////////////
   // GLOBAL SPHERE //
   ///////////////////
   {
      std::string testIdentifier("Test: global sphere inside block with no slip mapping");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(positionInsideBlock);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);

      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      sphereMappingChecker(testIdentifier,positionInsideBlock,false);
      sphereMappingChecker.checkGhostLayer(testIdentifier,positionInsideBlock,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere inside block with no slip mapping ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionInsideBlock ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionInsideBlock);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::FixedParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      sphereMappingChecker(testIdentifier,positionInsideBlock,false);
      sphereMappingChecker.checkGhostLayer(testIdentifier,positionInsideBlock,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   //////////////////////////////////
   // TEST SPHERE AT BLOCK BOARDER //
   //////////////////////////////////

   ////////////////////
   // REGULAR SPHERE //
   ////////////////////
   {
      std::string testIdentifier("Test: regular sphere at block boarder with no slip mapping ");
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
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      sphereMappingChecker(testIdentifier,positionAtBlockBoarder,false);
      sphereMappingChecker.checkGhostLayer(testIdentifier,positionAtBlockBoarder,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   ///////////////////
   // GLOBAL SPHERE //
   ///////////////////

   {
      std::string testIdentifier("Test: global sphere at block boarder with no slip mapping ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(positionAtBlockBoarder);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);

      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      sphereMappingChecker(testIdentifier,positionAtBlockBoarder,false);
      sphereMappingChecker.checkGhostLayer(testIdentifier,positionAtBlockBoarder,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }


   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere at block boarder with no slip mapping ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionAtBlockBoarder ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionAtBlockBoarder);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::FixedParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      sphereMappingChecker(testIdentifier,positionAtBlockBoarder,false);
      sphereMappingChecker.checkGhostLayer(testIdentifier,positionAtBlockBoarder,false);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }


   /////////////////////////////////////
   // TEST SPHERE AT PERIODIC BOARDER //
   /////////////////////////////////////

   ////////////////////
   // REGULAR SPHERE //
   ////////////////////
   {
      std::string testIdentifier("Test: regular sphere at periodic boarder with no slip mapping ");
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
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      sphereMappingChecker(testIdentifier,positionAtPeriodicBoarder,true);
      sphereMappingChecker.checkGhostLayer(testIdentifier,positionAtPeriodicBoarder,true);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   ///////////////////
   // GLOBAL SPHERE //
   ///////////////////

   // left out since periodicity does not affect global particles

   //////////////////
   // FIXED SPHERE //
   //////////////////
   {
      std::string testIdentifier("Test: fixed sphere at periodic boarder with no slip mapping ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create sphere
      if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), positionAtPeriodicBoarder ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(positionAtPeriodicBoarder);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::FixedParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      sphereMappingChecker(testIdentifier,positionAtPeriodicBoarder,true);
      sphereMappingChecker.checkGhostLayer(testIdentifier,positionAtPeriodicBoarder,true);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   ///////////////
   // HALFSPACE //
   ///////////////

   auto halfSpacePosition = Vector3<real_t>(22.5,0.,3.5);
   auto halfSpaceNormal = Vector3<real_t>(0.,0.,1.);
   auto halfSpaceShape = ss->create<mesa_pd::data::HalfSpace>( halfSpaceNormal );

   HalfSpaceMappingChecker halfSpaceMappingChecker(blocks, boundaryHandlingID);
   {
      std::string testIdentifier("Test: half space at bottom of domain with no slip mapping ");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      // create half space
      {
         mesa_pd::data::Particle&& p = *ps->create(true);
         p.setPosition(halfSpacePosition);
         p.setInteractionRadius(std::numeric_limits<real_t>::infinity());
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(halfSpaceShape);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
      }

      syncNextNeighborFunc(*ps, domain, overlap);

      // map
      ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), accessor, particleMappingKernel, accessor, NoSlip_Flag );

      if( writeVTK ) flagFieldVTK->write();

      // check mapping
      halfSpaceMappingChecker(testIdentifier,halfSpacePosition,halfSpaceNormal);
      halfSpaceMappingChecker.checkGhostLayer(testIdentifier,halfSpacePosition,halfSpaceNormal);
      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
      syncNextNeighborFunc(*ps, domain, overlap);
   }

   ///////////////////////
   // ROTATED HALFSPACE //
   ///////////////////////

   // removed because rotation is not supported by half space, see half space docu


   return 0;

}

} //namespace particle_mapping

int main( int argc, char **argv ){
   particle_mapping::main(argc, argv);
}
