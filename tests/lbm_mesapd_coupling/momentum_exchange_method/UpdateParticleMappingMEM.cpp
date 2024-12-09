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
//! \file UpdateParticleMapping.cpp
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

#include "field/AddToStorage.h"

#include "lbm/boundary/all.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
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

namespace update_particle_mapping
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
const FlagUID NoSlip_Flag( "no slip" );


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

template <typename ParticleAccessor>
class MyBoundaryHandling
{
public:

   using NoSlip_T = lbm::NoSlip< LatticeModel_T, flag_t >;
   using MO_T = lbm_mesapd_coupling::SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor >;
   using Type = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, MO_T >;

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
                                  NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
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
                  const BlockDataID & boundaryHandlingID) :
         blocks_( blocks ), boundaryHandlingID_( boundaryHandlingID )
   { }

   // check the mapping and check mapped volume against real sphere volume located at position pos
   void operator()(std::string testIdentifier, const Vector3<real_t> & pos, real_t sphereRadius )
   {
      real_t sphereVolume = math::pi * real_t(4) / real_t(3) * sphereRadius * sphereRadius * sphereRadius;
      uint_t cellCounter( uint_t(0) );
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         auto * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         auto *        flagField = boundaryHandling->getFlagField();

         auto xyzSize = flagField->xyzSize();

         for (auto cellIt : xyzSize) {
            Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter(*blockIt, cellIt);
            real_t distance = (cellCenter - pos).length();

            if (distance < sphereRadius) {
               WALBERLA_CHECK(boundaryHandling->isBoundary(cellIt),
                              testIdentifier << " Invalid mapping in cell " << cellIt
                                             << " with center at " << cellCenter
                                             << " - expected boundary cell. Distance cell center - particle center = "
                                             << distance << ".");
               ++cellCounter;
            }
         }
      }
      // mapped volume has to be - approximately - the same as the real volume
      real_t mappedVolume = real_c(cellCounter); // dx=1
      WALBERLA_CHECK(std::fabs( mappedVolume - sphereVolume ) / sphereVolume <= real_t(0.1),
                     testIdentifier << " Mapped volume " << mappedVolume << " does not fit to real sphere volume " << sphereVolume << ".");
   }

   // check the mapping of a plane for a given AABB
   void operator()(std::string testIdentifier, const math::AABB & planeAABB )
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         auto * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID_ );
         auto *        flagField = boundaryHandling->getFlagField();

         auto xyzSize = flagField->xyzSize();

         for (auto cellIt : xyzSize) {
            Vector3<real_t> cellCenter = blocks_->getBlockLocalCellCenter(*blockIt, cellIt);

            if (planeAABB.contains(cellCenter)) {
               WALBERLA_CHECK(boundaryHandling->isBoundary(cellIt),
                              testIdentifier << " Invalid mapping in cell " << cellIt
                                             << " with center at " << cellCenter
                                             << " - expected boundary cell since inside AABB " << planeAABB);
            }
         }
      }
   }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID boundaryHandlingID_;

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



/*!\brief Test case for the update particle mapping functionality
 *
 * The following scenarios are tested:
 *  - two spheres in contact (i.e. overlapping the same cell), then one moves away
 *  - sphere-plane in contact, then sphere moves away
 *  - sphere in contact with a no-slip boundary condition, then moves away
 *
 *  The expected behavior in all cases is that the mapping afterwards is still valid and existing boundary conditions are not overwritten.
 */


//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   bool writeVTK = true;
   const real_t omega  = real_t(1);
   const real_t dx     = real_t(1);
   const real_t radius = real_t(5);

   ///////////////////////////
   // DATA STRUCTURES SETUP //
   ///////////////////////////

   Vector3<uint_t> blocksPerDirection(uint_t(1), uint_t(1), uint_t(1));
   Vector3<uint_t> cellsPerBlock(uint_t(30), uint_t(30), uint_t(30));
   Vector3<bool> periodicity(false, false, false);

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

   //init data structures
   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   shared_ptr<ParticleAccessor_T> accessor = make_shared<ParticleAccessor_T>(ps, ss);
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );
   auto planeShape = ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(real_t(1), real_t(0), real_t(0)) );

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::fzyx, FieldGhostLayers );

   // add boundary handling
   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor), "boundary handling" );
   auto bhBlockSweep = BoundaryHandling_T::getBlockSweep( boundaryHandlingID );

   // reconstructor
   auto reconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, false);

   // vtk
   auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field" );
   flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );

   // testing utilities
   MappingChecker<BoundaryHandling_T> mappingChecker(blocks, boundaryHandlingID);
   MappingResetter<BoundaryHandling_T> mappingResetter(blocks, boundaryHandlingID, particleFieldID, accessor->getInvalidUid());

   // mapping functors
   auto regularParticleMapper = lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag, lbm_mesapd_coupling::RegularParticlesSelector(), true );
   auto globalParticleMapper = lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag, lbm_mesapd_coupling::GlobalParticlesSelector(), true );

   /////////////////////////////
   // SPHERE - SPHERE MAPPING //
   /////////////////////////////
   {
      std::string testIdentifier("Test: two spheres single moving");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      Vector3<real_t> pos1(10, 10, 9.5);
      Vector3<real_t> pos2(19, 10, 9.5);

      walberla::id_t sphere2Uid;

      // create sphere 1
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos1);
         p.setInteractionRadius(radius);
         p.setShapeID(sphereShape);
      }

      // create sphere 2
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos2);
         p.setInteractionRadius(radius);
         p.setShapeID(sphereShape);
         sphere2Uid = p.getUid();
      }

      // map
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (regularParticleMapper)(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check initial mapping
      mappingChecker(testIdentifier + " mapping check 1, sphere 1", pos1, radius);
      mappingChecker(testIdentifier + " mapping check 1, sphere 2", pos2, radius);

      // carry out boundary handling (sets forces)
      for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) bhBlockSweep(&(*blockIt));

      // check force on particles (should be zero)
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, [](const size_t idx, const ParticleAccessor_T& ac){ WALBERLA_CHECK_EQUAL(ac.getHydrodynamicForce(idx), Vector3<real_t>(real_t(0)));}, *accessor);

      // update position
      auto updatedPos2 = pos2 + Vector3<real_t>(real_t(1), real_t(0), real_t(0));
      accessor->setPosition(accessor->uidToIdx(sphere2Uid),updatedPos2 );

      // update mapping
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (regularParticleMapper)(&(*blockIt));

      // PDF reconstruction
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (*reconstructionManager)(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check updated mapping
      mappingChecker(testIdentifier + " mapping check 2, sphere 1", pos1, radius);
      mappingChecker(testIdentifier + " mapping check 2, sphere 2", updatedPos2, radius);

      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
   }


   /////////////////////////////
   // SPHERE - SPHERE MAPPING //
   /////////////////////////////
   {
      std::string testIdentifier("Test: two spheres both moving");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      Vector3<real_t> pos1(10, 10, 9.5);
      Vector3<real_t> pos2(19, 10, 9.5);

      walberla::id_t sphere1Uid;
      walberla::id_t sphere2Uid;

      // create sphere 1
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos1);
         p.setInteractionRadius(radius);
         p.setShapeID(sphereShape);
         sphere1Uid = p.getUid();
      }

      // create sphere 2
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos2);
         p.setInteractionRadius(radius);
         p.setShapeID(sphereShape);
         sphere2Uid = p.getUid();
      }

      // map
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (regularParticleMapper)(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check initial mapping
      mappingChecker(testIdentifier + " mapping check 1, sphere 1", pos1, radius);
      mappingChecker(testIdentifier + " mapping check 1, sphere 2", pos2, radius);

      // carry out boundary handling (sets forces)
      for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) bhBlockSweep(&(*blockIt));

      // check force on particles (should be zero)
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, [](const size_t idx, const ParticleAccessor_T& ac){ WALBERLA_CHECK_EQUAL(ac.getHydrodynamicForce(idx), Vector3<real_t>(real_t(0)));}, *accessor);

      // update position
      auto updatedPos1 = pos1 + Vector3<real_t>(-real_t(1), real_t(0), real_t(0));
      accessor->setPosition(accessor->uidToIdx(sphere1Uid),updatedPos1 );
      auto updatedPos2 = pos2 + Vector3<real_t>(real_t(1), real_t(0), real_t(0));
      accessor->setPosition(accessor->uidToIdx(sphere2Uid),updatedPos2 );

      // update mapping
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (regularParticleMapper)(&(*blockIt));

      // PDF reconstruction
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (*reconstructionManager)(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check updated mapping
      mappingChecker(testIdentifier + " mapping check 2, sphere 1", updatedPos1, radius);
      mappingChecker(testIdentifier + " mapping check 2, sphere 2", updatedPos2, radius);

      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
   }


   /////////////////////////////
   // SPHERE - PLANE MAPPING //
   /////////////////////////////
   {
      std::string testIdentifier("Test: sphere and plane");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      Vector3<real_t> pos1(1, 10, 9.5);
      Vector3<real_t> pos2(5, 10, 9.5);

      walberla::id_t sphere2Uid;

      // create plane
      {
         mesa_pd::data::Particle&& p = *ps->create(true);
         p.setPosition(pos1);
         p.setShapeID(planeShape);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
      }

      // create sphere 2
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos2);
         p.setInteractionRadius(radius);
         p.setShapeID(sphereShape);
         sphere2Uid = p.getUid();
      }

      // map
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (regularParticleMapper)(&(*blockIt));
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (globalParticleMapper)(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check initial mapping
      math::AABB planeAABB(real_t(0), real_t(0), real_t(0), pos1[0], real_c(cellsPerBlock[1]), real_c(cellsPerBlock[2]));
      mappingChecker(testIdentifier + " mapping check 1, plane", planeAABB);
      mappingChecker(testIdentifier + " mapping check 1, sphere", pos2, radius);

      // carry out boundary handling (sets forces)
      for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) bhBlockSweep(&(*blockIt));

      // check force on particles (should be zero)
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, [](const size_t idx, const ParticleAccessor_T& ac){ WALBERLA_CHECK_EQUAL(ac.getHydrodynamicForce(idx), Vector3<real_t>(real_t(0)));}, *accessor);

      // update position
      auto updatedPos2 = pos2 + Vector3<real_t>(real_t(1), real_t(0), real_t(0));
      accessor->setPosition(accessor->uidToIdx(sphere2Uid),updatedPos2 );

      // update mapping
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (regularParticleMapper)(&(*blockIt));

      // PDF reconstruction
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (*reconstructionManager)(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check updated mapping
      mappingChecker(testIdentifier + " mapping check 2, plane", planeAABB);
      mappingChecker(testIdentifier + " mapping check 2, sphere", updatedPos2, radius);

      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
   }

   ///////////////////////////////
   // SPHERE - BOUNDARY MAPPING //
   ///////////////////////////////
   // here, a plane is used to map the no slip boundary condition into the domain
   {
      std::string testIdentifier("Test: sphere and boundary");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      Vector3<real_t> pos1(1, 10, 9.5);
      Vector3<real_t> pos2(5, 10, 9.5);

      walberla::id_t sphere2Uid;

      // create plane
      {
         mesa_pd::data::Particle&& p = *ps->create(true);
         p.setPosition(pos1);
         p.setShapeID(planeShape);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
         mesa_pd::data::particle_flags::set(p.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
      }

      // create sphere 2
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos2);
         p.setInteractionRadius(radius);
         p.setShapeID(sphereShape);
         sphere2Uid = p.getUid();
      }

      // map
      lbm_mesapd_coupling::ParticleMappingKernel<BoundaryHandling_T> particleMappingKernel(blocks, boundaryHandlingID);
      ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), *accessor, particleMappingKernel, *accessor, NoSlip_Flag);
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (regularParticleMapper)(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check initial mapping
      math::AABB planeAABB(real_t(0), real_t(0), real_t(0), pos1[0], real_c(cellsPerBlock[1]), real_c(cellsPerBlock[2]));
      mappingChecker(testIdentifier + " mapping check 1, boundary", planeAABB);
      mappingChecker(testIdentifier + " mapping check 1, sphere", pos2, radius);

      // carry out boundary handling (sets forces)
      for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) bhBlockSweep(&(*blockIt));

      // check force on particles (should be zero)
      ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, [](const size_t idx, const ParticleAccessor_T& ac){ WALBERLA_CHECK_EQUAL(ac.getHydrodynamicForce(idx), Vector3<real_t>(real_t(0)));}, *accessor);

      // update position
      auto updatedPos2 = pos2 + Vector3<real_t>(real_t(1), real_t(0), real_t(0));
      accessor->setPosition(accessor->uidToIdx(sphere2Uid),updatedPos2 );

      // update mapping
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (regularParticleMapper)(&(*blockIt));

      // PDF reconstruction
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) (*reconstructionManager)(&(*blockIt));

      if( writeVTK ) flagFieldVTK->write();

      // check updated mapping
      mappingChecker(testIdentifier + " mapping check 2, boundary", planeAABB);
      mappingChecker(testIdentifier + " mapping check 2, sphere", updatedPos2, radius);

      mappingResetter();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

      // clean up
      ps->clear();
   }

   return 0;

}

} //namespace update_particle_mapping

int main( int argc, char **argv ){
   update_particle_mapping::main(argc, argv);
}
