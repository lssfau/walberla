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
//! \file PDFReconstruction.cpp
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
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/ExtrapolationDirectionFinder.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/DataTypes.h"

#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/SingleCast.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"

#include "vtk/all.h"
#include "field/vtk/all.h"

#include "vector"
#include "tuple"

namespace pdf_reconstruction
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
   using Type = BoundaryHandling< FlagField_T, Stencil_T, MO_T>;

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

template< typename PdfField_T, typename BoundaryHandling_T>
class ReconstructionChecker
{
public:
   ReconstructionChecker(const shared_ptr< StructuredBlockStorage > & blocks,const BlockDataID & pdfFieldID,
                         const BlockDataID & boundaryHandlingID, Vector3<real_t> expectedVelocity, real_t expectedDensity,
                         const FlagUID formerObstacle) :
   blocks_( blocks ), pdfFieldID_( pdfFieldID ), boundaryHandlingID_ (boundaryHandlingID ), expectedVelocity_( expectedVelocity ), expectedDensity_( expectedDensity ), formerObstacle_(formerObstacle)
   { }

   void operator()(std::string testIdentifier)
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt ) {
         auto *pdfField         = blockIt->getData<PdfField_T>(pdfFieldID_);
         auto *boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID_);
         auto *flagField        = boundaryHandling->getFlagField();

         auto xyzSize = pdfField->xyzSize();

         const flag_t formerObstacle = flagField->getFlag( formerObstacle_ );

         for (auto cellIt : xyzSize) {
            WALBERLA_CHECK(!isFlagSet(flagField->get(cellIt), formerObstacle),testIdentifier << "Found incorrect formerObstacle flag in cell " << cellIt << ".")

            if(boundaryHandling->isDomain(cellIt))
            {
               auto velocity = pdfField->getVelocity(cellIt);
               WALBERLA_CHECK_FLOAT_EQUAL(velocity[0], expectedVelocity_[0],testIdentifier << "Invalid velocity X in cell " << cellIt << ".");
               WALBERLA_CHECK_FLOAT_EQUAL(velocity[1], expectedVelocity_[1],testIdentifier << "Invalid velocity Y in cell " << cellIt << ".");
               WALBERLA_CHECK_FLOAT_EQUAL(velocity[2], expectedVelocity_[2],testIdentifier << "Invalid velocity Z in cell " << cellIt << ".");

               auto density = pdfField->getDensity(cellIt);
               WALBERLA_CHECK_FLOAT_EQUAL(density, expectedDensity_,testIdentifier << "Invalid density in cell " << cellIt << ".");
            }
         }
      }
   }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;

   Vector3<real_t> expectedVelocity_;
   real_t expectedDensity_;

   const FlagUID formerObstacle_;

};



template <typename BoundaryHandling_T>
class MappingChecker
{
public:
   MappingChecker(const shared_ptr< StructuredBlockStorage > & blocks,
                  const BlockDataID & boundaryHandlingID, real_t sphereRadius) :
         blocks_( blocks ), boundaryHandlingID_( boundaryHandlingID ),
         sphereRadius_( sphereRadius ), sphereVolume_( math::pi * real_t(4) / real_t(3) * sphereRadius * sphereRadius * sphereRadius )
   { }

   // check the mapping in the inner domain of the block and check mapped volume against real sphere volume
   void operator()(std::string testIdentifier, const Vector3<real_t> & pos, bool periodic )
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
   void checkGhostLayer(std::string testIdentifier, const Vector3<real_t> & pos, bool periodic )
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
            // reset particle field
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


template<typename BoundaryHandling_T, typename ParticleAccessor_T, typename ExtrapolationDirectionFinder_T>
void checkExtrapolationDirectionFinder( std::string identifier,
                                        shared_ptr< StructuredBlockForest > & blocks, shared_ptr<mesa_pd::data::ParticleStorage> & ps, shared_ptr<mesa_pd::data::ShapeStorage> ss,
                                        shared_ptr<ParticleAccessor_T> & accessor, BlockDataID pdfFieldID, BlockDataID boundaryHandlingID, BlockDataID particleFieldID,
                                        ExtrapolationDirectionFinder_T directionFinder, real_t radius, bool conserveMomentum)
{
   mesa_pd::domain::BlockForestDomain domain(blocks->getBlockForestPointer());
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );
   mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
   const real_t overlap = real_t( 1.5 );

   MappingResetter<BoundaryHandling_T> mappingResetter(blocks, boundaryHandlingID, particleFieldID, accessor->getInvalidUid());
   auto regularParticleMapper = lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag, mesa_pd::kernel::SelectAll(), conserveMomentum);

   std::string testIdentifier(identifier + " Test:");
   WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

   Vector3<real_t> spherePosition(real_t(9.5), real_t(9.5), real_t(9.5));

   std::array<Vector3<real_t>,3> referenceNormals = {{Vector3<real_t>(1,0,0), Vector3<real_t>(0,-1,0), Vector3<real_t>(-1,0,1)/sqrt(2)}};
   std::array<Vector3<real_t>,3> evaluationPoints = {{spherePosition + referenceNormals[0] * radius, spherePosition + referenceNormals[1] * radius, spherePosition + referenceNormals[2] * radius}};

   // create sphere
   walberla::id_t sphereUid = 0;
   if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), spherePosition ))
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(spherePosition);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);
      sphereUid = p.getUid();
   }

   mpi::allReduceInplace(sphereUid, mpi::SUM);

   syncNextNeighborFunc(*ps, domain, overlap);

   // map
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularParticleMapper(&(*blockIt));

   // check normal at evaluation points
   // since the normal must be aligned with some lattice direction, the normal could be different from the real one
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      for( uint_t i = 0; i < evaluationPoints.size(); ++i)
      {
         auto evalPoint = evaluationPoints[i];
         if( blocks->getAABB((*blockIt).getId()).contains(evalPoint))
         {
            auto cell = blocks->getBlockLocalCell(*blockIt,evalPoint);

            size_t idx = accessor->uidToIdx(sphereUid);
            if( idx != accessor->getInvalidIdx())
            {
               auto normal = directionFinder(&(*blockIt), cell.x(), cell.y(), cell.z(), idx, *accessor).getNormalized();
               for(uint_t comp = 0; comp < 3; ++comp) WALBERLA_CHECK_FLOAT_EQUAL(normal[comp], referenceNormals[i][comp], "Non-matching extrapolation direction for normal component " << comp << " in evaluation point " << evalPoint);
            }
         }
      }
   }

   WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

   // clean up
   mappingResetter();
   ps->clear();
   syncNextNeighborFunc(*ps, domain, overlap);


}


template<typename BoundaryHandling_T, typename ParticleAccessor_T, typename ReconstructionManager_T>
void checkReconstruction( std::string testIdentifier, Vector3<real_t> spherePosition, bool periodicCase,
                                   shared_ptr< StructuredBlockForest > & blocks, shared_ptr<mesa_pd::data::ParticleStorage> & ps, shared_ptr<mesa_pd::data::ShapeStorage> ss,
                                   shared_ptr<ParticleAccessor_T> & accessor, BlockDataID pdfFieldID, BlockDataID boundaryHandlingID, BlockDataID particleFieldID,
                                   ReconstructionManager_T reconstructionManager,
                                   real_t radius, Vector3<real_t> velocity, real_t density, bool conserveMomentum )
{


   mesa_pd::domain::BlockForestDomain domain(blocks->getBlockForestPointer());
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );
   mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
   const real_t overlap = real_t( 1.5 );

   MappingChecker<BoundaryHandling_T> mappingChecker(blocks, boundaryHandlingID, radius);
   MappingResetter<BoundaryHandling_T> mappingResetter(blocks, boundaryHandlingID, particleFieldID, accessor->getInvalidUid());
   ReconstructionChecker<PdfField_T, BoundaryHandling_T> reconstructionChecker(blocks, pdfFieldID, boundaryHandlingID, velocity, density, FormerMO_Flag );

   auto regularParticleMapper = lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag, mesa_pd::kernel::SelectAll(), conserveMomentum);

   WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

   // create sphere
   walberla::id_t sphereUid = 0;
   if (domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), spherePosition ))
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(spherePosition);
      p.setLinearVelocity(velocity);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);
      sphereUid = p.getUid();
   }
   syncNextNeighborFunc(*ps, domain, overlap);

   mpi::allReduceInplace(sphereUid, mpi::SUM);

   // map
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularParticleMapper(&(*blockIt));

   auto currentSpherePosition = spherePosition;
   // run several "timesteps" and advance the sphere in total one cell
   for( uint_t t = 0; t < 10; ++t )
   {
      currentSpherePosition = currentSpherePosition + velocity;
      blocks->mapToPeriodicDomain(currentSpherePosition);

      // advance sphere
      size_t idx = accessor->uidToIdx(sphereUid);
      if( idx != accessor->getInvalidIdx())
      {
         accessor->setPosition(idx, currentSpherePosition);
      }
      syncNextNeighborFunc(*ps, domain, overlap);

      // update mapping
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) regularParticleMapper(&(*blockIt));

      // carry out restore step
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) reconstructionManager(&(*blockIt));


      // check mapping
      mappingChecker(testIdentifier + " (t=" + std::to_string(t)+") ", currentSpherePosition, periodicCase);
      mappingChecker.checkGhostLayer(testIdentifier + " (t=" + std::to_string(t)+") ", currentSpherePosition, periodicCase);

      // check reconstruction
      reconstructionChecker(testIdentifier + " (t=" + std::to_string(t)+") ");
   }

   WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended successfully");

   // clean up
   mappingResetter();
   ps->clear();
   syncNextNeighborFunc(*ps, domain, overlap);
}


/*!\brief Test case for the Pdf reconstruction functionality
 *
 * It tests these different reconstructors:
 *  - EquilibriumReconstructor
 *  - EquilibriumAndNonEquilibriumReconstructor
 *  - ExtrapolationReconstructor
 *  - GradsMomentApproximationReconstructor
 *
 * The following features of the PdfReconstructionManager are tested as well:
 *  - reconstruction at block boarder
 *  - version for small obstacle fraction
 *  - periodicity
 *
 * Additionally, the following extrapolation direction finders are tested:
 *  - SphereNormalExtrapolationFinder
 *  - FlagFieldExtrapolationDirectionFinder
 *
 * ----------------------------------------------------------------------
 * #                      |                      |                      #
 * #                      |                      |                      #
 * #                      |                      |                      #
 * #                      |                      |                      #
 * # 5        1          2|  3                   |                     4#
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

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   const real_t omega  = real_t(1);
   const real_t dx     = real_t(1);
   const real_t radius = real_t(5);
   const Vector3<real_t> velocity( real_t(0.1), real_t(0), real_t(0) );
   const real_t density = real_t(1);

   bool conserveMomentum = true;

   ///////////////////////////
   // DATA STRUCTURES SETUP //
   ///////////////////////////

   Vector3<uint_t> blocksPerDirection(uint_t(3), uint_t(1), uint_t(1));
   Vector3<uint_t> cellsPerBlock(uint_t(20), uint_t(20), uint_t(20));
   Vector3<bool> periodicity(true, false, false);

   auto blocks = blockforest::createUniformBlockGrid( blocksPerDirection[0], blocksPerDirection[1], blocksPerDirection[2],
                                                      cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
                                                      dx,
                                                      0, false, false,
                                                      periodicity[0], periodicity[1], periodicity[2],
                                                      false );

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( omega );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         velocity, density,
                                                                         FieldGhostLayers, field::fzyx );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // rpd setup

   //init data structures
   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = make_shared<ParticleAccessor_T>(ps, ss);

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::fzyx, FieldGhostLayers );

   // add boundary handling
   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor), "boundary handling" );

   // vtk
   //auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field" );
   //flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );

   // test setups -> tuple of (setupName, spherePosition, periodicityTested)
   std::vector<std::tuple<std::string, Vector3<real_t>, bool> > testSetups;
   testSetups.push_back( std::make_tuple( "sphere inside block",          Vector3<real_t>(real_t(10), real_t(10), real_t(10)),        false) );
   testSetups.push_back( std::make_tuple( "sphere on block boarder",      Vector3<real_t>(real_t(19.5), real_t(10), real_t(10)),      false) );
   testSetups.push_back( std::make_tuple( "sphere on block boarder 2",    Vector3<real_t>(real_t(20)+radius, real_t(10), real_t(10)), false) );
   testSetups.push_back( std::make_tuple( "sphere on periodic boarder",   Vector3<real_t>(real_t(59.5), real_t(10), real_t(10)),      true) );
   testSetups.push_back( std::make_tuple( "sphere on periodic boarder 2", Vector3<real_t>(radius, real_t(10), real_t(10)),            true) );


   /////////////////////////////////////
   // EQUILIBRIUM RECONSTRUCTOR TESTS //
   /////////////////////////////////////

   auto equilibriumReconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, conserveMomentum);
   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Equilibrium Reconstructor Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *equilibriumReconstructionManager, radius, velocity, density, conserveMomentum );
   }

   /////////////////////////////////////
   // EQUILIBRIUM RECONSTRUCTOR TESTS //
   //  for small obstacle fractions   //
   /////////////////////////////////////

   auto equilibriumReconstructionSmallObstacleFractionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, conserveMomentum, true);
   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Equilibrium Reconstructor Small Obstacle Fraction Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *equilibriumReconstructionSmallObstacleFractionManager, radius, velocity, density, conserveMomentum );
   }

   /////////////////////////////////////////////////////////
   // EQUILIBRIUM AND NON-EQUILIBRIUM RECONSTRUCTOR TESTS //
   //       with Sphere Normal Extrapolator               //
   /////////////////////////////////////////////////////////

   auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(blocks);

   checkExtrapolationDirectionFinder<BoundaryHandling_T>( "Sphere Normal Extrapolation Direction Finder", blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                                          *sphereNormalExtrapolationDirectionFinder, radius, conserveMomentum );


   auto equilibriumAndNonEquilibriumSphereNormalReconstructor = lbm_mesapd_coupling::makeEquilibriumAndNonEquilibriumReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, 3);
   auto equilibriumAndNonEquilibriumSphereNormalReconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, equilibriumAndNonEquilibriumSphereNormalReconstructor, conserveMomentum);

   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Equilibrium and NonEquilibrium Reconstructor with Sphere Normal Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *equilibriumAndNonEquilibriumSphereNormalReconstructionManager, radius, velocity, density, conserveMomentum );
   }


   /////////////////////////////////////////////////////////
   // EQUILIBRIUM AND NON-EQUILIBRIUM RECONSTRUCTOR TESTS //
   //       with FlagField Normal Extrapolator            //
   /////////////////////////////////////////////////////////

   auto flagFieldNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::FlagFieldNormalExtrapolationDirectionFinder<BoundaryHandling_T> >(blocks,boundaryHandlingID);

   checkExtrapolationDirectionFinder<BoundaryHandling_T>( "FlagField Normal Extrapolation Direction Finder", blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                                          *flagFieldNormalExtrapolationDirectionFinder, radius, conserveMomentum );

   auto equilibriumAndNonEquilibriumFlagFieldNormalReconstructor = lbm_mesapd_coupling::makeEquilibriumAndNonEquilibriumReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, flagFieldNormalExtrapolationDirectionFinder, 2);
   auto equilibriumAndNonEquilibriumFlagFieldNormalReconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, equilibriumAndNonEquilibriumFlagFieldNormalReconstructor, conserveMomentum);

   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Equilibrium and NonEquilibrium Reconstructor with FlagField Normal Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *equilibriumAndNonEquilibriumFlagFieldNormalReconstructionManager, radius, velocity, density, conserveMomentum );
   }

   /////////////////////////////////////////////////////////
   // EQUILIBRIUM AND NON-EQUILIBRIUM RECONSTRUCTOR TESTS //
   //       with only one extrapolation cell              //
   /////////////////////////////////////////////////////////

   auto equilibriumAndNonEquilibriumSphereNormal1Reconstructor = lbm_mesapd_coupling::makeEquilibriumAndNonEquilibriumReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, 1);
   auto equilibriumAndNonEquilibriumSphereNormal1ReconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, equilibriumAndNonEquilibriumSphereNormal1Reconstructor, conserveMomentum);

   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Equilibrium and NonEquilibrium Reconstructor with one ext cell Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *equilibriumAndNonEquilibriumSphereNormal1ReconstructionManager, radius, velocity, density, conserveMomentum );
   }

   ///////////////////////////////////////
   // EXTRAPOLATION RECONSTRUCTOR TESTS //
   ///////////////////////////////////////

   auto extrapolationReconstructor = lbm_mesapd_coupling::makeExtrapolationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder);
   auto extrapolationReconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, extrapolationReconstructor, conserveMomentum);

   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Extrapolation Reconstructor Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *extrapolationReconstructionManager, radius, velocity, density, conserveMomentum );
   }

   ///////////////////////////////////////
   // EXTRAPOLATION RECONSTRUCTOR TESTS //
   //   with activated constrain        //
   ///////////////////////////////////////

   auto extrapolationConstrainReconstructor = lbm_mesapd_coupling::makeExtrapolationReconstructor<BoundaryHandling_T, lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder, true>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder);
   auto extrapolationConstrainReconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, extrapolationReconstructor, conserveMomentum);

   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Extrapolation Reconstructor with Constrain Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *extrapolationConstrainReconstructionManager, radius, velocity, density, conserveMomentum );
   }

   //////////////////////////////////////////////
   // GRADs MOMENT APPROXIMATION RECONSTRUCTOR //
   //////////////////////////////////////////////

   auto gradReconstructor = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, omega, false);
   auto gradReconstructorManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, gradReconstructor, conserveMomentum);

   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Grads Moment Approximation Reconstructor Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *gradReconstructorManager, radius, velocity, density, conserveMomentum );
   }

   //////////////////////////////////////////////
   // GRADs MOMENT APPROXIMATION RECONSTRUCTOR //
   //        with density recomputation        //
   //////////////////////////////////////////////

   auto gradReconstructorWithRecomp = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, omega, true);
   auto gradReconstructorWithRecompManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, gradReconstructorWithRecomp, conserveMomentum);

   for(auto & testSetup : testSetups)
   {
      checkReconstruction<BoundaryHandling_T>( "Grads Moment Approximation Reconstructor with Density Recomputation Test: " + std::get<0>(testSetup), std::get<1>(testSetup), std::get<2>(testSetup),
                                               blocks, ps, ss, accessor, pdfFieldID, boundaryHandlingID, particleFieldID,
                                               *gradReconstructorWithRecompManager, radius, velocity, density, conserveMomentum );
   }

   return 0;

}

} //namespace pdf_reconstruction

int main( int argc, char **argv ){
   pdf_reconstruction::main(argc, argv);
}
