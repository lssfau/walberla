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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/math/Random.h>
#include <core/mpi/Reduce.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>

#include <memory>
#include <string>

namespace walberla {
namespace mesa_pd {

class ParticleAccessorWithShape : public data::ParticleAccessor
{
public:
   ParticleAccessorWithShape(std::shared_ptr<data::ParticleStorage>& ps, std::shared_ptr<data::ShapeStorage>& ss)
         : ParticleAccessor(ps)
         , ss_(ss)
   {}

   const auto& getInvMass(const size_t p_idx) const {return ss_->shapes[ps_->getShapeID(p_idx)]->getInvMass();}

   const auto& getInvInertiaBF(const size_t p_idx) const {return ss_->shapes[ps_->getShapeID(p_idx)]->getInvInertiaBF();}

   data::BaseShape* getShape(const size_t p_idx) const {return ss_->shapes[ps_->getShapeID(p_idx)].get();}
private:
   std::shared_ptr<data::ShapeStorage> ss_;
};

/*
 * Tests linked cells
 * Fully periodic, simple cubic sphere array is created
 * Together with linked cells of size 2 * sphere spacing.
 * This allows to analytically calculate how many particles should reside in each linked cell and how many contacts should be found and treated.
 */
int main( const int particlesPerAxisPerProcess = 2 )
{
   using namespace walberla::timing;

   walberla::mpi::MPIManager::instance()->resetMPI();

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   //   logging::Logging::instance()->setFileLogLevel(logging::Logging::DETAIL);
   //   logging::Logging::instance()->includeLoggingToFile("CollisionDetection");

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * walberla::mpi::MPIManager::instance()->worldRank()) );

   const real_t radius = real_t(0.9);
   const real_t generationSpacing = 1_r;
   const int linkedCellMultipleOfGenerationSpacing = 2;
   WALBERLA_CHECK_EQUAL(particlesPerAxisPerProcess % linkedCellMultipleOfGenerationSpacing, 0, "Only a multiple of " << linkedCellMultipleOfGenerationSpacing << " is allowed");

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   const int centerParticles  = particlesPerAxisPerProcess * particlesPerAxisPerProcess * particlesPerAxisPerProcess;
   const int faceParticles    = particlesPerAxisPerProcess * particlesPerAxisPerProcess * 6;
   const int edgeParticles    = particlesPerAxisPerProcess * 4 * 3;
   const int cornerParticles  = 8;

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particlesPerAxisPerProcess);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(centerParticles);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(faceParticles);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(edgeParticles);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cornerParticles);

   // create forest
   const real_t eps = real_c(0.01); //shift to make contact points unambigous
   Vector3<uint_t>blocksPerDirection (2,2,2);
   int numBlocks = int(blocksPerDirection[0] * blocksPerDirection[1] * blocksPerDirection[2]);
   Vector3<bool> periodicity(true, true, true);
   auto forest = blockforest::createBlockForest( math::AABB(real_t(0),
                                                            real_t(0),
                                                            real_t(0),
                                                            real_c(particlesPerAxisPerProcess) * generationSpacing * real_c(blocksPerDirection[0]),
                                                            real_c(particlesPerAxisPerProcess) * generationSpacing * real_c(blocksPerDirection[1]),
                                                            real_c(particlesPerAxisPerProcess) * generationSpacing * real_c(blocksPerDirection[2])).getTranslated(Vec3(eps,eps,eps)),
                                                 blocksPerDirection, periodicity);
   domain::BlockForestDomain domain(forest);

   auto localDomain = forest->begin()->getAABB();
   WALBERLA_CHECK_EQUAL(forest->size(), 1, "please run with 8 processes -> 1 process per block");

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape accessor(ps, ss);
   const real_t linkedCellSize = real_c(linkedCellMultipleOfGenerationSpacing) * generationSpacing;
   data::LinkedCells lc(localDomain.getExtended(linkedCellSize), linkedCellSize );

   auto  smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));
   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(), Vector3<real_t>(generationSpacing, generationSpacing, generationSpacing) * real_c(0.5), generationSpacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));

         auto p                       = ps->create();
         p->setPosition( pt );
         p->setInteractionRadius( radius );
         p->setShapeID( smallSphere );
         p->setOwner( walberla::mpi::MPIManager::instance()->rank() );
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::allReduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   // Init kernels
   kernel::InsertParticleIntoLinkedCells ipilc;
   mpi::SyncNextNeighbors                SNN;

   // initial sync
   WALBERLA_CHECK_EQUAL(ps->size(), centerParticles);
   SNN(*ps, domain);
   int numKnownParticles = centerParticles + faceParticles + edgeParticles + cornerParticles;
   WALBERLA_CHECK_EQUAL(ps->size(), numKnownParticles);

   lc.clear();
   ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, lc);

   int interiorCellsPerDirection = particlesPerAxisPerProcess / linkedCellMultipleOfGenerationSpacing;

   int numLinkedCellsPerProcess = (interiorCellsPerDirection+2) * (interiorCellsPerDirection+2) * (interiorCellsPerDirection+2);
   WALBERLA_CHECK_EQUAL(lc.cells_.size(), numLinkedCellsPerProcess);

   const int numLinkedCellsCenter = interiorCellsPerDirection * interiorCellsPerDirection * interiorCellsPerDirection;
   const int numLinkedCellsFace = interiorCellsPerDirection * interiorCellsPerDirection * 6;
   const int numLinkedCellsEdge = interiorCellsPerDirection * 12;
   const int numLinkedCellsCorner = 8;

   const int numParticlesPerCenterCell = centerParticles / numLinkedCellsCenter;
   const int numParticlesPerFaceCell = faceParticles / numLinkedCellsFace;
   const int numParticlesPerEdgeCell = edgeParticles / numLinkedCellsEdge;
   const int numParticlesPerCornerCell = cornerParticles / numLinkedCellsCorner;

   int cellsCenter = 0;
   int cellsFace = 0;
   int cellsEdge = 0;
   int cellsCorner = 0;

   for (const auto& idx : lc.cells_)
   {
      int particleCounter = 0;
      int p_idx = idx;
      while (p_idx != -1)
      {
         ++particleCounter;
         p_idx = ps->getNextParticle(uint_c(p_idx));
      }

      if (particleCounter == numParticlesPerCenterCell) ++cellsCenter;
      else if (particleCounter == numParticlesPerFaceCell) ++cellsFace;
      else if (particleCounter == numParticlesPerEdgeCell) ++cellsEdge;
      else if (particleCounter == numParticlesPerCornerCell) ++cellsCorner;
      else WALBERLA_CHECK(false, "Linked cell found with unexpected number of particles " << particleCounter);

   }
   WALBERLA_CHECK_EQUAL(cellsCenter, numLinkedCellsCenter);
   WALBERLA_CHECK_EQUAL(cellsFace, numLinkedCellsFace);
   WALBERLA_CHECK_EQUAL(cellsEdge, numLinkedCellsEdge);
   WALBERLA_CHECK_EQUAL(cellsCorner, numLinkedCellsCorner);

   std::atomic<int64_t> contactsChecked (0);
   std::atomic<int64_t> contactsDetected(0);
   std::atomic<int64_t> contactsTreated (0);
   lc.forEachParticlePairHalf(true,
                              kernel::SelectAll(),
                              accessor,
                              [&](const size_t idx1, const size_t idx2, auto& ac)
   {
      collision_detection::AnalyticContactDetection acd;
      kernel::DoubleCast double_cast;
      mpi::ContactFilter contact_filter;
      ++contactsChecked;
      if (double_cast(idx1, idx2, ac, acd, ac ))
      {
         ++contactsDetected;
         if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
         {
            ++contactsTreated;
         }
      }
   },
   accessor );

   WALBERLA_LOG_DEVEL_ON_ROOT( "contacts checked/detected/treated: " << contactsChecked << " / " << contactsDetected << " / " << contactsTreated );

   WALBERLA_CHECK_LESS_EQUAL(contactsChecked, numKnownParticles * numKnownParticles / 2 ); // has to be better than n^2 checking
   WALBERLA_CHECK_EQUAL(contactsDetected, (centerParticles*26 + faceParticles*17 + edgeParticles*11 + cornerParticles*7) / 2 );
   const int totalNumberPairWiseContactsInDomain = 26 * int(numParticles) / 2; //fully-periodic simple cubic grid
   WALBERLA_CHECK_EQUAL(contactsTreated, totalNumberPairWiseContactsInDomain / numBlocks);

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mesa_pd::main( 2 );
   walberla::mesa_pd::main( 4 );
   return EXIT_SUCCESS;
}
