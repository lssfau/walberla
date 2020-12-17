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
//! \file SyncEquivalence.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"


#include "vtk/VTKOutput.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "pe/basic.h"
#include "pe/synchronization/SyncNextNeighbors.h"
#include "CheckVitalParameters.h"

#include "core/grid_generator/HCPIterator.h"
#include "core/logging/Logging.h"
#include "core/timing/TimingTree.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"

namespace walberla {
using namespace walberla::pe;
using namespace walberla::timing;

using BodyTuple = std::tuple<Sphere> ;

struct BodyData
{
   BodyData(const walberla::id_t uid, const Vec3& pos)
      : uid_(uid)
      , pos_(pos)
   {}

   walberla::id_t uid_;
   Vec3 pos_;
};

bool comp(const BodyData& a, const BodyData& b)
{
   return a.uid_ < b.uid_;
}

bool compareOwner(const Owner& a, const Owner& b)
{
   return a.blockID_ < b.blockID_;
}
bool compareContact(const Contact& a, const Contact& b)
{
   if (a.getBody1()->getID() != b.getBody1()->getID())
   {
      return a.getBody1()->getID() < b.getBody1()->getID();
   } else
   {
      return a.getBody2()->getID() < b.getBody2()->getID();
   }
}

struct SimInfo
{
   shared_ptr<BodyStorage>             globalBodyStorage;
   shared_ptr<StructuredBlockForest>   forest;
   BlockDataID                         storageID;
   BlockDataID                         ccdID;
   BlockDataID                         fcdID;
   shared_ptr<cr::ICR>                 cr;
};

void createSimulation(math::AABB& simulationDomain,
                      Vector3<uint_t> blocks,
                      Vector3<uint_t> cellsPerBlock,
                      real_t spacing,
                      real_t radius,
                      real_t vMax,
                      SimInfo& info)
{
    std::mt19937 generator;
    generator.seed(1337);
    //math::seedRandomGenerator( 1337 ); //static_cast<std::mt19937::result_type>(1337 * mpi::MPIManager::instance()->worldRank()) );

    info.globalBodyStorage = make_shared<BodyStorage>();

    // create forest
    info.forest = blockforest::createUniformBlockGrid(
           simulationDomain,
           blocks[0], blocks[1], blocks[2],                        // number of blocks in x,y,z direction
           cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],   // how many cells per block (x,y,z)
           0,                                                      // max blocks per process
           false, false,                                           // include metis / force metis
           true, true, true );                                     // full periodicity

    // add block data
    info.storageID           = info.forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
    info.ccdID               = info.forest->addBlockData(ccd::createHashGridsDataHandling( info.globalBodyStorage, info.storageID ), "CCD");
    info.fcdID               = info.forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

    info.cr = std::make_shared<cr::HCSITS>(info.globalBodyStorage, info.forest->getBlockStoragePointer(), info.storageID, info.ccdID, info.fcdID );

    int numParticles = int_c(0);

    for (auto blkIt = info.forest->begin(); blkIt != info.forest->end(); ++blkIt)
    {
        IBlock & currentBlock = *blkIt;
        for (auto it = grid_generator::HCPIterator(currentBlock.getAABB(), Vector3<real_t>(-5,-5,-5), spacing); it != grid_generator::HCPIterator(); ++it)
        {
            SphereID sp = pe::createSphere( *(info.globalBodyStorage.get()), info.forest->getBlockStorage(), info.storageID, static_cast<walberla::id_t>(mpi::MPIManager::instance()->worldRank() * 1000000 + numParticles), *it, radius);
            Vec3 rndVel(math::realRandom<real_t>(-vMax, vMax, generator), math::realRandom<real_t>(-vMax, vMax, generator), math::realRandom<real_t>(-vMax, vMax, generator));
            if (sp != nullptr) sp->setLinearVel(rndVel);
            if (sp != nullptr) ++numParticles;
        }
    }
    WALBERLA_LOG_INFO_ON_ROOT(numParticles << " particles created on root");
    WALBERLA_LOG_PROGRESS_ON_ROOT("*** SETUP - END ***");
}

int main( int argc, char ** argv )
{
    mpi::Environment mpiEnv(argc, argv);
    WALBERLA_UNUSED(mpiEnv);

    logging::Logging::instance()->setStreamLogLevel(logging::Logging::WARNING);

    math::AABB simulationDomain   = math::AABB( Vec3(0,0,0), Vec3(10, 10, 10));
    Vector3<uint_t> blocks        = Vector3<uint_t>(2, 2, 2);
    Vector3<uint_t> cellsPerBlock = Vector3<uint_t>(1, 1, 1);

    real_t spacing = real_c(1.0);

    real_t radius = real_c(0.4);

    real_t vMax = real_c(1.0);

    int simulationSteps = 100;
    real_t dt = real_c(0.01);

    // initialize body type ids
    SetBodyTypeIDs<BodyTuple>::execute();

    SimInfo sim1;
    SimInfo sim2;

    createSimulation(simulationDomain, blocks, cellsPerBlock, spacing, radius, vMax, sim1);
    createSimulation(simulationDomain, blocks, cellsPerBlock, spacing, radius, vMax, sim2);

    WALBERLA_CHECK_EQUAL(sim1.forest->size(), sim2.forest->size());
    WALBERLA_CHECK_EQUAL(sim1.forest->numberOfBlockDataItems(), sim2.forest->numberOfBlockDataItems());

    std::cout << std::setprecision(10);

     // synchronize particles
     syncNextNeighbors<BodyTuple>( sim1.forest->getBlockForest(), sim1.storageID);
     syncNextNeighbors<BodyTuple>( sim1.forest->getBlockForest(), sim1.storageID);
     syncShadowOwners<BodyTuple>( sim2.forest->getBlockForest(), sim2.storageID);
     syncShadowOwners<BodyTuple>( sim2.forest->getBlockForest(), sim2.storageID);

     WALBERLA_LOG_PROGRESS_ON_ROOT("*** SIMULATION - START ***");
     for (int i=0; i < simulationSteps; ++i)
     {
        WALBERLA_LOG_WARNING_ON_ROOT(i);

        sim1.cr->timestep( real_c(dt) );
        syncNextNeighbors<BodyTuple>( sim1.forest->getBlockForest(), sim1.storageID);

        sim2.cr->timestep( real_c(dt) );
        syncShadowOwners<BodyTuple>( sim2.forest->getBlockForest(), sim2.storageID);

        auto blkIt2 = sim2.forest->begin();
        for (auto blkIt1 = sim1.forest->begin(); blkIt1 != sim1.forest->end(); ++blkIt1, ++blkIt2)
        {
            Storage * storage1 = blkIt1->getData< Storage >( sim1.storageID );
            Storage * storage2 = blkIt2->getData< Storage >( sim2.storageID );
            BodyStorage& localStorage1 = (*storage1)[0];
            BodyStorage& localStorage2 = (*storage2)[0];
            WALBERLA_CHECK_EQUAL((*storage1)[0].size(), (*storage2)[0].size());
            WALBERLA_CHECK_EQUAL((*storage1)[1].size(), (*storage2)[1].size());

            fcd::IFCD * fcd1 = blkIt1->getData< fcd::IFCD >( sim1.fcdID );
            fcd::IFCD * fcd2 = blkIt2->getData< fcd::IFCD >( sim2.fcdID );
            WALBERLA_CHECK_EQUAL(fcd1->getContacts().size(), fcd2->getContacts().size());
            WALBERLA_LOG_WARNING_ON_ROOT( "Collisions: " << fcd1->getContacts().size() );
//            std::sort(bodyIt1->MPITrait.beginShadowOwners(), bodyIt1->MPITrait.endShadowOwners(), compareOwner);
//            std::sort(bodyIt2->MPITrait.beginShadowOwners(), bodyIt2->MPITrait.endShadowOwners(), compareOwner);

//            auto contact2It = fcd2->getContacts().begin();
//            for (auto contact1It = fcd1->getContacts().begin(); contact1It != fcd1->getContacts().end(); ++contact1It, ++ contact2It)
//            {
//               WALBERLA_CHECK_EQUAL(contact1It->getID(), contact2It->getID());
//               WALBERLA_CHECK_EQUAL(contact1It->getBody1()->getID(), contact2It->getBody1()->getID());
//               WALBERLA_CHECK_EQUAL(contact1It->getBody2()->getID(), contact2It->getBody2()->getID());
//               WALBERLA_CHECK_IDENTICAL(contact1It->getPosition(), contact2It->getPosition());
//               WALBERLA_CHECK_IDENTICAL(contact1It->getNormal(), contact2It->getNormal());
//            }

            auto bodyIt2 = localStorage2.begin();
            for (auto bodyIt1 = localStorage1.begin(); bodyIt1 != localStorage1.end(); ++bodyIt1, ++bodyIt2)
            {
               WALBERLA_CHECK_EQUAL( bodyIt1->getID(), bodyIt2->getID() );

               WALBERLA_CHECK_EQUAL(bodyIt1->MPITrait.sizeShadowOwners(), bodyIt2->MPITrait.sizeShadowOwners());
               std::sort(bodyIt1->MPITrait.beginShadowOwners(), bodyIt1->MPITrait.endShadowOwners(), compareOwner);
               std::sort(bodyIt2->MPITrait.beginShadowOwners(), bodyIt2->MPITrait.endShadowOwners(), compareOwner);
               auto shadowOwnersIt2 = bodyIt2->MPITrait.beginShadowOwners();
               for (auto shadowOwnersIt1 = bodyIt1->MPITrait.beginShadowOwners(); shadowOwnersIt1 != bodyIt1->MPITrait.endShadowOwners(); ++shadowOwnersIt1, ++shadowOwnersIt2)
               {
                  WALBERLA_CHECK_EQUAL(shadowOwnersIt1->rank_, shadowOwnersIt2->rank_);
                  WALBERLA_CHECK_EQUAL(shadowOwnersIt1->blockID_, shadowOwnersIt2->blockID_);
               }

               checkVitalParameters( static_cast<SphereID>(bodyIt1.getBodyID()), static_cast<SphereID>(bodyIt2.getBodyID()) );

            }
        }
     }

    return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}