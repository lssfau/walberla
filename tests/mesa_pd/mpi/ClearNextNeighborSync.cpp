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
//! \file ClearNextNeighborSync.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/mpi/ClearNextNeighborSync.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <blockforest/Initialization.h>
#include <core/debug/TestSubsystem.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/mpi/Reduce.h>
#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {

int main(int argc, char **argv)
{
   walberla::debug::enterTestMode();
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   const real_t spacing = real_c(1);
   const real_t radius = real_c(0.5);

   WALBERLA_LOG_INFO_ON_ROOT("*** MESA_PD ***");
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ac = data::ParticleAccessor(ps);
   auto ss = std::make_shared<data::ShapeStorage>();

   auto smallSphere = ss->create<data::Sphere>(radius);
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   auto forest = blockforest::createBlockForest(math::AABB(real_t(0),
                                                           real_t(0),
                                                           real_t(0),
                                                           real_t(6),
                                                           real_t(6),
                                                           real_t(6)),
                                                Vector3<uint_t>(2, 2, 2),
                                                Vector3<bool>(true, true, true));
   domain::BlockForestDomain domain(forest);

   WALBERLA_CHECK_EQUAL(forest->size(), 1, "please run with 8 processes -> 1 process per block");

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");
   for (auto &iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.2),
                                            spacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));

         auto p = ps->create();
         p->setPosition(pt);
         p->setInteractionRadius(radius);
         p->setShapeID(smallSphere);
         p->setOwner(walberla::mpi::MPIManager::instance()->rank());
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   mesa_pd::mpi::ClearNextNeighborSync CSNN;
   mesa_pd::mpi::SyncNextNeighbors SNN;
   SNN(*ps, domain);
   SNN(*ps, domain);
   SNN(*ps, domain);
   CSNN(ac);

   for(auto p : *ps)
   {
      using namespace mesa_pd::data::particle_flags;
      WALBERLA_CHECK(!isSet(p.getFlags(), GHOST));
      WALBERLA_CHECK(p.getGhostOwners().empty());
   }

   return EXIT_SUCCESS;
}
} //namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::mesa_pd::main( argc, argv );
}
