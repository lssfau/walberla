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
//
//======================================================================================================================

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>
#include <core/mpi/Reduce.h>

#include <iostream>
#include <memory>

namespace walberla {
namespace mesa_pd {

void periodicSync()
{
   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(0,0,0,2,2,2), // simulation domain
                                                 Vector3<uint_t>(2,2,2), // blocks in each direction
                                                 Vector3<bool>(true, true, true) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);

   //init data structures
   data::ParticleStorage ps(100);

   //initialize particles
   real_t spacing = real_c(1);
   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));

         auto p                       = ps.create();
         p->getPositionRef()          = pt;
         p->getInteractionRadiusRef() = spacing * real_t(0.5);
         p->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
         p->getTypeRef()              = 0;
         p->getLinearVelocityRef()    = (forest->getDomain().center() - p->getPosition()).getNormalized() * real_t(0.1);
      }
   }

   //init kernels
   mpi::SyncNextNeighbors SNN;
   SNN(ps, domain);

   for (auto i = 0; i < 100; ++i)
   {
      for (auto p : ps)
      {
         p.getPositionRef() += p.getLinearVelocity();
      }
      SNN(ps, domain);

      for (auto p : ps)
      {
         using namespace data::particle_flags;
         if (isSet(p->getFlags(), GHOST))
         {
            WALBERLA_CHECK_UNEQUAL(p.getOwner(), walberla::mpi::MPIManager::instance()->rank(), p);
         } else
         {
            WALBERLA_CHECK_EQUAL(p.getOwner(), walberla::mpi::MPIManager::instance()->rank(), p);
         }
      }
   }

   WALBERLA_CHECK_UNEQUAL(ps.size(), 0);
}

} //namespace mesa_pd
} //namespace walberla

/**
 * Particles on a regular grid with fly through the domain
 *
 * The domain is a 2x2x2 grid. Particles get initialized to cross the domain.
 * Correct migration is checked.
 * Also fewer processes than blocks is checked.
 */
int main( int argc, char ** argv )
{
   using namespace walberla;
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   walberla::mesa_pd::periodicSync();

   return EXIT_SUCCESS;
}
