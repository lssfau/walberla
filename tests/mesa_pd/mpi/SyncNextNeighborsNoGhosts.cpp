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
//! \file   SyncNextNeighborsNoGhosts.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/mpi/SyncNextNeighborsNoGhosts.h>

#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

namespace walberla {
namespace mesa_pd {

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   //init domain partitioning
   auto forest = blockforest::createBlockForest( math::AABB(0,0,0,10,10,10), // simulation domain
                                                 Vector3<uint_t>(2,1,1), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);

   data::ParticleStorage ps(100);

   Vec3 pt(2.5, 2.5, 2.5);
   if (forest->begin()->getAABB().contains(pt))
   {
      auto pIt = ps.create();
      pIt->setPosition(pt);
      pIt->setInteractionRadius(real_t(0));
      pIt->setOwner(walberla::mpi::MPIManager::instance()->rank());
   }

   if (forest->begin()->getAABB().contains(pt))
   {
      WALBERLA_CHECK_EQUAL(ps.size(), 1);
   } else
   {
      WALBERLA_CHECK_EQUAL(ps.size(), 0);
   }

   for (auto p : ps)
   {
      p.setPosition(Vec3(7.5, 2.5, 2.5));
   }

   mpi::SyncNextNeighborsNoGhosts SNN;

   SNN(ps, domain);

   if (forest->begin()->getAABB().contains(pt))
   {
      WALBERLA_CHECK_EQUAL(ps.size(), 0);
   } else
   {
      WALBERLA_CHECK_EQUAL(ps.size(), 1);
   }

   return 0;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
