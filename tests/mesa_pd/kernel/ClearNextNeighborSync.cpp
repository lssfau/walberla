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
//! \file   ClearNextNeighborSync.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/mpi/ClearNextNeighborSync.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/mpi/Reduce.h>

#include <iostream>
#include <memory>

namespace walberla {
namespace mesa_pd {

const real_t radius = real_t(1);

void createSphere(data::ParticleStorage& ps, domain::IDomain& domain, const Vec3& pos)
{
   auto owned = domain.isContainedInProcessSubdomain( uint_c(walberla::mpi::MPIManager::instance()->rank()), pos );
   if (owned)
   {
      data::Particle&& p          = *ps.create();
      p.getPositionRef()          = pos;
      p.getInteractionRadiusRef() = radius;
      p.getRotationRef()          = Rot3(Quat());
      p.getLinearVelocityRef()    = Vec3(1,2,3);
      p.getAngularVelocityRef()   = Vec3(4,5,6);
      p.getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   }
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(0,0,0,10,5,5), // simulation domain
                                                 Vector3<uint_t>(2,1,1), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);
   std::array< bool, 3 > periodic;
   periodic[0] = forest->isPeriodic(0);
   periodic[1] = forest->isPeriodic(1);
   periodic[2] = forest->isPeriodic(2);

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage> (100);
   data::ParticleAccessor ac(ps);

   //initialize particle
   createSphere(*ps, domain, Vec3(real_t(4.5),3,3));
   createSphere(*ps, domain, Vec3(real_t(5.5),3,3));

   //init kernels
   mpi::ClearNextNeighborSync CNNS;
   mpi::SyncNextNeighbors     SNN;

   //init
   WALBERLA_CHECK_EQUAL(ps->size(), 1);
   WALBERLA_CHECK_EQUAL(ps->getGhostOwners(0).size(), 0);

   //sync
   SNN(*ps, domain);

   WALBERLA_CHECK_EQUAL(ps->size(), 2);
   WALBERLA_CHECK_EQUAL(ps->getGhostOwners(0).size(), 1);
   WALBERLA_CHECK(data::particle_flags::isSet(ps->getFlags(1), data::particle_flags::GHOST));

   //clear
   CNNS(ac);

   WALBERLA_CHECK_EQUAL(ps->size(), 1);
   WALBERLA_CHECK_EQUAL(ps->getGhostOwners(0).size(), 0);

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
