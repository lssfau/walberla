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
//! \file   SyncGhostOwnersLarge.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/mpi/SyncGhostOwners.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/mpi/Reduce.h>

#include <iostream>
#include <memory>

namespace walberla {
namespace mesa_pd {

const real_t radius = real_t(20);

walberla::id_t createSphere(data::ParticleStorage& ps, domain::IDomain& domain)
{
   walberla::id_t uid = 0;
   auto owned = domain.isContainedInProcessSubdomain( uint_c(walberla::mpi::MPIManager::instance()->rank()), Vec3(0,0,0) );
   if (owned)
   {
      data::Particle&& p          = *ps.create();
      p.getPositionRef()          = Vec3(0,0,0);
      p.getInteractionRadiusRef() = radius;
      p.getRotationRef()          = Rot3(Quat());
      p.getLinearVelocityRef()    = Vec3(1,2,3);
      p.getAngularVelocityRef()   = Vec3(4,5,6);
      p.getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
      uid = p.getUid();
      WALBERLA_LOG_DETAIL("SPHERE CREATED");
   }

   walberla::mpi::allReduceInplace(uid, walberla::mpi::SUM);
   return uid;
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   //logging::Logging::instance()->setStreamLogLevel(logging::Logging::DETAIL);
   //logging::Logging::instance()->includeLoggingToFile("MESA_PD_Kernel_SyncGhostOwnersLarge");
   //logging::Logging::instance()->setFileLogLevel(logging::Logging::DETAIL);

   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(-5,-5,-5,25,25,25), // simulation domain
                                                 Vector3<uint_t>(3,3,3), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);
   std::array< bool, 3 > periodic;
   periodic[0] = forest->isPeriodic(0);
   periodic[1] = forest->isPeriodic(1);
   periodic[2] = forest->isPeriodic(2);

   //init data structures
   data::ParticleStorage ps(100);

   //initialize particle
   auto uid = createSphere(ps, domain);
   WALBERLA_LOG_DEVEL_ON_ROOT("uid: " << uid);

   //init kernels
   mpi::SyncGhostOwners SNN;

   SNN(ps, domain);
   SNN(ps, domain);
   SNN(ps, domain);
   SNN(ps, domain);

   std::vector<real_t> deltas {
      real_t(0),
            real_t(4.9),
            real_t(5.1),
            real_t(10),
            real_t(14.9),
            real_t(15.1),
            real_t(20),
            real_t(24.9)};

   for (auto delta : deltas)
   {
      WALBERLA_LOG_DEVEL(delta);
      auto pos = Vec3(1,0,0) * delta;
      WALBERLA_LOG_DETAIL("checking position: " << pos);
      // owner moves particle to new position
      auto pIt = ps.find(uid);
      if (pIt != ps.end())
      {
         if (!data::particle_flags::isSet(pIt->getFlags(), data::particle_flags::GHOST))
         {
            pIt->setPosition(pos);
         }
      }

      //sync
      SNN(ps, domain);

      //check
      if (sqDistancePointToAABB(pos, forest->begin()->getAABB()) <= radius * radius)
      {
         WALBERLA_CHECK_EQUAL(ps.size(), 1);
         if (forest->begin()->getAABB().contains(pos))
         {
            WALBERLA_CHECK(!data::particle_flags::isSet(ps.begin()->getFlags(), data::particle_flags::GHOST));
         } else
         {
            WALBERLA_CHECK(data::particle_flags::isSet(ps.begin()->getFlags(), data::particle_flags::GHOST));
         }
      } else
      {
         WALBERLA_CHECK_EQUAL(ps.size(), 0);
      }
   }


   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
