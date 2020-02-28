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
//! \file   ReduceProperty.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>
#include <mesa_pd/mpi/notifications/HeatFluxNotification.h>
#include <mesa_pd/mpi/ReduceProperty.h>
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

const real_t radius = real_t(3);

walberla::id_t createSphere(data::ParticleStorage& ps, domain::IDomain& domain)
{
   walberla::id_t uid = 0;
   auto owned = domain.isContainedInProcessSubdomain( uint_c(walberla::mpi::MPIManager::instance()->rank()), Vec3(9,9,9) );
   if (owned)
   {
      data::Particle&& p          = *ps.create();
      p.getPositionRef()          = Vec3(9,9,9);
      p.getInteractionRadiusRef() = radius;
      p.getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
      uid = p.getUid();
   }

   walberla::mpi::allReduceInplace(uid, walberla::mpi::SUM);
   return uid;
}

void main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   //init domain partitioning
   auto forest = blockforest::createBlockForest( AABB(0,0,0,20,20,20), // simulation domain
                                                 Vector3<uint_t>(2,2,2), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);

   //init data structures
   data::ParticleStorage ps(100);

   //initialize particle
   auto uid = createSphere(ps, domain);
   WALBERLA_LOG_DEVEL_ON_ROOT("uid: " << uid);

   //init kernels
   mpi::ReduceProperty    RP;
   mpi::SyncNextNeighbors SNN;

   //sync
   SNN(ps, domain);

   auto pIt = ps.find(uid);
   if (pIt != ps.end())
   {
      pIt->getForceRef() += Vec3(real_c(walberla::mpi::MPIManager::instance()->rank()));
      pIt->getTorqueRef() += Vec3(real_c(walberla::mpi::MPIManager::instance()->rank()));
      pIt->getHeatFluxRef() += real_c(walberla::mpi::MPIManager::instance()->rank());
   }

   RP.operator()<ForceTorqueNotification>(ps);
   RP.operator()<HeatFluxNotification>(ps);

   if (walberla::mpi::MPIManager::instance()->rank() == 0)
   {
      WALBERLA_CHECK_FLOAT_EQUAL( pIt->getForce(), Vec3(real_t(28)) );
      WALBERLA_CHECK_FLOAT_EQUAL( pIt->getTorque(), Vec3(real_t(28)) );
      WALBERLA_CHECK_FLOAT_EQUAL( pIt->getHeatFlux(), real_t(28) );
   } else
   {
      WALBERLA_CHECK_FLOAT_EQUAL( pIt->getForce(), Vec3(0) );
      WALBERLA_CHECK_FLOAT_EQUAL( pIt->getTorque(), Vec3(0) );
      WALBERLA_CHECK_FLOAT_EQUAL( pIt->getHeatFlux(), real_t(0) );
   }
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   walberla::mesa_pd::main(argc, argv);
   return EXIT_SUCCESS;
}
