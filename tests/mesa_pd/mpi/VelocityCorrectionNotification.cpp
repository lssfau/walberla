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


#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <mesa_pd/domain/BlockForestDomain.h>

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>

#include <mesa_pd/mpi/SyncNextNeighbors.h>
#include <mesa_pd/mpi/BroadcastProperty.h>
#include <mesa_pd/mpi/ReduceProperty.h>

#include <mesa_pd/mpi/notifications/VelocityCorrectionNotification.h>
#include <mesa_pd/mpi/notifications/VelocityUpdateNotification.h>

#include <iostream>

namespace walberla {

int main( int argc, char ** argv )
{

   using namespace mesa_pd;

   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   int np = walberla::mpi::MPIManager::instance()->numProcesses();

   // create blocks
   shared_ptr< blockforest::BlockForest > forest = blockforest::createBlockForest(
           math::AABB(-10,-10,-10, 10, 10, 10),
           Vector3<uint_t>(2,2,2),
           Vector3<bool>(false, false, false));

   domain::BlockForestDomain domain(forest);

   //init data structures
   data::ParticleStorage ps(10);

   //initialize particle
   const auto linVel = Vec3(1,2,3);
   const auto angVel = Vec3(-1,-2,-3);

   Vec3 pt(real_t(0.1), real_t(0.1),real_t(0.1));

   for (auto& iBlk : *forest)
   {
      if(iBlk.getAABB().contains(pt)) {
         auto p                       = ps.create();
         p->getPositionRef()          = pt;
         p->getInteractionRadiusRef() = real_t(1.0);
         p->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
         p->getTypeRef()              = 0;
         p->getLinearVelocityRef()    = linVel;
         p->getAngularVelocityRef()    = angVel;
         WALBERLA_LOG_INFO("Particle created.");
      }
   }

   // Sync (create ghost particles on the other processes)
   mesa_pd::mpi::SyncNextNeighbors snn;
   snn(ps, domain);

   auto relax_param = real_t(0.8);
   VelocityUpdateNotification::Parameters::relaxationParam = relax_param;
   mesa_pd::mpi::ReduceProperty reductionKernel;
   mesa_pd::mpi::BroadcastProperty broadcastKernel;

   WALBERLA_LOG_INFO("Particle-Storage size: " << ps.size());

   // Reduce dv
   // dv per process
   auto dvprocess = Vec3(real_t(0.1),real_t(0.1),real_t(0.1));
   auto dwprocess = Vec3(real_t(0.05),real_t(0.1),real_t(0.15));
   ps.setDv(0, dvprocess);
   ps.setDw(0, dwprocess);
   reductionKernel.operator()<VelocityCorrectionNotification>(ps);

   for (auto& iBlk : *forest)
   {
      if(iBlk.getAABB().contains(pt)) {
         WALBERLA_CHECK_FLOAT_EQUAL(ps.getLinearVelocity(0), linVel);
         WALBERLA_CHECK_FLOAT_EQUAL(ps.getAngularVelocity(0), angVel);
         WALBERLA_CHECK_FLOAT_EQUAL(ps.getDv(0), real_t(np) * dvprocess);
         WALBERLA_CHECK_FLOAT_EQUAL(ps.getDw(0), real_t(np) * dwprocess);
      }
   }

   broadcastKernel.operator()<VelocityUpdateNotification>(ps);


   // Broadcast v
   reductionKernel.operator()<VelocityCorrectionNotification>(ps);
   if(np > 1){
      WALBERLA_CHECK_FLOAT_EQUAL(ps.getLinearVelocity(0), linVel + relax_param * real_t(np) * dvprocess);
      WALBERLA_CHECK_FLOAT_EQUAL(ps.getAngularVelocity(0), angVel + relax_param * real_t(np) * dwprocess);
      WALBERLA_CHECK_FLOAT_EQUAL(ps.getDv(0), Vec3());
      WALBERLA_CHECK_FLOAT_EQUAL(ps.getDw(0), Vec3());
   }else{
      WALBERLA_CHECK_FLOAT_EQUAL(ps.getLinearVelocity(0) + ps.getDv(0), linVel + real_t(np) * dvprocess);
      WALBERLA_CHECK_FLOAT_EQUAL(ps.getAngularVelocity(0) + ps.getDw(0), angVel + real_t(np) * dwprocess);
   }

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
