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
//! \file   BlockForestDomain.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/domain/BlockForestDomain.h>

#include <blockforest/loadbalancing/StaticCurve.h>
#include <blockforest/Initialization.h>
#include <blockforest/SetupBlockForest.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

void main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();
   auto rank = mpi::MPIManager::instance()->rank();

   //   logging::Logging::instance()->setStreamLogLevel(logging::Logging::DETAIL);
   //   logging::Logging::instance()->includeLoggingToFile("MESA_PD_Kernel_SyncNextNeighbor");
   //   logging::Logging::instance()->setFileLogLevel(logging::Logging::DETAIL);

   //init domain partitioning
   auto sforest = std::make_unique<blockforest::SetupBlockForest>( );
   sforest->addWorkloadMemorySUIDAssignmentFunction( blockforest::uniformWorkloadAndMemoryAssignment );
   sforest->init( math::AABB(0,0,0,10,10,10), 2, 1, 1 , true, true, true);

   WALBERLA_LOG_INFO_ON_ROOT( "Balancing " << sforest->getNumberOfBlocks() << " blocks for " << 2 << " processes...");

   sforest->balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), 2, real_t(0), memory_t(0), false, true );

   auto forest = std::make_shared< BlockForest >( uint_c( walberla::MPIManager::instance()->rank() ), *sforest, false );

   domain::BlockForestDomain domain(forest);

   auto neighbors = domain.getNeighborProcesses();
   WALBERLA_CHECK_EQUAL(neighbors.size(), 1);
   WALBERLA_CHECK_EQUAL(neighbors[0], rank == 0 ? 1 : 0);

   WALBERLA_CHECK_EQUAL(domain.findContainingProcessRank(Vec3(2,2,2)), 0);
   WALBERLA_CHECK_EQUAL(domain.findContainingProcessRank(Vec3(7,7,7)), 1);

   WALBERLA_CHECK(domain.isContainedInProcessSubdomain(0, Vec3(2,2,2)));
   WALBERLA_CHECK(domain.isContainedInProcessSubdomain(1, Vec3(7,7,7)));

   WALBERLA_CHECK_EQUAL(domain.isContainedInProcessSubdomain(Vec3(2,2,2), real_t(1)), rank == 0 ? true : false);
   WALBERLA_CHECK_EQUAL(domain.isContainedInProcessSubdomain(Vec3(7,7,7), real_t(1)), rank == 0 ? false : true);

   WALBERLA_CHECK_EQUAL(domain.isContainedInProcessSubdomain(Vec3(real_t(4.5),2,2), real_t(1)), rank == 0 ? false : false);
   WALBERLA_CHECK_EQUAL(domain.isContainedInProcessSubdomain(Vec3(real_t(5.5),7,7), real_t(1)), rank == 0 ? false : false);

   if (rank != 0)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(0, Vec3(real_t(2),2,2), real_t(1)), true);
   if (rank != 1)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(1, Vec3(real_t(2),2,2), real_t(1)), false);

   if (rank != 0)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(0, Vec3(real_t(7),2,2), real_t(1)), false);
   if (rank != 1)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(1, Vec3(real_t(7),2,2), real_t(1)), true);

   if (rank != 0)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(0, Vec3(real_t(5.5),2,2), real_t(1)), true);
   if (rank != 1)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(1, Vec3(real_t(5.5),2,2), real_t(1)), true);

   if (rank != 0)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(0, Vec3(real_t(4.5),2,2), real_t(1)), true);
   if (rank != 1)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(1, Vec3(real_t(4.5),2,2), real_t(1)), true);

   if (rank != 0)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(0, Vec3(real_t(9.5),2,2), real_t(1)), true);
   if (rank != 1)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(1, Vec3(real_t(9.5),2,2), real_t(1)), true);

   if (rank != 0)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(0, Vec3(real_t(0.5),2,2), real_t(1)), true);
   if (rank != 1)
      WALBERLA_CHECK_EQUAL(domain.intersectsWithProcessSubdomain(1, Vec3(real_t(0.5),2,2), real_t(1)), true);
}

}

int main( int argc, char ** argv )
{
   walberla::main(argc, argv);
   return EXIT_SUCCESS;
}
