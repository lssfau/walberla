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
//! \file   InfiniteDomain.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/domain/InfiniteDomain.h>

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

   auto randomPoint = Vec3(1.23_r,2.34_r,3.45_r);
   domain::InfiniteDomain domain;
   WALBERLA_CHECK(domain.isContainedInProcessSubdomain(static_cast<uint_t>(rank), randomPoint));
   WALBERLA_CHECK(!domain.isContainedInProcessSubdomain(static_cast<uint_t>(rank + 1), randomPoint));
   WALBERLA_CHECK_EQUAL(domain.findContainingProcessRank(randomPoint), rank);
   auto pt = randomPoint;
   domain.periodicallyMapToDomain(pt);
   WALBERLA_CHECK_IDENTICAL(pt, randomPoint);
   WALBERLA_CHECK_EQUAL(domain.getNeighborProcesses().size(), 0);
   WALBERLA_CHECK(domain.intersectsWithProcessSubdomain(static_cast<uint_t>(rank), randomPoint, 1_r));
   WALBERLA_CHECK(!domain.intersectsWithProcessSubdomain(static_cast<uint_t>(rank + 1), randomPoint, 1_r));
   domain.correctParticlePosition(pt);
   WALBERLA_CHECK_IDENTICAL(pt, randomPoint);
}

}

int main( int argc, char ** argv )
{
   walberla::main(argc, argv);
   return EXIT_SUCCESS;
}
