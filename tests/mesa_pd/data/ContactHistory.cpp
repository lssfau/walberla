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
//! \file   ContactHistory.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ContactHistory.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <algorithm>
#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

void basic_test()
{
   //init data structures
   data::ContactHistory cs;

   cs.setImpactVelocityMagnitude(1.23456_r);
   cs.setIsSticking(true);
   cs.setTangentialSpringDisplacement(Vec3(1.23_r,2.345_r,3.56_r));

   mpi::SendBuffer sb;
   sb << cs;
   mpi::RecvBuffer rb(sb);

   data::ContactHistory cs_recv;
   rb >> cs_recv;

   WALBERLA_CHECK_IDENTICAL(cs.getImpactVelocityMagnitude(), cs_recv.getImpactVelocityMagnitude());
   WALBERLA_CHECK_IDENTICAL(cs.getIsSticking(), cs_recv.getIsSticking());
   WALBERLA_CHECK_IDENTICAL(cs.getTangentialSpringDisplacement(), cs_recv.getTangentialSpringDisplacement());

   WALBERLA_LOG_DEVEL( cs );
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   basic_test();

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
