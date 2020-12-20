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
#include <mesa_pd/data/ShapeStorage.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   //init data structures
   data::ParticleStorage ps(100);

   //initialize particle
   const auto linVel = Vec3(1,2,3);
   const auto angVel = Vec3(1,2,3);

   const auto force  = Vec3(1,2,3);
   const auto torque = Vec3(1,2,3);

   data::Particle&& p        = *ps.create();
   p.getPositionRef()        = Vec3(0,0,0);
   p.getRotationRef()        = Rot3(Quat());
   p.getLinearVelocityRef()  = linVel;
   p.getAngularVelocityRef() = angVel;
   p.getForceRef()           = force;
   p.getTorqueRef()          = torque;

   //init kernels
   WALBERLA_LOG_DEVEL(p);

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
