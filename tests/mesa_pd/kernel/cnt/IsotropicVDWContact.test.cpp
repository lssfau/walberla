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
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <mesa_pd/kernel/cnt/IsotropicVDWContact.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {
namespace mesa_pd {

int main(int argc, char **argv)
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   if (std::is_same<real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(2);
   data::Particle &&p1 = *ps->create();
   data::Particle &&p2 = *ps->create();

   data::ParticleAccessor accessor(ps);

   //init kernels
   kernel::cnt::IsotropicVDWContact vdWContact;

   auto calcForce = [&](const Vec3 pos1, const Vec3 pos2) {
      p1.setPosition(pos1);
      p2.setPosition(pos2);
      clear(p1.getForceRef());
      clear(p2.getForceRef());
      vdWContact(0, 1, accessor);
      return p1.getForce();
   };

   const Vec3 randomNormal = Vec3(1_r, 2_r, 3_r).getNormalized();

   WALBERLA_LOG_INFO("checking repulsion - equilibrium - attraction");
   vdWContact(0, 1, accessor);
   WALBERLA_CHECK_LESS       (dot(randomNormal, calcForce(Vec3(0), randomNormal * (vdWContact.equilibriumDistance() - 1_r))),
                              0_r);
   WALBERLA_CHECK_FLOAT_EQUAL(dot(randomNormal, calcForce(Vec3(0), randomNormal * vdWContact.equilibriumDistance())),
                              0_r);
   WALBERLA_CHECK_GREATER    (dot(randomNormal, calcForce(Vec3(0), randomNormal * (vdWContact.equilibriumDistance() + 1_r))),
                              0_r);

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main(int argc, char **argv)
{
   return walberla::mesa_pd::main(argc, argv);
}
