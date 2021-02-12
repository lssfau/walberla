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

#include <mesa_pd/kernel/cnt/WallContact.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {
namespace mesa_pd {

int main( int argc, char ** argv )
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
   data::SingleParticleAccessor p;
   p.setPosition(0, Vec3(0,0,0));

   //init kernels
   kernel::cnt::WallContact wallContact(0_r);

   auto calcForce = [&](const real_t zPos)
   {
      clear(p.getForceRef(0));
      p.getPositionRef(0)[2] = zPos;
      wallContact(0, p);
      WALBERLA_CHECK_FLOAT_EQUAL(wallContact.getLastForce(), p.getForce(0)[2]);
      return p.getForce(0)[2];
   };

   WALBERLA_LOG_INFO("checking that force is always repulsive");
   WALBERLA_CHECK_GREATER_EQUAL( calcForce(wallContact.z1 * 0.5_r + wallContact.r),
                                 0_r );
   WALBERLA_CHECK_GREATER_EQUAL( calcForce((wallContact.z1 + wallContact.z2) * 0.5_r + wallContact.r),
                                 0_r );
   WALBERLA_CHECK_GREATER_EQUAL( calcForce(wallContact.z2 * 2_r + wallContact.r),
                                 0_r );

   WALBERLA_CHECK_LESS_EQUAL( calcForce(-wallContact.z1 * 0.5_r - wallContact.r),
                              0_r );
   WALBERLA_CHECK_LESS_EQUAL( calcForce(-(wallContact.z1 + wallContact.z2) * 0.5_r - wallContact.r),
                              0_r );
   WALBERLA_CHECK_LESS_EQUAL( calcForce(-wallContact.z2 * 2_r - wallContact.r),
                              0_r );

   WALBERLA_LOG_INFO("checking that force is symmetric");
   WALBERLA_CHECK_FLOAT_EQUAL(calcForce(wallContact.z1 * 0.5_r + wallContact.r),
                              -calcForce(-wallContact.z1 * 0.5_r - wallContact.r) );
   WALBERLA_CHECK_FLOAT_EQUAL(calcForce((wallContact.z1 + wallContact.z2) * 0.5_r + wallContact.r),
                              -calcForce(-(wallContact.z1 + wallContact.z2) * 0.5_r - wallContact.r) );
   WALBERLA_CHECK_FLOAT_EQUAL(calcForce(wallContact.z2 * 2_r + wallContact.r),
                              -calcForce(-wallContact.z2 * 2_r - wallContact.r) );

   WALBERLA_LOG_INFO("checking smoothness of the force");
   WALBERLA_CHECK_FLOAT_EQUAL( calcForce(std::nextafter(wallContact.z1 + wallContact.r, 0_r)),
                               calcForce(std::nextafter(wallContact.z1 + wallContact.r, std::numeric_limits<real_t>::infinity())) );
   WALBERLA_CHECK_FLOAT_EQUAL( calcForce(std::nextafter(wallContact.z2 + wallContact.r, 0_r)),
                               calcForce(std::nextafter(wallContact.z2 + wallContact.r, std::numeric_limits<real_t>::infinity())) );

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
