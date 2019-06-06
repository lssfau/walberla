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
//! \file   VelocityVerlet.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/VelocityVerlet.h>

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

   const real_t k = real_t(0.1);

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   data::ParticleAccessor accessor(ps);

   const Vec3 startingVelocity(real_t(1),real_t(2),real_t(3));

   //initialize particle
   data::Particle&& p        = *ps->create();
   p.getPositionRef()        = Vec3(0,0,0);
   p.getInvMassRef()         = real_t(1.1);
   p.getLinearVelocityRef()  = startingVelocity;
   p.getForceRef()           = Vec3(0,0,0);
   p.getOldForceRef()        = Vec3(0,0,0);

   const real_t dt = real_t(0.1);
   const real_t w = std::sqrt(k*accessor.getInvMass(0));
   const Vec3 A = accessor.getLinearVelocity(0) / w;
   WALBERLA_LOG_DEVEL_VAR(w);
   WALBERLA_LOG_DEVEL_VAR(A);

   auto analyticPosition = [A, w](const real_t t){return A * std::sin(w*t);};
   auto analyticVelocity = [A, w](const real_t t){return A * std::cos(w*t) * w;};

   //init kernels
   kernel::VelocityVerletPreForceUpdate  preForce( dt );
   kernel::VelocityVerletPostForceUpdate postForce( dt );

   p.getPositionRef()        = analyticPosition(-dt);
   p.getLinearVelocityRef()  = analyticVelocity(-dt);
   p.getOldForceRef()        = - k * p.getPosition();
   ps->forEachParticle(false, kernel::SelectAll(), accessor, preForce, accessor);
   p.getLinearVelocityRef()  = analyticVelocity(real_t(0));
   p.getForceRef()           = - k * p.getPosition();
   ps->forEachParticle(false, kernel::SelectAll(), accessor, postForce, accessor);

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getPosition(), Vec3(0), real_t(1e-2));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getLinearVelocity(), startingVelocity, real_t(1e-2));

   for (auto i = 1; i < 500; ++i)
   {
      ps->forEachParticle(false, kernel::SelectAll(), accessor, preForce, accessor);
      p.getForceRef()           = - k * p.getPosition();
      ps->forEachParticle(false, kernel::SelectAll(), accessor, postForce, accessor);

      //check force
      WALBERLA_CHECK_FLOAT_EQUAL(p.getForce(), Vec3(0), p);

      //check velocity
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getLinearVelocity(),
                                         analyticVelocity(real_c(i) * dt),
                                         real_t(1e-2),
                                         "iteration: " << i << "\n" <<
                                         "t: " << real_c(i)*dt << "\n" <<
                                         p);

      //check position
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getPosition(),
                                         analyticPosition(real_c(i) * dt),
                                         real_t(1e-2),
                                         "iteration: " << i << "\n" <<
                                         "t: " << real_c(i)*dt << "\n" <<
                                         p);
//      WALBERLA_LOG_DEVEL_VAR(p.getPosition());
//      WALBERLA_LOG_DEVEL_VAR(analyticPosition(real_c(i) * dt));
   }

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
