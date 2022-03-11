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

#include <mesa_pd/data/ParticleAccessorWithShape.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Sphere.h>

#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/VelocityVerlet.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <cmath>
#include <fstream>
#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

class ParticleAccessorWithShape : public data::ParticleAccessor
{
public:
   ParticleAccessorWithShape(std::shared_ptr<data::ParticleStorage>& ps)
      : ParticleAccessor(ps)
   {
      sp.updateMassAndInertia(real_t(1234));
   }

   const auto& getInvMass(const size_t /*p_idx*/) const {return sp.getInvMass();}

   const auto& getInvInertiaBF(const size_t /*p_idx*/) const {return sp.getInvInertiaBF();}
   const auto& getInertiaBF(const size_t /*p_idx*/) const {return sp.getInertiaBF();}

   data::BaseShape* getShape(const size_t /*p_idx*/) {return &sp;}
private:
   data::Sphere sp = data::Sphere(real_t(0.6));
};

int main()
{
   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   ParticleAccessorWithShape accessor(ps);

   const Vec3 startingVelocity(real_t(1),real_t(2),real_t(3));

   //initialize particle
   data::Particle&& p        = *ps->create();
   p.getPositionRef()        = Vec3(0,0,0);
   p.getLinearVelocityRef()  = startingVelocity;
   p.getAngularVelocityRef() = startingVelocity;
   p.getForceRef()           = Vec3(0,0,0);
   p.getOldForceRef()        = Vec3(0,0,0);

   const real_t dt = real_t(0.001);

   const real_t k = real_t(10000);
   const real_t w = std::sqrt(k*accessor.getInvMass(0));
   const Vec3   A = accessor.getLinearVelocity(0) / w;
   WALBERLA_LOG_DEVEL_VAR(w);
   WALBERLA_LOG_DEVEL_VAR(A);

   auto analyticPosition = [A, w](const real_t t){return A * std::sin(w*t);};
   auto analyticVelocity = [A, w](const real_t t){return A * std::cos(w*t) * w;};

   const real_t kappa = real_t(1000);
   const real_t omega = std::sqrt(kappa*accessor.getInvInertiaBF(0)[0]);
   const Vec3   Alpha = accessor.getAngularVelocity(0) / omega;
   WALBERLA_LOG_DEVEL_VAR(omega);
   WALBERLA_LOG_DEVEL_VAR(Alpha);

   auto analyticRotation = [Alpha, omega](const real_t t){return Alpha * std::sin(omega*t);};
   auto analyticAngVel   = [Alpha, omega](const real_t t){return Alpha * std::cos(omega*t) * omega;};

   //init kernels
   kernel::VelocityVerletPreForceUpdate  preForce( dt );
   kernel::VelocityVerletPostForceUpdate postForce( dt );

   p.getPositionRef()        = analyticPosition(-dt);
   p.getLinearVelocityRef()  = analyticVelocity(-dt);
   p.getOldForceRef()        = - k * p.getPosition();

   p.getRotationRef()        = Rot3(Quat(analyticRotation(-dt), analyticRotation(-dt).length()));
   p.getAngularVelocityRef() = analyticAngVel(-dt);
   p.getOldTorqueRef()       = - kappa * p.getRotation().getQuaternion().getAxis() * p.getRotation().getQuaternion().getAngle();

   ps->forEachParticle(false, kernel::SelectAll(), accessor, preForce, accessor);

   p.getLinearVelocityRef()  = analyticVelocity(real_t(0));
   p.getForceRef()           = - k * p.getPosition();

   p.getAngularVelocityRef() = analyticAngVel(real_t(0));
   p.getTorqueRef()          = - kappa * p.getRotation().getQuaternion().getAxis() * p.getRotation().getQuaternion().getAngle();
   ps->forEachParticle(false, kernel::SelectAll(), accessor, postForce, accessor);

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getPosition(), Vec3(0), real_t(1e-2));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getLinearVelocity(), startingVelocity, real_t(1e-2));

   std::fstream fout("VelocityVerlet.txt", std::ofstream::out);
   for (auto i = 1; i < 3000; ++i)
   {
      ps->forEachParticle(false, kernel::SelectAll(), accessor, preForce, accessor);
      p.getForceRef()           = - k * p.getPosition();
      p.getTorqueRef()          = - kappa * p.getRotation().getQuaternion().getAxis() * p.getRotation().getQuaternion().getAngle();
      ps->forEachParticle(false, kernel::SelectAll(), accessor, postForce, accessor);

      //check force&torque
      WALBERLA_CHECK_FLOAT_EQUAL(p.getForce(), Vec3(0), p);
      WALBERLA_CHECK_FLOAT_EQUAL(p.getTorque(), Vec3(0), p);

      //check velocity
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getLinearVelocity(),
                                         analyticVelocity(real_c(i) * dt),
                                         real_t(1e-4),
                                         "iteration: " << i << "\n" <<
                                         "t: " << real_c(i)*dt << "\n" <<
                                         p);

      //check angular velocity
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getAngularVelocity(),
                                         analyticAngVel(real_c(i) * dt),
                                         real_t(1e-4),
                                         "iteration: " << i << "\n" <<
                                         "t: " << real_c(i)*dt << "\n" <<
                                         p);
      //            WALBERLA_LOG_DEVEL_VAR(p.getAngularVelocity());
      //            WALBERLA_LOG_DEVEL_VAR(analyticAngVel(real_c(i) * dt));

      //check position
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getPosition(),
                                         analyticPosition(real_c(i) * dt),
                                         real_t(1e-4),
                                         "iteration: " << i << "\n" <<
                                         "t: " << real_c(i)*dt << "\n" <<
                                         p);

      //check rotation
      auto angSol = Rot3(Quat(analyticRotation(real_c(i) * dt), analyticRotation(real_c(i) * dt).length()));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(p.getRotation().getQuaternion().getAxis() * p.getRotation().getQuaternion().getAngle(),
                                         angSol.getQuaternion().getAxis() * angSol.getQuaternion().getAngle(),
                                         real_t(1e-4),
                                         "iteration: " << i << "\n" <<
                                         "t: " << real_c(i)*dt << "\n" <<
                                         p);
      //      WALBERLA_LOG_DEVEL_VAR(p.getPosition());
      //      WALBERLA_LOG_DEVEL_VAR(analyticPosition(real_c(i) * dt));

      //            WALBERLA_LOG_DEVEL_VAR(/*p.getRotation().getQuaternion().getAxis() * */p.getRotation().getQuaternion().getAngle());
      //            WALBERLA_LOG_DEVEL_VAR(/*angSol.getQuaternion().getAxis() * */angSol.getQuaternion().getAngle());
      fout << p.getPosition().length() << " " <<
              analyticPosition(real_c(i) * dt).length() << " " <<
              p.getLinearVelocity().length() << " " <<
              analyticVelocity(real_c(i) * dt).length() << " " <<
              p.getRotation().getQuaternion().getAngle() << " " <<
              angSol.getQuaternion().getAngle() << " " <<
              p.getAngularVelocity().length() << " " <<
              analyticAngVel(real_c(i) * dt).length() << std::endl;
   }
   fout.close();

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   if (std::is_same<walberla::real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   return walberla::main();
}
