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
//! \file   ExplicitEulerWithShape.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleAccessor.h>

#include <mesa_pd/kernel/ExplicitEulerWithShape.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

class SingleParticleAccessor : public data::SingleParticleAccessor
{
public:
   const walberla::real_t& getInvMass(const size_t /*p_idx*/) const {return invMass_;}
   void setInvMass(const size_t /*p_idx*/, const walberla::real_t& v) { invMass_ = v;}
   const walberla::mesa_pd::Mat3& getInvInertiaBF(const size_t /*p_idx*/) const {return invInertiaBF_;}
   void setInvInertiaBF(const size_t /*p_idx*/, const walberla::mesa_pd::Mat3& v) { invInertiaBF_ = v;}

   walberla::real_t        invMass_;
   walberla::mesa_pd::Mat3 invInertiaBF_;
};

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   //init data structures
   SingleParticleAccessor accessor;

   //initialize particle
   const auto linVel = Vec3(1,2,3);
   const auto angVel = Vec3(1,2,3);

   const auto force  = Vec3(1,2,3);
   const auto torque = Vec3(1,2,3);

   accessor.setPosition(        0, Vec3(0,0,0));
   accessor.setRotation(        0, Rot3(Quat()));
   accessor.setLinearVelocity(  0, linVel);
   accessor.setAngularVelocity( 0, angVel);
   accessor.setForce(           0, force);
   accessor.setTorque(          0, torque);
   accessor.setInvMass(         0, real_t(1.23456));
   accessor.setInvInertiaBF(    0, Mat3(real_t(1.23456), real_t(0), real_t(0), real_t(0), real_t(1.23456), real_t(0), real_t(0), real_t(0), real_t(1.23456)));

   //init kernels
   const real_t dt = real_t(1);
   kernel::ExplicitEulerWithShape integrator( dt );

   integrator(0, accessor);

   const auto& R = accessor.getRotation(0).getMatrix();
   const auto wdot = R * accessor.getInvInertiaBF(0) * R.getTranspose() * torque;

   //check force
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getForce(0), Vec3(0));
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getTorque(0), Vec3(0));

   //check velocity
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getLinearVelocity(0), force * accessor.getInvMass(0) * dt + linVel);
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getAngularVelocity(0), wdot * dt + angVel);

   //check position
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getPosition(0), linVel * dt + force * accessor.getInvMass(0) * dt * dt);
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getRotation(0).getQuaternion(), Quat( (wdot * dt + angVel).getNormalized(), (wdot * dt + angVel).length() * dt ));

   accessor.setPosition(        0, Vec3(0,0,0));
   accessor.setRotation(        0, Rot3(Quat()));
   accessor.setLinearVelocity(  0, linVel);
   accessor.setAngularVelocity( 0, angVel);
   accessor.setForce(           0, force);
   accessor.setTorque(          0, torque);
   accessor.setInvMass(         0, real_t(1.23456));
   accessor.setInvInertiaBF(    0, Mat3(real_t(1.23456), real_t(0), real_t(0), real_t(0), real_t(1.23456), real_t(0), real_t(0), real_t(0), real_t(1.23456)));
   data::particle_flags::set( accessor.getFlagsRef(0), data::particle_flags::FIXED );

   integrator(0, accessor);

   //check force
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getForce(0), Vec3(0));
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getTorque(0), Vec3(0));

   //check velocity
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getLinearVelocity(0), linVel);
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getAngularVelocity(0), angVel);

   //check position
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getPosition(0), Vec3(0,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getRotation(0).getQuaternion(), Quat());

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
