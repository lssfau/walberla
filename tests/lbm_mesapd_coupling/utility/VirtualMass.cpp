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
//! \file VirtualMass.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Lukas Werner <lks.werner@fau.de>
//! \brief Evaluates the algorithms behind the virtual mass kernels by comparing with reference values.
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"

#include "lbm_mesapd_coupling/utility/virtualmass/AddVirtualForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/virtualmass/InitializeVirtualMassKernel.h"
#include "lbm_mesapd_coupling/utility/virtualmass/ParticleAccessorWithShapeVirtualMassWrapper.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"

namespace virtual_mass_test
{
using namespace walberla;

int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor = mesa_pd::data::ParticleAccessorWithShape;
   ParticleAccessor accessor(ps, ss);

   real_t sphereRadius = real_t(1);
   real_t sphereDensity = real_t(1);
   auto sphereShape = ss->create<mesa_pd::data::Sphere>(sphereRadius);
   ss->shapes[sphereShape]->updateMassAndInertia(sphereDensity);

   mesa_pd::data::Particle&& p1 = *ps->create();
   p1.setShapeID(sphereShape);
   auto idx = p1.getIdx();

   using VirtualMass_ParticleAccessor_T = lbm_mesapd_coupling::ParticleAccessorWithShapeVirtualMassWrapper<ParticleAccessor>;
   auto virtualMassAccessor = walberla::make_shared<VirtualMass_ParticleAccessor_T>(ps, ss);

   auto C_v = real_t(0.5);
   auto C_v_omega = real_t(0.5);
   auto fluidDensity = real_t(0.5);

   lbm_mesapd_coupling::InitializeVirtualMassKernel virtualMassInit;
   virtualMassInit(idx, accessor, C_v, C_v_omega, fluidDensity);

   const real_t sphereVirtualMass = C_v * fluidDensity * accessor.getVolume(idx);

   WALBERLA_CHECK_FLOAT_EQUAL(sphereVirtualMass, accessor.getVirtualMass(idx));
   WALBERLA_CHECK_FLOAT_EQUAL(real_t(1.) / (accessor.getMass(idx) + sphereVirtualMass),
                              virtualMassAccessor->getInvMass(idx));
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getMass(idx) + sphereVirtualMass,
                              virtualMassAccessor->getMass(idx));

   const real_t sphereAngularVirtualMass = C_v_omega * fluidDensity * accessor.getVolume(idx);
   const mesa_pd::Mat3 sphereVirtualInertiaBF = accessor.getInertiaBF(idx) * accessor.getInvMass(idx) * sphereAngularVirtualMass;

   WALBERLA_CHECK_FLOAT_EQUAL(sphereVirtualInertiaBF, accessor.getVirtualInertiaBF(idx));
   WALBERLA_CHECK_FLOAT_EQUAL((accessor.getInertiaBF(idx) + sphereVirtualInertiaBF).getInverse(),
                              virtualMassAccessor->getInvInertiaBF(idx));
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getInertiaBF(idx) + sphereVirtualInertiaBF,
                              virtualMassAccessor->getInertiaBF(idx));

   lbm_mesapd_coupling::AddVirtualForceTorqueKernel addVirtualForceTorque(ps);

   accessor.setForce(idx, Vector3(real_t(1)));
   auto oldForce = accessor.getForce(idx);
   addVirtualForceTorque(idx, accessor);
   WALBERLA_CHECK_FLOAT_EQUAL(oldForce, accessor.getForce(idx));
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getOldLinearAcceleration(idx),
                              Vector3<real_t>(real_t(1) / (accessor.getMass(idx) + accessor.getVirtualMass(idx))));

   Vector3 oldAcceleration(real_t(2));
   accessor.setOldLinearAcceleration(idx, oldAcceleration);
   accessor.setOldAngularAcceleration(idx, oldAcceleration);
   oldForce = Vector3(real_t(2));
   accessor.setForce(idx, oldForce);
   Vector3 oldTorque(real_t(2));
   accessor.setTorque(idx, oldTorque);
   addVirtualForceTorque(idx, accessor);
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getForce(idx), oldForce + sphereVirtualMass*oldAcceleration);
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getForce(idx),
                              accessor.getOldLinearAcceleration(idx)*(accessor.getMass(idx)+accessor.getVirtualMass(idx)));
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getTorque(idx), oldTorque + sphereVirtualInertiaBF*oldAcceleration);
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getTorque(idx),
                              (accessor.getInertiaBF(idx)+accessor.getVirtualInertiaBF(idx))*accessor.getOldAngularAcceleration(idx));

   return 0;
}

} // namespace virtual_mass_test

int main( int argc, char **argv ){
   virtual_mass_test::main(argc, argv);
}
