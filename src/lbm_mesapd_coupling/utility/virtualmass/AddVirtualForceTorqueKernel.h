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
//! \file AddVirtualForceTorqueKernel.h
//! \ingroup lbm_mesapd_coupling
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#pragma once

#include <utility>

#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/IAccessor.h"
#include "mesa_pd/common/ParticleFunctions.h"
#include "mesa_pd/data/ParticleStorage.h"

namespace walberla {
namespace lbm_mesapd_coupling {

/**
 * Kernel that sets a virtual force and torque on particles. It accesses a virtual mass and inertia, which have to be
 * set beforehand using InitializeVirtualMassKernel.
 * During the calculation of virtual force and torque linear and angular accelerations are used, for which acceleration
 * estimators need to be supplied.
 *
 * This kernel requires the following particle attributes:
 *  ps.addProperty("virtualMass",       "walberla::real_t",        defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")
 *  ps.addProperty("virtualInertiaBF",  "walberla::mesa_pd::Mat3", defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")
 *  ps.add_property("oldLinearAcceleration",  "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")
 *  ps.add_property("oldAngularAcceleration", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")
 */
class AddVirtualForceTorqueKernel
{
public:
   explicit AddVirtualForceTorqueKernel(shared_ptr<mesa_pd::data::ParticleStorage> ps) : ps_(std::move(ps)){};

   template <typename Accessor_T>
   void operator()(const size_t idx, Accessor_T& ac) {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, Accessor_T>::value, "please provide a valid accessor");

      WALBERLA_CHECK_FLOAT_UNEQUAL(ac.getVirtualMass(idx), real_t(0.),
                                   "No virtual mass set on body. Was the virtualMass kernel not called before?");

      //// force
      // obtain the acceleration estimation
      Vector3<real_t> approximatedCurrentLinearAcceleration = ac.getOldLinearAcceleration(idx);

      // prepare the next acceleration estimation for the following time step
      const Vector3<real_t> virtualForce = ac.getVirtualMass(idx) * approximatedCurrentLinearAcceleration;
      const Vector3<real_t> currentLinearAcceleration = (ac.getForce(idx) + virtualForce) /
                                                          (ac.getVirtualMass(idx) + ac.getMass(idx));
      ac.setOldLinearAcceleration(idx, currentLinearAcceleration);

      mesa_pd::addForceAtomic(idx, ac, virtualForce);

      //// torque
      // obtain the acceleration estimation
      Vector3<real_t> approximatedCurrentAngularAcceleration = ac.getOldAngularAcceleration(idx);

      // prepare the next acceleration estimation for the following time step
      const Vector3<real_t> virtualTorque = ac.getVirtualInertiaBF(idx) * approximatedCurrentAngularAcceleration;
      const Vector3<real_t> angularAcceleration = math::transformMatrixRART(ac.getRotation(idx).getMatrix(),
        ac.getInvInertiaBFIncludingVirtual(idx)) * (ac.getTorque(idx) + virtualTorque);
      ac.setOldAngularAcceleration(idx, angularAcceleration);

      mesa_pd::addTorqueAtomic(idx, ac, virtualTorque);
   }

private:
   const shared_ptr<mesa_pd::data::ParticleStorage> ps_;
};

}
}