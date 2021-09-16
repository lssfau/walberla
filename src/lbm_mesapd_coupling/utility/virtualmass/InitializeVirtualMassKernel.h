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
//! \file InitializeVirtualMassKernel.h
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

namespace walberla {
namespace lbm_mesapd_coupling {

/**
 * This kernel calculates virtual mass and inertia and sets them as attributes on a particle.
 *
 * It requires the following particle attributes:
 *  ps.addProperty("virtualMass",                  "walberla::real_t",        defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")
 *  ps.addProperty("invMassIncludingVirtual",      "walberla::real_t",        defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")
 *  ps.addProperty("virtualInertiaBF",             "walberla::mesa_pd::Mat3", defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")
 *  ps.addProperty("invInertiaBFIncludingVirtual", "walberla::mesa_pd::Mat3", defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")
 */
class InitializeVirtualMassKernel {
public:
   InitializeVirtualMassKernel() = default;

   template <typename Accessor>
   void operator()(size_t i, Accessor& ac, real_t C_v, real_t C_v_omega, real_t fluidDensity) const;
};


template <typename Accessor>
inline void InitializeVirtualMassKernel::operator()(const size_t i, Accessor& ac, const real_t C_v, const real_t C_v_omega,
        const real_t fluidDensity) const {
   static_assert(std::is_base_of<mesa_pd::data::IAccessor, Accessor>::value, "please provide a valid accessor");

   const real_t virtualMass = C_v * fluidDensity * ac.getVolume(i);
   ac.setVirtualMass(i, virtualMass);
   ac.setInvMassIncludingVirtual(i, real_t(1.) / (ac.getMass(i) + virtualMass));

   const real_t angularVirtualMass = C_v_omega * fluidDensity * ac.getVolume(i);
   const mesa_pd::Mat3 virtualInertiaBF = ac.getInertiaBF(i) * ac.getInvMass(i) * angularVirtualMass;
   const mesa_pd::Mat3 inertiaBFIncludingVirtual = ac.getInertiaBF(i) + virtualInertiaBF;
   ac.setVirtualInertiaBF(i, virtualInertiaBF);
   ac.setInvInertiaBFIncludingVirtual(i, inertiaBFIncludingVirtual.getInverse());
}

} //namespace lbm_mesapd_coupling
} //namespace walberla