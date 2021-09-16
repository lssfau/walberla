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
//! \file ParticleAccessorWithShapeVirtualMassWrapper.h
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>

namespace walberla {
namespace lbm_mesapd_coupling {

/**
 * !\brief Wraps a ParticleAccessor to change its getInvMass and getInvInertiaBF methods to include virtual mass.
 * @tparam T A ParticleAccessor providing functions getInvMass and getInvInertiaBF, and their respective getInvMassIncludingVirtual and getInvInertiaBFIncludingVirtual.
 */
template <typename T>
class ParticleAccessorWithShapeVirtualMassWrapper : public T {
public:
  ParticleAccessorWithShapeVirtualMassWrapper(std::shared_ptr<mesa_pd::data::ParticleStorage>& ps, std::shared_ptr<mesa_pd::data::ShapeStorage>& ss)
      : T(ps, ss) {
   }

   auto getInvMass(const size_t p_idx) const {
      return T::getInvMassIncludingVirtual(p_idx);
   }

   const auto getMass(const size_t p_idx) const {
      return T::getMass(p_idx) + T::getVirtualMass(p_idx);
   }

   auto getInvInertiaBF(const size_t p_idx) const {
      return T::getInvInertiaBFIncludingVirtual(p_idx);
   }

   const auto getInertiaBF(const size_t p_idx) const {
      return T::getInertiaBF(p_idx) + T::getVirtualInertiaBF(p_idx);
   }
};

}
}