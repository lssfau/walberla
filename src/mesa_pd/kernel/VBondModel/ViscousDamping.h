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

#pragma once

#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

#include <core/math/Constants.h>
#include <core/logging/Logging.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {
namespace VBondModel {

class ViscousDamping
{
public:
   ViscousDamping(real_t forceDampingFactor, real_t torqueDampingFactor)
   : forceDampingFactor_(forceDampingFactor)
   , torqueDampingFactor_(torqueDampingFactor)
   {}

   template<typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor &ac) const;

   auto getForceDampingFactor() const {return forceDampingFactor_;}
   auto getTorqueDampingFactor() const {return torqueDampingFactor_;}
private:
   const real_t forceDampingFactor_;
   const real_t torqueDampingFactor_;
};

template<typename Accessor>
inline void ViscousDamping::operator()(const size_t p_idx1,
                                       const size_t p_idx2,
                                       Accessor &ac) const
{
   Vec3 velDampingForce = (ac.getLinearVelocity(p_idx1) - ac.getLinearVelocity(p_idx2)) * forceDampingFactor_;
   addForceAtomic(p_idx1, ac, -velDampingForce);
   addForceAtomic(p_idx2, ac,  velDampingForce);

   Vec3 angDampingTorque = (ac.getAngularVelocity(p_idx1) - ac.getAngularVelocity(p_idx2)) * torqueDampingFactor_;
   addTorqueAtomic(p_idx1, ac, -angDampingTorque);
   addTorqueAtomic(p_idx2, ac,  angDampingTorque);

}

} //namespace VBondModel
} //namespace kernel
} //namespace mesa_pd
} //namespace walberla