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
//! \file AddAccelerationOnParticlesKernel.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "mesa_pd/common/ParticleFunctions.h"
#include "mesa_pd/data/IAccessor.h"

#include <functional>

namespace walberla {
namespace lbm_mesapd_coupling {

/*
 * Kernel that adds a constant acceleration via F = m * a on particles
 *
 * Note that setting forces on ghost particles results in duplicated forces on these particles.
 *
 * If the force is constant for all particles, it is better to use AddForceOnParticlesKernel.
 *
 */
class AddAccelerationOnParticlesKernel
{

public:

   explicit AddAccelerationOnParticlesKernel( const Vector3<real_t> & acceleration )
   : acceleration_( acceleration )
     { }


   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac) const
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      auto force = acceleration_ / ac.getInvMass(idx);
      
      mesa_pd::addForceAtomic(idx, ac, force);
   }

   void setAcceleration( const Vector3<real_t> & newAcceleration )
   {
      acceleration_ = newAcceleration;
   }

   Vector3<real_t> getAcceleration() const
   {
      return acceleration_;
   }

private:
   Vector3<real_t> acceleration_;
};

} // namespace lbm_mesapd_coupling
} // namespace walberla
