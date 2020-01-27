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
//! \file AddForceOnParticlesKernel.h
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
 * Kernel that adds a constant force on particles
 *
 * Note that setting forces on ghost particles results in duplicated forces on these particles.
 *
 */
class AddForceOnParticlesKernel
{

public:

   explicit AddForceOnParticlesKernel( const Vector3<real_t> & force )
   : force_( force )
     { }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac) const
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      mesa_pd::addForceAtomic(idx, ac, force_);
   }

   void setForce( const Vector3<real_t> & newForce )
   {
      force_ = newForce;
   }

   Vector3<real_t> getForce() const
   {
      return force_;
   }

private:
   Vector3<real_t> force_;
};

} // namespace lbm_mesapd_coupling
} // namespace walberla
