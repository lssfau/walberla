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
//! \file InitializeHydrodynamicForceTorqueForAveragingKernel.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesa_pd/data/IAccessor.h"

namespace walberla {
namespace lbm_mesapd_coupling {

/*
 * Kernel that initializes the old hydrodynamic force/torque property of a particle with the currently set one.
 * This should be used when starting the simulation (from anew or from checkpoint) and after load balancing.
 * Only then, the following averaging kernel (AverageHydrodynamicForceTorqueKernel) applies the correct amount of force.
 */
class InitializeHydrodynamicForceTorqueForAveragingKernel
{

public:

   InitializeHydrodynamicForceTorqueForAveragingKernel( ) = default;

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac) const
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");
      
      ac.setOldHydrodynamicForce( idx, ac.getHydrodynamicForce(idx) );
      ac.setOldHydrodynamicTorque( idx, ac.getHydrodynamicTorque(idx) );
   }
};

} // namespace lbm_mesapd_coupling
} // namespace walberla
