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
//! \file AverageHydrodynamicForceTorqueKernel.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesa_pd/data/IAccessor.h"

namespace walberla {
namespace lbm_mesapd_coupling {

/*
 * Kernel that averages the force/torque over two time steps in a rolling average fashion.
 * This is said to damp oscillations of the interaction force/torque.
 * See Ladd - "Numerical simulations of particulate suspensions via a discretized Boltzmann equation. Part 1. Theoretical foundation", 1994, p.
 *
 * When reducing hyd. force/torque before, this should usually be carried out only on local particles. (recommended)
 * If not, it must be carried out on local and ghost particles.
 */
class AverageHydrodynamicForceTorqueKernel
{

public:

   AverageHydrodynamicForceTorqueKernel( ) = default;

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac) const
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      auto hydForce = ac.getHydrodynamicForce(idx);
      auto hydTorque = ac.getHydrodynamicTorque(idx);

      auto averageForce = real_t(0.5) * ( hydForce + ac.getOldHydrodynamicForce(idx) );
      auto averageTorque = real_t(0.5) * ( hydTorque + ac.getOldHydrodynamicTorque(idx) );

      // swap old values with new ones
      ac.setOldHydrodynamicForce( idx, hydForce );
      ac.setOldHydrodynamicTorque( idx, hydTorque );

      // set values to averaged ones
      ac.setHydrodynamicForce ( idx, averageForce );
      ac.setHydrodynamicTorque ( idx, averageTorque );
   }
};

} // namespace lbm_mesapd_coupling
} // namespace walberla
