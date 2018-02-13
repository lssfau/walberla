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
//! \file AddedMassForceCorrelations.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include <cmath>

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Correlation functions for the added mass force.
 *
 * To be compatible with the interface of AddedMassForceEvaluator, all functions use the following signature:
 * Vector3<real_t> ( const Vector3<real_t> & timeDerivativeFluidVel, const Vector3<real_t> & timeDerivativeBodyVel,
 *                   const real_t & bodyVolume, const real_t & fluidDensity )
 */

Vector3<real_t> addedMassForceFinn( const Vector3<real_t> & timeDerivativeFluidVel, const Vector3<real_t> & timeDerivativeBodyVel,
                                    const real_t & bodyVolume, const real_t & fluidDensity )
{
   // formula from Finn et al(2016)
   const real_t Coeffam = real_t(0.5);
   return bodyVolume * fluidDensity * Coeffam * ( timeDerivativeFluidVel - timeDerivativeBodyVel );
}

Vector3<real_t> noAddedMassForce( const Vector3<real_t> &, const Vector3<real_t> &, const real_t &, const real_t & )
{
   return Vector3<real_t>(real_t(0));
}

} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
