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
//! \file LiftForceCorrelations.h
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

/*!\brief Correlation functions for the lift force.
 *
 * To be compatible with the interface of LiftForceEvaluator, all functions use the following signature:
 * Vector3<real_t> ( const Vector3<real_t> & fluidVel, const Vector3<real_t> & curlFluidVel, const Vector3<real_t> & particleVel,
 *                   real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
 */

// Saffman lift force
Vector3<real_t> liftForceSaffman ( const Vector3<real_t> & fluidVel, const Vector3<real_t> & curlFluidVel, const Vector3<real_t> & particleVel,
                                   real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   const real_t absCurlVel = curlFluidVel.length();
   if( absCurlVel < real_t(1e-10) ) return Vector3<real_t>(real_t(0));

   // Finn et al (2016) for spheres
   const real_t Cl = real_t(1.61) * std::sqrt( ( fluidDynamicViscosity * fluidDensity) / absCurlVel );
   return Cl * diameter * diameter * ( ( fluidVel - particleVel ) % curlFluidVel );

   // Sun, Xiao (2016)
   //const real_t Cl = real_t(1.6);
   //return Cl * fluidDensity * std::sqrt( fluidDynamicViscosity / fluidDensity ) * diameter * diameter * ( ( fluidVel - particleVel ) % curlFluidVel );

}

Vector3<real_t> noLiftForce ( const Vector3<real_t> &, const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t, real_t )
{
   return Vector3<real_t>(real_t(0));
}

} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
