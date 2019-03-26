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
//! \file DragForceCorrelations.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Constants.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {


/*!\brief Various correlation functions for the drag force.
 *
 * These functions calculate the drag force for fluid-particle interactions based on different empirical correlations
 * from literature.
 * Always be aware that those empirical formulas were obtained by different setups and are thus not generally applicable!
 *
 * To be compatible with the interface of InteractionForceEvaluator, all functions use the following signature:
 * Vector3<real_t> ( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel, real_t solidVolumeFraction,
 *                   real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
 *
 */


// helper functions -> often used in more advanced drag correlations

// equation to calculate the drag coefficient on isolated spherical particle
// Schiller, L., Naumann, A., 1935. A drag coefficient correlation. Vdi Zeitung 77, 318-320.
real_t dragCoeffSchillerNaumann( real_t reynoldsNumber )
{
   WALBERLA_ASSERT_GREATER_EQUAL( reynoldsNumber, real_t(0) );

   return ( reynoldsNumber < real_t(1000) ) ? real_t(24) * ( real_t(1) + real_t(0.15) * std::pow(reynoldsNumber, real_t(0.687) ) ) / reynoldsNumber
                                            : real_t(0.44);
}

// Coefficient from Stokes' law for drag, only valid for Stokes regime (low Reynolds numbers)
// = 3 * M_PI * mu * D * fluidVolumeFraction
real_t dragCoeffStokes ( real_t fluidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity )
{
   return real_t(3) * math::M_PI * diameter * fluidDynamicViscosity * fluidVolumeFraction;
}

// threshold value for absolute relative velocity
// if it is below this value, a drag force of 0 is set, to avoid instabilities stemming from divisions by this small value
const real_t thresholdAbsoluteVelocityDifference = real_t(1e-10);


//////////////////////
//                  //
//   CORRELATIONS   //
//                  //
//////////////////////

// Stokes drag law
Vector3<real_t> dragForceStokes( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                                 real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t /*fluidDensity*/ )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, real_t(0) );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, real_t(1) );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(real_t(0));

   real_t fluidVolumeFraction = real_t(1) - solidVolumeFraction;

   return dragCoeffStokes( fluidVolumeFraction, diameter, fluidDynamicViscosity ) * velDiff;
}


// S. Ergun, Fluid flow through packed columns. Chemical Engineering Progress 48 (1952), 89-94.
// Y. C. Wen, Y.H. Yu, Mechanics of fluidization. Chemical Engineering Progress Symposium Series 62 (1966), 100-111.
// see also Beetstra, van der Hoef, Kuipers, "Drag Force of Intermediate Reynolds Number Flow Past Mono- and Bidisperse Arrays of Spheres" (2007)
Vector3<real_t>  dragForceErgunWenYu( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                                      real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, real_t(0) );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, real_t(1) );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(real_t(0));

   real_t fluidVolumeFraction = real_t(1) - solidVolumeFraction;

   if( fluidVolumeFraction < real_t(0.8) )
   {
      // Ergun relation
      real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;
      real_t fDrag = real_t(150) * solidVolumeFraction / ( real_t(18) * fluidVolumeFraction * fluidVolumeFraction ) +
                     real_t(1.75) / ( real_t(18) * fluidVolumeFraction * fluidVolumeFraction ) * reynoldsNumber;
      return fDrag * dragCoeffStokes( fluidVolumeFraction, diameter, fluidDynamicViscosity ) * velDiff;
   } else
   {
      // Wen & Yu correlation
      real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;
      real_t fDrag = dragCoeffSchillerNaumann( reynoldsNumber ) * reynoldsNumber / real_t(24) * std::pow( fluidVolumeFraction, real_t(-3.7) );
      return fDrag * dragCoeffStokes( fluidVolumeFraction, diameter, fluidDynamicViscosity ) * velDiff;
   }
}

// drag correlation proposed by Tang et al. - "A New Drag Correlation from Fully Resolved Simulations of Flow Past
// Monodisperse Static Arrays of Spheres", AiChE, 2014
Vector3<real_t> dragForceTang( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                               real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, real_t(0) );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, real_t(1) );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(real_t(0));

   real_t fluidVolumeFraction = real_t(1) - solidVolumeFraction;
   real_t fluidVolumeFractionP2 = fluidVolumeFraction * fluidVolumeFraction;
   real_t inv_fluidVolumeFractionP4 = real_t(1) / (fluidVolumeFractionP2 * fluidVolumeFractionP2);
   real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;

   // Eq.21 from the paper
   real_t fDrag = real_t(10) * solidVolumeFraction / fluidVolumeFractionP2 + fluidVolumeFractionP2 * ( real_t(1) + real_t(1.5) * std::sqrt(solidVolumeFraction) )
                + ( real_t(0.11) * solidVolumeFraction * ( real_t(1) + solidVolumeFraction ) - real_t(0.00456) * inv_fluidVolumeFractionP4
                + ( real_t(0.169) * fluidVolumeFraction + real_t(0.0644) * inv_fluidVolumeFractionP4 ) * std::pow( reynoldsNumber, -real_t(0.343) ) ) * reynoldsNumber;

   return fDrag * dragCoeffStokes( fluidVolumeFraction, diameter, fluidDynamicViscosity ) * velDiff;

}


// drag correlation based on findings from Felice (1994)
// used e.g. in Kafui et al (2002)
Vector3<real_t> dragForceFelice( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                                 real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, real_t(0) );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, real_t(1) );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(real_t(0));

   real_t fluidVolumeFraction = real_t(1) - solidVolumeFraction;

   real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;

   real_t temp1 = ( real_t(0.63) + real_t(4.8) / std::sqrt( reynoldsNumber ) );
   real_t dragCoeff = temp1 * temp1;

   real_t temp2 = real_t(1.5) - std::log10( reynoldsNumber );
   real_t chi = real_t(3.7) - std::pow( real_t(0.65), (- real_t(0.5) * temp2 * temp2 ) );

   return real_t(0.125) * dragCoeff * fluidDensity * math::M_PI * diameter * diameter * absVelDiff *
          std::pow( fluidVolumeFraction, real_t(2) - chi) * velDiff;

}

// drag correlation based on findings from Tenneti, Garg, Subramaniam (2011)
// used e.g. in Finn, Li, Apte - Particle based modelling and simulation of natural sand dynamics in the wave bottom boundary layer (2016)
// could be generalized also for non-spherical particles, see Finn et al (2016)
Vector3<real_t> dragForceTenneti( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                                  real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, real_t(0) );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, real_t(1) );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   const real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(real_t(0));

   const real_t fluidVolumeFraction = real_t(1) - solidVolumeFraction;

   const real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;

   const real_t fvfCubed = fluidVolumeFraction * fluidVolumeFraction * fluidVolumeFraction;
   const real_t A = real_t(5.81) * solidVolumeFraction / fvfCubed + real_t(0.48) * std::cbrt( solidVolumeFraction ) / ( fvfCubed * fluidVolumeFraction );

   const real_t svfCubed = solidVolumeFraction * solidVolumeFraction * solidVolumeFraction;
   const real_t B = svfCubed * reynoldsNumber * ( real_t(0.95) + real_t(0.61) * svfCubed / ( fluidVolumeFraction * fluidVolumeFraction ) );

   // version from Finn et al.
   const real_t CdRe0Sphere = real_t(1) + real_t(0.15) *  std::pow( reynoldsNumber, real_t(0.687) );

   const real_t CdRe = fluidVolumeFraction * ( CdRe0Sphere / fvfCubed + A + B );

   return real_t(3) * math::M_PI * diameter * fluidDynamicViscosity * fluidVolumeFraction * CdRe * velDiff;

}


Vector3<real_t> noDragForce( const Vector3<real_t> & /*fluidVel*/, const Vector3<real_t> & /*particleVel*/,
                             real_t /*solidVolumeFraction*/, real_t /*diameter*/, real_t /*fluidDynamicViscosity*/, real_t /*fluidDensity*/ )
{
   return Vector3<real_t>(real_t(0));
}

} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
