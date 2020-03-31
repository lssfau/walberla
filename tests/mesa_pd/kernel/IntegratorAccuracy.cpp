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
//! \file   IntegratorAccuracy.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"

#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/kernel/ExplicitEuler.h"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/all.h"

#include <iostream>

namespace integrator_accuracy
{

using namespace walberla;

using namespace walberla::mesa_pd;


Vec3 getForce(const Vec3& pos, real_t k)
{
   return -k * pos;
}

real_t analyticalTrajectory(real_t amplitude, real_t timeStep, real_t omega, real_t phase)
{
   return amplitude * std::cos(omega * timeStep + phase);
}

real_t analyticalVelocity(real_t amplitude, real_t timeStep, real_t omega, real_t phase)
{
   return -amplitude * omega * std::sin(omega * timeStep + phase);
}


/*
 * Simulates a harmonic oscillator to test the accuracy of the integrators.
 * The error of the maximum position and the maximum velocity is compared against the analytical values.
 * Via command line arguments, the simulation can be adapted.
 * Currently tested integrators:
 *  - explicit Euler (default)
 *  - velocity verlet (-useVV)
 */
int main( int argc, char ** argv )
{
   mpi::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   real_t amplitude       = real_t(1.5);
   real_t k               = real_t(0.1);
   real_t mass            = real_t(0.9);
   real_t dt              = real_t(0.2);
   bool useVelocityVerlet = false;
   real_t phaseFraction   = real_t(0);
   real_t periods         = real_t(1.1);

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--useVV" )         == 0 ) { useVelocityVerlet = true; continue; }
      if( std::strcmp( argv[i], "--dt" )            == 0 ) { dt = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--amplitude" )     == 0 ) { amplitude = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--k" )             == 0 ) { k = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--mass" )          == 0 ) { mass = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--phaseFraction" ) == 0 ) { phaseFraction = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--periods" )       == 0 ) { periods = real_c(std::atof( argv[++i] )); continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }


   real_t phase = phaseFraction * math::pi;
   real_t omega = std::sqrt(k / mass);
   real_t durationOnePeriod = real_t(2) * math::pi / omega;
   uint_t timeSteps = uint_c(periods * durationOnePeriod / dt);

   WALBERLA_LOG_INFO("omega = " << omega << ", T = " << real_t(2) * math::pi / omega << ", time steps = " << timeSteps << ", phase = " << phase << ", periods = " << periods);

   //initialize particle
   const auto pos = Vec3(0,0,analyticalTrajectory(amplitude, real_t(0), omega, phase));
   const auto linVel = Vec3(0,0,analyticalVelocity(amplitude, real_t(0), omega, phase));

   WALBERLA_LOG_INFO("Initial pos = " << pos << ", vel = " << linVel);

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(1);
   auto ss = std::make_shared<data::ShapeStorage>();
   auto  smallSphere = ss->create<data::Sphere>( real_t(1) );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));
   data::ParticleAccessorWithShape accessor(ps, ss);
   
   data::Particle&& p = *ps->create();
   p.setPosition(pos);
   p.setLinearVelocity(linVel);
   p.setForce(getForce(pos, k));
   p.setOldForce(getForce(pos, k));
   p.setInvMass(real_t(1) / mass);
   p.setShapeID(smallSphere);

   // velocity verlet
   kernel::VelocityVerletPreForceUpdate  preForce( dt );
   kernel::VelocityVerletPostForceUpdate postForce( dt );

   // explicit euler
   kernel::ExplicitEuler explEuler( dt );

   real_t maxVel = 0.;
   real_t maxRefVel = 0.;
   real_t maxPos = 0.;
   real_t maxRefPos = 0.;
   
   for (auto i = uint_t(1); i <= timeSteps; ++i)
   {

      if( useVelocityVerlet ) ps->forEachParticle(false, kernel::SelectAll(), accessor, preForce, accessor);
      p.setForce( getForce(p.getPosition(), k) );
      auto force = p.getForce();

      if( useVelocityVerlet ) ps->forEachParticle(false, kernel::SelectAll(), accessor, postForce, accessor);
      else  ps->forEachParticle(false, kernel::SelectAll(), accessor, explEuler, accessor);

      real_t refPos = analyticalTrajectory(amplitude, real_c(i) * dt, omega, phase);
      real_t refVel = analyticalVelocity(amplitude, real_c(i) * dt, omega, phase);

      WALBERLA_LOG_INFO(i << ": pos = " << p.getPosition()[2] << " " << refPos
                          << " || vel = " << p.getLinearVelocity()[2] << " " << refVel
                          << " || force = " << force[2] );

      maxPos = std::max(maxPos, std::abs(p.getPosition()[2]));
      maxRefPos = std::max(maxRefPos, std::abs(refPos));

      maxVel = std::max(maxVel, std::abs(p.getLinearVelocity()[2]));
      maxRefVel = std::max(maxRefVel, std::abs(refVel));
   }

   real_t relativePositionError = ( maxPos - maxRefPos ) / maxRefPos;
   WALBERLA_LOG_INFO("error in position = " << relativePositionError * 100. << "%");

   real_t relativeVelocityError = ( maxVel - maxRefVel ) / maxRefVel;
   WALBERLA_LOG_INFO("error in velocity = " << relativeVelocityError * 100. << "%");

   if( useVelocityVerlet )
   {
      WALBERLA_CHECK_LESS(relativePositionError, real_t(0.01), "Error in position too large!");
      WALBERLA_CHECK_LESS(relativeVelocityError, real_t(0.01), "Error in velocity too large!");
   }
   else
   {
      WALBERLA_CHECK_LESS(relativePositionError, real_t(0.11), "Error in position too large!");
      WALBERLA_CHECK_LESS(relativeVelocityError, real_t(0.10), "Error in velocity too large!");
   }


   return EXIT_SUCCESS;
}

} //namespace integrator_accuracy

int main( int argc, char ** argv )
{
   return integrator_accuracy::main(argc, argv);
}

