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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "mesa_pd/data/ParticleAccessor.h"

#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/SemiImplicitEuler.h"
#include "mesa_pd/kernel/VelocityVerlet.h"

#include "core/Environment.h"
#include "core/math/all.h"

#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

class SingleParticleAccessorWithShape : public data::SingleParticleAccessor
{
public:
   void setInvMass(const size_t /*p_idx*/, const real_t &val)
   { invMass_ = val; }

   const auto &getInvMass(const size_t /*p_idx*/) const
   { return invMass_; }

   void setInvInertiaBF(const size_t /*p_idx*/, const Mat3 &val)
   { invInertiaBF_ = val; }

   const auto &getInvInertiaBF(const size_t /*p_idx*/) const
   { return invInertiaBF_; }

private:
   real_t invMass_;
   Mat3 invInertiaBF_;
};

struct Oscillator
{
   real_t amplitude = 1.5_r;
   real_t k = 0.1_r;
   real_t damping = 0_r;
   real_t mass = 0.9_r;
   real_t dt = 0.2_r;
   real_t phaseFraction = 0_r;
   real_t periods = 10_r;
   real_t phase = phaseFraction * math::pi;
   real_t dampingRatio = damping / (2_r * std::sqrt(mass * k));
   real_t omega = std::sqrt(k / mass) * std::sqrt(1_r - dampingRatio * dampingRatio);
   real_t decay = std::sqrt(k / mass) * dampingRatio;
   real_t durationOnePeriod = 2_r * math::pi / omega;
   uint_t timeSteps = uint_c(periods * durationOnePeriod / dt);

   Oscillator(int argc, char **argv)
   {
      for (int i = 1; i < argc; ++i)
      {
         if (std::strcmp(argv[i], "--dt") == 0)
         {
            dt = real_c(std::atof(argv[++i]));
            continue;
         }
         if (std::strcmp(argv[i], "--amplitude") == 0)
         {
            amplitude = real_c(std::atof(argv[++i]));
            continue;
         }
         if (std::strcmp(argv[i], "--k") == 0)
         {
            k = real_c(std::atof(argv[++i]));
            continue;
         }
         if (std::strcmp(argv[i], "--damping") == 0)
         {
            damping = real_c(std::atof(argv[++i]));
            continue;
         }
         if (std::strcmp(argv[i], "--mass") == 0)
         {
            mass = real_c(std::atof(argv[++i]));
            continue;
         }
         if (std::strcmp(argv[i], "--phaseFraction") == 0)
         {
            phaseFraction = real_c(std::atof(argv[++i]));
            continue;
         }
         if (std::strcmp(argv[i], "--periods") == 0)
         {
            periods = real_c(std::atof(argv[++i]));
            continue;
         }
         WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
      }

      update();
   }

   void update()
   {
      phase = phaseFraction * math::pi;
      dampingRatio = damping / (2_r * std::sqrt(mass * k));
      omega = std::sqrt(k / mass) * std::sqrt(1_r - dampingRatio * dampingRatio);
      decay = std::sqrt(k / mass) * dampingRatio;
      durationOnePeriod = 2_r * math::pi / omega;
      timeSteps = uint_c(periods * durationOnePeriod / dt);
   }

   Vec3 getForce(const Vec3 &pos, const Vec3 &vel) const
   {
      return -k * pos - damping * vel;
   }

   real_t getEnergy(const real_t &pos, const real_t &vel) const
   {
      return 0.5_r * mass * vel * vel + 0.5_r * k * pos * pos;
   }

   real_t analyticalPos(real_t t) const
   {
      return amplitude * std::exp(-decay * t) * std::cos(omega * t + phase);
   }

   real_t analyticalVel(real_t t) const
   {
      return -decay * amplitude * std::exp(-decay * t) * std::cos(omega * t + phase)
             -amplitude * std::exp(-decay * t) * omega * std::sin(omega * t + phase);
   }
};

struct ExplicitEuler
{
   ExplicitEuler(real_t dt) : integrator(dt) {}
   void operator()(SingleParticleAccessorWithShape& particle,
                   const Oscillator& osc)
   {
      particle.setForce(0, osc.getForce(particle.getPosition(0),
                                           particle.getLinearVelocity(0)));
      integrator(0, particle);
   }
   kernel::ExplicitEuler integrator;
};

struct SemiImplicitEuler
{
   SemiImplicitEuler(real_t dt) : integrator(dt) {}
   void operator()(SingleParticleAccessorWithShape& particle,
                   const Oscillator& osc)
   {
      particle.setForce(0, osc.getForce(particle.getPosition(0),
                                           particle.getLinearVelocity(0)));
      integrator(0, particle);
   }
   kernel::SemiImplicitEuler integrator;
};

struct VelocityVerlet
{
   VelocityVerlet(real_t dt) : preVV(dt), postVV(dt) {}
   void operator()(SingleParticleAccessorWithShape& particle,
                   const Oscillator& osc)
   {
      preVV(0, particle);
      particle.setForce(0, osc.getForce(particle.getPosition(0),
                                           particle.getLinearVelocity(0)));
      postVV(0, particle);
   }
   kernel::VelocityVerletPreForceUpdate preVV;
   kernel::VelocityVerletPostForceUpdate postVV;
};

struct AccuracyResult
{
   real_t maxPosDeviation;
   real_t maxVelDeviation;
   real_t maxEneDeviation;
};

template <typename Integrator>
AccuracyResult checkIntegrator(const Oscillator& osc)
{
   //init data structures
   SingleParticleAccessorWithShape particle;

   //first dummy argument is needed to fulfill accessor interface
   particle.setPosition(0, Vec3(0, 0, osc.analyticalPos(0_r)));
   particle.setLinearVelocity(0, Vec3(0, 0, osc.analyticalVel(0_r)));
   particle.setInvMass(0, 1_r / osc.mass);
   particle.setForce(0, osc.getForce(Vec3(0, 0, osc.analyticalPos(real_t(0))),
                                        Vec3(0, 0, osc.analyticalVel(real_t(0)))));
   particle.setOldForce(0, osc.getForce(Vec3(0, 0, osc.analyticalPos(-osc.dt)),
                                           Vec3(0, 0, osc.analyticalVel(-osc.dt))));

   // explicit euler
   Integrator integrator(osc.dt);

   real_t maxPosDeviation = 0_r;
   real_t maxVelDeviation = 0_r;
   real_t maxEneDeviation = 0_r;

   for (auto i = uint_t(0); i <= osc.timeSteps; ++i)
   {
      real_t refPos = osc.analyticalPos(real_c(i) * osc.dt);
      real_t refVel = osc.analyticalVel(real_c(i) * osc.dt);
      real_t refEne = osc.getEnergy(refPos, refVel);

      maxPosDeviation = std::max(maxPosDeviation, std::abs(particle.getPosition(0)[2] - refPos));
      maxVelDeviation = std::max(maxVelDeviation, std::abs(particle.getLinearVelocity(0)[2] - refVel));
      maxEneDeviation = std::max(maxEneDeviation, std::abs(osc.getEnergy(particle.getPosition(0)[2], particle.getLinearVelocity(0)[2]) - refEne));

      std::cout << real_t(i) * osc.dt << " "
                << refPos << " "
                << refVel << " "
                << refEne << " "
                << particle.getPosition(0)[2] << " "
                << particle.getLinearVelocity(0)[2] << " "
                << osc.getEnergy(particle.getPosition(0)[2], particle.getLinearVelocity(0)[2]) << " "
                << maxPosDeviation << " "
                << maxVelDeviation << std::endl;

      integrator(particle, osc);
   }

   return {maxPosDeviation, maxVelDeviation, maxEneDeviation};
}

int main(int argc, char **argv)
{
   mpi::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   Oscillator osc(argc, argv);

   AccuracyResult res;
//   res = walberla::checkIntegrator<walberla::ExplicitEuler>(osc);
//   res = walberla::checkIntegrator<walberla::SemiImplicitEuler>(osc);
   res = walberla::checkIntegrator<walberla::VelocityVerlet>(osc);

   return EXIT_SUCCESS;
}

} //namespace walberla


/*
 * Simulates a harmonic oscillator to test the accuracy of the integrators.
 * Playground for integrator analysis. The corresponding unit test is located at
 * tests/mesa_pd/kernel/IntegratorAccuracy.cpp
 */
int main(int argc, char **argv)
{
   return walberla::main(argc, argv);
}

