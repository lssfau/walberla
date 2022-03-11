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

   const auto getInvInertiaBF(const size_t /*p_idx*/) const // dummy
   { return Mat3(real_t(0)); }

   const auto getInertiaBF(const size_t /*p_idx*/) const // dummy
   { return Mat3(real_t(0)); }

private:
   real_t invMass_;
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

      integrator(particle, osc);
   }

   return {maxPosDeviation, maxVelDeviation, maxEneDeviation};
}

int main(int argc, char **argv)
{
   mpi::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   if (std::is_same<real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   Oscillator osc;

   AccuracyResult res;
   osc.dt = 0.1_r;
   osc.update();
   res = walberla::checkIntegrator<walberla::ExplicitEuler>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 1.03068993874562120e+00_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 3.33225581358688350e-01_r);
   res = walberla::checkIntegrator<walberla::SemiImplicitEuler>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 2.92576360173544339e-02_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 1.45120298258364505e-03_r);
   res = walberla::checkIntegrator<walberla::VelocityVerlet>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 4.24116598691555435e-03_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 1.45059589844465098e-03_r);

   osc.dt = 0.2_r;
   osc.update();
   res = walberla::checkIntegrator<walberla::ExplicitEuler>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 2.76387000209972467e+00_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 8.88433661185269896e-01_r);
   res = walberla::checkIntegrator<walberla::SemiImplicitEuler>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 6.70464626869100577e-02_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 5.81009940233766925e-03_r);
   res = walberla::checkIntegrator<walberla::VelocityVerlet>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 1.69147419522671719e-02_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 5.78432113979018836e-03_r);

   osc.dt = 0.4_r;
   osc.update();
   res = walberla::checkIntegrator<walberla::ExplicitEuler>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 1.04680753378045406e+01_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 3.34470580215144420e+00_r);
   res = walberla::checkIntegrator<walberla::SemiImplicitEuler>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 1.68291142727994780e-01_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 2.33193930919295134e-02_r);
   res = walberla::checkIntegrator<walberla::VelocityVerlet>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxPosDeviation, 6.71909796751584687e-02_r);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxVelDeviation, 2.29906986019860066e-02_r);

   //check energy conservation
   osc.dt = 0.4_r;
   osc.periods = 1000;
   osc.update();
   //res = walberla::checkIntegrator<walberla::ExplicitEuler>(osc);
   //explicit euler is not symplectic
   res = walberla::checkIntegrator<walberla::SemiImplicitEuler>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxEneDeviation, 8.03571427989904774e-03_r);
   res = walberla::checkIntegrator<walberla::VelocityVerlet>(osc);
   WALBERLA_CHECK_FLOAT_EQUAL( res.maxEneDeviation, 4.99960610503419334e-04_r);

   return EXIT_SUCCESS;
}

} //namespace walberla


/*
 * Simulates a harmonic oscillator to test the accuracy of the integrators.
 */
int main(int argc, char **argv)
{
   return walberla::main(argc, argv);
}

