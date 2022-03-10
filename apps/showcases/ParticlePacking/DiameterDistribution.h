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
//! \file   DiameterDistribution.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include <random>

#include "Utility.h"

namespace walberla {
namespace mesa_pd {


class DiameterGenerator
{
public:
   virtual real_t get() = 0;
   virtual ~DiameterGenerator() = default;
};

class LogNormal: public DiameterGenerator
{
public:
   LogNormal(real_t mu, real_t var, uint_t seed = 42)
         : gen_(seed)
   {
      auto m = std::log(mu*mu/std::sqrt(var+mu*mu));
      auto s = std::sqrt(log(var/(mu*mu) + 1_r));
      dist_ = std::lognormal_distribution<> (m, s);
   }

   real_t get() override
   {
      return real_c(dist_(gen_));
   }
private:
   std::mt19937 gen_;
   std::lognormal_distribution<> dist_;
};

class Uniform: public DiameterGenerator
{
public:
   Uniform(real_t diameter)
         : diameter_(diameter)
   { }

   real_t get() override
   {
      return diameter_;
   }
private:
   real_t diameter_;
};

class DiscreteSieving: public DiameterGenerator
{
public:
   DiscreteSieving(std::vector<real_t> diameters, std::vector<real_t> massFractions, uint_t seed,
                   real_t normalParticleVolume, real_t totalParticleMass = 1_r, real_t particleDensity = 1_r)
         : diameters_(diameters), gen_(seed)
   {
      WALBERLA_CHECK_EQUAL(diameters.size(), massFractions.size(), "Number of entries in diameter and mass-fraction array has to be the same!");
      WALBERLA_CHECK_FLOAT_EQUAL(real_t(1), std::accumulate(massFractions.begin(), massFractions.end(), real_t(0)), "Sum of mass fractions has to be 1!");
      auto particleNumbers = transferMassFractionsToParticleNumbers(massFractions, diameters, normalParticleVolume, totalParticleMass, particleDensity);
      std::string outString = "Discrete Sieving: Expected particle numbers per diameter: | ";
      bool issueWarning = false;
      for(const auto & p : particleNumbers){
         issueWarning |= (p > 0_r && p < 1_r);
         outString += std::to_string(int(p)) + " | ";
      }
      auto totalParticles = std::accumulate(particleNumbers.begin(), particleNumbers.end(), real_t(0));
      outString += " -> total = " + std::to_string(int(totalParticles));
      WALBERLA_LOG_INFO_ON_ROOT(outString);
      if(issueWarning) WALBERLA_LOG_WARNING_ON_ROOT("Attention: The simulated particle mass is not enough to have at least a single particle of all size classes!");
      dist_ = std::discrete_distribution<uint_t>(particleNumbers.begin(), particleNumbers.end());
   }

   real_t get() override
   {
      return diameters_[dist_(gen_)];
   }
private:
   std::vector<real_t> diameters_;
   std::mt19937 gen_;
   std::discrete_distribution<uint_t> dist_;
};


class ContinuousSieving: public DiameterGenerator
{
public:
   ContinuousSieving(std::vector<real_t> sieveSizes, std::vector<real_t> massFractions, uint_t seed,
                     real_t normalParticleVolume, real_t totalParticleMass = 1_r, real_t particleDensity = 1_r)
         :  gen_(seed)
   {
      WALBERLA_CHECK_EQUAL(sieveSizes.size(), massFractions.size()+1, "Number of entries in sieves has to be one larger than the mass-fraction array!");
      WALBERLA_CHECK_FLOAT_EQUAL(real_t(1), std::accumulate(massFractions.begin(), massFractions.end(), real_t(0)), "Sum of mass fractions has to be 1!");

      auto meanDiameters = getMeanDiametersFromSieveSizes(sieveSizes);
      auto particleNumbers = transferMassFractionsToParticleNumbers(massFractions, meanDiameters, normalParticleVolume, totalParticleMass, particleDensity);
      std::string outString = "Continuous Sieving: Expected particle numbers per diameter: | ";
      bool issueWarning = false;
      for(const auto & p : particleNumbers){
         issueWarning |= (p > 0_r && p < 1_r);
         outString += std::to_string(int(p)) + " | ";
      }
      auto totalParticles = std::accumulate(particleNumbers.begin(), particleNumbers.end(), real_t(0));
      outString += " -> total = " + std::to_string(int(totalParticles));
      WALBERLA_LOG_INFO_ON_ROOT(outString);
      if(issueWarning) WALBERLA_LOG_WARNING_ON_ROOT("Attention: The simulated particle mass is not enough to have at least a single particle of all size classes!");
      std::vector<real_t> logSieveSizes(sieveSizes.size(),0_r); // =  phi-scale -> uniformly distributed
      for(uint_t i = 0; i < logSieveSizes.size(); ++i)
      {
         logSieveSizes[i] = - std::log2(sieveSizes[i]);
      }

      dist_ = std::piecewise_constant_distribution<real_t>(logSieveSizes.begin(), logSieveSizes.end(), particleNumbers.begin());
   }

   real_t get() override
   {
      return std::pow(2_r,-dist_(gen_)); //re-transfer from phi-scale to actual diameter
   }
private:
   std::mt19937 gen_;
   std::piecewise_constant_distribution<real_t> dist_;
};


} // namespace mesa_pd
} // namespace walberla
