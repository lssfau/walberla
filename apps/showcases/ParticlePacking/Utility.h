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
//! \file   Utility.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/mpi/MPITextFile.h"

#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/CylindricalBoundary.h"

#include <iterator>
#include <algorithm>
#include <functional>

namespace walberla {
namespace mesa_pd {

template< typename T>
std::vector<T> parseStringToVector(std::string inputStr)
{
   std::istringstream iss(inputStr);
   return std::vector<T>{std::istream_iterator<T>(iss),std::istream_iterator<T>()};
}

real_t radiusFromSphereVolume(real_t volume)
{
   return std::cbrt(real_t(3) / ( real_t(4) * math::pi) * volume);
}

real_t diameterFromSphereVolume(real_t volume)
{
   return real_t(2) * radiusFromSphereVolume(volume);
}

uint_t getIndexOfSecondSemiAxis(Vector3<real_t> semiAxes)
{
   std::set<uint_t> indices = {0,1,2};
   indices.erase(semiAxes.indexOfMin());
   indices.erase(semiAxes.indexOfMax());
   return *indices.begin();
}

void sortVector(Vector3<real_t> & v)
{
   std::sort(v.data(),v.data()+3);
}

real_t sizeFromSemiAxes(Vec3 semiAxes)
{
   sortVector(semiAxes);
   return 2_r * std::sqrt((semiAxes[0] * semiAxes[0] + semiAxes[1] * semiAxes[1]) / 2_r); // factor of 2 because those are SEMI axis
}

// assuming a mass-equivalent ellipsoid
Vec3 semiAxesFromInertiaTensor(const Matrix3<real_t> & inertiaTensor, real_t mass)
{
   Vec3 semiAxes;
   semiAxes[0] = std::sqrt( (- inertiaTensor(0,0) + inertiaTensor(1,1) + inertiaTensor(2,2)) / (2_r * mass / 5_r) );
   semiAxes[1] = std::sqrt( (- inertiaTensor(1,1) + inertiaTensor(2,2) + inertiaTensor(0,0)) / (2_r * mass / 5_r) );
   semiAxes[2] = std::sqrt( (- inertiaTensor(2,2) + inertiaTensor(0,0) + inertiaTensor(1,1)) / (2_r * mass / 5_r) );
   return semiAxes;
}


std::vector<real_t> getMeanDiametersFromSieveSizes(std::vector<real_t> sieveSizes)
{
   //since grain sizes are logarithmically distributed, it is practice to use the geometric mean
   std::vector<real_t> meanDiameters(sieveSizes.size()-1,0_r);
   for(uint_t i = 0; i < meanDiameters.size(); ++i)
   {
      meanDiameters[i] = std::sqrt(sieveSizes[i] * sieveSizes[i+1]);
   }
   return meanDiameters;
}

// if totalMass and density are given actual numbers, the resulting particle numbers are a good estimate for the actual numbers
// else, the resulting numbers are directly proportional to the actual ones, which is sufficient to define the distributions
std::vector<real_t> transferMassFractionsToParticleNumbersFromAvgVolumes(std::vector<real_t> massFractions, std::vector<real_t> avgVolumePerSizeFraction, real_t totalMass = real_t(1), real_t density = real_t(1) )
{
   WALBERLA_CHECK_EQUAL(avgVolumePerSizeFraction.size(), massFractions.size(), "Number of entries in volume and mass-fraction array has to be the same!");
   std::vector<real_t> particleNumbers(massFractions.size(), real_t(0));
   for(uint_t n = 0; n < massFractions.size(); ++n )
   {
      if(avgVolumePerSizeFraction[n] > 0_r) particleNumbers[n] = totalMass * massFractions[n] / (density * avgVolumePerSizeFraction[n]);
   }
   return particleNumbers;
}

// note: normalVolume is the volume of a typical particle with diameter = 1. For sphere: PI / 6
std::vector<real_t> transferMassFractionsToParticleNumbers(std::vector<real_t> massFractions, std::vector<real_t> diameters, real_t normalVolume = math::pi / real_t(6), real_t totalMass = real_t(1), real_t density = real_t(1) )
{
   WALBERLA_CHECK_EQUAL(diameters.size(), massFractions.size(), "Number of entries in diameter and mass-fraction array has to be the same!");
   std::vector<real_t> avgVolumePerSizeFraction(massFractions.size(), real_t(0));
   for(uint_t n = 0; n < massFractions.size(); ++n )
   {
      avgVolumePerSizeFraction[n] = normalVolume * diameters[n] * diameters[n] * diameters[n];
   }
   return transferMassFractionsToParticleNumbersFromAvgVolumes(massFractions, avgVolumePerSizeFraction, totalMass, density);
}

real_t computePercentileFromSieveDistribution(std::vector<real_t> diameters, std::vector<real_t> massFractions, real_t percentile)
{
   WALBERLA_CHECK(percentile > 0_r && percentile < 100_r, "Percentile is a value between 0 and 100");
   if(diameters[0] > diameters[1])
   {
      // reverse order to have it ascending
      std::reverse(diameters.begin(), diameters.end());
      std::reverse(massFractions.begin(), massFractions.end());
   }

   std::vector<real_t> cdf(massFractions.size(),0_r);
   std::partial_sum(massFractions.begin(), massFractions.end(), cdf.begin());

   for(uint_t i = 0; i < cdf.size()-1; ++i)
   {
      if(cdf[i] <= percentile/100_r && percentile/100_r <= cdf[i+1] )
      {
         real_t f_small = cdf[i];
         real_t f_large = cdf[i+1];
         if(f_small <= 0.000001_r && f_large>= 0.999999_r)
         {
            // special case of uniform distribution -> percentile is this value
            return diameters[i+1];
         }
         real_t phi_small = - std::log2(diameters[i]);
         real_t phi_large = - std::log2(diameters[i+1]);
         // logarithmic interpolation of diameter value
         real_t phi_percentile = phi_small + (phi_large - phi_small) / (f_large - f_small) * (percentile / 100_r - f_small);
         return std::pow(2_r,-phi_percentile);
      }
   }
   return diameters[0];
}


auto createPlane( std::shared_ptr<data::ParticleStorage> particleStorage,
                  const Vec3& pos,
                  const Vec3& normal)
{
   auto p = particleStorage->create(true);
   p->setPosition( pos );
   p->setBaseShape( std::make_shared<data::HalfSpace>( normal ) );
   p->getBaseShapeRef()->updateMassAndInertia(real_t(1));
   p->setOwner( walberla::mpi::MPIManager::instance()->rank() );
   p->setType( 0 );
   p->setInteractionRadius(std::numeric_limits<real_t>::infinity());
   data::particle_flags::set(p->getFlagsRef(), data::particle_flags::GLOBAL);
   data::particle_flags::set(p->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(p->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   return p;
}

auto createCylindricalBoundary( std::shared_ptr<data::ParticleStorage> particleStorage,
                                const Vec3& pos,
                                const Vec3& axis, real_t radius)
{
   auto p = particleStorage->create(true);
   p->setPosition( pos );
   p->setBaseShape( std::make_shared<data::CylindricalBoundary>( radius, axis ) );
   p->getBaseShapeRef()->updateMassAndInertia(real_t(1));
   p->setOwner( walberla::mpi::MPIManager::instance()->rank() );
   p->setType( 0 );
   p->setInteractionRadius(std::numeric_limits<real_t>::infinity());
   data::particle_flags::set(p->getFlagsRef(), data::particle_flags::GLOBAL);
   data::particle_flags::set(p->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(p->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   return p;
}


void writeParticleInformationToFile(const std::string& filename, const std::string& particleInfoStr, bool logToProcessLocalFiles)
{

   if(logToProcessLocalFiles)
   {
      std::ofstream file;
      file.open( filename.c_str());
      file << particleInfoStr;
      file.close();
   }else{
      walberla::mpi::writeMPITextFile( filename, particleInfoStr );
   }

}



} // namespace mesa_pd
} // namespace walberla
