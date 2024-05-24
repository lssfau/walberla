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
//! \file   ParticleUtility.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/mpi/Broadcast.h"
#include "core/mpi/MPITextFile.h"
#include "core/mpi/Reduce.h"

#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Sphere.h"

namespace walberla
{
namespace piping
{

// Some functions in this file (as the one below) are based on showcases/Antidunes/Utility.cpp

void writeSphereInformationToFile(const std::string& filename, walberla::mesa_pd::data::ParticleStorage& ps,
                                  const Vector3< real_t >& domainSize, const int precision = 12)
{
   std::ostringstream ossData;
   ossData << std::setprecision(precision);

   WALBERLA_ROOT_SECTION() { ossData << domainSize[0] << " " << domainSize[1] << " " << domainSize[2] << "\n"; }

   for (auto pIt : ps)
   {
      using namespace walberla::mesa_pd::data;
      if (pIt->getBaseShape()->getShapeType() != Sphere::SHAPE_TYPE) continue;
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(pIt->getFlags(), GHOST)) continue;
      auto sp = static_cast< Sphere* >(pIt->getBaseShape().get());

      auto position = pIt->getPosition();

      ossData << position[0] << " " << position[1] << " " << position[2] << " " << sp->getRadius() << '\n';
   }

   walberla::mpi::writeMPITextFile(filename, ossData.str());
}

bool sphereBoxOverlap(const mesa_pd::Vec3& spherePosition, const real_t sphereRadius, const mesa_pd::Vec3& boxPosition,
                      const mesa_pd::Vec3& boxEdgeLength)
{
   if ((spherePosition[0] + sphereRadius < boxPosition[0] - boxEdgeLength[0] / real_t(2)) ||
       (spherePosition[1] + sphereRadius < boxPosition[1] - boxEdgeLength[1] / real_t(2)) ||
       (spherePosition[2] + sphereRadius < boxPosition[2] - boxEdgeLength[2] / real_t(2)) ||
       (spherePosition[0] - sphereRadius > boxPosition[0] + boxEdgeLength[0] / real_t(2)) ||
       (spherePosition[1] - sphereRadius > boxPosition[1] + boxEdgeLength[1] / real_t(2)) ||
       (spherePosition[2] - sphereRadius > boxPosition[2] + boxEdgeLength[2] / real_t(2)))
   {
      return false;
   }
   return true;
}

void initSpheresFromFile(const std::string& fileName, walberla::mesa_pd::data::ParticleStorage& ps,
                         const walberla::mesa_pd::domain::IDomain& domain, walberla::real_t particleDensity,
                         math::AABB& simulationDomain, const Vector3< uint_t >& domainSize,
                         const mesa_pd::Vec3& boxPosition, const mesa_pd::Vec3& boxEdgeLength, real_t& maxDiameter)
{
   using namespace walberla::mesa_pd::data;

   auto rank = walberla::mpi::MPIManager::instance()->rank();

   std::string textFile;

   WALBERLA_ROOT_SECTION()
   {
      std::ifstream t(fileName.c_str());
      if (!t) { WALBERLA_ABORT("Invalid input file " << fileName << "\n"); }
      std::stringstream buffer;
      buffer << t.rdbuf();
      textFile = buffer.str();
   }

   walberla::mpi::broadcastObject(textFile);

   std::istringstream fileIss(textFile);
   std::string line;

   // first line contains generation domain sizes
   std::getline(fileIss, line);
   Vector3< real_t > generationDomainSize_SI(0_r);
   std::istringstream firstLine(line);
   firstLine >> generationDomainSize_SI[0] >> generationDomainSize_SI[1] >> generationDomainSize_SI[2];
   real_t scalingFactor = real_t(domainSize[0]) / generationDomainSize_SI[0];
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(generationDomainSize_SI)
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(scalingFactor)

   real_t minParticleDiameter = std::numeric_limits< real_t >::max();
   real_t maxParticleDiameter = real_t(0);

   while (std::getline(fileIss, line))
   {
      std::istringstream iss(line);

      ParticleStorage::position_type position;
      real_t radius;
      iss >> position[0] >> position[1] >> position[2] >> radius;
      position *= scalingFactor;
      radius *= scalingFactor;

      WALBERLA_CHECK(simulationDomain.contains(position),
                     "Particle read from file is not contained in simulation domain");

      if (!domain.isContainedInProcessSubdomain(uint_c(rank), position)) continue;
      if (sphereBoxOverlap(position, radius, boxPosition, boxEdgeLength)) continue;

      auto pIt = ps.create();
      pIt->setPosition(position);
      pIt->getBaseShapeRef() = std::make_shared< data::Sphere >(radius);
      pIt->getBaseShapeRef()->updateMassAndInertia(particleDensity);
      pIt->setInteractionRadius(radius);
      pIt->setOwner(rank);
      pIt->setType(0);

      minParticleDiameter = std::min(real_t(2) * radius, minParticleDiameter);
      maxParticleDiameter = std::max(real_t(2) * radius, maxParticleDiameter);

      WALBERLA_CHECK_EQUAL(iss.tellg(), -1);
   }

   WALBERLA_MPI_SECTION() { walberla::mpi::allReduceInplace(minParticleDiameter, walberla::mpi::MIN); }
   WALBERLA_MPI_SECTION() { walberla::mpi::allReduceInplace(maxParticleDiameter, walberla::mpi::MAX); }
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(minParticleDiameter)
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(maxParticleDiameter)
   // Maximum particle diameter is used for the size of the linked cells
   maxDiameter = maxParticleDiameter;
}

template< typename ParticleAccessor_T >
void getParticleVelocities(const ParticleAccessor_T& ac, uint_t& numParticles, real_t& maxVelocity,
                           real_t& averageVelocity)
{
   maxVelocity     = real_t(0);
   averageVelocity = real_t(0);
   numParticles    = uint_t(0);

   for (uint_t i = 0; i < ac.size(); ++i)
   {
      if (isSet(ac.getFlags(i), walberla::mesa_pd::data::particle_flags::GHOST)) continue;
      if (isSet(ac.getFlags(i), walberla::mesa_pd::data::particle_flags::GLOBAL)) continue;

      ++numParticles;
      real_t velMagnitude = ac.getLinearVelocity(i).length();
      maxVelocity         = std::max(maxVelocity, velMagnitude);
      averageVelocity += velMagnitude;
   }

   WALBERLA_MPI_SECTION()
   {
      walberla::mpi::allReduceInplace(maxVelocity, walberla::mpi::MAX);
      walberla::mpi::allReduceInplace(averageVelocity, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(numParticles, walberla::mpi::SUM);
   }

   averageVelocity /= real_t(numParticles);
}

auto createPlane(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& pos, const mesa_pd::Vec3& normal)
{
   auto p0 = ps.create(true);
   p0->setPosition(pos);
   p0->setBaseShape(std::make_shared< mesa_pd::data::HalfSpace >(normal));
   // Mass is set to infinity internally for HalfSpace (independent of the density that is set here)
   p0->getBaseShapeRef()->updateMassAndInertia(real_t(1));
   p0->setOwner(walberla::mpi::MPIManager::instance()->rank());
   p0->setType(1);
   p0->setInteractionRadius(std::numeric_limits< real_t >::infinity());
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::GLOBAL);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::NON_COMMUNICATING);
   return p0;
}

auto createBox(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& pos, const mesa_pd::Vec3& edgeLength,
               const bool movingBucket)
{
   auto p0 = ps.create(true);
   p0->setPosition(pos);
   p0->setBaseShape(std::make_shared< mesa_pd::data::Box >(edgeLength));
   if (movingBucket)
   {
      // TODO: replace the density of 2.0
      p0->getBaseShapeRef()->updateMassAndInertia(real_t(2.0));
   }
   else
   {
      // If the bucket is fixed, its collision behaviour should be the same as for the bounding planes
      p0->getBaseShapeRef()->updateMassAndInertia(std::numeric_limits< real_t >::infinity());
   }
   p0->setOwner(walberla::mpi::MPIManager::instance()->rank());
   p0->setType(1);
   p0->setInteractionRadius(std::numeric_limits< real_t >::infinity());
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::GLOBAL);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   return p0->getUid();
}

template< typename ParticleAccessor_T, typename Sync_T, typename CollisionResponse_T >
void settleParticles(const uint_t numTimeSteps, const shared_ptr< ParticleAccessor_T >& accessor,
                     const shared_ptr< mesa_pd::data::ParticleStorage >& ps,
                     const walberla::mesa_pd::domain::IDomain& domain, mesa_pd::data::LinkedCells& linkedCells,
                     Sync_T& syncNextNeighborFunc, CollisionResponse_T& collisionResponse,
                     const real_t& particleDensityRatio, const real_t& gravitationalAcceleration, const bool& useOpenMP)
{
   // Increase the settling speed
   const real_t timeStepSizeParticles = real_t(10);
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce(timeStepSizeParticles);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(timeStepSizeParticles);
   mesa_pd::mpi::ReduceProperty reduceProperty;
   mesa_pd::mpi::ReduceContactHistory reduceAndSwapContactHistory;
   mesa_pd::kernel::InsertParticleIntoLinkedCells ipilc;

   WALBERLA_LOG_INFO_ON_ROOT("Starting initial particle settling...")

   for (uint_t t = uint_t(0); t < numTimeSteps; ++t)
   {
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
      syncNextNeighborFunc(*ps, domain);

      linkedCells.clear();
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, ipilc, *accessor, linkedCells);

      // collision
      linkedCells.forEachParticlePairHalf(
         useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
         [&collisionResponse, &domain, &timeStepSizeParticles](const size_t idx1, const size_t idx2, auto& ac) {
            mesa_pd::collision_detection::AnalyticContactDetection acd;
            mesa_pd::kernel::DoubleCast double_cast;
            mesa_pd::mpi::ContactFilter contact_filter;
            if (double_cast(idx1, idx2, ac, acd, ac))
            {
               if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
               {
                  collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(),
                                    acd.getPenetrationDepth(), timeStepSizeParticles);
               }
            }
         },
         *accessor);
      reduceAndSwapContactHistory(*ps);

      // gravity - buoyancy
      ps->forEachParticle(
         useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
         [particleDensityRatio, gravitationalAcceleration](const size_t idx, auto& ac) {
            mesa_pd::addForceAtomic(
               idx, ac,
               Vector3< real_t >(real_t(0), real_t(0),
                                 -(particleDensityRatio - real_c(1)) * ac.getVolume(idx) * gravitationalAcceleration));
         },
         *accessor);

      reduceProperty.operator()< mesa_pd::ForceTorqueNotification >(*ps);

      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
      syncNextNeighborFunc(*ps, domain);

      if (t % (numTimeSteps / uint_t(10)) == 0)
      {
         real_t maxVelocity;
         real_t averageVelocity;
         uint_t numAveragedParticles;

         getParticleVelocities(*accessor, numAveragedParticles, maxVelocity, averageVelocity);
         WALBERLA_LOG_INFO_ON_ROOT("Timestep "
                                   << t << " / " << numTimeSteps << ", average velocity = " << averageVelocity
                                   << ", max velocity = " << maxVelocity << ", #particles = " << numAveragedParticles);
      }
   }

   // Velocities should be 0 after settling such that the simulation starts from rest
   ps->forEachParticle(
      useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
      [](const size_t idx, auto& ac) {
         ac.setLinearVelocity(idx, Vector3(real_t(0)));
         ac.setAngularVelocity(idx, Vector3(real_t(0)));
      },
      *accessor);
   syncNextNeighborFunc(*ps, domain);
}

} // namespace piping
} // namespace walberla
