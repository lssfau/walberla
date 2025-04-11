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
//! \file PSMWrapperSweepsCPU.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

#include "lbm/sweeps/StreamPull.h"
#include "lbm/sweeps/SweepBase.h"

#include "lbm_mesapd_coupling/DataTypesCodegen.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/common/ParticleFunctions.h"

#include "timeloop/SweepTimeloop.h"

#include <cassert>

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{
namespace gpu
{

template< typename ParticleAccessor_T, typename ParticleSelector_T, int Weighting_T >
class SetParticleTemperaturesSweep
{
 public:
   SetParticleTemperaturesSweep(const shared_ptr< StructuredBlockStorage >& bs,
                                const shared_ptr< ParticleAccessor_T >& ac,
                                const ParticleSelector_T& mappingParticleSelector,
                                ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA)
      : bs_(bs), ac_(ac), mappingParticleSelector_(mappingParticleSelector),
        particleAndVolumeFractionSoA_(particleAndVolumeFractionSoA)
   {}
   void operator()(IBlock* block)
   {
      // Check that uids of the particles have not changed since the last mapping to avoid incorrect indices
      std::vector< walberla::id_t > currentUIDs;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_)) { currentUIDs.push_back(ac_->getUid(idx)); }
      }
      WALBERLA_ASSERT(particleAndVolumeFractionSoA_.mappingUIDs == currentUIDs);

      size_t numMappedParticles = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_)) { numMappedParticles++; }
      }

      if (numMappedParticles == uint_t(0)) return;

      size_t arraySizes = numMappedParticles * sizeof(real_t) * 1;

      // Allocate unified memory for the particle information required for computing the temperature (used in
      // the solid collision operator for temperature field)
      real_t* temperatures = (real_t*) malloc(arraySizes);
      memset(temperatures, 0, arraySizes);

      // Store particle information inside memory
      size_t idxMapped = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            temperatures[idxMapped] = ac_->getTemperature(idx);
            WALBERLA_LOG_INFO("temperature from accessor is  " << ac_->getTemperature(idx));
            idxMapped++;
         }
      }

      auto nOverlappingParticlesField =
         block->getData< nOverlappingParticlesField_T >(particleAndVolumeFractionSoA_.nOverlappingParticlesFieldID);
      auto idxField = block->getData< idxField_T >(particleAndVolumeFractionSoA_.idxFieldID);
      auto particleTemperaturesField =
         block->getData< particleTemperaturesField_T >(particleAndVolumeFractionSoA_.particleTemperaturesFieldID);
      WALBERLA_FOR_ALL_CELLS_XYZ(
         particleTemperaturesField, for (uint_t p = 0; p < nOverlappingParticlesField->get(x, y, z); p++) {
            particleTemperaturesField->get(x, y, z, p) = temperatures[idxField->get(x, y, z, p)];
         })
      free(temperatures);
   }

 private:
   shared_ptr< StructuredBlockStorage > bs_;
   const shared_ptr< ParticleAccessor_T > ac_;
   const ParticleSelector_T& mappingParticleSelector_;
   ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA_;
};

template< typename ParticleAccessor_T, typename ParticleSelector_T, int Weighting_T >
class SetParticleVelocitiesSweep
{
 public:
   SetParticleVelocitiesSweep(const shared_ptr< StructuredBlockStorage >& bs,
                              const shared_ptr< ParticleAccessor_T >& ac,
                              const ParticleSelector_T& mappingParticleSelector,
                              ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA)
      : bs_(bs), ac_(ac), mappingParticleSelector_(mappingParticleSelector),
        particleAndVolumeFractionSoA_(particleAndVolumeFractionSoA)
   {}
   void operator()(IBlock* block)
   {
      // Check that uids of the particles have not changed since the last mapping to avoid incorrect indices
      std::vector< walberla::id_t > currentUIDs;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_)) { currentUIDs.push_back(ac_->getUid(idx)); }
      }
      WALBERLA_ASSERT(particleAndVolumeFractionSoA_.mappingUIDs == currentUIDs);

      size_t numMappedParticles = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_)) { numMappedParticles++; }
      }

      if (numMappedParticles == uint_t(0)) return;

      size_t arraySizes = numMappedParticles * sizeof(real_t) * 3;

      // Allocate unified memory for the particle information required for computing the velocity at a WF point (used in
      // the solid collision operator)
      real_t* linearVelocities = (real_t*) malloc(arraySizes);
      memset(linearVelocities, 0, arraySizes);
      real_t* angularVelocities = (real_t*) malloc(arraySizes);
      memset(angularVelocities, 0, arraySizes);

      // Store particle information inside memory
      size_t idxMapped = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            for (size_t d = 0; d < 3; ++d)
            {
               linearVelocities[idxMapped * 3 + d]  = ac_->getLinearVelocity(idx)[d];
               angularVelocities[idxMapped * 3 + d] = ac_->getAngularVelocity(idx)[d];
            }
            idxMapped++;
         }
      }

      auto nOverlappingParticlesField =
         block->getData< nOverlappingParticlesField_T >(particleAndVolumeFractionSoA_.nOverlappingParticlesFieldID);
      auto idxField = block->getData< idxField_T >(particleAndVolumeFractionSoA_.idxFieldID);
      auto particleVelocitiesField =
         block->getData< particleVelocitiesField_T >(particleAndVolumeFractionSoA_.particleVelocitiesFieldID);

      // For every cell, compute the particle velocities of the overlapping particles evaluated at the cell center
      const real_t dx = block->getAABB().xSize() / real_t(nOverlappingParticlesField->xSize());
      WALBERLA_FOR_ALL_CELLS_XYZ(
         particleVelocitiesField, const Vector3< real_t > cellCenter =
                                     Vector3< real_t >(real_t(x) + real_t(0.5) * dx, real_t(y) + real_t(0.5) * dx,
                                                       real_t(z) + real_t(0.5) * dx) +
                                     block->getAABB().minCorner();
         for (uint_t p = 0; p < nOverlappingParticlesField->get(x, y, z); p++) {
            Vector3< real_t > particleVelocityAtWFPoint =
               Vector3< real_t >(linearVelocities[idxField->get(x, y, z, p) * 3 + 0],
                                 linearVelocities[idxField->get(x, y, z, p) * 3 + 1],
                                 linearVelocities[idxField->get(x, y, z, p) * 3 + 2]) +
               cross(Vector3< real_t >(angularVelocities[idxField->get(x, y, z, p) * 3 + 0],
                                       angularVelocities[idxField->get(x, y, z, p) * 3 + 1],
                                       angularVelocities[idxField->get(x, y, z, p) * 3 + 2]),
                     Vector3< real_t >(
                        cellCenter[0] - particleAndVolumeFractionSoA_.positions[idxField->get(x, y, z, p) * 3 + 0],
                        cellCenter[1] - particleAndVolumeFractionSoA_.positions[idxField->get(x, y, z, p) * 3 + 1],
                        cellCenter[2] - particleAndVolumeFractionSoA_.positions[idxField->get(x, y, z, p) * 3 + 2]));
            particleVelocitiesField->get(x, y, z, p * 3 + 0) = particleVelocityAtWFPoint[0];
            particleVelocitiesField->get(x, y, z, p * 3 + 1) = particleVelocityAtWFPoint[1];
            particleVelocitiesField->get(x, y, z, p * 3 + 2) = particleVelocityAtWFPoint[2];
         })

      free(linearVelocities);
      free(angularVelocities);
   }

 private:
   shared_ptr< StructuredBlockStorage > bs_;
   const shared_ptr< ParticleAccessor_T > ac_;
   const ParticleSelector_T& mappingParticleSelector_;
   ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA_;
};

template< typename ParticleAccessor_T, typename ParticleSelector_T, int Weighting_T >
class ReduceParticleForcesSweep
{
 public:
   ReduceParticleForcesSweep(const shared_ptr< StructuredBlockStorage >& bs, const shared_ptr< ParticleAccessor_T >& ac,
                             const ParticleSelector_T& mappingParticleSelector,
                             const ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA)
      : bs_(bs), ac_(ac), mappingParticleSelector_(mappingParticleSelector),
        particleAndVolumeFractionSoA_(particleAndVolumeFractionSoA)
   {}
   void operator()(IBlock* block)
   {
      // Check that uids of the particles have not changed since the last mapping to avoid incorrect indices
      std::vector< walberla::id_t > currentUIDs;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_)) { currentUIDs.push_back(ac_->getUid(idx)); }
      }
      WALBERLA_ASSERT(particleAndVolumeFractionSoA_.mappingUIDs == currentUIDs);

      const real_t dxCurrentLevel      = bs_->dx(bs_->getLevel(*block));
      const real_t lengthScalingFactor = dxCurrentLevel;
      const real_t forceScalingFactor  = lengthScalingFactor * lengthScalingFactor;

      size_t numMappedParticles = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_)) { numMappedParticles++; }
      }

      if (numMappedParticles == uint_t(0)) return;

      size_t arraySizes = numMappedParticles * sizeof(real_t) * 3;

      // Allocate memory for the reduction of the particle forces and torques
      real_t* hydrodynamicForces = (real_t*) malloc(arraySizes);
      memset(hydrodynamicForces, 0, arraySizes);
      real_t* hydrodynamicTorques = (real_t*) malloc(arraySizes);
      memset(hydrodynamicTorques, 0, arraySizes);

      auto nOverlappingParticlesField =
         block->getData< nOverlappingParticlesField_T >(particleAndVolumeFractionSoA_.nOverlappingParticlesFieldID);
      auto idxField = block->getData< idxField_T >(particleAndVolumeFractionSoA_.idxFieldID);
      auto particleForcesField =
         block->getData< particleForcesField_T >(particleAndVolumeFractionSoA_.particleForcesFieldID);

      // For every cell, reduce the hydrodynamic forces and torques of the overlapping particles
      const real_t dx = block->getAABB().xSize() / real_t(nOverlappingParticlesField->xSize());
      // Do not use WALBERLA_FOR_ALL_CELLS_XYZ to avoid race condition if waLBerla is built with OpenMP, i.e., this loop
      // is not OpenMP parallel
#ifdef WALBERLA_BUILD_WITH_OPENMP
#   pragma omp parallel for schedule(static)
#endif
      for (cell_idx_t z = cell_idx_t(0); z < cell_idx_t(particleForcesField->zSize()); ++z)
      {
         for (cell_idx_t y = cell_idx_t(0); y < cell_idx_t(particleForcesField->ySize()); ++y)
         {
            for (cell_idx_t x = cell_idx_t(0); x < cell_idx_t(particleForcesField->xSize()); ++x)
            {
               const Vector3< real_t > cellCenter =
                  Vector3< real_t >(real_t(x) + real_t(0.5) * dx, real_t(y) + real_t(0.5) * dx,
                                    real_t(z) + real_t(0.5) * dx) +
                  block->getAABB().minCorner();
               for (uint_t p = 0; p < nOverlappingParticlesField->get(x, y, z); p++)
               {
                  Vector3< real_t > forceOnParticle(particleForcesField->get(x, y, z, p * 3 + 0),
                                                    particleForcesField->get(x, y, z, p * 3 + 1),
                                                    particleForcesField->get(x, y, z, p * 3 + 2));
                  forceOnParticle[0] *= forceScalingFactor;
                  forceOnParticle[1] *= forceScalingFactor;
                  forceOnParticle[2] *= forceScalingFactor;
#ifdef WALBERLA_BUILD_WITH_OPENMP
#   pragma omp critical
#endif
                  {
                     hydrodynamicForces[idxField->get(x, y, z, p) * 3 + 0] += forceOnParticle[0];
                     hydrodynamicForces[idxField->get(x, y, z, p) * 3 + 1] += forceOnParticle[1];
                     hydrodynamicForces[idxField->get(x, y, z, p) * 3 + 2] += forceOnParticle[2];
                  }
                  Vector3< real_t > torqueOnParticle = cross(
                     Vector3< real_t >(
                        cellCenter[0] - particleAndVolumeFractionSoA_.positions[idxField->get(x, y, z, p) * 3 + 0],
                        cellCenter[1] - particleAndVolumeFractionSoA_.positions[idxField->get(x, y, z, p) * 3 + 1],
                        cellCenter[2] - particleAndVolumeFractionSoA_.positions[idxField->get(x, y, z, p) * 3 + 2]),
                     forceOnParticle);

#ifdef WALBERLA_BUILD_WITH_OPENMP
#   pragma omp critical
#endif
                  {
                     hydrodynamicTorques[idxField->get(x, y, z, p) * 3 + 0] += torqueOnParticle[0];
                     hydrodynamicTorques[idxField->get(x, y, z, p) * 3 + 1] += torqueOnParticle[1];
                     hydrodynamicTorques[idxField->get(x, y, z, p) * 3 + 2] += torqueOnParticle[2];
                  }
               }
            }
         }
      }

      // Copy forces and torques of particles
      size_t idxMapped = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            for (size_t d = 0; d < 3; ++d)
            {
               ac_->getHydrodynamicForceRef(idx)[d] += hydrodynamicForces[idxMapped * 3 + d];
               ac_->getHydrodynamicTorqueRef(idx)[d] += hydrodynamicTorques[idxMapped * 3 + d];
            }
            idxMapped++;
         }
      }

      free(hydrodynamicForces);
      free(hydrodynamicTorques);
   }

 private:
   shared_ptr< StructuredBlockStorage > bs_;
   const shared_ptr< ParticleAccessor_T > ac_;
   const ParticleSelector_T& mappingParticleSelector_;
   const ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA_;
};

} // namespace gpu
} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
