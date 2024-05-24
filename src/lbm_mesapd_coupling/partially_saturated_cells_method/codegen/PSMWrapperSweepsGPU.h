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
//! \file PSMWrapperSweepsGPU.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

#include "gpu/FieldIndexing.h"
#include "gpu/GPUField.h"
#include "gpu/Kernel.h"
#include "gpu/sweeps/GPUSweepBase.h"

#include "lbm/sweeps/StreamPull.h"
#include "lbm/sweeps/SweepBase.h"

#include "lbm_mesapd_coupling/DataTypesCodegen.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/common/ParticleFunctions.h"

#include "timeloop/SweepTimeloop.h"

#include <cassert>

#include "PSMWrapperKernels.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{
namespace gpu
{

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
      real_t* linearVelocities;
      WALBERLA_GPU_CHECK(gpuMallocManaged(&linearVelocities, arraySizes));
      WALBERLA_GPU_CHECK(gpuMemset(linearVelocities, 0, arraySizes));
      real_t* angularVelocities;
      WALBERLA_GPU_CHECK(gpuMallocManaged(&angularVelocities, arraySizes));
      WALBERLA_GPU_CHECK(gpuMemset(angularVelocities, 0, arraySizes));

      // Store particle information inside unified memory to communicate information to the GPU
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
         block->getData< nOverlappingParticlesFieldGPU_T >(particleAndVolumeFractionSoA_.nOverlappingParticlesFieldID);
      auto idxField = block->getData< idxFieldGPU_T >(particleAndVolumeFractionSoA_.idxFieldID);
      auto particleVelocitiesField =
         block->getData< particleVelocitiesFieldGPU_T >(particleAndVolumeFractionSoA_.particleVelocitiesFieldID);

      // For every cell, compute the particle velocities of the overlapping particles evaluated at the cell center
      auto velocitiesKernel = walberla::gpu::make_kernel(&(SetParticleVelocities));
      velocitiesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< uint_t >::xyz(*nOverlappingParticlesField));
      velocitiesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< id_t >::xyz(*idxField));
      velocitiesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< real_t >::xyz(*particleVelocitiesField));
      velocitiesKernel.addParam(linearVelocities);
      velocitiesKernel.addParam(angularVelocities);
      velocitiesKernel.addParam(particleAndVolumeFractionSoA_.positions);
      const double3 blockStart = { block->getAABB().minCorner()[0], block->getAABB().minCorner()[1],
                                   block->getAABB().minCorner()[2] };
      velocitiesKernel.addParam(blockStart);
      velocitiesKernel.addParam(block->getAABB().xSize() / real_t(nOverlappingParticlesField->xSize()));
      velocitiesKernel();

      WALBERLA_GPU_CHECK(gpuFree(linearVelocities));
      WALBERLA_GPU_CHECK(gpuFree(angularVelocities));
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

      // TODO: for multiple blocks per process, this data is transferred multiple times per time step (unnecessarily)
      // Allocate unified memory for the reduction of the particle forces and torques on the GPU
      real_t* hydrodynamicForces;
      WALBERLA_GPU_CHECK(gpuMallocManaged(&hydrodynamicForces, arraySizes));
      // Using unsafeAtomicAdd() required coarse grained memory, see:
      // https://fs.hlrs.de/projects/par/events/2023/GPU-AMD/day3/11.%20AMD_Node_Memory_Model.pdf
#ifdef WALBERLA_BUILD_WITH_HIP
      int deviceId = -1;
      WALBERLA_GPU_CHECK(hipGetDevice(&deviceId));
      WALBERLA_GPU_CHECK(hipMemAdvise(hydrodynamicForces, arraySizes, hipMemAdviseSetCoarseGrain, deviceId));
#endif
      WALBERLA_GPU_CHECK(gpuMemset(hydrodynamicForces, 0, arraySizes));
      real_t* hydrodynamicTorques;
      WALBERLA_GPU_CHECK(gpuMallocManaged(&hydrodynamicTorques, arraySizes));
#ifdef WALBERLA_BUILD_WITH_HIP
      WALBERLA_GPU_CHECK(hipMemAdvise(hydrodynamicTorques, arraySizes, hipMemAdviseSetCoarseGrain, deviceId));
#endif
      WALBERLA_GPU_CHECK(gpuMemset(hydrodynamicTorques, 0, arraySizes));

      auto nOverlappingParticlesField =
         block->getData< nOverlappingParticlesFieldGPU_T >(particleAndVolumeFractionSoA_.nOverlappingParticlesFieldID);
      auto idxField = block->getData< idxFieldGPU_T >(particleAndVolumeFractionSoA_.idxFieldID);
      auto particleForcesField =
         block->getData< particleForcesFieldGPU_T >(particleAndVolumeFractionSoA_.particleForcesFieldID);

      const double3 blockStart = { block->getAABB().minCorner()[0], block->getAABB().minCorner()[1],
                                   block->getAABB().minCorner()[2] };

      // For every cell, reduce the hydrodynamic forces and torques of the overlapping particles
      auto forcesKernel = walberla::gpu::make_kernel(&(ReduceParticleForces));
      forcesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< uint_t >::xyz(*nOverlappingParticlesField));
      forcesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< id_t >::xyz(*idxField));
      forcesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< real_t >::xyz(*particleForcesField));
      forcesKernel.addParam(hydrodynamicForces);
      forcesKernel.addParam(hydrodynamicTorques);
      forcesKernel.addParam(particleAndVolumeFractionSoA_.positions);
      forcesKernel.addParam(blockStart);
      forcesKernel.addParam(block->getAABB().xSize() / real_t(nOverlappingParticlesField->xSize()));
      forcesKernel.addParam(forceScalingFactor);
      forcesKernel();

      WALBERLA_GPU_CHECK(gpuDeviceSynchronize());

      // Copy forces and torques of particles from GPU to CPU
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

      WALBERLA_GPU_CHECK(gpuFree(hydrodynamicForces));
      WALBERLA_GPU_CHECK(gpuFree(hydrodynamicTorques));
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
