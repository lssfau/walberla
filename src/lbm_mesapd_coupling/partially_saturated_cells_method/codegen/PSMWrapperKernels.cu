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
//! \file PSMWrapperKernels.cu
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \brief Provide two kernels that need to be called before and after the PSM kernel
//
//======================================================================================================================

#include "PSMUtilityGPU.h"
#include "PSMWrapperKernels.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{
namespace gpu
{

__global__ void SetParticleVelocities(walberla::gpu::FieldAccessor< uint_t > nOverlappingParticlesField,
                                      walberla::gpu::FieldAccessor< uint_t > idxField,
                                      walberla::gpu::FieldAccessor< real_t > particleVelocitiesField,
                                      real_t* __restrict__ const linearVelocities,
                                      real_t* __restrict__ const angularVelocities,
                                      real_t* __restrict__ const positions, const double3 blockStart, const real_t dx)
{
   const uint3 blockIdx_uint3  = make_uint3(blockIdx.x, blockIdx.y, blockIdx.z);
   const uint3 threadIdx_uint3 = make_uint3(threadIdx.x, threadIdx.y, threadIdx.z);

   nOverlappingParticlesField.set(blockIdx_uint3, threadIdx_uint3);
   idxField.set(blockIdx_uint3, threadIdx_uint3);
   particleVelocitiesField.set(blockIdx_uint3, threadIdx_uint3);

   // Cell center is needed in order to compute the particle velocity at this WF point
   const real_t cellCenter[] = { real_t(blockStart.x + (threadIdx.x + real_t(0.5)) * dx), // NOLINT(*-avoid-c-arrays)
                                 real_t(blockStart.y + (blockIdx.x + real_t(0.5)) * dx),
                                 real_t(blockStart.z + (blockIdx.y + real_t(0.5)) * dx) };

   // Compute the particle velocity at the cell center for all overlapping particles
   for (uint_t p = 0; p < nOverlappingParticlesField.get(); p++)
   {
      real_t particleVelocityAtWFPoint[] = { real_t(0.0), real_t(0.0), real_t(0.0) }; // NOLINT(*-avoid-c-arrays)
      getVelocityAtWFPoint(particleVelocityAtWFPoint, &linearVelocities[idxField.get(p) * 3],
                           &angularVelocities[idxField.get(p) * 3], &positions[idxField.get(p) * 3], cellCenter);
      particleVelocitiesField.get(p * 3 + 0) = particleVelocityAtWFPoint[0];
      particleVelocitiesField.get(p * 3 + 1) = particleVelocityAtWFPoint[1];
      particleVelocitiesField.get(p * 3 + 2) = particleVelocityAtWFPoint[2];
   }
}

__global__ void ReduceParticleForces(walberla::gpu::FieldAccessor< uint_t > nOverlappingParticlesField,
                                     walberla::gpu::FieldAccessor< id_t > idxField,
                                     walberla::gpu::FieldAccessor< real_t > particleForcesField,
                                     real_t* __restrict__ const hydrodynamicForces,
                                     real_t* __restrict__ const hydrodynamicTorques,
                                     real_t* __restrict__ const positions, const double3 blockStart, const real_t dx,
                                     const real_t forceScalingFactor)
{
   const uint3 blockIdx_uint3  = make_uint3(blockIdx.x, blockIdx.y, blockIdx.z);
   const uint3 threadIdx_uint3 = make_uint3(threadIdx.x, threadIdx.y, threadIdx.z);

   nOverlappingParticlesField.set(blockIdx_uint3, threadIdx_uint3);
   idxField.set(blockIdx_uint3, threadIdx_uint3);
   particleForcesField.set(blockIdx_uint3, threadIdx_uint3);

   // Cell center is needed in order to compute the particle velocity at this WF point
   const real_t cellCenter[] = { real_t(blockStart.x + (threadIdx.x + real_t(0.5)) * dx), // NOLINT(*-avoid-c-arrays)
                                 real_t(blockStart.y + (blockIdx.x + real_t(0.5)) * dx),
                                 real_t(blockStart.z + (blockIdx.y + real_t(0.5)) * dx) };

   // Reduce the forces for all overlapping particles
   for (uint_t p = 0; p < nOverlappingParticlesField.get(); p++)
   {
      real_t forceOnParticle[] = { particleForcesField.get(p * 3 + 0), // NOLINT(*-avoid-c-arrays)
                                   particleForcesField.get(p * 3 + 1), particleForcesField.get(p * 3 + 2) };
      forceOnParticle[0] *= forceScalingFactor;
      forceOnParticle[1] *= forceScalingFactor;
      forceOnParticle[2] *= forceScalingFactor;
      addHydrodynamicForceTorqueAtWFPosAtomic(&hydrodynamicForces[idxField.get(p) * 3],
                                              &hydrodynamicTorques[idxField.get(p) * 3], forceOnParticle,
                                              &positions[idxField.get(p) * 3], cellCenter);
   }
}

} // namespace gpu
} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla