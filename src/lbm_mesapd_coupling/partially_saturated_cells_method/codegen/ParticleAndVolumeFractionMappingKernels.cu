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
//! \file ParticleAndVolumeFractionMappingKernels.cu
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#include "lbm_mesapd_coupling/DataTypesCodegen.h"

#include <cassert>

#include "ParticleAndVolumeFractionMappingKernels.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{
namespace gpu
{

// Functions to calculate Bs
template< int Weighting_T >
__device__ void calculateWeighting(real_t* __restrict__ const  /*weighting*/, const real_t& /*epsilon*/,
                                   const real_t& /*tau*/)
{
   WALBERLA_STATIC_ASSERT(Weighting_T == 1 || Weighting_T == 2);
}
template<>
__device__ void calculateWeighting< 1 >(real_t* __restrict__ const weighting, const real_t& epsilon,
                                        const real_t& /*tau*/)
{
   *weighting = epsilon;
}
template<>
__device__ void calculateWeighting< 2 >(real_t* __restrict__ const weighting, const real_t& epsilon, const real_t& tau)
{
   *weighting = epsilon * (tau - real_t(0.5)) / ((real_t(1) - epsilon) + (tau - real_t(0.5)));
}

template< int Weighting_T >
__global__ void superSampling(walberla::gpu::FieldAccessor< uint_t > nOverlappingParticlesField,
                              walberla::gpu::FieldAccessor< real_t > BsField,
                              walberla::gpu::FieldAccessor< id_t > idxField,
                              walberla::gpu::FieldAccessor< real_t > BField, const real_t omega,
                              const real_t* __restrict__ const spherePositions,
                              const real_t* __restrict__ const sphereRadii, const double3 blockStart, const real_t dx,
                              const int3 nSamples, const size_t* __restrict__ const numParticlesSubBlocks,
                              const size_t* __restrict__ const particleIDsSubBlocks, const uint3 subBlocksPerDim)
{
   const uint3 blockIdx_uint3  = make_uint3(blockIdx.x, blockIdx.y, blockIdx.z);
   const uint3 threadIdx_uint3 = make_uint3(threadIdx.x, threadIdx.y, threadIdx.z);

   nOverlappingParticlesField.set(blockIdx_uint3, threadIdx_uint3);
   BsField.set(blockIdx_uint3, threadIdx_uint3);
   idxField.set(blockIdx_uint3, threadIdx_uint3);
   BField.set(blockIdx_uint3, threadIdx_uint3);

   // Clear the fields
   for (uint i = 0; i < MaxParticlesPerCell; i++)
   {
      BsField.get(i)  = real_t(0.0);
      idxField.get(i) = size_t(0);
   }
   nOverlappingParticlesField.get() = uint_t(0);
   BField.get()                     = real_t(0.0);

   double3 sampleDistance = { 1.0 / (nSamples.x + 1) * dx, 1.0 / (nSamples.y + 1) * dx, 1.0 / (nSamples.z + 1) * dx };
   double3 startSamplingPoint = { (blockStart.x + threadIdx.x * dx + sampleDistance.x),
                                  (blockStart.y + blockIdx.x * dx + sampleDistance.y),
                                  (blockStart.z + blockIdx.y * dx + sampleDistance.z) };
   const ulong3 subBlockIndex = { size_t(real_t(threadIdx.x) / blockDim.x * real_t(subBlocksPerDim.x)),
                                  size_t(real_t(blockIdx.x) / gridDim.x * real_t(subBlocksPerDim.y)),
                                  size_t(real_t(blockIdx.y) / gridDim.y * real_t(subBlocksPerDim.z)) };
   size_t linearizedSubBlockIndex =
      subBlockIndex.z * subBlocksPerDim.x * subBlocksPerDim.y + subBlockIndex.y * subBlocksPerDim.x + subBlockIndex.x;

   for (uint i = 0; i < numParticlesSubBlocks[linearizedSubBlockIndex]; i++)
   {
      // SoA
      size_t idxMapped =
         particleIDsSubBlocks[linearizedSubBlockIndex + i * subBlocksPerDim.x * subBlocksPerDim.y * subBlocksPerDim.z];
      double3 currentSamplingPoint = startSamplingPoint;

      double3 minCornerSphere = { spherePositions[idxMapped * 3] - sphereRadii[idxMapped],
                                  spherePositions[idxMapped * 3 + 1] - sphereRadii[idxMapped],
                                  spherePositions[idxMapped * 3 + 2] - sphereRadii[idxMapped] };
      double3 maxCornerSphere = { spherePositions[idxMapped * 3] + sphereRadii[idxMapped],
                                  spherePositions[idxMapped * 3 + 1] + sphereRadii[idxMapped],
                                  spherePositions[idxMapped * 3 + 2] + sphereRadii[idxMapped] };

      double overlapFraction = 0.0;

      if (startSamplingPoint.x + dx > minCornerSphere.x && startSamplingPoint.x < maxCornerSphere.x &&
          startSamplingPoint.y + dx > minCornerSphere.y && startSamplingPoint.y < maxCornerSphere.y &&
          startSamplingPoint.z + dx > minCornerSphere.z && startSamplingPoint.z < maxCornerSphere.z)
      {
         for (uint_t z = 0; z < nSamples.z; z++)
         {
            currentSamplingPoint.y = startSamplingPoint.y;
            for (uint_t y = 0; y < nSamples.y; y++)
            {
               currentSamplingPoint.x = startSamplingPoint.x;
               for (uint_t x = 0; x < nSamples.x; x++)
               {
                  if ((currentSamplingPoint.x - spherePositions[idxMapped * 3]) *
                            (currentSamplingPoint.x - spherePositions[idxMapped * 3]) +
                         (currentSamplingPoint.y - spherePositions[idxMapped * 3 + 1]) *
                            (currentSamplingPoint.y - spherePositions[idxMapped * 3 + 1]) +
                         (currentSamplingPoint.z - spherePositions[idxMapped * 3 + 2]) *
                            (currentSamplingPoint.z - spherePositions[idxMapped * 3 + 2]) <=
                      sphereRadii[idxMapped] * sphereRadii[idxMapped])
                  {
                     overlapFraction += 1.0;
                  }
                  currentSamplingPoint.x += sampleDistance.x;
               }
               currentSamplingPoint.y += sampleDistance.y;
            }
            currentSamplingPoint.z += sampleDistance.z;
         }

         // store overlap fraction only if there is an intersection
         if (overlapFraction > 0.0)
         {
            assert(nOverlappingParticlesField.get() < MaxParticlesPerCell);
            BsField.get(nOverlappingParticlesField.get()) = overlapFraction;
            BsField.get(nOverlappingParticlesField.get()) *= 1.0 / (nSamples.x * nSamples.y * nSamples.z);
            calculateWeighting< Weighting_T >(&BsField.get(nOverlappingParticlesField.get()),
                                              BsField.get(nOverlappingParticlesField.get()), real_t(1.0) / omega);
            idxField.get(nOverlappingParticlesField.get()) = idxMapped;
            BField.get() += BsField.get(nOverlappingParticlesField.get());
            nOverlappingParticlesField.get() += 1;
         }
      }
   }

   // Normalize fraction field (Bs) if sum over all fractions (B) > 1
   if (BField.get() > 1)
   {
      for (uint i = 0; i < nOverlappingParticlesField.get(); i++)
      {
         BsField.get(i) /= BField.get();
      }
      BField.get() = 1.0;
   }
}

// Based on the following paper: https://doi.org/10.1108/EC-02-2016-0052
template< int Weighting_T >
__global__ void
   linearApproximation(walberla::gpu::FieldAccessor< uint_t > nOverlappingParticlesField,
                       walberla::gpu::FieldAccessor< real_t > BsField, walberla::gpu::FieldAccessor< id_t > idxField,
                       walberla::gpu::FieldAccessor< real_t > BField, const real_t omega,
                       const real_t* __restrict__ const spherePositions, const real_t* __restrict__ const sphereRadii,
                       const real_t* __restrict__ const f_rs, const double3 blockStart, const real_t dx,
                       const size_t* __restrict__ const numParticlesSubBlocks,
                       const size_t* __restrict__ const particleIDsSubBlocks, const uint3 subBlocksPerDim)
{
   const uint3 blockIdx_uint3  = make_uint3(blockIdx.x, blockIdx.y, blockIdx.z);
   const uint3 threadIdx_uint3 = make_uint3(threadIdx.x, threadIdx.y, threadIdx.z);

   nOverlappingParticlesField.set(blockIdx_uint3, threadIdx_uint3);
   BsField.set(blockIdx_uint3, threadIdx_uint3);
   idxField.set(blockIdx_uint3, threadIdx_uint3);
   BField.set(blockIdx_uint3, threadIdx_uint3);

   // Clear the fields
   for (uint i = 0; i < MaxParticlesPerCell; i++)
   {
      BsField.get(i)  = real_t(0.0);
      idxField.get(i) = size_t(0);
   }
   nOverlappingParticlesField.get() = uint_t(0);
   BField.get()                     = real_t(0.0);

   const double3 cellCenter   = { (blockStart.x + (threadIdx.x + 0.5) * dx), (blockStart.y + (blockIdx.x + 0.5) * dx),
                                  (blockStart.z + (blockIdx.y + 0.5) * dx) };
   const ulong3 subBlockIndex = { size_t(real_t(threadIdx.x) / blockDim.x * real_t(subBlocksPerDim.x)),
                                  size_t(real_t(blockIdx.x) / gridDim.x * real_t(subBlocksPerDim.y)),
                                  size_t(real_t(blockIdx.y) / gridDim.y * real_t(subBlocksPerDim.z)) };
   size_t linearizedSubBlockIndex =
      subBlockIndex.z * subBlocksPerDim.x * subBlocksPerDim.y + subBlockIndex.y * subBlocksPerDim.x + subBlockIndex.x;

   for (uint i = 0; i < numParticlesSubBlocks[linearizedSubBlockIndex]; i++)
   {
      size_t idxMapped =
         particleIDsSubBlocks[linearizedSubBlockIndex + i * subBlocksPerDim.x * subBlocksPerDim.y * subBlocksPerDim.z];
      double3 minCornerSphere = { spherePositions[idxMapped * 3] - sphereRadii[idxMapped],
                                  spherePositions[idxMapped * 3 + 1] - sphereRadii[idxMapped],
                                  spherePositions[idxMapped * 3 + 2] - sphereRadii[idxMapped] };
      double3 maxCornerSphere = { spherePositions[idxMapped * 3] + sphereRadii[idxMapped],
                                  spherePositions[idxMapped * 3 + 1] + sphereRadii[idxMapped],
                                  spherePositions[idxMapped * 3 + 2] + sphereRadii[idxMapped] };
      if (cellCenter.x + dx > minCornerSphere.x && cellCenter.x - dx < maxCornerSphere.x &&
          cellCenter.y + dx > minCornerSphere.y && cellCenter.y - dx < maxCornerSphere.y &&
          cellCenter.z + dx > minCornerSphere.z && cellCenter.z - dx < maxCornerSphere.z)
      {
         const double3 cellSphereVector = { spherePositions[idxMapped * 3] - cellCenter.x,
                                            spherePositions[idxMapped * 3 + 1] - cellCenter.y,
                                            spherePositions[idxMapped * 3 + 2] - cellCenter.z };

         const real_t D = sqrt(cellSphereVector.x * cellSphereVector.x + cellSphereVector.y * cellSphereVector.y +
                               cellSphereVector.z * cellSphereVector.z) -
                          sphereRadii[idxMapped];

         real_t epsilon = -D + f_rs[idxMapped];
         epsilon        = max(epsilon, 0.0);
         epsilon        = min(epsilon, 1.0);

         // Store overlap fraction only if there is an intersection
         if (epsilon > 0.0)
         {
            // Check that the maximum number of overlapping particles has not yet been reached
            assert(nOverlappingParticlesField.get() < MaxParticlesPerCell);
            BsField.get(nOverlappingParticlesField.get()) = epsilon;
            calculateWeighting< Weighting_T >(&BsField.get(nOverlappingParticlesField.get()),
                                              BsField.get(nOverlappingParticlesField.get()), real_t(1.0) / omega);
            idxField.get(nOverlappingParticlesField.get()) = idxMapped;
            BField.get() += BsField.get(nOverlappingParticlesField.get());
            nOverlappingParticlesField.get() += 1;
         }
      }
   }

   // Normalize fraction field (Bs) if sum over all fractions (B) > 1
   if (BField.get() > 1)
   {
      for (uint i = 0; i < nOverlappingParticlesField.get(); i++)
      {
         BsField.get(i) /= BField.get();
      }
      BField.get() = 1.0;
   }
}

template< int Weighting_T >
__global__ void boxMapping(walberla::gpu::FieldAccessor< uint_t > nOverlappingParticlesField,
                           walberla::gpu::FieldAccessor< real_t > BsField,
                           walberla::gpu::FieldAccessor< id_t > idxField, walberla::gpu::FieldAccessor< real_t > BField,
                           const real_t omega, const double3 boxPositionMin, const double3 boxPositionMax,
                           const double3 blockStart, const real_t dx, const id_t idxMapped)
{
   const uint3 blockIdx_uint3  = make_uint3(blockIdx.x, blockIdx.y, blockIdx.z);
   const uint3 threadIdx_uint3 = make_uint3(threadIdx.x, threadIdx.y, threadIdx.z);

   nOverlappingParticlesField.set(blockIdx_uint3, threadIdx_uint3);
   BsField.set(blockIdx_uint3, threadIdx_uint3);
   idxField.set(blockIdx_uint3, threadIdx_uint3);
   BField.set(blockIdx_uint3, threadIdx_uint3);

   const double3 cellCenter = { (blockStart.x + (threadIdx.x + 0.5) * dx), (blockStart.y + (blockIdx.x + 0.5) * dx),
                                (blockStart.z + (blockIdx.y + 0.5) * dx) };
   const double3 cellMin    = { cellCenter.x - dx * real_t(0.5), cellCenter.y - dx * real_t(0.5),
                                cellCenter.z - dx * real_t(0.5) };
   const double3 cellMax    = { cellCenter.x + dx * real_t(0.5), cellCenter.y + dx * real_t(0.5),
                                cellCenter.z + dx * real_t(0.5) };

   const real_t xOverlap        = max(real_t(0), min(boxPositionMax.x, cellMax.x) - max(boxPositionMin.x, cellMin.x));
   const real_t yOverlap        = max(real_t(0), min(boxPositionMax.y, cellMax.y) - max(boxPositionMin.y, cellMin.y));
   const real_t zOverlap        = max(real_t(0), min(boxPositionMax.z, cellMax.z) - max(boxPositionMin.z, cellMin.z));
   const real_t overlapFraction = xOverlap * yOverlap * zOverlap / (dx * dx * dx);

   if (overlapFraction > real_t(0))
   {
      assert(nOverlappingParticlesField.get() < MaxParticlesPerCell);

      BsField.get(nOverlappingParticlesField.get()) = overlapFraction;
      calculateWeighting< Weighting_T >(&BsField.get(nOverlappingParticlesField.get()),
                                        BsField.get(nOverlappingParticlesField.get()), real_t(1.0) / omega);
      idxField.get(nOverlappingParticlesField.get()) = idxMapped;
      BField.get() += BsField.get(nOverlappingParticlesField.get());
      nOverlappingParticlesField.get() += 1;

      // TODO: it can happen that the BsField for spheres is normalized twice, one here and in the sphere mapping
      // Normalize fraction field (Bs) if sum over all fractions (B) > 1
      if (BField.get() > 1)
      {
         for (uint i = 0; i < nOverlappingParticlesField.get(); i++)
         {
            BsField.get(i) /= BField.get();
         }
         BField.get() = 1.0;
      }
   }
}

auto instance0_with_weighting_1 = superSampling< 1 >;
auto instance1_with_weighting_2 = superSampling< 2 >;
auto instance2_with_weighting_1 = linearApproximation< 1 >;
auto instance3_with_weighting_2 = linearApproximation< 2 >;
auto instance4_with_weighting_1 = boxMapping< 1 >;
auto instance5_with_weighting_2 = boxMapping< 2 >;

} // namespace gpu
} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
