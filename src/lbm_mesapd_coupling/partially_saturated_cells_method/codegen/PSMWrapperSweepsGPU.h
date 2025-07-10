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
class SetParticleTemperaturesSweep
{
 public:
   SetParticleTemperaturesSweep(const shared_ptr< StructuredBlockStorage >& bs,
                                const shared_ptr< ParticleAccessor_T >& ac,
                                const ParticleSelector_T& mappingParticleSelector,
                                ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA,
                                BlockDataID & densityConcentrationFieldCPUGPUID,bool uniformParticleTemperature = false)
      : bs_(bs), ac_(ac), mappingParticleSelector_(mappingParticleSelector),
        particleAndVolumeFractionSoA_(particleAndVolumeFractionSoA), densityConcentrationFieldCPUGPUID_(densityConcentrationFieldCPUGPUID),
        uniformParticleTemperature_(uniformParticleTemperature)
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
      real_t* temperatures_h = (real_t*) malloc(arraySizes);
      memset(temperatures_h, 0, arraySizes);

      // Store particle information inside memory
      size_t idxMapped = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            temperatures_h[idxMapped] = ac_->getTemperature(idx);
            idxMapped++;
         }
      }

      real_t* temperatures;

      WALBERLA_GPU_CHECK(gpuMalloc(&temperatures, arraySizes));
      WALBERLA_GPU_CHECK(gpuMemcpy(temperatures, temperatures_h, arraySizes, gpuMemcpyHostToDevice));
      auto nOverlappingParticlesField =
         block->getData< nOverlappingParticlesFieldGPU_T >(particleAndVolumeFractionSoA_.nOverlappingParticlesFieldID);
      auto idxField = block->getData< idxFieldGPU_T >(particleAndVolumeFractionSoA_.idxFieldID);
      auto particleTemperaturesField =
         block->getData< particleTemperaturesFieldGPU_T >(particleAndVolumeFractionSoA_.particleTemperaturesFieldID);
      WALBERLA_LOG_INFO_ON_ROOT("set temperatures reached till here on GPU");
      auto BsField =
         block->getData< BsFieldGPU_T >(particleAndVolumeFractionSoA_.BsFieldID);

      // For every cell, compute the particle velocities of the overlapping particles evaluated at the cell center
      auto temperaturesKernel = walberla::gpu::make_kernel(&(SetParticleTemperatures));
      temperaturesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< uint_t >::xyz(*nOverlappingParticlesField));
      temperaturesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< id_t >::xyz(*idxField));
      temperaturesKernel.addFieldIndexingParam(walberla::gpu::FieldIndexing< real_t >::xyz(*particleTemperaturesField));
      temperaturesKernel.addParam(temperatures);

      const double3 blockStart = { block->getAABB().minCorner()[0], block->getAABB().minCorner()[1],
                                   block->getAABB().minCorner()[2] };
      //temperaturesKernel.addParam(blockStart);
      //temperaturesKernel.addParam(block->getAABB().xSize() / real_t(nOverlappingParticlesField->xSize()));
      temperaturesKernel();
      WALBERLA_GPU_CHECK(gpuFree(temperatures));
      free(temperatures_h);;

   }

 private:
   shared_ptr< StructuredBlockStorage > bs_;
   const shared_ptr< ParticleAccessor_T > ac_;
   const ParticleSelector_T& mappingParticleSelector_;
   ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA_;
   const BlockDataID &densityConcentrationFieldCPUGPUID_;
   const bool uniformParticleTemperature_;
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

      // Allocate memory for the particle information required for computing the velocity at a WF point (used in
      // the solid collision operator)
      real_t* linearVelocities_h  = (real_t*) malloc(arraySizes);
      real_t* angularVelocities_h = (real_t*) malloc(arraySizes);

      // Store particle information inside memory to communicate information to the GPU
      size_t idxMapped = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            for (size_t d = 0; d < 3; ++d)
            {
               linearVelocities_h[idxMapped * 3 + d]  = ac_->getLinearVelocity(idx)[d];
               angularVelocities_h[idxMapped * 3 + d] = ac_->getAngularVelocity(idx)[d];
            }
            idxMapped++;
         }
      }

      real_t* linearVelocities;
      WALBERLA_GPU_CHECK(gpuMalloc(&linearVelocities, arraySizes));
      WALBERLA_GPU_CHECK(gpuMemcpy(linearVelocities, linearVelocities_h, arraySizes, gpuMemcpyHostToDevice));
      real_t* angularVelocities;
      WALBERLA_GPU_CHECK(gpuMalloc(&angularVelocities, arraySizes));
      WALBERLA_GPU_CHECK(gpuMemcpy(angularVelocities, angularVelocities_h, arraySizes, gpuMemcpyHostToDevice));

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
      free(linearVelocities_h);

      WALBERLA_GPU_CHECK(gpuFree(angularVelocities));
      free(angularVelocities_h);
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

      // Allocate memory for the reduction of the particle forces and torques on the GPU
      real_t* hydrodynamicForces;
      WALBERLA_GPU_CHECK(gpuMalloc(&hydrodynamicForces, arraySizes));
      WALBERLA_GPU_CHECK(gpuMemset(hydrodynamicForces, 0, arraySizes));

      real_t* hydrodynamicTorques;
      WALBERLA_GPU_CHECK(gpuMalloc(&hydrodynamicTorques, arraySizes));
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

      real_t* hydrodynamicForces_h = (real_t*) malloc(arraySizes);
      WALBERLA_GPU_CHECK(gpuMemcpy(hydrodynamicForces_h, hydrodynamicForces, arraySizes, gpuMemcpyDeviceToHost));

      real_t* hydrodynamicTorques_h = (real_t*) malloc(arraySizes);
      WALBERLA_GPU_CHECK(gpuMemcpy(hydrodynamicTorques_h, hydrodynamicTorques, arraySizes, gpuMemcpyDeviceToHost));

      // Copy forces and torques of particles from GPU to CPU
      size_t idxMapped = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            for (size_t d = 0; d < 3; ++d)
            {
               ac_->getHydrodynamicForceRef(idx)[d] += hydrodynamicForces_h[idxMapped * 3 + d];
               ac_->getHydrodynamicTorqueRef(idx)[d] += hydrodynamicTorques_h[idxMapped * 3 + d];
            }
            idxMapped++;
         }
      }

      WALBERLA_GPU_CHECK(gpuFree(hydrodynamicForces));
      free(hydrodynamicForces_h);

      WALBERLA_GPU_CHECK(gpuFree(hydrodynamicTorques));
      free(hydrodynamicTorques_h);
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
