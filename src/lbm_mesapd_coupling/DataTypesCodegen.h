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
//! \file DataTypesGPU.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
#   include "gpu/AddGPUFieldToStorage.h"
#   include "gpu/GPUField.h"
#endif

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{
namespace gpu
{

const uint_t MaxParticlesPerCell = MAX_PARTICLES_PER_CELL; // MAX_PARTICLES_PER_CELL comes from CMake

// nOverlappingParticlesField is used to store the amount of overlapping particles per cell
// B denotes the local weighting factor and is calculated by taking the sum of all local particle
// weighting factor Bs. The naming of the variables is based on the following paper:
// https://doi.org/10.1016/j.compfluid.2017.05.033
// idxField is used to store the indices of the overlapping particles
// particleVelocitiesField is used to store the velocities of the overlapping particles evaluated at the cell center
// particleForcesField is used to store the hydrodynamic forces of the cell acting on the overlapping particles

using nOverlappingParticlesField_T = GhostLayerField< uint_t, 1 >;
using BsField_T                    = GhostLayerField< real_t, MaxParticlesPerCell >;
using idxField_T                   = GhostLayerField< size_t, MaxParticlesPerCell >;
using BField_T                     = GhostLayerField< real_t, 1 >;
using particleVelocitiesField_T    = GhostLayerField< real_t, MaxParticlesPerCell * 3 >;
using particleForcesField_T        = GhostLayerField< real_t, MaxParticlesPerCell * 3 >;
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
using nOverlappingParticlesFieldGPU_T = walberla::gpu::GPUField< uint_t >;
using BsFieldGPU_T                    = walberla::gpu::GPUField< real_t >;
using idxFieldGPU_T                   = walberla::gpu::GPUField< size_t >;
using BFieldGPU_T                     = walberla::gpu::GPUField< real_t >;
using particleVelocitiesFieldGPU_T    = walberla::gpu::GPUField< real_t >;
using particleForcesFieldGPU_T        = walberla::gpu::GPUField< real_t >;
#endif

// The ParticleAndVolumeFractionSoA encapsulates the data needed by the routines involved in the coupling
template< int Weighting_T >
struct ParticleAndVolumeFractionSoA_T
{
   BlockDataID nOverlappingParticlesFieldID;
   BlockDataID BsFieldID;
   BlockDataID idxFieldID;
   BlockDataID BFieldID;
   BlockDataID particleVelocitiesFieldID;
   BlockDataID particleForcesFieldID;
   // relaxation rate omega is used for Weighting_T != 1
   real_t omega_;
   // UIDs of the particles are stored during mapping, and it is checked that they are the same during the PSM kernel.
   // This prevents running into troubles due to changed indices
   std::vector< walberla::id_t > mappingUIDs;
   // Store positions globally to avoid copying them from CPU to GPU in multiple sweeps
   real_t* positions = nullptr;

   // nrOfGhostLayers is also 1 for the fields that do not need a ghost layer since the generated sweeps can only handle
   // fields with the same number of ghost layerserated kernels)
   ParticleAndVolumeFractionSoA_T(const shared_ptr< StructuredBlockStorage >& bs, const real_t omega)
   {
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      nOverlappingParticlesFieldID = walberla::gpu::addGPUFieldToStorage< nOverlappingParticlesFieldGPU_T >(
         bs, "number of overlapping particles field GPU", uint_t(1), field::fzyx, uint_t(1), true);
      BsFieldID  = walberla::gpu::addGPUFieldToStorage< BsFieldGPU_T >(bs, "Bs field GPU", MaxParticlesPerCell,
                                                                       field::fzyx, uint_t(1), true);
      idxFieldID = walberla::gpu::addGPUFieldToStorage< idxFieldGPU_T >(bs, "idx field GPU", MaxParticlesPerCell,
                                                                        field::fzyx, uint_t(1), true);
      BFieldID = walberla::gpu::addGPUFieldToStorage< BFieldGPU_T >(bs, "B field GPU", 1, field::fzyx, uint_t(1), true);
      particleVelocitiesFieldID = walberla::gpu::addGPUFieldToStorage< particleVelocitiesFieldGPU_T >(
         bs, "particle velocities field GPU", MaxParticlesPerCell * 3, field::fzyx, uint_t(1), true);
      particleForcesFieldID = walberla::gpu::addGPUFieldToStorage< particleForcesFieldGPU_T >(
         bs, "particle forces field GPU", MaxParticlesPerCell * 3, field::fzyx, uint_t(1), true);
#else
      nOverlappingParticlesFieldID = field::addToStorage< nOverlappingParticlesField_T >(
         bs, "number of overlapping particles field CPU", uint_t(0), field::fzyx, uint_t(1), true);
      BsFieldID  = field::addToStorage< BsField_T >(bs, "Bs field CPU", real_t(0), field::fzyx, uint_t(1), true);
      idxFieldID = field::addToStorage< idxField_T >(bs, "idx field CPU", uint_t(0), field::fzyx, uint_t(1), true);
      BFieldID   = field::addToStorage< BField_T >(bs, "B field CPU", real_t(0), field::fzyx, uint_t(1), true);
      particleVelocitiesFieldID = field::addToStorage< particleVelocitiesField_T >(
         bs, "particle velocities field CPU", real_t(0), field::fzyx, uint_t(1), true);
      particleForcesFieldID = field::addToStorage< particleForcesField_T >(bs, "particle forces field CPU", real_t(0),
                                                                           field::fzyx, uint_t(1), true);
#endif
      omega_ = omega;
   }

   ~ParticleAndVolumeFractionSoA_T()
   {
      if (positions != nullptr)
      {
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         WALBERLA_GPU_CHECK(gpuFree(positions));
#else
         free(positions);
#endif
      }
   }
};

} // namespace gpu
} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
