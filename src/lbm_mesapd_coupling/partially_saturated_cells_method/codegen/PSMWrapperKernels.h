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
//! \file PSMWrapperKernels.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#include "gpu/FieldAccessor.h"

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
                                      real_t* __restrict__ const positions, const double3 blockStart,
                                      const real_t dx);

__global__ void ReduceParticleForces(walberla::gpu::FieldAccessor< uint_t > nOverlappingParticlesField,
                                     walberla::gpu::FieldAccessor< id_t > idxField,
                                     walberla::gpu::FieldAccessor< real_t > particleForcesField,
                                     real_t* __restrict__ const hydrodynamicForces,
                                     real_t* __restrict__ const hydrodynamicTorques,
                                     real_t* __restrict__ const positions, const double3 blockStart, const real_t dx,
                                     const real_t forceScalingFactor);

} // namespace gpu
} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
