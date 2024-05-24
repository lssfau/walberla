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
//! \file PSMUtilityGPU.cuh
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

__device__ void cross(real_t* __restrict__ const crossResult, const real_t* __restrict__ const lhs,
                      const real_t* __restrict__ const rhs)
{
   crossResult[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
   crossResult[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
   crossResult[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
}

__device__ void getVelocityAtWFPoint(real_t* __restrict__ const velocityAtWFPoint,
                                     const real_t* __restrict__ const linearVelocity,
                                     const real_t* __restrict__ const angularVelocity,
                                     const real_t* __restrict__ const position, const real_t* __restrict__ const wf_pt)
{
   real_t crossResult[3];
   real_t rhs[] = { wf_pt[0] - position[0], wf_pt[1] - position[1], wf_pt[2] - position[2] };
   cross(crossResult, angularVelocity, rhs);
   velocityAtWFPoint[0] = linearVelocity[0] + crossResult[0];
   velocityAtWFPoint[1] = linearVelocity[1] + crossResult[1];
   velocityAtWFPoint[2] = linearVelocity[2] + crossResult[2];
}

__device__ void addHydrodynamicForceAtWFPosAtomic(real_t* __restrict__ const particleForce,
                                                  real_t* __restrict__ const particleTorque,
                                                  const real_t* __restrict__ const f,
                                                  const real_t* __restrict__ const pos,
                                                  const real_t* __restrict__ const wf_pt)
{
#if defined(WALBERLA_BUILD_WITH_CUDA)
   atomicAdd(&(particleForce[0]), f[0]);
   atomicAdd(&(particleForce[1]), f[1]);
   atomicAdd(&(particleForce[2]), f[2]);
#endif

   // Using unsafeAtomicAdd ensures that HW FP Atomics are used instead of CAS loops, see:
   // https://fs.hlrs.de/projects/par/events/2023/GPU-AMD/day3/11.%20AMD_Node_Memory_Model.pdf
#ifdef WALBERLA_BUILD_WITH_HIP
   unsafeAtomicAdd(&(particleForce[0]), f[0]);
   unsafeAtomicAdd(&(particleForce[1]), f[1]);
   unsafeAtomicAdd(&(particleForce[2]), f[2]);
#endif

   real_t torque[] = { 0.0, 0.0, 0.0 };
   real_t lhs[]    = { wf_pt[0] - pos[0], wf_pt[1] - pos[1], wf_pt[2] - pos[2] };
   cross(torque, lhs, f);

#if defined(WALBERLA_BUILD_WITH_CUDA)
   atomicAdd(&(particleTorque[0]), torque[0]);
   atomicAdd(&(particleTorque[1]), torque[1]);
   atomicAdd(&(particleTorque[2]), torque[2]);
#endif

#ifdef WALBERLA_BUILD_WITH_HIP
   unsafeAtomicAdd(&(particleTorque[0]), torque[0]);
   unsafeAtomicAdd(&(particleTorque[1]), torque[1]);
   unsafeAtomicAdd(&(particleTorque[2]), torque[2]);
#endif
}

} // namespace gpu
} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
