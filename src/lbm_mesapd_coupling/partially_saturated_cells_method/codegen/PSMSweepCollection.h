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
//! \file PSMSweepCollection.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
#   include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMWrapperSweepsGPU.h"
#   include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/ParticleAndVolumeFractionMappingSweepsGPU.h"
#else
#   include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMWrapperSweepsCPU.h"
#   include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/ParticleAndVolumeFractionMappingSweepsCPU.h"
#endif

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{
namespace gpu
{

// The deviceSyncWrapper can be used so that the timeloop measures the correct device runtime
inline auto deviceSyncWrapper = [](std::function< void(IBlock*) > sweep) {
   return [sweep](IBlock* b) {
      sweep(b);
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      WALBERLA_GPU_CHECK(gpuDeviceSynchronize());
#endif
   };
};

template< typename ParticleAccessor_T, typename ParticleSelector_T, int Weighting_T >
class PSMSweepCollection
{
 public:
   PSMSweepCollection(const shared_ptr< StructuredBlockStorage >& bs, const shared_ptr< ParticleAccessor_T >& ac,
                      const ParticleSelector_T& ps,
                      ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA,
                      const Vector3< uint_t > particleSubBlockSize = Vector3< uint_t >(10))
      : particleMappingSweep(SphereFractionMappingSweep< ParticleAccessor_T, ParticleSelector_T, Weighting_T >(
           bs, ac, ps, particleAndVolumeFractionSoA, particleSubBlockSize)),
        setParticleVelocitiesSweep(SetParticleVelocitiesSweep< ParticleAccessor_T, ParticleSelector_T, Weighting_T >(
           bs, ac, ps, particleAndVolumeFractionSoA)),
        reduceParticleForcesSweep(ReduceParticleForcesSweep< ParticleAccessor_T, ParticleSelector_T, Weighting_T >(
           bs, ac, ps, particleAndVolumeFractionSoA))
   {}
   SphereFractionMappingSweep< ParticleAccessor_T, ParticleSelector_T, Weighting_T > particleMappingSweep;
   SetParticleVelocitiesSweep< ParticleAccessor_T, ParticleSelector_T, Weighting_T > setParticleVelocitiesSweep;
   ReduceParticleForcesSweep< ParticleAccessor_T, ParticleSelector_T, Weighting_T > reduceParticleForcesSweep;
};

template< typename SweepCollection, typename PSMSweep >
void addPSMSweepsToTimeloop(SweepTimeloop& timeloop, SweepCollection& psmSweepCollection, PSMSweep& psmSweep,
                            bool synchronize = true)
{
   if (synchronize)
   {
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.particleMappingSweep), "Particle mapping");
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.setParticleVelocitiesSweep),
                              "Set particle velocities");
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweep), "PSM sweep");
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.reduceParticleForcesSweep),
                              "Reduce particle forces");
   }
   else
   {
      timeloop.add() << Sweep(psmSweepCollection.particleMappingSweep, "Particle mapping");
      timeloop.add() << Sweep(psmSweepCollection.setParticleVelocitiesSweep, "Set particle velocities");
      timeloop.add() << Sweep(psmSweep, "PSM sweep");
      timeloop.add() << Sweep(psmSweepCollection.reduceParticleForcesSweep, "Reduce particle forces");
   };
}

template< typename SweepCollection, typename PSMSweep, typename Communication >
void addPSMSweepsToTimeloops(SweepTimeloop& commTimeloop, SweepTimeloop& timeloop, Communication& comm,
                             SweepCollection& psmSweepCollection, PSMSweep& psmSweep, bool synchronize = true)
{
   if (synchronize)
   {
      commTimeloop.add() << BeforeFunction([&]() { comm.startCommunication(); })
                         << Sweep(deviceSyncWrapper(psmSweepCollection.particleMappingSweep), "Particle mapping");
      commTimeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.setParticleVelocitiesSweep),
                                  "Set particle velocities");
      commTimeloop.add() << Sweep(deviceSyncWrapper(psmSweep.getInnerSweep()), "PSM inner sweep")
                         << AfterFunction([&]() { comm.wait(); }, "LBM Communication (wait)");
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweep.getOuterSweep()), "PSM outer sweep");
      timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.reduceParticleForcesSweep),
                              "Reduce particle forces");
   }
   else
   {
      commTimeloop.add() << BeforeFunction([&]() { comm.startCommunication(); })
                         << Sweep(psmSweepCollection.particleMappingSweep, "Particle mapping");
      commTimeloop.add() << Sweep(psmSweepCollection.setParticleVelocitiesSweep, "Set particle velocities");
      commTimeloop.add() << Sweep(psmSweep.getInnerSweep(), "PSM inner sweep")
                         << AfterFunction([&]() { comm.wait(); }, "LBM Communication (wait)");
      timeloop.add() << Sweep(psmSweep.getOuterSweep(), "PSM outer sweep");
      timeloop.add() << Sweep(psmSweepCollection.reduceParticleForcesSweep, "Reduce particle forces");
   };
}

} // namespace gpu
} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
