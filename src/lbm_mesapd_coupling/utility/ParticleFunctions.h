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
//! \file ParticleFunctions.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/AABB.h"
#include "core/math/Limits.h"

#include "mesa_pd/data/Flags.h"
#include "mesa_pd/data/IAccessor.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/kernel/SingleCast.h"

namespace walberla {
namespace lbm_mesapd_coupling {


/**
 * Force is applied at the center of mass.
 */
template <typename ParticleAccessor_T>
inline void addHydrodynamicForceAtomic(const size_t p_idx, ParticleAccessor_T& ac, const Vector3<real_t>& f)
{
   // Increasing the force and torque on this particle
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicForceRef(p_idx)[0]  += f[0];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicForceRef(p_idx)[1]  += f[1];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicForceRef(p_idx)[2]  += f[2];
}

template <typename ParticleAccessor_T>
inline void addHydrodynamicForceAtWFPosAtomic(const size_t p_idx, ParticleAccessor_T& ac, const Vector3<real_t>& f, const Vector3<real_t>& wf_pt)
{
   // Increasing the force and torque on this particle
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicForceRef(p_idx)[0]  += f[0];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicForceRef(p_idx)[1]  += f[1];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicForceRef(p_idx)[2]  += f[2];

   const auto t = cross(( wf_pt - ac.getPosition(p_idx) ), f);

#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicTorqueRef(p_idx)[0] += t[0];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicTorqueRef(p_idx)[1] += t[1];
#ifdef _OPENMP
#pragma omp atomic
#endif
   ac.getHydrodynamicTorqueRef(p_idx)[2] += t[2];
}


} // namespace lbm_mesapd_coupling
} // namespace walberla
