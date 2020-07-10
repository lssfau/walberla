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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/domain/IDomain.h>

#include <core/logging/Logging.h>
#include <core/mpi/MPIManager.h>

namespace walberla {
namespace mesa_pd {

/**
 * Returns the "surface" velocity at a certain point given in world frame coordinates.
 */
template <typename Accessor>
inline Vec3 getVelocityAtWFPoint(const size_t p_idx, Accessor& ac, const Vec3& wf_pt)
{
   return ac.getLinearVelocity(p_idx) + cross(ac.getAngularVelocity(p_idx), ( wf_pt - ac.getPosition(p_idx) ));
}

/**
 * Force is applied at the center of mass.
 */
template <typename Accessor>
inline void addForceAtomic(const size_t p_idx, Accessor& ac, const Vec3& f)
{
   // Increasing the force and torque on this particle
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getForceRef(p_idx)[0]  += f[0];
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getForceRef(p_idx)[1]  += f[1];
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getForceRef(p_idx)[2]  += f[2];
}

template <typename Accessor>
inline void addForceAtWFPosAtomic(const size_t p_idx, Accessor& ac, const Vec3& f, const Vec3& wf_pt)
{
   // Increasing the force and torque on this particle
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getForceRef(p_idx)[0]  += f[0];
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getForceRef(p_idx)[1]  += f[1];
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getForceRef(p_idx)[2]  += f[2];

   const auto t = cross(( wf_pt - ac.getPosition(p_idx) ), f);

   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getTorqueRef(p_idx)[0] += t[0];
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getTorqueRef(p_idx)[1] += t[1];
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getTorqueRef(p_idx)[2] += t[2];
}

/**
 * Torque is directly applied on the particle.
 */
template <typename Accessor>
inline void addTorqueAtomic(const size_t p_idx, Accessor& ac, const Vec3& t)
{
   // Increasing the torque on this particle
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getTorqueRef(p_idx)[0]  += t[0];
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getTorqueRef(p_idx)[1]  += t[1];
   {%- if module.enableOpenMP %}
   #ifdef _OPENMP
   #pragma omp atomic
   #endif
   {%- endif %};
   ac.getTorqueRef(p_idx)[2]  += t[2];
}


} //namespace mesa_pd
} //namespace walberla
