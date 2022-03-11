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
//! \file
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
 * Transformations between world frame (WF) and body frame (BF) coordinates
 */
template <typename Accessor>
inline Vec3 transformPositionFromWFtoBF(const size_t p_idx, Accessor& ac, const Vec3& positionWF)
{
   return ac.getRotation(p_idx).getMatrix().getTranspose() * ( positionWF - ac.getPosition(p_idx)  );
}

template <typename Accessor>
inline Vec3 transformVectorFromWFtoBF(const size_t p_idx, Accessor& ac, const Vec3& vectorWF)
{
   return ac.getRotation(p_idx).getMatrix().getTranspose() * vectorWF;
}

template <typename Accessor>
inline Vec3 transformPositionFromBFtoWF(const size_t p_idx, Accessor& ac, const Vec3& positionBF)
{
   return ac.getPosition(p_idx) + ac.getRotation(p_idx).getMatrix() * positionBF;
}

template <typename Accessor>
inline Vec3 transformVectorFromBFtoWF(const size_t p_idx, Accessor& ac, const Vec3& vectorBF)
{
   return ac.getRotation(p_idx).getMatrix() * vectorBF;
}

/**
 * Transform (inverse) particle's moment of inertia from body frame coordinates (as stored by shape) to world frame.
 */
template <typename Accessor>
inline Mat3 getInvInertia(const size_t p_idx, Accessor& ac)
{
   return math::transformMatrixRART(ac.getRotation(p_idx).getMatrix(), ac.getInvInertiaBF(p_idx));
}

template <typename Accessor>
inline Mat3 getInertia(const size_t p_idx, Accessor& ac)
{
   return math::transformMatrixRART(ac.getRotation(p_idx).getMatrix(), ac.getInertiaBF(p_idx));
}

/**
 * Force is applied at the center of mass.
 */
template <typename Accessor>
inline void addForceAtomic(const size_t p_idx, Accessor& ac, const Vec3& f)
{
   // Increasing the force and torque on this particle;
   ac.getForceRef(p_idx)[0]  += f[0];;
   ac.getForceRef(p_idx)[1]  += f[1];;
   ac.getForceRef(p_idx)[2]  += f[2];
}

template <typename Accessor>
inline void addForceAtWFPosAtomic(const size_t p_idx, Accessor& ac, const Vec3& f, const Vec3& wf_pt)
{
   // Increasing the force and torque on this particle;
   ac.getForceRef(p_idx)[0]  += f[0];;
   ac.getForceRef(p_idx)[1]  += f[1];;
   ac.getForceRef(p_idx)[2]  += f[2];

   const auto t = cross(( wf_pt - ac.getPosition(p_idx) ), f);;
   ac.getTorqueRef(p_idx)[0] += t[0];;
   ac.getTorqueRef(p_idx)[1] += t[1];;
   ac.getTorqueRef(p_idx)[2] += t[2];
}

/**
 * Torque is directly applied on the particle.
 */
template <typename Accessor>
inline void addTorqueAtomic(const size_t p_idx, Accessor& ac, const Vec3& t)
{
   // Increasing the torque on this particle;
   ac.getTorqueRef(p_idx)[0]  += t[0];;
   ac.getTorqueRef(p_idx)[1]  += t[1];;
   ac.getTorqueRef(p_idx)[2]  += t[2];
}


} //namespace mesa_pd
} //namespace walberla