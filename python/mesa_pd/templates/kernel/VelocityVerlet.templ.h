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
//! \file VelocityVerlet.h
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
#include <mesa_pd/data/IAccessor.h>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Velocity verlet integration for all particles.
 *
 * Velocit verlet integration is a two part kernel. preForceUpdate has to be
 * called before the force calculation and postFroceUpdate afterwards. The
 * integration is only complete when both functions are called. The integration
 * is symplectic.
 *
 * Wachs, A. Particle-scale computational approaches to model dry and saturated granular flows of non-Brownian, non-cohesive, and non-spherical rigid bodies. Acta Mech 230, 1919–1980 (2019). https://doi.org/10.1007/s00707-019-02389-9
 *
 * This kernel requires the following particle accessor interface
 * \code
   {%- for prop in interface %}
   {%- if 'g' in prop.access %}
 * const {{prop.type}}& get{{prop.name | capFirst}}(const size_t p_idx) const;
   {%- endif %}
   {%- if 's' in prop.access %}
 * void set{{prop.name | capFirst}}(const size_t p_idx, const {{prop.type}}& v);
   {%- endif %}
   {%- if 'r' in prop.access %}
 * {{prop.type}}& get{{prop.name | capFirst}}Ref(const size_t p_idx);
   {%- endif %}
 *
   {%- endfor %}
 * \endcode
 * \ingroup mesa_pd_kernel
 */
class VelocityVerletPreForceUpdate
{
public:
   VelocityVerletPreForceUpdate(const real_t dt) : dt_(dt) {}

   template <typename Accessor>
   void operator()(const size_t i, Accessor& ac) const;

   real_t dt_;
};

/// \see VelocityVerletPreForceUpdate
class VelocityVerletPostForceUpdate
{
public:
   VelocityVerletPostForceUpdate(const real_t dt) : dt_(dt) {}

   template <typename Accessor>
   void operator()(const size_t i, Accessor& ac) const;

   real_t dt_;
};

template <typename Accessor>
inline void VelocityVerletPreForceUpdate::operator()(const size_t p_idx, Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (!data::particle_flags::isSet( ac.getFlags(p_idx), data::particle_flags::FIXED))
   {
      ac.setPosition(p_idx, ac.getPosition(p_idx) +
                            ac.getLinearVelocity(p_idx) * dt_ +
                            real_t(0.5) * ac.getInvMass(p_idx) * ac.getOldForce(p_idx) * dt_ * dt_);

      {%- if bIntegrateRotation %}
      const Vec3 wdot = math::transformMatrixRART(ac.getRotation(p_idx).getMatrix(),
                                                  ac.getInvInertiaBF(p_idx)) * ac.getOldTorque(p_idx);

      // Calculating the rotation angle
      const Vec3 phi( ac.getAngularVelocity(p_idx) * dt_ + real_t(0.5) * wdot * dt_ * dt_);

      // Calculating the new orientation
      auto rotation = ac.getRotation(p_idx);
      rotation.rotate( phi );
      ac.setRotation(p_idx, rotation);
      {%- endif %}
   }
}

template <typename Accessor>
inline void VelocityVerletPostForceUpdate::operator()(const size_t p_idx, Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (!data::particle_flags::isSet( ac.getFlags(p_idx), data::particle_flags::FIXED))
   {
      ac.setLinearVelocity(p_idx, ac.getLinearVelocity(p_idx) +
                                  real_t(0.5) * ac.getInvMass(p_idx) * (ac.getOldForce(p_idx) + ac.getForce(p_idx)) * dt_);

      {%- if bIntegrateRotation %}
      const auto torque = ac.getOldTorque(p_idx) + ac.getTorque(p_idx);
      const Vec3 wdot = math::transformMatrixRART(ac.getRotation(p_idx).getMatrix(),
                                                  ac.getInvInertiaBF(p_idx)) * torque;

      ac.setAngularVelocity(p_idx, ac.getAngularVelocity(p_idx) +
                                   real_t(0.5) * wdot * dt_ );
      {%- endif %}
   }

   ac.setOldForce(p_idx,       ac.getForce(p_idx));
   ac.setForce(p_idx,          Vec3(real_t(0), real_t(0), real_t(0)));

   {%- if bIntegrateRotation %}
   ac.setOldTorque(p_idx,      ac.getTorque(p_idx));
   ac.setTorque(p_idx,         Vec3(real_t(0), real_t(0), real_t(0)));
   {%- endif %}
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
