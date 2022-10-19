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
#include <mesa_pd/data/IAccessor.h>
{%- if bIntegrateRotation %}
#include <mesa_pd/common/ParticleFunctions.h>
{%- endif %}

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Velocity verlet integration for all particles.
 *
 * Velocity verlet integration is a two part kernel. preForceUpdate has to be
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
 * const {{prop.type}}& get{{prop.name | capFirst}}(const size_t idx) const;
   {%- endif %}
   {%- if 's' in prop.access %}
 * void set{{prop.name | capFirst}}(const size_t idx, const {{prop.type}}& v);
   {%- endif %}
   {%- if 'r' in prop.access %}
 * {{prop.type}}& get{{prop.name | capFirst}}Ref(const size_t idx);
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
inline void VelocityVerletPreForceUpdate::operator()(const size_t idx, Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (!data::particle_flags::isSet( ac.getFlags(idx), data::particle_flags::FIXED))
   {
      ac.setPosition(idx, ac.getPosition(idx) +
                            ac.getLinearVelocity(idx) * dt_ +
                            real_t(0.5) * ac.getInvMass(idx) * ac.getOldForce(idx) * dt_ * dt_);

      {%- if bIntegrateRotation %}
      {%- if bUseFullAngularMomentumEquation %}
      // computation done in body frame: d(omega)/ dt = J^-1 ((J*omega) x omega + T), update in world frame
      // see Wachs, 2019, doi:10.1007/s00707-019-02389-9, Eq. 27
      // note: this implementation (pre and post) is experimental as it is in principle unclear in which order
      //       angular velocities and rotations (affect again the transformations WF - BF) have to be carried out
      const auto omegaBF = transformVectorFromWFtoBF(idx, ac, ac.getAngularVelocity(idx));
      const auto torqueBF = transformVectorFromWFtoBF(idx, ac, ac.getOldTorque(idx));
      const Vec3 wdotBF = ac.getInvInertiaBF(idx) * ( ( ac.getInertiaBF(idx) * omegaBF ) % omegaBF + torqueBF );
      const Vec3 wdot = transformVectorFromBFtoWF(idx, ac, wdotBF);
      {%- else %}
      // note: contribution (J*omega) x omega is ignored here -> see template for other variant
      const Vec3 wdot = math::transformMatrixRART(ac.getRotation(idx).getMatrix(),
                                                  ac.getInvInertiaBF(idx)) * ac.getOldTorque(idx);
      {%- endif %}


      // Calculating the rotation angle
      const Vec3 phi( ac.getAngularVelocity(idx) * dt_ + real_t(0.5) * wdot * dt_ * dt_);

      // Calculating the new orientation
      auto rotation = ac.getRotation(idx);
      rotation.rotate( phi );
      ac.setRotation(idx, rotation);
      {%- endif %}
   }
}

template <typename Accessor>
inline void VelocityVerletPostForceUpdate::operator()(const size_t idx, Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (!data::particle_flags::isSet( ac.getFlags(idx), data::particle_flags::FIXED))
   {
      ac.setLinearVelocity(idx, ac.getLinearVelocity(idx) +
                                  real_t(0.5) * ac.getInvMass(idx) * (ac.getOldForce(idx) + ac.getForce(idx)) * dt_);

      {%- if bIntegrateRotation %}
      {%- if bUseFullAngularMomentumEquation %}
      const auto omegaBF = transformVectorFromWFtoBF(idx, ac, ac.getAngularVelocity(idx));
      const auto torqueBF = transformVectorFromWFtoBF(idx, ac, 0.5_r * (ac.getOldTorque(idx) + ac.getTorque(idx)));
      const Vec3 wdotBF = ac.getInvInertiaBF(idx) * ( ( ac.getInertiaBF(idx) * omegaBF ) % omegaBF + torqueBF );
      const Vec3 wdot = transformVectorFromBFtoWF(idx, ac, wdotBF);
      {%- else %}
      const auto torque = 0.5_r * (ac.getOldTorque(idx) + ac.getTorque(idx));
      const Vec3 wdot = math::transformMatrixRART(ac.getRotation(idx).getMatrix(),
                                                  ac.getInvInertiaBF(idx)) * torque;
      {%- endif %}

      ac.setAngularVelocity(idx, ac.getAngularVelocity(idx) +
                                   wdot * dt_ );
      {%- endif %}
   }

   ac.setOldForce(idx,       ac.getForce(idx));
   ac.setForce(idx,          Vec3(real_t(0), real_t(0), real_t(0)));

   {%- if bIntegrateRotation %}
   ac.setOldTorque(idx,      ac.getTorque(idx));
   ac.setTorque(idx,         Vec3(real_t(0), real_t(0), real_t(0)));
   {%- endif %}
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
