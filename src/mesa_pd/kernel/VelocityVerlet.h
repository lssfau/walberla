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
#include <mesa_pd/common/ParticleFunctions.h>

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
 * const walberla::mesa_pd::Vec3& getPosition(const size_t idx) const;
 * void setPosition(const size_t idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::Vec3& getLinearVelocity(const size_t idx) const;
 * void setLinearVelocity(const size_t idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::real_t& getInvMass(const size_t idx) const;
 *
 * const walberla::mesa_pd::Vec3& getForce(const size_t idx) const;
 * void setForce(const size_t idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::Vec3& getOldForce(const size_t idx) const;
 * void setOldForce(const size_t idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::Rot3& getRotation(const size_t idx) const;
 * void setRotation(const size_t idx, const walberla::mesa_pd::Rot3& v);
 *
 * const walberla::mesa_pd::Vec3& getAngularVelocity(const size_t idx) const;
 * void setAngularVelocity(const size_t idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::Mat3& getInvInertiaBF(const size_t idx) const;
 *
 * const walberla::mesa_pd::Mat3& getInertiaBF(const size_t idx) const;
 *
 * const walberla::mesa_pd::Vec3& getTorque(const size_t idx) const;
 * void setTorque(const size_t idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::Vec3& getOldTorque(const size_t idx) const;
 * void setOldTorque(const size_t idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::data::particle_flags::FlagT& getFlags(const size_t idx) const;
 *
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
      // note: contribution (J*omega) x omega is ignored here -> see template for other variant
      const Vec3 wdot = math::transformMatrixRART(ac.getRotation(idx).getMatrix(),
                                                  ac.getInvInertiaBF(idx)) * ac.getOldTorque(idx);


      // Calculating the rotation angle
      const Vec3 phi( ac.getAngularVelocity(idx) * dt_ + real_t(0.5) * wdot * dt_ * dt_);

      // Calculating the new orientation
      auto rotation = ac.getRotation(idx);
      rotation.rotate( phi );
      ac.setRotation(idx, rotation);
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
      const auto torque = 0.5_r * (ac.getOldTorque(idx) + ac.getTorque(idx));
      const Vec3 wdot = math::transformMatrixRART(ac.getRotation(idx).getMatrix(),
                                                  ac.getInvInertiaBF(idx)) * torque;

      ac.setAngularVelocity(idx, ac.getAngularVelocity(idx) +
                                   wdot * dt_ );
   }

   ac.setOldForce(idx,       ac.getForce(idx));
   ac.setForce(idx,          Vec3(real_t(0), real_t(0), real_t(0)));
   ac.setOldTorque(idx,      ac.getTorque(idx));
   ac.setTorque(idx,         Vec3(real_t(0), real_t(0), real_t(0)));
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla