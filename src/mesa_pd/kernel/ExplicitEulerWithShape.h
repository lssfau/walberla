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
//! \file ExplicitEuler.h
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
 * Kernel which explicitly integrates all particles in time.
 *
 * This kernel requires the following particle accessor interface
 * \code
 * const walberla::mesa_pd::Vec3& getPosition(const size_t p_idx) const;
 * void setPosition(const size_t p_idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::Vec3& getLinearVelocity(const size_t p_idx) const;
 * void setLinearVelocity(const size_t p_idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::real_t& getInvMass(const size_t p_idx) const;
 *
 * const walberla::mesa_pd::Vec3& getForce(const size_t p_idx) const;
 * void setForce(const size_t p_idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::data::particle_flags::FlagT& getFlags(const size_t p_idx) const;
 *
 * const walberla::mesa_pd::Rot3& getRotation(const size_t p_idx) const;
 * void setRotation(const size_t p_idx, const walberla::mesa_pd::Rot3& v);
 *
 * const walberla::mesa_pd::Vec3& getAngularVelocity(const size_t p_idx) const;
 * void setAngularVelocity(const size_t p_idx, const walberla::mesa_pd::Vec3& v);
 *
 * const walberla::mesa_pd::Mat3& getInvInertiaBF(const size_t p_idx) const;
 *
 * const walberla::mesa_pd::Vec3& getTorque(const size_t p_idx) const;
 * void setTorque(const size_t p_idx, const walberla::mesa_pd::Vec3& v);
 *
 * \endcode
 *
 * \pre  All forces and torques acting on the particles have to be set.
 * \post All forces and torques are reset to 0.
 * \ingroup mesa_pd_kernel
 */
class ExplicitEuler
{
public:
   explicit ExplicitEuler(const real_t dt) : dt_(dt) {}

   template <typename Accessor>
   void operator()(const size_t i, Accessor& ac) const;
private:
   real_t dt_ = real_t(0.0);
};

template <typename Accessor>
inline void ExplicitEuler::operator()(const size_t idx,
                                               Accessor& ac) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (!data::particle_flags::isSet( ac.getFlags(idx), data::particle_flags::FIXED))
   {
      ac.setPosition      (idx, ac.getInvMass(idx) * ac.getForce(idx) * dt_ * dt_ + ac.getLinearVelocity(idx) * dt_ + ac.getPosition(idx));
      ac.setLinearVelocity(idx, ac.getInvMass(idx) * ac.getForce(idx) * dt_ + ac.getLinearVelocity(idx));
      const Vec3 wdot = math::transformMatrixRART(ac.getRotation(idx).getMatrix(),
                                                  ac.getInvInertiaBF(idx)) * ac.getTorque(idx);

      // Calculating the rotation angle
      const Vec3 phi( ac.getAngularVelocity(idx) * dt_ + wdot * dt_ * dt_);

      // Calculating the new orientation
      auto rotation = ac.getRotation(idx);
      rotation.rotate( phi );
      ac.setRotation(idx, rotation);

      ac.setAngularVelocity(idx, wdot * dt_ + ac.getAngularVelocity(idx));
   }

   ac.setForce (idx, Vec3(real_t(0), real_t(0), real_t(0)));
   ac.setTorque(idx, Vec3(real_t(0), real_t(0), real_t(0)));
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla