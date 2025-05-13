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
//! \author Tobias Leemann <tobias.leemann@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/data/ContactStorage.h>
#include <mesa_pd/data/Flags.h>
#include <core/logging/Logging.h>
#include <core/math/Constants.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Kernel performing integration of the particles after the HCSITS iteration is done.
 * Call this kernel on all particles j to integrate them by a timestep of size dt.
 * The Speed Limiter limits the number of body radii, a particle can travel in one timestep.
 * The speed limit factor defines a number of radii that are allowed in a timestep, e.g. a factor of 1 means, that
 * a particle can only travel one time its radius in each timestep.
 *
 * \ingroup mesa_pd_kernel
 */
class IntegrateParticlesHCSITS
{
   public:
   // Default constructor sets the default values
   IntegrateParticlesHCSITS() :
   {%- for prop in properties %}
   {{prop.name}}_( {{prop.defValue}} ){{ "," if not loop.last}}
   {%- endfor %}
   {}

   // Getter and Setter Functions
   {%- for prop in properties %}
   inline const {{prop.type}}& get{{prop.name | capFirst}}() const {return {{prop.name}}_;}
   inline void set{{prop.name | capFirst}}({{prop.type}} v) { {{prop.name}}_ = v;}
   {% endfor %}

   template <typename PAccessor>
   void operator()(size_t j, PAccessor& ac, real_t dt);

   template <typename PAccessor>
   void integratePositions(PAccessor& ac, size_t body, Vec3 v, Vec3 w, real_t dt ) const;


   private:
   // List properties
   {% for prop in properties %}
   {{prop.type}} {{prop.name }}_;
   {%- endfor %}
};


template <typename PAccessor>
inline void IntegrateParticlesHCSITS::operator()(size_t j, PAccessor& ac, real_t dt)
{
   using namespace data::particle_flags;
   static_assert(std::is_base_of_v<data::IAccessor, PAccessor>, "please provide a valid accessor");

   auto body_flags = ac.getFlagsRef(j);
   if (isSet(body_flags, FIXED)){
      integratePositions(ac, j, ac.getLinearVelocity(j), ac.getAngularVelocity(j), dt );
   }else{
      integratePositions(ac, j, ac.getLinearVelocity(j) + ac.getDv(j), ac.getAngularVelocity(j) + ac.getDw(j), dt );
   }
}




//*************************************************************************************************
/*!\brief Time integration of the position and orientation of a given body.
 *
 * \param body The body whose position and orientation to time integrate
 * \param v The linear velocity to use for time integration of the position.
 * \param w The angular velocity to use for time integration of the orientation.
 * \param dt The time step size.
 * \return void
 *
 * Performs an Euler time integration of the positions of the given body. Velocities are damped if
 * indicated by the settings and stored back in the body properties. The bounding box is
 * recalculated and it is redetermined whether the body is awake or not. Also the data
 * structure tracking the contacts attached to the body are cleared and
 */
template <typename PAccessor>
inline void IntegrateParticlesHCSITS::integratePositions(PAccessor& ac, size_t body, Vec3 v, Vec3 w, real_t dt ) const {
   using namespace data::particle_flags;
   auto body_flags = ac.getFlags(body);
   if (isSet(body_flags, FIXED))
      return;


   if (getSpeedLimiterActive()) {
      const auto speed = v.length();
      const auto edge = ac.getInteractionRadius(body);
      if (speed * dt > edge * getSpeedLimitFactor()) {
         v = v * (edge * getSpeedLimitFactor() / dt / speed);
      }

      const real_t maxPhi = real_t(2) * math::pi * getSpeedLimitFactor();
      const real_t phi = w.length() * dt;
      if (phi > maxPhi) {
         w = w / phi * maxPhi;
      }
   }

   // Calculating the translational displacement
   ac.getPositionRef(body) = (ac.getPosition(body) + v * dt);

   // Calculating the rotation angle
   const Vec3 phi(w * dt);

   // Calculating the new orientation
   WALBERLA_ASSERT(!math::isnan(phi), " phi: " << phi << " w: " << w << " dt: " << dt );
   ac.getRotationRef(body).rotate(phi, phi.length());

   // Storing the velocities back in the body properties
   ac.getLinearVelocityRef(body) = v;
   ac.getAngularVelocityRef(body) = w;
}
//*************************************************************************************************



} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
