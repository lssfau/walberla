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
#include <mesa_pd/data/Flags.h>
#include <core/logging/Logging.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {
/**
 * Init the datastructures for the particles for later use of the HCSITS-Solver.
 * Call this kernel on all particles that will be treated with HCSITS before performing any relaxation timesteps.
 * Use setGlobalAcceleration() to set an acceleration action uniformly across all particles (e.g. gravity)
 * \ingroup mesa_pd_kernel
 */
class InitParticlesForHCSITS
{
public:
   // Default constructor sets the default values
   InitParticlesForHCSITS() :
   {%- for prop in properties %}
   {{prop.name}}_( {{prop.defValue}} ){{ "," if not loop.last}}
   {%- endfor %}
   {}

   // Getter and Setter Functions for properties
   {%- for prop in properties %}
   inline const {{prop.type}}& get{{prop.name | capFirst}}() const {return {{prop.name}}_;}
   inline void set{{prop.name | capFirst}}({{prop.type}} v) { {{prop.name}}_ = v;}
   {% endfor %}

   template <typename Accessor>
   void operator()(size_t j, Accessor& ac, real_t dt);

   template <typename Accessor>
   void initializeVelocityCorrections(Accessor& ac, size_t body, Vec3& dv, Vec3& dw, real_t dt ) const;

private:
   // List properties
   {% for prop in properties %}
   {{prop.type}} {{prop.name }}_;
   {%- endfor %}

};


template <typename Accessor>
inline void InitParticlesForHCSITS::operator()(size_t j, Accessor& ac, real_t dt)
{
   using namespace data::particle_flags;
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");
   auto particle_flags = ac.getFlagsRef(j);
   if (isSet(particle_flags, GLOBAL)){
      WALBERLA_CHECK( isSet(particle_flags, FIXED), "Global bodies need to have infinite mass as they are not communicated!" );
      initializeVelocityCorrections( ac, j, ac.getDvRef(j), ac.getDwRef(j), dt ); // use applied external forces to calculate starting velocity
   }else{
      initializeVelocityCorrections( ac, j, ac.getDvRef(j), ac.getDwRef(j), dt ); // use applied external forces to calculate starting velocity
      if(!isSet(particle_flags, FIXED)){ // Update velocities with global acceleration and angular velocity with euler eqn
         ac.getLinearVelocityRef(j) = ac.getLinearVelocity(j) + getGlobalAcceleration() * dt;
         const auto omegaBF = transformVectorFromWFtoBF(j, ac, ac.getAngularVelocity(j));
         ac.getAngularVelocityRef(j) = ac.getAngularVelocity(j) + dt * transformVectorFromBFtoWF(j, ac, ( ac.getInvInertiaBF(j) * ( ( ac.getInertiaBF(j) * omegaBF ) % omegaBF ) ) );
      }
   }
}


//*************************************************************************************************
/*!\brief Calculates the initial velocity corrections of a given body.
 * \param ac The particle accessor
 * \param body The body whose velocities to time integrate
 * \param dv On return the initial linear velocity correction.
 * \param w On return the initial angular velocity correction.
 * \param dt The time step size.
 * \return void
 *
 * Calculates the velocity corrections effected by external forces and torques in an explicit Euler
 * time integration of the velocities of the given body. For fixed objects the velocity corrections
 * are set to zero. External forces and torques are reset if indicated by the settings.
 */
template <typename Accessor>
inline void InitParticlesForHCSITS::initializeVelocityCorrections(Accessor& ac, size_t body, Vec3& dv, Vec3& dw, real_t dt ) const
{
   dv = ( ac.getInvMass(body) * dt ) * ac.getForce(body);
   dw = dt * ( getInvInertia(body, ac) * ac.getTorque(body) );

   ac.getForceRef(body) = Vec3();
   ac.getTorqueRef(body) = Vec3();
}
//*************************************************************************************************


} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
