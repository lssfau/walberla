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
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/mpi/notifications/NotificationType.h>
#include <mesa_pd/mpi/notifications/reset.h>
#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace mesa_pd {

/**
 * Adds up LinearVelocity + Parameters::relaxationParam * velocity_correction in a particle and transmits the result to
 * all its ghost particles.
 * After adding up and sending the result, the velocity corrections of the main particle are reset to 0.
 * After receiving new values, the velocity_corrections of the ghost_particles are reset to 0.
 * This notification is used during relaxation with the HCSITS solver.
 * Notification for use with BroadcastProperty only.
 *
 */
class VelocityUpdateNotification
{
public:

   struct Parameters
   {
      static real_t relaxationParam;
      id_t uid_;
      Vec3 v_; /* Linear velocity */
      Vec3 w_; /* Angular velocity */
   };

   inline explicit VelocityUpdateNotification( data::Particle& p ) : p_(p)  {}

   data::Particle& p_;
};

real_t VelocityUpdateNotification::Parameters::relaxationParam = real_t(0.8);

// Update method for broadcast
void update(data::Particle&& p, const VelocityUpdateNotification::Parameters& objparam) {
   // Reset the velocity corrections dv/dw of ghost particle
   p.getDvRef() = Vec3();
   p.getDwRef() = Vec3();
   p.getLinearVelocityRef() = objparam.v_;
   p.getAngularVelocityRef() = objparam.w_;
}

template<>
void reset<VelocityUpdateNotification>(data::Particle& p )
{
   p.setDv( Vec3(real_t(0)) );
   p.setDw( Vec3(real_t(0)) );
   p.setLinearVelocity( Vec3(real_t(0)) );
   p.setAngularVelocity( Vec3(real_t(0)) );
}

}  // namespace mesa_pd
}  // namespace walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
          typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::VelocityUpdateNotification& obj )
{
   // Perform the addition
   obj.p_.getLinearVelocityRef() = obj.p_.getLinearVelocity() + mesa_pd::VelocityUpdateNotification::Parameters::relaxationParam * obj.p_.getDv();
   obj.p_.getAngularVelocityRef() = obj.p_.getAngularVelocity() + mesa_pd::VelocityUpdateNotification::Parameters::relaxationParam * obj.p_.getDw();
   // Reset the velocity corrections dv/dw of main particle
   obj.p_.getDvRef() = mesa_pd::Vec3();
   obj.p_.getDwRef() = mesa_pd::Vec3();
   buf.addDebugMarker( "ft" );
   buf << obj.p_.getUid();
   buf << obj.p_.getLinearVelocity();
   buf << obj.p_.getAngularVelocity();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::VelocityUpdateNotification::Parameters& objparam )
{
   buf.readDebugMarker( "ft" );
   buf >> objparam.uid_;
   buf >> objparam.v_;
   buf >> objparam.w_;
   return buf;
}

template< >
struct BufferSizeTrait< mesa_pd::VelocityUpdateNotification > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<id_t>::size +
                              BufferSizeTrait<mesa_pd::Vec3>::size +
                              BufferSizeTrait<mesa_pd::Vec3>::size +
                              mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla
