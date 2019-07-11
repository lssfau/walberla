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
//! \file NewGhostParticleNotification.h
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

#include <core/mpi/Datatype.h>
#include <core/mpi/MPIWrapper.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace mesa_pd {

/**
 * This notification is send to the owner of a particle
 * to signal that a new ghost particle exists and the ghost particle list should be updated.
 */
class NewGhostParticleNotification
{
public:
   struct Parameters
   {
      id_t    uid_;
      walberla::mpi::MPIRank newOwner_;
   };

   inline explicit NewGhostParticleNotification( const data::Particle& particle, const walberla::mpi::MPIRank newOwner )
      : particle_(particle)
      , newOwner_(newOwner)
   {}
   const data::Particle& particle_;
   walberla::mpi::MPIRank newOwner_;
};

template<>
struct NotificationTrait<NewGhostParticleNotification>
{
   static const NotificationType id = NEW_GHOST_PARTICLE_NOTIFICATION;
};

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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::NewGhostParticleNotification& obj )
{
   buf.addDebugMarker( "cn" );
   buf << obj.particle_.getUid();
   buf << obj.newOwner_;
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::NewGhostParticleNotification::Parameters& objparam )
{
   buf.readDebugMarker( "cn" );
   buf >> objparam.uid_;
   buf >> objparam.newOwner_;
   return buf;
}

} // mpi
} // walberla