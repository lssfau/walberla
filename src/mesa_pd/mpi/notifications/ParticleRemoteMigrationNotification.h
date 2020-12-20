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

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/mpi/notifications/NotificationType.h>

#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace mesa_pd {

/**
 * The ownership for one of the ghost particles has changed.
 */
class ParticleRemoteMigrationNotification {
public:
   struct Parameters {
      id_t uid_;
      int newOwner_;
   };

   inline ParticleRemoteMigrationNotification( const data::Particle& particle, const int& newOwner )
      : particle_(particle), newOwner_(newOwner) {}
   const data::Particle& particle_;
   const int newOwner_;
};

template<>
struct NotificationTrait<ParticleRemoteMigrationNotification>
{
   static const NotificationType id = PARTICLE_REMOTE_MIGRATION_NOTIFICATION;
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ParticleRemoteMigrationNotification& obj )
{
   buf.addDebugMarker( "rm" );
   buf << obj.particle_.getUid();
   buf << obj.newOwner_;
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ParticleRemoteMigrationNotification::Parameters& objparam )
{
   buf.readDebugMarker( "rm" );
   buf >> objparam.uid_;
   buf >> objparam.newOwner_;
   return buf;
}

template <>
struct BufferSizeTrait< mesa_pd::ParticleRemoteMigrationNotification > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<id_t>::size + BufferSizeTrait<int>::size + mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla