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
 * The specified particle should be removed from the process.
 */
class ParticleRemovalNotification {
public:
   struct Parameters {
      id_t uid_;
   };

   inline explicit ParticleRemovalNotification( const data::Particle& particle ) : particle_(particle) {}
   const data::Particle& particle_;
};

template<>
struct NotificationTrait<ParticleRemovalNotification>
{
   static const NotificationType id = PARTICLE_REMOVAL_NOTIFICATION;
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ParticleRemovalNotification& obj )
{
   buf.addDebugMarker( "un" );
   buf << obj.particle_.getUid();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ParticleRemovalNotification::Parameters& objparam )
{
   buf.readDebugMarker( "un" );
   buf >> objparam.uid_;
   return buf;
}

template<>
struct BufferSizeTrait< mesa_pd::ParticleRemovalNotification > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<id_t>::size + mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla
