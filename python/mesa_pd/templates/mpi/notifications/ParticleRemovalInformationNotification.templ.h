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

#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace mesa_pd {

/**
 * The ParticleRemovalInformationNotification class is used to signal other processes that a
 * shadow copy was destroyed.
 */
class ParticleRemovalInformationNotification
{
public:
   struct Parameters
   {
      id_t    uid_;
      walberla::mpi::MPIRank owner_;
   };

   inline explicit ParticleRemovalInformationNotification( const data::Particle& particle )
      : particle_(particle)
   {}
   const data::Particle& particle_;
};

template<>
struct NotificationTrait<ParticleRemovalInformationNotification>
{
   static const NotificationType id = PARTICLE_REMOVAL_INFORMATION_NOTIFICATION;
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ParticleRemovalInformationNotification& obj )
{
   buf.addDebugMarker( "ri" );
   buf << obj.particle_.getUid();
   buf << static_cast<walberla::mpi::MPIRank>(obj.particle_.getOwner());
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ParticleRemovalInformationNotification::Parameters& objparam )
{
   buf.readDebugMarker( "ri" );
   buf >> objparam.uid_;
   buf >> objparam.owner_;
   return buf;
}

} // mpi
} // walberla
