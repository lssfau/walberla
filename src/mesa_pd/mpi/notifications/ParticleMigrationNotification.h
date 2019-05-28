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
//! \file ParticleMigrationNotification.h
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

#include <core/mpi/BufferDataTypeExtensions.h>
#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>


namespace walberla {
namespace mesa_pd {

/**
 * Migrate the particle to this process. Making the receiver the new owner.
 */
class ParticleMigrationNotification {
public:
   struct Parameters {
      id_t uid_;
      std::vector<int> ghostOwners_ {};
      walberla::mesa_pd::Vec3 oldForce_ {real_t(0)};
      walberla::mesa_pd::Vec3 oldTorque_ {real_t(0)};
   };

   inline explicit ParticleMigrationNotification( const data::Particle& particle ) : particle_(particle) {}
   const data::Particle& particle_;
};

template<>
struct NotificationTrait<ParticleMigrationNotification>
{
   static const NotificationType id = PARTICLE_MIGRATION_NOTIFICATION;
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ParticleMigrationNotification& obj )
{
   buf.addDebugMarker( "mn" );
   buf << obj.particle_.getUid();
   buf << obj.particle_.getGhostOwners();
   buf << obj.particle_.getOldForce();
   buf << obj.particle_.getOldTorque();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ParticleMigrationNotification::Parameters& objparam )
{
   buf.readDebugMarker( "mn" );
   buf >> objparam.uid_;
   buf >> objparam.ghostOwners_;
   buf >> objparam.oldForce_;
   buf >> objparam.oldTorque_;
   return buf;
}

template<>
struct BufferSizeTrait< mesa_pd::ParticleMigrationNotification > {
   static const bool constantSize = false;
};

} // mpi
} // walberla