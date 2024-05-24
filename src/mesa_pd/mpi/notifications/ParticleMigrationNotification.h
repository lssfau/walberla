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
      std::unordered_set<walberla::mpi::MPIRank> ghostOwners_ {};
      walberla::mesa_pd::Vec3 oldForce_ {real_t(0)};
      walberla::mesa_pd::Vec3 oldTorque_ {real_t(0)};
      walberla::mesa_pd::Vec3 hydrodynamicForce_ {real_t(0)};
      walberla::mesa_pd::Vec3 hydrodynamicTorque_ {real_t(0)};
      walberla::mesa_pd::Vec3 oldHydrodynamicForce_ {real_t(0)};
      walberla::mesa_pd::Vec3 oldHydrodynamicTorque_ {real_t(0)};
      walberla::real_t totalDisplacement_ {real_t(0)};
      walberla::real_t collisionForceNorm_ {real_t(0)};
      walberla::real_t virtualMass_ {real_t(0)};
      walberla::real_t invMassIncludingVirtual_ {real_t(0)};
      walberla::mesa_pd::Vec3 oldLinearAcceleration_ {real_t(0)};
      walberla::mesa_pd::Mat3 invInertiaBF_ {real_t(0)};
      walberla::mesa_pd::Mat3 virtualInertiaBF_ {real_t(0)};
      walberla::mesa_pd::Mat3 invInertiaBFIncludingVirtual_ {real_t(0)};
      walberla::mesa_pd::Vec3 oldAngularAcceleration_ {real_t(0)};
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
   buf << obj.particle_.getHydrodynamicForce();
   buf << obj.particle_.getHydrodynamicTorque();
   buf << obj.particle_.getOldHydrodynamicForce();
   buf << obj.particle_.getOldHydrodynamicTorque();
   buf << obj.particle_.getTotalDisplacement();
   buf << obj.particle_.getCollisionForceNorm();
   buf << obj.particle_.getVirtualMass();
   buf << obj.particle_.getInvMassIncludingVirtual();
   buf << obj.particle_.getOldLinearAcceleration();
   buf << obj.particle_.getInvInertiaBF();
   buf << obj.particle_.getVirtualInertiaBF();
   buf << obj.particle_.getInvInertiaBFIncludingVirtual();
   buf << obj.particle_.getOldAngularAcceleration();
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
   buf >> objparam.hydrodynamicForce_;
   buf >> objparam.hydrodynamicTorque_;
   buf >> objparam.oldHydrodynamicForce_;
   buf >> objparam.oldHydrodynamicTorque_;
   buf >> objparam.totalDisplacement_;
   buf >> objparam.collisionForceNorm_;
   buf >> objparam.virtualMass_;
   buf >> objparam.invMassIncludingVirtual_;
   buf >> objparam.oldLinearAcceleration_;
   buf >> objparam.invInertiaBF_;
   buf >> objparam.virtualInertiaBF_;
   buf >> objparam.invInertiaBFIncludingVirtual_;
   buf >> objparam.oldAngularAcceleration_;
   return buf;
}

template<>
struct BufferSizeTrait< mesa_pd::ParticleMigrationNotification > {
   static const bool constantSize = false;
};

} // mpi
} // walberla