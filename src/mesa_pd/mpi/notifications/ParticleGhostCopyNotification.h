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
#include <mesa_pd/mpi/ShapePackUnpack.h>
#include <mesa_pd/mpi/notifications/NotificationType.h>

#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace mesa_pd {

/**
 * A complete particle copy for a new ghost particle.
 *
 * Copies all properties marked ON_GHOST_CREATION or ALWAYS.
 */
class ParticleGhostCopyNotification
{
public:
   struct Parameters
   {
      walberla::id_t uid {UniqueID<data::Particle>::invalidID()};
      walberla::mesa_pd::Vec3 position {real_t(0)};
      walberla::real_t interactionRadius {real_t(0)};
      walberla::mesa_pd::data::particle_flags::FlagT flags {};
      int owner {-1};
      walberla::mesa_pd::Vec3 linearVelocity {real_t(0)};
      walberla::real_t invMass {real_t(1)};
      size_t shapeID {};
      walberla::mesa_pd::Rot3 rotation {};
      walberla::mesa_pd::Vec3 angularVelocity {real_t(0)};
      walberla::real_t radiusAtTemperature {real_t(0)};
      uint_t type {0};
      std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> oldContactHistory {};
      walberla::real_t temperature {real_t(0)};
      int64_t clusterID {-1};
      int64_t segmentID {-1};
   };

   inline explicit ParticleGhostCopyNotification( const data::Particle& particle ) : particle_(particle) {}
   const data::Particle& particle_;
};

inline data::ParticleStorage::iterator createNewParticle(data::ParticleStorage& ps, const ParticleGhostCopyNotification::Parameters& data)
{
   WALBERLA_ASSERT_EQUAL(ps.find(data.uid), ps.end(), "Particle with same uid already existent!");

   auto pIt = ps.create(data.uid);
   pIt->setUid(data.uid);
   pIt->setPosition(data.position);
   pIt->setInteractionRadius(data.interactionRadius);
   pIt->setFlags(data.flags);
   pIt->setOwner(data.owner);
   pIt->setLinearVelocity(data.linearVelocity);
   pIt->setInvMass(data.invMass);
   pIt->setShapeID(data.shapeID);
   pIt->setRotation(data.rotation);
   pIt->setAngularVelocity(data.angularVelocity);
   pIt->setRadiusAtTemperature(data.radiusAtTemperature);
   pIt->setType(data.type);
   pIt->setOldContactHistory(data.oldContactHistory);
   pIt->setTemperature(data.temperature);
   pIt->setClusterID(data.clusterID);
   pIt->setSegmentID(data.segmentID);
   return pIt;
}

template<>
struct NotificationTrait<ParticleGhostCopyNotification>
{
   static const NotificationType id = PARTICLE_GHOST_COPY_NOTIFICATION;
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ParticleGhostCopyNotification& obj )
{
   buf.addDebugMarker( "cn" );
   buf << obj.particle_.getUid();
   buf << obj.particle_.getPosition();
   buf << obj.particle_.getInteractionRadius();
   buf << obj.particle_.getFlags();
   buf << obj.particle_.getOwner();
   buf << obj.particle_.getLinearVelocity();
   buf << obj.particle_.getInvMass();
   buf << obj.particle_.getShapeID();
   buf << obj.particle_.getRotation();
   buf << obj.particle_.getAngularVelocity();
   buf << obj.particle_.getRadiusAtTemperature();
   buf << obj.particle_.getType();
   buf << obj.particle_.getOldContactHistory();
   buf << obj.particle_.getTemperature();
   buf << obj.particle_.getClusterID();
   buf << obj.particle_.getSegmentID();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ParticleGhostCopyNotification::Parameters& objparam )
{
   buf.readDebugMarker( "cn" );
   buf >> objparam.uid;
   buf >> objparam.position;
   buf >> objparam.interactionRadius;
   buf >> objparam.flags;
   buf >> objparam.owner;
   buf >> objparam.linearVelocity;
   buf >> objparam.invMass;
   buf >> objparam.shapeID;
   buf >> objparam.rotation;
   buf >> objparam.angularVelocity;
   buf >> objparam.radiusAtTemperature;
   buf >> objparam.type;
   buf >> objparam.oldContactHistory;
   buf >> objparam.temperature;
   buf >> objparam.clusterID;
   buf >> objparam.segmentID;
   return buf;
}

} // mpi
} // walberla