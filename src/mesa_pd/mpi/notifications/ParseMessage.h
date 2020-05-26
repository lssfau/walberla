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
//! \file ParseMessage.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Parsing of messages
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/domain/IDomain.h>
#include <mesa_pd/mpi/notifications/NewGhostParticleNotification.h>
#include <mesa_pd/mpi/notifications/NotificationType.h>
#include <mesa_pd/mpi/notifications/ParticleGhostCopyNotification.h>
#include <mesa_pd/mpi/notifications/ParticleMigrationNotification.h>
#include <mesa_pd/mpi/notifications/ParticleRemoteMigrationNotification.h>
#include <mesa_pd/mpi/notifications/ParticleRemovalNotification.h>
#include <mesa_pd/mpi/notifications/ParticleRemovalInformationNotification.h>
#include <mesa_pd/mpi/notifications/ParticleUpdateNotification.h>

#include <core/debug/Debug.h>
#include <core/logging/Logging.h>
#include <core/mpi/RecvBuffer.h>

namespace walberla {
namespace mesa_pd {

class ParseMessage
{
public:
   void operator()(int sender,
                   walberla::mpi::RecvBuffer& rb,
                   data::ParticleStorage& ps,
                   const domain::IDomain& domain);
   void allowMultipleGhostCopyNotifications(bool val) {allowMultipleGhostCopyNotifications_ = val;}
private:
   int receiver_ = int_c( walberla::mpi::MPIManager::instance()->rank() );

   bool allowMultipleGhostCopyNotifications_ = false;
};

inline
void ParseMessage::operator()(int sender,
                              walberla::mpi::RecvBuffer& rb,
                              data::ParticleStorage& ps,
                              const domain::IDomain& domain)
{
   NotificationType notificationType;
   rb >> notificationType;

   switch( notificationType ) {
   case PARTICLE_GHOST_COPY_NOTIFICATION: {
      typename ParticleGhostCopyNotification::Parameters objparam;
      rb >> objparam;

      WALBERLA_LOG_DETAIL( "Received PARTICLE_GHOST_COPY_NOTIFICATION for particle " << objparam.uid << "from neighboring process with rank " << sender );

      if ( ps.find(objparam.uid) == ps.end() )
      {
         auto pIt = createNewParticle(ps, objparam);

         domain.correctParticlePosition(pIt->getPositionRef());

         //WALBERLA_CHECK(!data::particle_flags::isSet(pIt->getFlags(), data::particle_flags::GHOST));
         data::particle_flags::set(pIt->getFlagsRef(), data::particle_flags::GHOST);
      } else
      {
         if (allowMultipleGhostCopyNotifications_)
         {
            //superfluous ghost creation messages might be send during ghost owner sync
            WALBERLA_LOG_DETAIL("Ghost particle with id " << objparam.uid << " already existend.");
         } else
         {
            WALBERLA_ABORT("Ghost particle with id " << objparam.uid << " already existend.");
         }
      }

      WALBERLA_LOG_DETAIL( "Processed PARTICLE_GHOST_COPY_NOTIFICATION for particle " << objparam.uid << "."  );

      break;
   }
   case PARTICLE_UPDATE_NOTIFICATION: {
      typename ParticleUpdateNotification::Parameters objparam;
      rb >> objparam;

      WALBERLA_LOG_DETAIL( "Received PARTICLE_UPDATE_NOTIFICATION for particle " << objparam.uid <<
                           " from neighboring process with rank " << sender );

      auto pIt = ps.find( objparam.uid );
      WALBERLA_CHECK_UNEQUAL( pIt, ps.end() );

      WALBERLA_CHECK_EQUAL( pIt->getOwner(), sender, "Update notifications must be sent by owner." );
      WALBERLA_CHECK(data::particle_flags::isSet(pIt->getFlags(), data::particle_flags::GHOST),
                     "Update notification must only concern shadow copies.");
      pIt->setUid(objparam.uid);
      pIt->setPosition(objparam.position);
      pIt->setLinearVelocity(objparam.linearVelocity);
      pIt->setRotation(objparam.rotation);
      pIt->setAngularVelocity(objparam.angularVelocity);
      pIt->setRadiusAtTemperature(objparam.radiusAtTemperature);
      pIt->setOldContactHistory(objparam.oldContactHistory);
      pIt->setTemperature(objparam.temperature);

      domain.correctParticlePosition(pIt->getPositionRef());

      WALBERLA_LOG_DETAIL( "Processed PARTICLE_UPDATE_NOTIFICATION." );

      break;
   }
   case PARTICLE_MIGRATION_NOTIFICATION: {
      ParticleMigrationNotification::Parameters objparam;
      rb >> objparam;

      WALBERLA_LOG_DETAIL( "Received PARTICLE_MIGRATION_NOTIFICATION for particle " << objparam.uid_ <<
                           " from neighboring process with rank " << sender );

      auto pIt = ps.find( objparam.uid_ );
      WALBERLA_CHECK_UNEQUAL( pIt, ps.end(),
                              "Object with id: " << objparam.uid_ << " not found! Cannot transfer ownership!" );
      WALBERLA_CHECK_EQUAL( sender, pIt->getOwner(), "Migration notifications must be sent by previous owner.");
      WALBERLA_CHECK(data::particle_flags::isSet(pIt->getFlags(), data::particle_flags::GHOST),
                     "Migration notification must only concern ghost particles");

      WALBERLA_CHECK( domain.isContainedInProcessSubdomain( uint_c(receiver_), pIt->getPosition() ),
                      "Receiving particle migration even though we do not own it." );

      pIt->setOwner(receiver_);
      data::particle_flags::unset(pIt->getFlagsRef(), data::particle_flags::GHOST);
      pIt->setGhostOwners(objparam.ghostOwners_);
      pIt->setOldForce(objparam.oldForce_);
      pIt->setOldTorque(objparam.oldTorque_);
      pIt->setHydrodynamicForce(objparam.hydrodynamicForce_);
      pIt->setHydrodynamicTorque(objparam.hydrodynamicTorque_);
      pIt->setOldHydrodynamicForce(objparam.oldHydrodynamicForce_);
      pIt->setOldHydrodynamicTorque(objparam.oldHydrodynamicTorque_);

      WALBERLA_LOG_DETAIL( "Processed PARTICLE_MIGRATION_NOTIFICATION." );

      break;
   }
   case PARTICLE_REMOTE_MIGRATION_NOTIFICATION: {
      ParticleRemoteMigrationNotification::Parameters objparam;
      rb >> objparam;

      WALBERLA_LOG_DETAIL( "Received PARTICLE_REMOTE_MIGRATION_NOTIFICATION for particle " << objparam.uid_ <<
                           " from neighboring process with rank " << sender <<
                           " (previous owner):\nnew owner = " << objparam.newOwner_ );

      auto pIt = ps.find( objparam.uid_ );
      WALBERLA_CHECK_UNEQUAL( pIt, ps.end() );

      WALBERLA_CHECK_EQUAL( sender, pIt->getOwner(),
                            "Remote migration notifications must be sent by previous owner." );
      WALBERLA_CHECK(data::particle_flags::isSet(pIt->getFlags(), data::particle_flags::GHOST),
                     "Particles in remote migration notifications must be available as ghost particles in local process.");
      WALBERLA_CHECK_UNEQUAL( objparam.newOwner_, receiver_,
                              "Particles in remote migration notifications may not migrate to local process." );

      pIt->setOwner(objparam.newOwner_);

      WALBERLA_LOG_DETAIL( "Processed PARTICLE_REMOTE_MIGRATION_NOTIFICATION." );

      break;
   }
   case PARTICLE_REMOVAL_NOTIFICATION: {
      ParticleRemovalNotification::Parameters objparam;
      rb >> objparam;

      WALBERLA_LOG_DETAIL( "Received PARTICLE_REMOVAL_NOTIFICATION for particle " << objparam.uid_ <<
                           " from neighboring process with rank " << sender << " (owner)." );

      // Remove ghost particle as prompted.
      auto pIt = ps.find( objparam.uid_ );
      WALBERLA_CHECK_UNEQUAL( pIt, ps.end() );

      WALBERLA_CHECK(data::particle_flags::isSet(pIt->getFlags(), data::particle_flags::GHOST),
                     "Only ghost particles should be removed by this message.");

      WALBERLA_CHECK_EQUAL( pIt->getOwner(), sender,
                            "Only owner is allowed to send removal notifications." );

      ps.erase(pIt);

      WALBERLA_LOG_DETAIL( "Processed PARTICLE_REMOVAL_NOTIFICATION" );

      break;
   }
   case NEW_GHOST_PARTICLE_NOTIFICATION: {
      NewGhostParticleNotification::Parameters objparam;
      rb >> objparam;

      WALBERLA_LOG_DETAIL( "Received new ghost particle notification for particle " <<
                           objparam.uid_ <<
                           " from neighboring process with rank " <<
                           sender <<
                           "." );

      auto pIt = ps.find( objparam.uid_ );
      WALBERLA_CHECK_UNEQUAL( pIt, ps.end() );

      pIt->getGhostOwnersRef().insert( objparam.newOwner_ );

      WALBERLA_LOG_DETAIL( "Processed new ghost particle notification" );

      break;
   }
   case PARTICLE_REMOVAL_INFORMATION_NOTIFICATION: {
      ParticleRemovalInformationNotification::Parameters objparam;
      rb >> objparam;

      WALBERLA_LOG_DETAIL( "Received particle removal information notification for particle " <<
                           objparam.uid_ <<
                           " from neighboring process with rank " <<
                           sender <<
                           "." );

      if (objparam.owner_ == receiver_)
      {
         using namespace walberla::mesa_pd::data::particle_flags;
         auto pIt = ps.find( objparam.uid_ );
         WALBERLA_CHECK_UNEQUAL( pIt, ps.end() );
         WALBERLA_CHECK(!isSet( pIt->getFlags(), GHOST));

         pIt->getGhostOwnersRef().erase( sender );
         pIt->getNeighborStateRef().erase( sender );
      } else
      {
         using namespace walberla::mesa_pd::data::particle_flags;
         auto pIt = ps.find( objparam.uid_ );
         if (pIt != ps.end() )
         {
            WALBERLA_CHECK(isSet( pIt->getFlags(), GHOST));
            pIt->getNeighborStateRef().erase( sender );
         }
      }

      WALBERLA_LOG_DETAIL( "Processed rigid body removal information notification" );

      break;
   }
   default:
      WALBERLA_ABORT( "Received invalid notification type: " << notificationType << " from sender: " << sender );
   }
}

}  // namespace mesa_pd
}  // namespace walberla