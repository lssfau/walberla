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
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Parsing of messages
//
//======================================================================================================================

#pragma once

#include "pe/rigidbody/BodyStorage.h"
#include "pe/communication/DynamicMarshalling.h"
#include "pe/communication/RigidBodyCopyNotification.h"
#include "pe/communication/RigidBodyDeletionNotification.h"
#include "pe/communication/RigidBodyForceNotification.h"
#include "pe/communication/RigidBodyMigrationNotification.h"
#include "pe/communication/RigidBodyNewShadowCopyNotification.h"
#include "pe/communication/RigidBodyRemoteMigrationNotification.h"
#include "pe/communication/RigidBodyRemovalInformationNotification.h"
#include "pe/communication/RigidBodyRemovalNotification.h"
#include "pe/communication/RigidBodyUpdateNotification.h"
#include "pe/communication/RigidBodyVelocityCorrectionNotification.h"
#include "pe/communication/RigidBodyVelocityUpdateNotification.h"
#include "pe/communication/rigidbody/Sphere.h"
#include "pe/BlockFunctions.h"

#include "core/mpi/RecvBuffer.h"
#include "blockforest/Block.h"
#include "domain_decomposition/BlockStorage.h"
#include "core/debug/Debug.h"

namespace walberla {
namespace pe {
namespace communication {

template <typename BodyTypeTuple>
void parseMessage(Owner sender, mpi::RecvBuffer& rb, const domain_decomposition::BlockStorage& blockStorage, IBlock& block, BodyStorage& localStorage, BodyStorage& shadowStorage)
{
   const Owner receiver( int_c( mpi::MPIManager::instance()->rank() ), block.getId().getID() );

   NotificationType notificationType;

   // receiving shadow copies [N], shadow copy updates [DN], (local) migrations [N], (remote) migrations [D], deletions [DN] and removal notifications [DN] from neighbors (N) and distant processes (D)
      unmarshal( rb, notificationType );

      switch( notificationType ) {
         case rigidBodyCopyNotification: {
            typename RigidBodyCopyNotification::Parameters objparam;
            unmarshal( rb, objparam );

            WALBERLA_LOG_DETAIL( "Received " << objparam.geomType_ << " copy notification from neighboring process with rank " << sender );

            BodyPtr obj = UnmarshalDynamically<BodyTypeTuple>::execute(rb, objparam.geomType_, blockStorage.getDomain(), block.getAABB());
            obj->setRemote( true );
            obj->MPITrait.setBlockState( sender.blockID_ );

            if (shadowStorage.find( obj->getSystemID() ) == shadowStorage.end())
            {
               WALBERLA_LOG_DETAIL( "Adding new shadow copy with id " << obj->getSystemID() << " to domain " << block.getAABB() << ".\n" << *obj);
               shadowStorage.add( std::move(obj) );
            } else
            {
               WALBERLA_LOG_DETAIL( "Shadow copy with id " << obj->getSystemID() << " already existend.");
            }

            WALBERLA_LOG_DETAIL( "Processed " << objparam.geomType_ << " copy notification."  );

            break;
         }
         case rigidBodyUpdateNotification: {
            typename RigidBodyUpdateNotification::Parameters objparam;
            unmarshal( rb, objparam );

            WALBERLA_LOG_DETAIL( "Received rigid body update notification for body " << objparam.sid_ << " from neighboring process with rank " << sender << ":\nv = " << objparam.v_ << "\nw = " << objparam.w_ << "\nposition = " << objparam.gpos_ << "\nquaternion = " << objparam.q_);

            auto bodyIt = shadowStorage.find( objparam.sid_ );
            WALBERLA_ASSERT_UNEQUAL( bodyIt, shadowStorage.end() );
            BodyID b( bodyIt.getBodyID() );

            WALBERLA_ASSERT( b->MPITrait.getOwner().blockID_ == sender.blockID_, "Update notifications must be sent by owner.\n" << b->MPITrait.getOwner().blockID_ << " != "<< sender.blockID_ );
            WALBERLA_ASSERT( b->isRemote(), "Update notification must only concern shadow copies." );

            correctBodyPosition(blockStorage.getDomain(), block.getAABB().center(), objparam.gpos_);
            b->setPosition   ( objparam.gpos_ );
            b->setOrientation( objparam.q_    );
            b->setLinearVel  ( objparam.v_    );
            b->setAngularVel ( objparam.w_    );

            WALBERLA_LOG_DETAIL( "Processed rigid body update notification.");//:\n" << b );

            break;
         }
         case rigidBodyMigrationNotification: {
            RigidBodyMigrationNotification::Parameters objparam;
            unmarshal( rb, objparam );

            WALBERLA_LOG_DETAIL( "Received rigid body migration notification for body " << objparam.sid_ << " from neighboring process with rank " << sender );
            WALBERLA_LOG_DETAIL( "New owner will be " << receiver );

            auto bodyIt = shadowStorage.find( objparam.sid_ );
            if ( bodyIt == shadowStorage.end() ) WALBERLA_ABORT( "Object with id: " << objparam.sid_ << " not found in shadowStorage! Cannot transfer ownership! \nlocal domain: " << block.getAABB() );
            BodyID b( bodyIt.getBodyID() );

            localStorage.add( shadowStorage.release( b ) );

            WALBERLA_ASSERT( sender.blockID_ == b->MPITrait.getOwner().blockID_, "Migration notifications must be sent by previous owner.\n" << b->MPITrait.getOwner().blockID_ << " != "<< sender.blockID_ );
            WALBERLA_ASSERT( b->isRemote(), "Bodies in migration notifications must be available as shadow copies in local process." );
            WALBERLA_ASSERT( block.getAABB().contains( b->getPosition() ), "Receiving body migration even though we do not own it." );

            b->MPITrait.setOwner( receiver );
            b->setRemote( false );

            WALBERLA_ASSERT_EQUAL(b->MPITrait.sizeShadowOwners(), 0);
            b->MPITrait.clearShadowOwners();
            // reconstruct list of shadow copy holders
            for( std::size_t i = 0; i < objparam.reglist_.size(); ++i ) {
                  b->MPITrait.registerShadowOwner( objparam.reglist_[i] );
            }

            WALBERLA_LOG_DETAIL( "Processed rigid body migration notification:\n" << *b );

            break;
         }
         case rigidBodyRemoteMigrationNotification: {
            RigidBodyRemoteMigrationNotification::Parameters objparam;
            unmarshal( rb, objparam );

            WALBERLA_LOG_DETAIL( "Received rigid body remote migration notification for body " << objparam.sid_ << " from neighboring process with rank " << sender << " (previous owner):\nnew owner = " << objparam.to_ );

            auto bodyIt = shadowStorage.find( objparam.sid_ );
            WALBERLA_ASSERT_UNEQUAL( bodyIt, shadowStorage.end() );
            BodyID b( bodyIt.getBodyID() );

            WALBERLA_ASSERT( sender.blockID_ == b->MPITrait.getOwner().blockID_, "Remote migration notifications must be sent by previous owner.\n" << sender.blockID_ << " != " << b->MPITrait.getOwner().blockID_  );
            WALBERLA_ASSERT( b->isRemote(), "Bodies in remote migration notifications must be available as shadow copies in local process." );
            WALBERLA_ASSERT( objparam.to_.blockID_ != receiver.blockID_, "Bodies in remote migration notifications may not migrate to local process." );

            b->MPITrait.setOwner( objparam.to_ );

            WALBERLA_LOG_DETAIL( "Processed rigid body remote migration notification:\n" << *b );

            break;
         }
         case rigidBodyRemovalNotification: {
            RigidBodyRemovalNotification::Parameters objparam;
            unmarshal( rb, objparam );

            WALBERLA_LOG_DETAIL( "Received rigid body removal notification for body " << objparam.sid_ << " from neighboring process with rank " << sender << " (owner)." );

            // Remove shadow copy as prompted.
            auto bodyIt = shadowStorage.find( objparam.sid_ );
            WALBERLA_ASSERT_UNEQUAL( bodyIt, shadowStorage.end() );
            BodyID b( bodyIt.getBodyID() );

            // TODO assert that we indeed do not need the shadow copy anymore
            WALBERLA_ASSERT( b->MPITrait.getOwner().blockID_ == sender.blockID_, "Only owner is allowed to send removal notifications.\n" << b->MPITrait.getOwner().blockID_ << " != "<< sender.blockID_ );

            shadowStorage.remove( b );

            WALBERLA_LOG_DETAIL( "Processed rigid body removal notification" );

            break;
         }
      case rigidBodyDeletionNotification: {
         RigidBodyDeletionNotification::Parameters objparam;
         unmarshal( rb, objparam );

         WALBERLA_LOG_DETAIL( "Received rigid body deletion notification for body " << objparam.sid_ << " from neighboring process with rank " << sender << " (owner)." );

         // Remove invalid shadow copy.
         auto bodyIt = shadowStorage.find( objparam.sid_ );
         WALBERLA_ASSERT_UNEQUAL( bodyIt, shadowStorage.end() );
         BodyID b( bodyIt.getBodyID() );

         WALBERLA_ASSERT( b->MPITrait.getOwner().blockID_ == sender.blockID_, "Only owner is allowed to send deletion notifications.\n" << b->MPITrait.getOwner().blockID_ << " != "<< sender.blockID_ );

         shadowStorage.remove( b );

         WALBERLA_LOG_DETAIL( "Processed rigid body deletion notification" );

         break;
      }
      case rigidBodyNewShadowCopyNotification: {
         RigidBodyNewShadowCopyNotification::Parameters objparam;
         unmarshal( rb, objparam );

         WALBERLA_LOG_DETAIL( "Received rigid body new shadow copy notification for body " << objparam.sid_ << " from neighboring process with rank " << sender << "." );

         // Remove invalid shadow copy.
         auto bodyIt = localStorage.find( objparam.sid_ );
         WALBERLA_ASSERT_UNEQUAL( bodyIt, localStorage.end() );
         BodyID b( bodyIt.getBodyID() );

         b->MPITrait.registerShadowOwner( objparam.newOwner_ );

         WALBERLA_LOG_DETAIL( "Processed rigid body new shadow copy notification" );

         break;
      }
      case rigidBodyRemovalInformationNotification: {
         RigidBodyRemovalInformationNotification::Parameters objparam;
         unmarshal( rb, objparam );

         WALBERLA_LOG_DETAIL( "Received rigid body removal information notification for body " << objparam.sid_ << " from neighboring process with rank " << sender << "." );

         if (objparam.owner_ == receiver)
         {
            auto bodyIt = localStorage.find( objparam.sid_ );
            WALBERLA_ASSERT_UNEQUAL( bodyIt, localStorage.end() );
            BodyID b( bodyIt.getBodyID() );
            b->MPITrait.deregisterShadowOwner( sender );
            b->MPITrait.unsetBlockState( sender.blockID_ );
         } else
         {
            auto bodyIt = shadowStorage.find( objparam.sid_ );
            if (bodyIt != shadowStorage.end() )
            {
               BodyID b( bodyIt.getBodyID() );
               b->MPITrait.unsetBlockState( sender.blockID_ );
            }
         }

         WALBERLA_LOG_DETAIL( "Processed rigid body removal information notification" );

         break;
      }
         default:
            throw std::runtime_error( "Received invalid notification type." );
      }
}

inline
void parseForceReduceMessage(Owner sender, mpi::RecvBuffer& rb, const domain_decomposition::BlockStorage& /*blockStorage*/, IBlock& /*block*/, BodyStorage& localStorage, BodyStorage& /*shadowStorage*/)
{
   WALBERLA_UNUSED(sender);

   NotificationType notificationType;

   // receiving shadow copies [N], shadow copy updates [DN], (local) migrations [N], (remote) migrations [D], deletions [DN] and removal notifications [DN] from neighbors (N) and distant processes (D)
   unmarshal( rb, notificationType );

   switch( notificationType ) {
   case rigidBodyForceNotification: {
      RigidBodyForceNotification::Parameters objparam;
      unmarshal( rb, objparam );

      auto bodyIt = localStorage.find( objparam.sid_ );
      WALBERLA_ASSERT_UNEQUAL( bodyIt, localStorage.end() );
      BodyID b( bodyIt.getBodyID() );

      WALBERLA_ASSERT( !b->isRemote(), "Update notification must only concern local bodies." );

      b->addForce( objparam.f_ );
      b->addTorque( objparam.tau_ );

      WALBERLA_LOG_DETAIL( "Received rigid body force contribution from neighboring block " << sender.blockID_ <<
                           ":\nf = " << objparam.f_ << "\ntau = " << objparam.tau_ << "\nf_total = " << b->getForce() <<
                           "\ntau_total = " << b->getTorque() << "\n" );
      break;
   }
      default:
         throw std::runtime_error( "Received invalid notification type." );
   }
}

inline
void parseForceDistributeMessage(Owner sender, mpi::RecvBuffer& rb, const domain_decomposition::BlockStorage& /*blockStorage*/, IBlock& /*block*/, BodyStorage& /*localStorage*/, BodyStorage& shadowStorage)
{
   WALBERLA_UNUSED(sender);

   NotificationType notificationType;
   unmarshal( rb, notificationType );

   switch( notificationType ) {
   case rigidBodyForceNotification: {
      RigidBodyForceNotification::Parameters objparam;
      unmarshal( rb, objparam );

      auto bodyIt = shadowStorage.find( objparam.sid_ );
      WALBERLA_ASSERT_UNEQUAL( bodyIt, shadowStorage.end() );
      BodyID b( bodyIt.getBodyID() );

      WALBERLA_ASSERT( b->isRemote(), "Update notification must only concern shadow bodies." );

      b->setForce( objparam.f_ );
      b->setTorque( objparam.tau_ );

      WALBERLA_LOG_DETAIL( "Received rigid body force from owning block " << sender.blockID_ <<
                           ":\nf = " << objparam.f_ << "\ntau = " << objparam.tau_ << "\nf_total = " << b->getForce() <<
                           "\ntau_total = " << b->getTorque() << "\n" );
      break;
   }
      default:
         throw std::runtime_error( "Received invalid notification type." );
   }
}

}  // namespace communication
}  // namespace pe
}  // namespace walberla
