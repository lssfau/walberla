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
//! \file Synchronization.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/BlockFunctions.h"
#include "pe/rigidbody/BodyStorage.h"

#include "pe/communication/ParseMessage.h"
#include "pe/communication/DynamicMarshalling.h"
#include "pe/communication/RigidBodyCopyNotification.h"
#include "pe/communication/RigidBodyDeletionNotification.h"
#include "pe/communication/RigidBodyForceNotification.h"
#include "pe/communication/RigidBodyMigrationNotification.h"
#include "pe/communication/RigidBodyRemoteMigrationNotification.h"
#include "pe/communication/RigidBodyRemovalNotification.h"
#include "pe/communication/RigidBodyUpdateNotification.h"
#include "pe/communication/RigidBodyVelocityCorrectionNotification.h"
#include "pe/communication/RigidBodyVelocityUpdateNotification.h"
#include "pe/communication/PackNotification.h"

#include "RemoveAndNotify.h"

#include "blockforest/BlockForest.h"
#include "core/mpi/BufferSystem.h"
#include "core/timing/TimingTree.h"
#include "domain_decomposition/MapPointToPeriodicDomain.h"

namespace walberla {
namespace pe {

template <typename BodyTypeTuple>
void generateSynchonizationMessages(mpi::BufferSystem& bs, const Block& block, BodyStorage& localStorage, BodyStorage& shadowStorage, const real_t dx, const bool syncNonCommunicatingBodies)
{
   using namespace walberla::pe::communication;

   const Owner me( int_c( mpi::MPIManager::instance()->rank() ), block.getId().getID() );
   const math::AABB domain = block.getBlockStorage().getDomain();

   WALBERLA_LOG_DETAIL( "Assembling of body synchronization message starts..." );

   // position update
   for( auto body = localStorage.begin(); body != localStorage.end(); )
   {
      //correct position to make sure body is always inside the domain!
      if (!body->isFixed())
      {
         auto pos = body->getPosition();
         block.getBlockStorage().mapToPeriodicDomain(pos);
         body->setPosition(pos);
      }

      bool isInsideDomain = domain.contains( body->getAABB(), -dx );

      WALBERLA_ASSERT( !body->isRemote(), "Central body storage contains remote bodies." );

      if( !body->isCommunicating() && !syncNonCommunicatingBodies ) {
         ++body;
         continue;
      }

      const Vec3 gpos( body->getPosition() );
      BodyID     b   ( body.getBodyID() );

      if (body->isMarkedForDeletion())
      {
         // delete it
         WALBERLA_LOG_DETAIL( "Sending deletion notifications for body " << body->getSystemID() << " due to manual deletion." );

         body = removeAndNotify( me, bs, localStorage, body );

         // Note: Attached shadow copies are not deleted here. Instead we rely on the deferred deletion since we no
         // longer need the shadow copy: The owner of an attached shadow copy will receive a deletion notification, detach
         // the attachable, delete the shadow copy of the deleted body and send us a removal notification of the body
         // of which we own a shadow copy in the next position update since (probably) we no longer require the body but
         // are still part of its registration list.
         continue;
      }

      // Note: At this point we know that the body was locally owned before the position update.

      if( b->MPITrait.getOwner().blockID_ != me.blockID_ )
      {
         WALBERLA_LOG_WARNING( "Found remote body with sid " << b->getSystemID() << " in central storage:\n" << *body );
      }

      WALBERLA_ASSERT_EQUAL( b->MPITrait.getOwner().blockID_, me.blockID_, "Owner field in local body storage does not contain the own process rank." );

      WALBERLA_LOG_DETAIL( "Processing local body " << b->getSystemID() );

      // Update (nearest) neighbor processes.
      for( uint_t nb = uint_t(0); nb < block.getNeighborhoodSize(); ++nb )
      {
         Owner nbProcess = Owner(int_c( block.getNeighborProcess(nb) ), block.getNeighborId(nb).getID() );
         if( (isInsideDomain ? block.getNeighborAABB(nb).intersects( b->getAABB(), dx ) : block.getBlockStorage().periodicIntersect(block.getNeighborAABB(nb), b->getAABB(), dx)) )
         {
            // The body is needed by the process.

            if( body->MPITrait.isShadowOwnerRegistered( nbProcess ) ) {
               mpi::SendBuffer& buffer( bs.sendBuffer(nbProcess.rank_) );

               WALBERLA_LOG_DETAIL( "Sending update notification for body " << b->getSystemID() << " to process " << (nbProcess) );

               packNotification(me, nbProcess, buffer, RigidBodyUpdateNotification( *b ));
            }
            else {
               mpi::SendBuffer& buffer( bs.sendBuffer(nbProcess.rank_) );

               WALBERLA_LOG_DETAIL( "Sending shadow copy notification for body " << b->getSystemID() << " to process " << (nbProcess) );

               packNotification(me, nbProcess, buffer, RigidBodyCopyNotification( *b ));
               MarshalDynamically<BodyTypeTuple>::execute( buffer, *b );

               b->MPITrait.registerShadowOwner( nbProcess );
            }
         }
         else
         {
            if( body->MPITrait.isShadowOwnerRegistered( nbProcess ) )
            {
               // In case the rigid body no longer intersects the remote process nor interacts with it but is registered,
               // send removal notification.
               mpi::SendBuffer& buffer( bs.sendBuffer(nbProcess.rank_) );

               WALBERLA_LOG_DETAIL( "Sending removal notification for body " << b->getSystemID() << " to process " << block.getNeighborProcess(nb) );

               packNotification(me, nbProcess, buffer, RigidBodyRemovalNotification( *b ));

               body->MPITrait.deregisterShadowOwner( nbProcess );
            }
         }
      }

      // Update remote processes (no intersections possible; (long-range) interactions only).
      // TODO iterate over all processes attached bodies are owned by (skipping nearest neighbors)
      // depending on registration send update or copy

      if( !block.getAABB().contains( gpos ) )
      {
         // Body is no longer locally owned (body is about to be migrated).
         Owner owner( findContainingProcess( block, gpos ) );

         WALBERLA_LOG_DETAIL( "Local body " << b->getSystemID() << " is no longer on process " << body->MPITrait.getOwner() << " but on process " << owner );

         if( owner.rank_ < 0 ) {
            // No owner found: Outflow condition.
            WALBERLA_LOG_DETAIL( "Sending deletion notifications for body " << body->getSystemID() << " due to outflow." );

            // Registered processes receive removal notification in the remove() routine.
            //todelete.push_back( body.getBodyID() );
            body = removeAndNotify( me, bs, localStorage, body );

            // Note: Attached shadow copies are not deleted here. Instead we rely on the deferred deletion since we no
            // longer need the shadow copy: The owner of an attached shadow copy will receive a deletion notification, detach
            // the attachable, delete the shadow copy of the deleted body and send us a removal notification of the body
            // of which we own a shadow copy in the next position update since (probably) we no longer require the body but
            // are still part of its registration list.
            continue;
         } else
         {
            WALBERLA_ASSERT_UNEQUAL( owner.blockID_, me.blockID_, "Position is " << gpos << "\nlocal Block is: " << block.getAABB() );

            // New owner found among neighbors.
            WALBERLA_ASSERT_UNEQUAL( owner.blockID_, block.getId().getID(), "Migration is restricted to neighboring blocks." );

            WALBERLA_LOG_DETAIL( "Sending migration notification for body " << b->getSystemID() << " to process " << owner << "." );
            WALBERLA_LOG_DETAIL_SECTION()
            {
               std::stringstream log;
               log << "Process registration list before migration: [" ;
               for( auto it = b->MPITrait.beginShadowOwners(); it != b->MPITrait.endShadowOwners(); ++it ) {
                  if( it != b->MPITrait.beginShadowOwners() )
                     log << ", ";
                  log << (*it);
               }
               log << "]\n";
               WALBERLA_LOG_DETAIL(log.str());
            }

            // Set new owner and transform to shadow copy
            b->MPITrait.setOwner( owner );
            b->setRemote( true );

            // Move body to shadow copy storage.
            {
               auto pos = b->getPosition();
               correctBodyPosition(domain, block.getAABB().center(), pos);
               b->setPosition(pos);
            }
            shadowStorage.add( localStorage.release( body ) );

            // Note: We cannot determine here whether we require the body since we do not have up to date shadow copies.
            // However, we will most likely have to keep the body since it typically still intersects the process subdomain.

            // Correct registration list (exclude new owner and us - the old owner) and
            // notify registered processes (except for new owner) of (remote) migration since they possess a shadow copy.
            b->MPITrait.deregisterShadowOwner( owner );

            for( auto it = b->MPITrait.beginShadowOwners(); it != b->MPITrait.endShadowOwners(); ++it ) {
               mpi::SendBuffer& buffer( bs.sendBuffer(it->rank_) );

               WALBERLA_LOG_DETAIL( "Sending remote migration notification for body " << b->getSystemID() << " to process " << (*it) );

               packNotification(me, *it, buffer, RigidBodyRemoteMigrationNotification( *b, owner ));
            }

            // Send migration notification to new owner
            {
               b->MPITrait.registerShadowOwner( me );

               mpi::SendBuffer& buffer( bs.sendBuffer(owner.rank_) );

               packNotification(me, owner, buffer, RigidBodyMigrationNotification( *b ));

               // Note: The new owner misses shadow copies of all attached bodies. Since we do not have an up to date view
               // we need to relay them later.
            }

            b->MPITrait.clearShadowOwners();

            continue;
         }
      } else
      {
         // Body still is locally owned after position update.
         WALBERLA_LOG_DETAIL( "Owner of body " << b->getSystemID() << " is still process " << body->MPITrait.getOwner() );
      }

      ++body;
   }

   WALBERLA_LOG_DETAIL( "Assembling of body synchronization message ended." );
}

template <typename BodyTypeTuple>
void syncNextNeighbors( BlockForest& forest, BlockDataID storageID, WcTimingTree* tt = nullptr, const real_t dx = real_t(0), const bool syncNonCommunicatingBodies = false )
{
   if (tt != nullptr) tt->start("Sync");
   if (tt != nullptr) tt->start("Assembling Body Synchronization");
   mpi::BufferSystem bs( mpi::MPIManager::instance()->comm() );

   for (auto it = forest.begin(); it != forest.end(); ++it)
   {
      Block * block = dynamic_cast< Block * >( &(*it) );
      Storage* storage  = block->getData< Storage >( storageID );
      BodyStorage * localStorage  = &(*storage)[0];
      BodyStorage * shadowStorage = &(*storage)[1];

      for( uint_t i = uint_t(0); i != block->getNeighborhoodSize(); ++i )
      {
         auto neighborRank = block->getNeighborProcess(i);
         if (bs.sendBuffer(neighborRank).isEmpty())
         {
            // fill empty buffers with a dummy byte to force transmission
            bs.sendBuffer(neighborRank) << walberla::uint8_c(0);
         }
      }
      generateSynchonizationMessages<BodyTypeTuple>(bs, *block, *localStorage, *shadowStorage, dx, syncNonCommunicatingBodies);
   }
   if (tt != nullptr) tt->stop("Assembling Body Synchronization");

   // size of buffer is unknown and changes with each send
   bs.setReceiverInfoFromSendBufferState(false, true);
   bs.sendAll();

   if (tt != nullptr) tt->start("Parsing Body Synchronization");
   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of body synchronization response starts..." );
   for( auto it = bs.begin(); it != bs.end(); ++it )
   {
      walberla::uint8_t tmp;
      it.buffer() >> tmp;
      while( !it.buffer().isEmpty() )
      {
         IBlockID::IDType sender;
         IBlockID::IDType receiver;
         it.buffer() >> sender;
         it.buffer() >> receiver;
         auto blk = forest.getBlock(receiver);
         WALBERLA_CHECK(blk != nullptr, receiver << " not on this process!");
         IBlock& block = *blk;
         Storage* storage  = block.getData< Storage >( storageID );
         BodyStorage& localStorage  = (*storage)[0];
         BodyStorage& shadowStorage = (*storage)[1];
         pe::communication::parseMessage<BodyTypeTuple>(Owner(it.rank(), sender), it.buffer(), forest, block, localStorage, shadowStorage);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of body synchronization response ended." );
   if (tt != nullptr) tt->stop("Parsing Body Synchronization");
   if (tt != nullptr) tt->stop("Sync");
}

}  // namespace pe
}  // namespace walberla
