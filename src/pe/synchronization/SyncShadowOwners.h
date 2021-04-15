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
//! \file SyncShadowOwners.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <blockforest/BlockForest.h>
#include <core/DataTypes.h>
#include <core/mpi/BufferSystem.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/Reduce.h>
#include <core/mpi/SendBuffer.h>

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

#include "core/timing/TimingTree.h"

namespace walberla {
namespace pe {

template <typename BodyTypeTuple>
void updateAndMigrate( BlockForest& forest, BlockDataID storageID, const bool syncNonCommunicatingBodies )
{
   using namespace walberla::pe::communication;
   //==========================================================
   // STEP1: Update & Migrate
   //==========================================================
   mpi::BufferSystem bs( mpi::MPIManager::instance()->comm(),  123);

   WALBERLA_LOG_DETAIL( "Assembling of Update&Migrate starts..." );
   std::set<mpi::MPIRank> recvRanks; // potential message senders
   for (auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt)
   {
      Block& block = *(dynamic_cast< Block * >( &(*blockIt) ));
      Storage* storage  = block.getData< Storage >( storageID );
      BodyStorage& localStorage  = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      const Owner me( int_c( mpi::MPIManager::instance()->rank() ), block.getId().getID() );
      const math::AABB& blkAABB = block.getAABB();

      // schedule receives
      for( auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); ++bodyIt)
      {
         if (bodyIt->isCommunicating() || syncNonCommunicatingBodies)
         {
            recvRanks.insert(bodyIt->MPITrait.getOwner().rank_);
         }
      }

      for( auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); )
      {
         BodyID b (bodyIt.getBodyID());

         //correct position to make sure body is always inside the domain!
         if (!b->isFixed())
         {
            auto pos = b->getPosition();
            block.getBlockStorage().mapToPeriodicDomain(pos);
            b->setPosition(pos);
         }

         if( !b->isCommunicating() && !syncNonCommunicatingBodies ) {
            ++bodyIt;
            continue;
         }

         if (bodyIt->isMarkedForDeletion())
         {
            // delete it
            WALBERLA_LOG_DETAIL( "Sending deletion notifications for body " << bodyIt->getSystemID() << " due to manual deletion." );

            bodyIt = removeAndNotify( me, bs, localStorage, bodyIt );

            // Note: Attached shadow copies are not deleted here. Instead we rely on the deferred deletion since we no
            // longer need the shadow copy: The owner of an attached shadow copy will receive a deletion notification, detach
            // the attachable, delete the shadow copy of the deleted body and send us a removal notification of the body
            // of which we own a shadow copy in the next position update since (probably) we no longer require the body but
            // are still part of its registration list.
            continue;
         }

         // Update
         for (auto it = b->MPITrait.beginShadowOwners(); it != b->MPITrait.endShadowOwners(); ++it)
         {
            WALBERLA_LOG_DETAIL( "Sending update notification for body " << b->getSystemID() << " to process " << (*it) );
            mpi::SendBuffer& sb = bs.sendBuffer(it->rank_);
            if (sb.isEmpty()) sb << walberla::uint8_c(0);
            packNotification(me, *it, sb, RigidBodyUpdateNotification( *b ));
         }
         if (!blkAABB.contains( b->getPosition() ))
         {
            Owner owner( findContainingProcess( block, b->getPosition() ) );
            if ( owner.rank_ < 0)
            {
               // No owner found: Outflow condition.
               WALBERLA_LOG_DETAIL( "Sending deletion notifications for body " << bodyIt->getSystemID() << " due to outflow." );

               //delete body
               bodyIt = removeAndNotify(me, bs, localStorage, bodyIt);

               continue;
            } else
            {
               // Migrate
               WALBERLA_ASSERT( owner.blockID_ != block.getId().getID(), "Migration is restricted to neighboring blocks." );

               // Set new owner and transform to shadow copy
               b->MPITrait.setOwner( owner );
               b->setRemote( true );

               // Move body to shadow copy storage.
               {
                  auto pos = b->getPosition();
                  correctBodyPosition(forest.getDomain(), block.getAABB().center(), pos);
                  b->setPosition(pos);
               }
               shadowStorage.add( localStorage.release( bodyIt ) );

               b->MPITrait.deregisterShadowOwner( owner );

               // Send remote migration notifications
               for( auto it = b->MPITrait.beginShadowOwners(); it != b->MPITrait.endShadowOwners(); ++it ) {
                  WALBERLA_LOG_DETAIL( "Sending remote migration notification for body " << b->getSystemID() << " to process " << (*it) );
                  mpi::SendBuffer& sb( bs.sendBuffer(it->rank_) );
                  if (sb.isEmpty()) sb << walberla::uint8_c(0);
                  packNotification(me, *it, sb, RigidBodyRemoteMigrationNotification( *b, owner ));
               }

               // Send migration notification to new owner
               {
                  b->MPITrait.registerShadowOwner( me );
                  WALBERLA_LOG_DETAIL( "Sending migration notification for body " << b->getSystemID() << " to process " << (owner) );
                  mpi::SendBuffer& sb( bs.sendBuffer(owner.rank_) );
                  if (sb.isEmpty()) sb << walberla::uint8_c(0);
                  packNotification(me, owner, sb, RigidBodyMigrationNotification( *b ));
               }
               b->MPITrait.clearShadowOwners();
               continue;
            }
         }
         ++bodyIt;
      }
   }
   WALBERLA_LOG_DETAIL( "Assembling of Update&Migrate ended." );

   WALBERLA_LOG_DETAIL( "UM: number of recv " << recvRanks.size());
   bs.setReceiverInfo(recvRanks, true);
   bs.sendAll();

   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of Update&Migrate starts..." );
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
         Block * block = dynamic_cast< Block * >( blk );
         Storage* storage  = block->getData< Storage >( storageID );
         BodyStorage& localStorage  = (*storage)[0];
         BodyStorage& shadowStorage = (*storage)[1];
         pe::communication::parseMessage<BodyTypeTuple>(Owner(it.rank(), sender), it.buffer(), forest, *block, localStorage, shadowStorage);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of Update&Migrate ended." );
}

template <typename BodyTypeTuple>
void checkAndResolveOverlap( BlockForest& forest, BlockDataID storageID, const real_t dx, const bool syncNonCommunicatingBodies )
{
   using namespace walberla::pe::communication;
   //==========================================================
   // STEP2: Check&Resolve
   //==========================================================
   mpi::BufferSystem bs( mpi::MPIManager::instance()->comm(), 124 );

   //init buffers
   const std::vector<uint_t>& nbProcesses = forest.getNeighboringProcesses();
   for( auto it = nbProcesses.begin(); it != nbProcesses.end(); ++it)
   {
      bs.sendBuffer( *it ) << walberla::uint8_c(0);
   }
   bs.sendBuffer( forest.getProcess() ) << walberla::uint8_c(0);

   WALBERLA_LOG_DETAIL( "Assembling of Check&Resolve starts..." );
   std::set<mpi::MPIRank> recvRanks; // potential message senders
   for (auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt)
   {
      Block& block = *(dynamic_cast< Block * >( &(*blockIt) ));
      Storage* storage  = block.getData< Storage >( storageID );
      BodyStorage& localStorage  = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      const Owner me( int_c(block.getProcess()), block.getId().getID() );
      //      const math::AABB& blkAABB = block->getAABB();

      for( auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt)
      {
         BodyID b (bodyIt.getBodyID());

         if( !b->isCommunicating() && !syncNonCommunicatingBodies ) continue;

         bool isInsideDomain = forest.getDomain().contains( b->getAABB(), -dx );

         mpi::SendBuffer& sbMaster = bs.sendBuffer(b->MPITrait.getOwner().rank_);
         if (sbMaster.isEmpty()) sbMaster << walberla::uint8_c(0);

         // Update (nearest) neighbor processes.
         for( uint_t nb = uint_t(0); nb < block.getNeighborhoodSize(); ++nb )
         {
            Owner nbProcess = Owner(int_c( block.getNeighborProcess(nb) ), block.getNeighborId(nb).getID() );

            mpi::SendBuffer& sb = bs.sendBuffer(nbProcess.rank_);
            if (sb.isEmpty()) sb << walberla::uint8_c(0);

            if (b->MPITrait.getOwner() == nbProcess) continue; // dont send to owner!!
            if (b->MPITrait.getBlockState( nbProcess.blockID_ )) continue; // only send to neighbor which do not know this body
            //            WALBERLA_LOG_DEVEL("neighbor aabb: " << block.getNeighborAABB(nb));
            //            WALBERLA_LOG_DEVEL("isInsideDomain: " << isInsideDomain);
            //            WALBERLA_LOG_DEVEL("body AABB: " << b->getAABB());
            //            WALBERLA_LOG_DEVEL("neighbor AABB: " << block.getNeighborAABB(nb));

            if( (isInsideDomain ? block.getNeighborAABB(nb).intersects( b->getAABB(), dx ) : block.getBlockStorage().periodicIntersect(block.getNeighborAABB(nb), b->getAABB(), dx)) )
            {
               WALBERLA_LOG_DETAIL( "Sending copy notification for body " << b->getSystemID() << " to process " << (nbProcess) << "\n master: " << b->MPITrait.getOwner());
               packNotification(me, nbProcess, sb, RigidBodyCopyNotification( *b ));
               MarshalDynamically<BodyTypeTuple>::execute( sb, *b );
               packNotification(me, b->MPITrait.getOwner(), sbMaster, RigidBodyNewShadowCopyNotification( *b, nbProcess ));
               b->MPITrait.setBlockState( nbProcess.blockID_ );
            }
         }
      }
      for( auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); )
      {
         BodyID b (bodyIt.getBodyID());
         WALBERLA_ASSERT(!b->isGlobal(), "Global body in ShadowStorage!");
         bool isInsideDomain = forest.getDomain().contains( b->getAABB(), -dx );

         mpi::SendBuffer& sbMaster = bs.sendBuffer(b->MPITrait.getOwner().rank_);
         if (sbMaster.isEmpty()) sbMaster << walberla::uint8_c(0);

         // Update (nearest) neighbor processes.
         for( uint_t nb = uint_t(0); nb < block.getNeighborhoodSize(); ++nb )
         {
            Owner nbProcess = Owner(int_c( block.getNeighborProcess(nb) ), block.getNeighborId(nb).getID() );

            mpi::SendBuffer& sb = bs.sendBuffer(nbProcess.rank_);
            if (sb.isEmpty()) sb << walberla::uint8_c(0);

            if (b->MPITrait.getOwner() == nbProcess) continue; // dont send to owner!!
            if (b->MPITrait.getBlockState( nbProcess.blockID_ )) continue; // only send to neighbor which do not know this body
            if( (isInsideDomain ? block.getNeighborAABB(nb).intersects( b->getAABB(), dx ) : block.getBlockStorage().periodicIntersect(block.getNeighborAABB(nb), b->getAABB(), dx)) )
            {
               WALBERLA_LOG_DETAIL( "Sending copy notification for body " << b->getSystemID() << " to process " << (nbProcess) );
               packNotification(me, nbProcess, sb, RigidBodyCopyNotification( *b ));
               MarshalDynamically<BodyTypeTuple>::execute( sb, *b );
               packNotification(me, b->MPITrait.getOwner(), sbMaster, RigidBodyNewShadowCopyNotification( *b, nbProcess ));
               b->MPITrait.setBlockState( nbProcess.blockID_ );
            }
         }

         if ( !block.getAABB().intersects(b->getAABB(), dx) )
         {
            // Delete
            // inform nearest neighbor processes.
            for( uint_t nb = uint_t(0); nb < block.getNeighborhoodSize(); ++nb )
            {
               Owner nbProcess = Owner(int_c( block.getNeighborProcess(nb) ), block.getNeighborId(nb).getID() );

               WALBERLA_LOG_DETAIL( "Sending removal information notification for body " << b->getSystemID() << " to process " << (nbProcess) );
               mpi::SendBuffer& sb = bs.sendBuffer(nbProcess.rank_);
               if (sb.isEmpty()) sb << walberla::uint8_c(0);
               packNotification(me, nbProcess, sb, RigidBodyRemovalInformationNotification( *b ));
            }

            //notify owner
            WALBERLA_LOG_DETAIL( "Sending removal information notification for body " << b->getSystemID() << " to process " << (b->MPITrait.getOwner()) );
            mpi::SendBuffer& sb = bs.sendBuffer(b->MPITrait.getOwner().rank_);
            if (sb.isEmpty()) sb << walberla::uint8_c(0);
            packNotification(me, b->MPITrait.getOwner(), sb, RigidBodyRemovalInformationNotification( *b ));

            bodyIt = shadowStorage.remove( bodyIt );

            continue;
         }
         ++bodyIt;
      }
      // schedule receives
      for( auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt)
      {
         if (!bodyIt->isCommunicating() && !syncNonCommunicatingBodies) continue;
         for (auto it = bodyIt->MPITrait.beginShadowOwners(); it != bodyIt->MPITrait.endShadowOwners(); ++it)
         {
            recvRanks.insert(it->rank_);
         }
      }
   }

   for( auto it = nbProcesses.begin(); it != nbProcesses.end(); ++it)
   {
      recvRanks.insert( int_c(*it) );
   }
   recvRanks.insert( int_c(forest.getProcess()) );
   WALBERLA_LOG_DETAIL( "Assembling of Check&Resolve ended." );

   // size of buffer is unknown and changes with each send
   WALBERLA_LOG_DETAIL( "CR: number of recv " << recvRanks.size());
   bs.setReceiverInfo(recvRanks, true);
   bs.sendAll();

   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of Check&Resolve starts..." );
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
         Block * block = dynamic_cast< Block * >( blk );
         Storage* storage  = block->getData< Storage >( storageID );
         BodyStorage& localStorage  = (*storage)[0];
         BodyStorage& shadowStorage = (*storage)[1];
         pe::communication::parseMessage<BodyTypeTuple>(Owner(it.rank(), sender), it.buffer(), forest, *block, localStorage, shadowStorage);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of Check&Resolve ended." );
}

template <typename BodyTypeTuple>
void syncShadowOwners( BlockForest& forest, BlockDataID storageID, WcTimingTree* tt = nullptr, const real_t dx = real_t(0), const bool syncNonCommunicatingBodies = false )
{
   if (tt != nullptr) tt->start("Sync");

   //==========================================================
   // STEP1: Update & Migrate
   //==========================================================
   if (tt != nullptr) tt->start("Update&Migrate");
   updateAndMigrate<BodyTypeTuple>( forest, storageID, syncNonCommunicatingBodies);
   if (tt != nullptr) tt->stop("Update&Migrate");

   //==========================================================
   // STEP2: Check & Resolve
   //==========================================================
   if (tt != nullptr) tt->start("Check&Resolve");
   checkAndResolveOverlap<BodyTypeTuple>( forest, storageID, dx, syncNonCommunicatingBodies);
   if (tt != nullptr) tt->stop("Check&Resolve");

   if (tt != nullptr) tt->stop("Sync");
}

}
}
