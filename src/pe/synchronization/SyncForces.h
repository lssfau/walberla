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
//! \file SynchronizeForces.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>
#include <core/mpi/BufferSystem.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/Reduce.h>
#include <core/mpi/SendBuffer.h>
#include <domain_decomposition/BlockStorage.h>

#include "pe/BlockFunctions.h"
#include "pe/rigidbody/BodyStorage.h"

#include "pe/communication/Marshalling.h"
#include "pe/communication/RigidBodyForceNotification.h"
#include "pe/communication/PackNotification.h"
#include "pe/communication/ParseMessage.h"

namespace walberla {
namespace pe {

inline
void reduceForces( BlockStorage& blocks, BlockDataID storageID )
{
   using namespace walberla::pe::communication;

   mpi::BufferSystem bs( mpi::MPIManager::instance()->comm(), 126 );
   std::set<mpi::MPIRank> recvRanks; // potential message senders

   // Sending local force contributions of shadow copies to owner.
   WALBERLA_LOG_DETAIL( "Assembling of force reduction message starts...");
   for (auto it = blocks.begin(); it != blocks.end(); ++it)
   {
      Block * block = dynamic_cast< Block * >( &(*it) );
      Storage* storage  = block->getData< Storage >( storageID );
      BodyStorage * localStorage  = &(*storage)[0];
      BodyStorage * shadowStorage = &(*storage)[1];

      const Owner me( int_c( mpi::MPIManager::instance()->rank() ), block->getId().getID() );

      // pack messages
      for( auto bodyIt = shadowStorage->begin(); bodyIt != shadowStorage->end(); ++bodyIt)
      {
         WALBERLA_LOG_DETAIL( "__Processing shadow copy of body " << bodyIt->getSystemID() << ".\n");

         auto neighborRank = bodyIt->MPITrait.getOwner().rank_;
         mpi::SendBuffer& sb = bs.sendBuffer(neighborRank);
         if (sb.isEmpty())
         {
            // fill empty buffers with a dummy byte to force transmission
            bs.sendBuffer(neighborRank) << walberla::uint8_c(0);
         }

         if( !bodyIt->hasForce() ) {
            // If we did not apply any forces do not send anything.
            continue;
         }

         WALBERLA_LOG_DETAIL( "__Sending force contribution " << bodyIt->getForce() << ", " << bodyIt->getTorque() << " of body " <<
                              bodyIt->getSystemID() << " to owner block " << bodyIt->MPITrait.getOwner().blockID_ << ".\n");
         packNotification(me, bodyIt->MPITrait.getOwner(), sb, RigidBodyForceNotification( *bodyIt ) );

      }

      // schedule receives
      for( auto bodyIt = localStorage->begin(); bodyIt != localStorage->end(); ++bodyIt)
      {
         for (auto sownerIt = bodyIt->MPITrait.beginShadowOwners(); sownerIt != bodyIt->MPITrait.endShadowOwners(); ++sownerIt)
         {
            recvRanks.insert(sownerIt->rank_);
         }
      }
   }

   WALBERLA_LOG_DETAIL( "Assembling of force reduction message ended." );

   WALBERLA_LOG_DETAIL( "FS: number of recv " << recvRanks.size());
   bs.setReceiverInfo(recvRanks, true);
   bs.sendAll();

   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of force reduction message starts..." );
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
         auto blk = blocks.getBlock(receiver);
         WALBERLA_CHECK(blk != nullptr, receiver << " not on this process!");
         IBlock& block = *blk;
         Storage* storage  = block.getData< Storage >( storageID );
         BodyStorage& localStorage  = (*storage)[0];
         BodyStorage& shadowStorage = (*storage)[1];
         pe::communication::parseForceReduceMessage(Owner(it.rank(), sender), it.buffer(), blocks, block, localStorage, shadowStorage);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of force reduction message ended." );
}

inline
void reduceForces( BlockStorage& blocks, BlockDataID storageID, BodyStorage& globalBodyStorage )
{
   using namespace walberla::pe::communication;

   reduceForces(blocks, storageID);

   WALBERLA_LOG_DETAIL( "Sync force on global bodies..." );
   {
      size_t i;
      std::vector<real_t> reductionBuffer( globalBodyStorage.size() * 6, real_c(0) );

      i = 0;
      for( auto it = globalBodyStorage.begin(); it != globalBodyStorage.end(); ++it )
      {
         if (it->hasInfiniteMass()) continue;

         const Vec3 f( it->getForce() ), tau( it->getTorque() );
         reductionBuffer[i++] = f[0];
         reductionBuffer[i++] = f[1];
         reductionBuffer[i++] = f[2];
         reductionBuffer[i++] = tau[0];
         reductionBuffer[i++] = tau[1];
         reductionBuffer[i++] = tau[2];
      }

      mpi::allReduceInplace(reductionBuffer, mpi::SUM);

      i = 0;
      for( auto it = globalBodyStorage.begin(); it != globalBodyStorage.end(); ++it )
      {
         if (it->hasInfiniteMass()) continue;
         it->setForce ( Vec3( reductionBuffer[i], reductionBuffer[i + 1], reductionBuffer[i + 2] ) );
         it->setTorque( Vec3( reductionBuffer[i + 3], reductionBuffer[i + 4], reductionBuffer[i + 5] ) );
      }
   }
   WALBERLA_LOG_DETAIL( "Sync force on global bodies finished." );
}

inline
void distributeForces( BlockStorage& blocks, BlockDataID storageID )
{
   using namespace walberla::pe::communication;

   mpi::BufferSystem bs( mpi::MPIManager::instance()->comm(), 127 );
   std::set<mpi::MPIRank> recvRanks; // potential message senders

   // distributing forces from owner to shadow copies
   WALBERLA_LOG_DETAIL( "Assembling of force distribution message starts...");
   for (auto it = blocks.begin(); it != blocks.end(); ++it)
   {
      Block * block = dynamic_cast< Block * >( &(*it) );
      Storage* storage  = block->getData< Storage >( storageID );
      BodyStorage * localStorage  = &(*storage)[0];
      BodyStorage * shadowStorage = &(*storage)[1];

      const Owner me( int_c( mpi::MPIManager::instance()->rank() ), block->getId().getID() );

      // pack messages
      for( auto bodyIt = localStorage->begin(); bodyIt != localStorage->end(); ++bodyIt)
      {
         WALBERLA_LOG_DETAIL( "__Processing local body " << bodyIt->getSystemID() << ".\n");

         for (auto sownerIt = bodyIt->MPITrait.beginShadowOwners(); sownerIt != bodyIt->MPITrait.endShadowOwners(); ++sownerIt)
         {
            auto neighborRank = sownerIt->rank_;
            mpi::SendBuffer& sb = bs.sendBuffer(neighborRank);
            if (sb.isEmpty())
            {
               // fill empty buffers with a dummy byte to force transmission
               bs.sendBuffer(neighborRank) << walberla::uint8_c(0);
            }

            if( !bodyIt->hasForce() ) {
               // If we did not apply any forces do not send anything.
               continue;
            }

            WALBERLA_LOG_DETAIL( "__Sending force contribution " << bodyIt->getForce() << ", " << bodyIt->getTorque() << " of body " <<
                                 bodyIt->getSystemID() << " to shadow owner " << sownerIt->blockID_ << ".\n");
            packNotification(me, *sownerIt, sb, RigidBodyForceNotification( *bodyIt ) );
         }
      }

      // schedule receives
      for( auto bodyIt = shadowStorage->begin(); bodyIt != shadowStorage->end(); ++bodyIt)
      {
         recvRanks.insert( bodyIt->MPITrait.getOwner().rank_ );
      }
   }

   WALBERLA_LOG_DETAIL( "Assembling of force distribution message ended." );

   WALBERLA_LOG_DETAIL( "DistributeForce: number of recv " << recvRanks.size());
   bs.setReceiverInfo(recvRanks, true);
   bs.sendAll();

   // Receiving the updates from the owner processes of shadow copies
   WALBERLA_LOG_DETAIL( "Parsing of force distribution message starts..." );
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
         auto blk = blocks.getBlock(receiver);
         WALBERLA_CHECK(blk != nullptr, receiver << " not on this process!");
         IBlock& block = *blk;
         Storage* storage  = block.getData< Storage >( storageID );
         BodyStorage& localStorage  = (*storage)[0];
         BodyStorage& shadowStorage = (*storage)[1];
         pe::communication::parseForceDistributeMessage(Owner(it.rank(), sender), it.buffer(), blocks, block, localStorage, shadowStorage);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of force synchronization message ended." );
}

}  // namespace pe
}  // namespace walberla
