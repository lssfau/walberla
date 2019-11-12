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
//! \file InfoCollection.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "InfoCollection.h"

#include "pe/fcd/IFCD.h"
#include "pe/rigidbody/BodyStorage.h"

#include "blockforest/BlockID.h"
#include "core/mpi/BufferSystem.h"

namespace walberla {
namespace pe {

void createWithNeighborhoodLocalShadow( const BlockForest& bf,
                                        const BlockDataID storageID,
                                        blockforest::InfoCollection& ic )
{
   using namespace walberla::blockforest;
   ic.clear();

   mpi::BufferSystem bs( MPIManager::instance()->comm(), 756 );

   for (auto blockIt = bf.begin(); blockIt != bf.end(); ++blockIt)
   {
      const Block* block   = static_cast<const Block*> (&(*blockIt));
      Storage const *     storage       = block->getData< Storage >( storageID );
      BodyStorage const & localStorage  = (*storage)[StorageType::LOCAL];
      BodyStorage const & shadowStorage = (*storage)[StorageType::SHADOW];
      ic.insert( InfoCollectionPair(block->getId(), BlockInfo(localStorage.size(), shadowStorage.size())) );
      for( uint_t nb = uint_t(0); nb < block->getNeighborhoodSize(); ++nb )
      {
         bs.sendBuffer( block->getNeighborProcess(nb) ) << InfoCollectionPair(block->getId(), BlockInfo(localStorage.size(), shadowStorage.size()));
      }

      for (uint_t branchID = 0; branchID < 8; ++branchID)
      {
         const auto childID   = BlockID(block->getId(), branchID);
         const auto childAABB = bf.getAABBFromBlockId(childID);
         uint_t local  = 0;
         for (auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt)
         {
            if (childAABB.contains(bodyIt->getPosition()))
               ++local;
         }
         uint_t shadow  = 0;
         for (auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); ++bodyIt)
         {
            if (childAABB.contains(bodyIt->getPosition()))
               ++shadow;
         }
         ic.insert( InfoCollectionPair(childID, BlockInfo(local, shadow)) );

         for( uint_t nb = uint_t(0); nb < block->getNeighborhoodSize(); ++nb )
         {
            bs.sendBuffer( block->getNeighborProcess(nb) ) << InfoCollectionPair(childID, BlockInfo(local, shadow));
         }
      }
   }

   // size of buffer is unknown and changes with each send
   bs.setReceiverInfoFromSendBufferState(false, true);
   bs.sendAll();

   for( auto recvIt = bs.begin(); recvIt != bs.end(); ++recvIt )
   {
      while( !recvIt.buffer().isEmpty() )
      {
         InfoCollectionPair val;
         recvIt.buffer() >> val;
         ic.insert(val);
      }
   }
}

void createWithNeighborhoodContactsShadow( BlockForest& bf,
                                           const BlockDataID storageID,
                                           const BlockDataID fcdID,
                                           blockforest::InfoCollection& ic )
{
   using namespace walberla::blockforest;

   ic.clear();

   mpi::BufferSystem bs( MPIManager::instance()->comm(), 756 );

   for (auto blockIt = bf.begin(); blockIt != bf.end(); ++blockIt)
   {
      Block* block         = static_cast<Block*> (&(*blockIt));
      Storage const *     storage       = block->getData< Storage >( storageID );
      BodyStorage const & shadowStorage = (*storage)[StorageType::SHADOW];
      fcd::IFCD *     fcd         = block->getData< fcd::IFCD >( fcdID );
      ic.insert( InfoCollectionPair(block->getId(), BlockInfo(fcd->getContacts().size(), shadowStorage.size())) );
      for( uint_t nb = uint_t(0); nb < block->getNeighborhoodSize(); ++nb )
      {
         bs.sendBuffer( block->getNeighborProcess(nb) ) << InfoCollectionPair(block->getId(), BlockInfo(fcd->getContacts().size(), shadowStorage.size()));
      }

      for (uint_t branchID = 0; branchID < 8; ++branchID)
      {
         const auto childID   = BlockID(block->getId(), branchID);
         const auto childAABB = bf.getAABBFromBlockId(childID);
         uint_t localContacts  = 0;
         auto& contacts = fcd->getContacts();
         for (auto cIt = contacts.begin(); cIt != contacts.end(); ++cIt)
         {
            if (childAABB.contains(cIt->getPosition()))
               ++localContacts;
         }
         uint_t shadow  = 0;
         for (auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); ++bodyIt)
         {
            if (childAABB.contains(bodyIt->getPosition()))
               ++shadow;
         }
         ic.insert( InfoCollectionPair(childID, BlockInfo(localContacts, shadow)) );

         for( uint_t nb = uint_t(0); nb < block->getNeighborhoodSize(); ++nb )
         {
            bs.sendBuffer( block->getNeighborProcess(nb) ) << InfoCollectionPair(childID, BlockInfo(localContacts, shadow));
         }
      }
   }

   // size of buffer is unknown and changes with each send
   bs.setReceiverInfoFromSendBufferState(false, true);
   bs.sendAll();

   for( auto recvIt = bs.begin(); recvIt != bs.end(); ++recvIt )
   {
      while( !recvIt.buffer().isEmpty() )
      {
         InfoCollectionPair val;
         recvIt.buffer() >> val;
         ic.insert(val);
      }
   }
}

}
}
