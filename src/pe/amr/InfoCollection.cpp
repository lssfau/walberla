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

#include "pe/rigidbody/BodyStorage.h"

#include "blockforest/BlockID.h"
#include "core/mpi/BufferSystem.h"

namespace walberla {
namespace pe {

void createWithNeighborhood(const BlockForest& bf, const BlockDataID storageID, InfoCollection& ic )
{
   ic.clear();

   mpi::BufferSystem bs( MPIManager::instance()->comm(), 756 );

   for (auto blockIt = bf.begin(); blockIt != bf.end(); ++blockIt)
   {
      const blockforest::Block* block   = static_cast<const blockforest::Block*> (&(*blockIt));
      Storage const *     storage       = block->getData< Storage >( storageID );
      BodyStorage const & localStorage  = (*storage)[StorageType::LOCAL];
      BodyStorage const & shadowStorage = (*storage)[StorageType::SHADOW];
      ic.insert( InfoCollection::value_type(block->getId(), BlockInfo(localStorage.size(), shadowStorage.size())) );

      for( uint_t nb = uint_t(0); nb < block->getNeighborhoodSize(); ++nb )
      {
         bs.sendBuffer( block->getNeighborProcess(nb) ) << InfoCollection::value_type(block->getId(), BlockInfo(localStorage.size(), shadowStorage.size()));
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
