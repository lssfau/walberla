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
//! \file InfoCollection.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockInfo.h"

#include "blockforest/BlockForest.h"
#include "blockforest/BlockID.h"
#include "core/mpi/BufferSystem.h"

#include "pe/ccd/ICCD.h"
#include "pe/fcd/IFCD.h"
#include "pe/rigidbody/BodyStorage.h"

#include <map>

namespace walberla {
namespace pe_coupling {

typedef std::map<blockforest::BlockID, BlockInfo>  InfoCollection;
typedef std::pair<blockforest::BlockID, BlockInfo> InfoCollectionPair;

template <typename BoundaryHandling_T>
void createWithNeighborhood(BlockForest& bf, const BlockDataID boundaryHandlingID,
                            const BlockDataID bodyStorageID, const BlockDataID ccdID, const BlockDataID fcdID,
                            const uint_t numberOfPeSubCycles,
                            InfoCollection& ic )
{
   ic.clear();

   mpi::BufferSystem bs( MPIManager::instance()->comm(), 856 );

   for (auto blockIt = bf.begin(); blockIt != bf.end(); ++blockIt)
   {
      auto * block = static_cast<blockforest::Block*> (&(*blockIt));

      // evaluate LBM quantities
      BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID );
      auto xyzSize = boundaryHandling->getFlagField()->xyzSize();
      const uint_t numCells = xyzSize.numCells();
      uint_t numFluidCells(0), numNearBoundaryCells(0);
      for( auto cellIt = xyzSize.begin(); cellIt != xyzSize.end(); ++cellIt)
      {
         if( boundaryHandling->isDomain(*cellIt) )
         {
            ++numFluidCells;
         }
         if( boundaryHandling->isNearBoundary(*cellIt))
         {
            ++numNearBoundaryCells;
         }
      }

      // evaluate PE quantities
      auto * bodyStorage = block->getData<pe::Storage>(bodyStorageID);
      pe::BodyStorage const & localStorage  = (*bodyStorage)[pe::StorageType::LOCAL];
      pe::BodyStorage const & shadowStorage = (*bodyStorage)[pe::StorageType::SHADOW];
      const uint_t numLocalParticles = localStorage.size();
      const uint_t numShadowParticles = shadowStorage.size();

      auto * ccd = block->getData<pe::ccd::ICCD>(ccdID);
      auto * fcd = block->getData<pe::fcd::IFCD>(fcdID);
      ccd->generatePossibleContacts();
      pe::Contacts& contacts = fcd->generateContacts( ccd->getPossibleContacts() );
      const uint_t numContacts = contacts.size();

      BlockInfo blockInfo(numCells, numFluidCells, numNearBoundaryCells, numLocalParticles, numShadowParticles, numContacts, numberOfPeSubCycles);
      InfoCollectionPair infoCollectionEntry(block->getId(), blockInfo);

      ic.insert( infoCollectionEntry );

      for( auto nb = uint_t(0); nb < block->getNeighborhoodSize(); ++nb )
      {
         bs.sendBuffer( block->getNeighborProcess(nb) ) << infoCollectionEntry;
      }

      //note: is it necessary to add child blocks already into the info collection?
      // here, we still have full geometrical information and can probably determine number of fluid and near boundary cells more easily
      // however, the interesting (and most costly) blocks are never refined and thus their child infos is never needed
      // see pe/amr/InfoCollection.cpp for an example

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

void getBlockInfoFromInfoCollection( const PhantomBlock * block, const shared_ptr<InfoCollection>& ic, BlockInfo & blockInfo );


} // namespace pe_coupling
} // namespace walberla
