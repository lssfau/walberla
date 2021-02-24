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

#include "mesa_pd/data/Flags.h"

#include <map>

namespace walberla {
namespace lbm_mesapd_coupling {
namespace amr {

using InfoCollection = std::map<blockforest::BlockID, BlockInfo>  ;
using InfoCollectionPair = std::pair<blockforest::BlockID, BlockInfo>;

template <typename BoundaryHandling_T, typename ParticleAccessor_T>
void updateAndSyncInfoCollection(BlockForest& bf, const BlockDataID boundaryHandlingID,
                                 const ParticleAccessor_T& accessor, const uint_t numberOfRPDSubCycles,
                                 InfoCollection& ic ) {
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

      // evaluate MESAPD quantities
      // count block local and (possible) ghost particles here
      uint_t numLocalParticles = 0, numGhostParticles = 0;
      for (size_t idx = 0; idx < accessor.size(); ++idx) {

         using namespace walberla::mesa_pd::data::particle_flags;
         if( isSet(accessor.getFlags(idx), GLOBAL) ) continue;

         if ( block->getAABB().contains(accessor.getPosition(idx)) ) {
            numLocalParticles++;
         } else if (block->getAABB().sqDistance(accessor.getPosition(idx)) < accessor.getInteractionRadius(idx)*accessor.getInteractionRadius(idx)) {
            numGhostParticles++;
         }
      }

      BlockInfo blockInfo(numCells, numFluidCells, numNearBoundaryCells, numLocalParticles, numGhostParticles, numberOfRPDSubCycles);
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

   // info collection has to be distirbuted to neighboring processes such that later on when coarsening was applied,
   // the weight of the coarsened block can be computed
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

void getBlockInfoFromInfoCollection( const PhantomBlock * block, const shared_ptr<InfoCollection>& ic,
                                     BlockInfo & blockInfo )
{
   WALBERLA_ASSERT_NOT_NULLPTR(block);

   if (block->sourceBlockIsLarger())
   {
      // block is a result of refinement -> BlockInfo object only available for the father block
      // there should be no particles on the block (otherwise it would not have been refined)
      // and refinement in LBM does not change the number of cells
      // we assume that the number of fluid and near boundary cells also stays the same
      // (ATTENTION: not true for blocks intersecting with a boundary!)
      // -> we can use the information of the father block for weight assignment

      auto infoIt = ic->find( block->getId().getFatherId() );
      WALBERLA_CHECK_UNEQUAL( infoIt, ic->end(), "Father block with ID " << block->getId().getFatherId() << " not found in info collection!" );

      // check the above mentioned assumptions
      WALBERLA_ASSERT_EQUAL(infoIt->second.numberOfLocalParticles, uint_t(0));

      blockInfo = infoIt->second;
   }
   else if (block->sourceBlockHasTheSameSize())
   {
      auto infoIt = ic->find( block->getId() );
      WALBERLA_CHECK_UNEQUAL( infoIt, ic->end(), "Block with ID " << block->getId() << " not found in info collection!" );
      blockInfo = infoIt->second;
   }
   else
   {
      // source block of block is smaller

      // block is a result of coarsening -> BlockInfo object is available on all 8 child blocks
      // there should be no particles on the block (otherwise it would not have been coarsened)
      // and refinement in LBM does not change the number of cells
      // we assume that the number of fluid and near boundary cells will be the average of all 8 child blocks
      // -> we can use the information of the child blocks for weight assignment

      blockforest::BlockID childIdForInit(block->getId(), 0);
      auto childForInitIt = ic->find( childIdForInit );
      WALBERLA_CHECK_UNEQUAL( childForInitIt, ic->end(), "Child block with ID " << childIdForInit << " not found in info collection!" );
      BlockInfo combinedInfo = childForInitIt->second;
      uint_t numFluidCells(0);
      uint_t numNearBoundaryCells(0);
      for (uint_t child = 0; child < 8; ++child)
      {
         blockforest::BlockID childId(block->getId(), child);
         auto childIt = ic->find( childId );
         WALBERLA_CHECK_UNEQUAL( childIt, ic->end(), "Child block with ID " << childId << " not found in info collection!" );
         numFluidCells += childIt->second.numberOfFluidCells;
         numNearBoundaryCells += childIt->second.numberOfNearBoundaryCells;

         // check above mentioned assumptions
         WALBERLA_ASSERT_EQUAL(childIt->second.numberOfLocalParticles, uint_t(0));
      }
      // total number of cells remains unchanged
      combinedInfo.numberOfFluidCells = uint_c(numFluidCells / uint_t(8)); //average
      combinedInfo.numberOfNearBoundaryCells = uint_c( numNearBoundaryCells / uint_t(8) ); //average
      combinedInfo.numberOfLocalParticles = uint_t(0);
      // number of rpd sub cycles stays the same

      blockInfo = combinedInfo;
   }
}

/*
 * Provides mapping from phantom block to weight evaluation via info collection
 */
class WeightEvaluationFunctor
{
public:
   WeightEvaluationFunctor(const shared_ptr<InfoCollection>& ic,
                           const std::function<real_t(const BlockInfo&)> & weightEvaluationFct) :
         ic_(ic), weightEvaluationFct_(weightEvaluationFct) {}

   real_t operator()(const PhantomBlock * block)
   {
      BlockInfo blockInfo;
      getBlockInfoFromInfoCollection(block, ic_, blockInfo);
      return weightEvaluationFct_(blockInfo);
   }

private:
   shared_ptr<InfoCollection> ic_;
   std::function<real_t(const BlockInfo&)> weightEvaluationFct_;
};


} // namepace amr
} // namespace lbm_mesapd_coupling
} // namespace walberla
