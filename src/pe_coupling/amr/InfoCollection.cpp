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


#include "InfoCollection.h"

namespace walberla {
namespace pe_coupling {


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
      WALBERLA_ASSERT_EQUAL(infoIt->second.numberOfLocalBodies, uint_t(0));
      WALBERLA_ASSERT_EQUAL(infoIt->second.numberOfShadowBodies, uint_t(0));
      WALBERLA_ASSERT_EQUAL(infoIt->second.numberOfContacts, uint_t(0));

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
      uint_t numFluidCells(0), numNearBoundaryCells(0);
      for (uint_t child = 0; child < 8; ++child)
      {
         blockforest::BlockID childId(block->getId(), child);
         auto childIt = ic->find( childId );
         WALBERLA_CHECK_UNEQUAL( childIt, ic->end(), "Child block with ID " << childId << " not found in info collection!" );
         numFluidCells += childIt->second.numberOfFluidCells;
         numNearBoundaryCells += childIt->second.numberOfNearBoundaryCells;

         // check above mentioned assumptions
         WALBERLA_ASSERT_EQUAL(childIt->second.numberOfLocalBodies, uint_t(0));
         WALBERLA_ASSERT_EQUAL(childIt->second.numberOfShadowBodies, uint_t(0));
         WALBERLA_ASSERT_EQUAL(childIt->second.numberOfContacts, uint_t(0));
      }
      // total number of cells remains unchanged
      combinedInfo.numberOfFluidCells = uint_c(numFluidCells / uint_t(8)); //average
      combinedInfo.numberOfNearBoundaryCells = uint_c( numNearBoundaryCells / uint_t(8) ); //average
      combinedInfo.numberOfLocalBodies = uint_t(0);
      combinedInfo.numberOfShadowBodies = uint_t(0);
      combinedInfo.numberOfContacts = uint_t(0); //sum
      // number of pe sub cycles stays the same

      blockInfo = combinedInfo;
   }

}

} // namespace pe_coupling
} // namespace walberla
