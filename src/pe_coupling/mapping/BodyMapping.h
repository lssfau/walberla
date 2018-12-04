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
//! \file BodyMapping.h
//! \ingroup pe_coupling
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "BodyBBMapping.h"

#include "core/debug/Debug.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/FlagField.h"

#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/BodyIterators.h"

#include "pe_coupling/utility/BodySelectorFunctions.h"

#include <functional>

namespace walberla {
namespace pe_coupling {

// general mapping functions for a given single body on a given single block
template< typename BoundaryHandling_T >
void mapBody( pe::BodyID body, IBlock & block, StructuredBlockStorage & blockStorage,
              const BlockDataID & boundaryHandlingID, const FlagUID & obstacle )
{
   WALBERLA_ASSERT_EQUAL( &block.getBlockStorage(), &(blockStorage.getBlockStorage()) );

   BoundaryHandling_T * boundaryHandling = block.getData< BoundaryHandling_T >( boundaryHandlingID );
   auto * flagField = boundaryHandling->getFlagField();

   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );
   WALBERLA_ASSERT( flagField->flagExists( obstacle ) );

   const auto obstacleFlag = flagField->getFlag( obstacle );

   CellInterval cellBB = getCellBB( body, block, blockStorage, flagField->nrOfGhostLayers() );

   if( cellBB.empty() ) return;

   Vector3<real_t> startCellCenter = blockStorage.getBlockLocalCellCenter( block, cellBB.min() );
   const real_t dx = blockStorage.dx( blockStorage.getLevel(block) );
   const real_t dy = blockStorage.dy( blockStorage.getLevel(block) );
   const real_t dz = blockStorage.dz( blockStorage.getLevel(block) );

   real_t cz = startCellCenter[2];
   for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
   {
      real_t cy = startCellCenter[1];
      for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
      {
         real_t cx = startCellCenter[0];
         for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
         {
            if( body->containsPoint(cx,cy,cz) )
            {
               boundaryHandling->forceBoundary( obstacleFlag, x, y, z );
            }
            cx += dx;
         }
         cy += dy;
      }
      cz += dz;
   }
}



/*!\brief Mapping function to map all bodies - with certain properties - onto all blocks with the 'obstacle' flag
 *
 * All bodies (from bodyStorage and globalBodyStorage) are iterated and mapped to all blocks.
 * Cells that are inside the bodies are set to 'obstacle'.
 *
 * Whether or not a body is mapped depends on the return value of the 'mappingBodySelectorFct'.
 */
template< typename BoundaryHandling_T >
void mapBodies( StructuredBlockStorage & blockStorage, const BlockDataID & boundaryHandlingID,
                const BlockDataID & bodyStorageID, pe::BodyStorage & globalBodyStorage,
                const FlagUID & obstacle,
                const std::function<bool(pe::BodyID)> & mappingBodySelectorFct = selectAllBodies )
{
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {

      for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( mappingBodySelectorFct(bodyIt.getBodyID()))
            mapBody<BoundaryHandling_T>( bodyIt.getBodyID(), *blockIt, blockStorage, boundaryHandlingID, obstacle );
      }

      for( auto bodyIt = globalBodyStorage.begin(); bodyIt != globalBodyStorage.end(); ++bodyIt )
      {
         if( mappingBodySelectorFct(bodyIt.getBodyID()))
            mapBody< BoundaryHandling_T >( bodyIt.getBodyID(), *blockIt, blockStorage, boundaryHandlingID, obstacle );
      }
   }
}


} // namespace pe_coupling
} // namespace walberla
