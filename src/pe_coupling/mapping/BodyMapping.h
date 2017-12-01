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
#include "field/iterators/IteratorMacros.h"

#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/BodyIterators.h"

namespace walberla {
namespace pe_coupling {

// general mapping functions for a given single body on a given single block
template< typename BoundaryHandling_T >
void mapBody( const pe::BodyID & body, IBlock & block, StructuredBlockStorage & blockStorage,
              BoundaryHandling_T * boundaryHandlingPtr, const FlagUID & obstacle )
{
   WALBERLA_ASSERT_EQUAL( &block.getBlockStorage(), &(blockStorage.getBlockStorage()) );
   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandlingPtr );

   auto * flagField = boundaryHandlingPtr->getFlagField();

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
               boundaryHandlingPtr->forceBoundary( obstacleFlag, x, y, z );
            cx += dx;
         }
         cy += dy;
      }
      cz += dz;
   }
}

template< typename BoundaryHandling_T >
void mapBody( const pe::BodyID & body, IBlock & block, StructuredBlockStorage & blockStorage,
              const BlockDataID & boundaryHandlingID, const FlagUID & obstacle )
{
   WALBERLA_ASSERT_EQUAL( &block.getBlockStorage(), &(blockStorage.getBlockStorage()) );

   BoundaryHandling_T * boundaryHandling = block.getData< BoundaryHandling_T >( boundaryHandlingID );

   mapBody(body, block, blockStorage, boundaryHandling, obstacle );
}



// mapping function to map all bodies from the body storage - with certain properties - onto all blocks
template< typename BoundaryHandling_T >
void mapBodies( StructuredBlockStorage & blockStorage, const BlockDataID & boundaryHandlingID,
                const BlockDataID & bodyStorageID, const FlagUID & obstacle,
                const bool fixedBodiesOnly = true, const bool moBodiesOnly = true )
{
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
      for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         WALBERLA_UNUSED(moBodiesOnly); // undo when other coupling algorithms are available

         if( fixedBodiesOnly && bodyIt->isFixed() )
            continue;

         mapBody<BoundaryHandling_T>( *bodyIt, *blockIt, blockStorage, boundaryHandlingID, obstacle );
      }
   }
}


// mapping function to map all global bodies - with certain properties - onto all blocks
template< typename BoundaryHandling_T >
void mapGlobalBodies( StructuredBlockStorage & blockStorage, const BlockDataID & boundaryHandlingID,
                      pe::BodyStorage & globalBodyStorage, const FlagUID & obstacle,
                      const bool fixedBodiesOnly = true, const bool moBodiesOnly = true )
{
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
      for( auto bodyIt = globalBodyStorage.begin(); bodyIt != globalBodyStorage.end(); ++bodyIt )
      {
         WALBERLA_UNUSED(moBodiesOnly); // undo when other coupling algorithms are available

         if( fixedBodiesOnly && bodyIt->isFixed() )
            continue;

         mapBody< BoundaryHandling_T >( *bodyIt, *blockIt, blockStorage, boundaryHandlingID, obstacle );
      }
   }
}

// mapping function to map a given single global body onto all blocks
template< typename BoundaryHandling_T >
void mapGlobalBody( const id_t globalBodySystemID,
                    StructuredBlockStorage & blockStorage, const BlockDataID & boundaryHandlingID,
                    pe::BodyStorage & globalBodyStorage, const FlagUID & obstacle )
{
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
      auto bodyIt = globalBodyStorage.find( globalBodySystemID );
      if( bodyIt != globalBodyStorage.end() )
      {
         mapBody< BoundaryHandling_T >( *bodyIt, *blockIt, blockStorage, boundaryHandlingID, obstacle );
      }
   }
}


// mapping function to map all global bodies - with certain properties - onto a given single block
template< typename BoundaryHandling_T >
void mapGlobalBodiesOnBlock( IBlock & block,
                             StructuredBlockStorage & blockStorage, BoundaryHandling_T * boundaryHandlingPtr,
                             pe::BodyStorage & globalBodyStorage, const FlagUID & obstacle,
                             const bool fixedBodiesOnly = true, const bool moBodiesOnly = true )
{
   for( auto bodyIt = globalBodyStorage.begin(); bodyIt != globalBodyStorage.end(); ++bodyIt)
   {
      WALBERLA_UNUSED(moBodiesOnly); // undo when other coupling algorithms are available

      if( fixedBodiesOnly && bodyIt->isFixed() )
         continue;

      mapBody< BoundaryHandling_T >( *bodyIt, block, blockStorage, boundaryHandlingPtr, obstacle );
   }
}


// mapping function to map a given single global body onto a given single block
template< typename BoundaryHandling_T >
void mapGlobalBodyOnBlock( const id_t globalBodySystemID, IBlock & block,
                           StructuredBlockStorage & blockStorage, BoundaryHandling_T * boundaryHandlingPtr,
                           pe::BodyStorage & globalBodyStorage, const FlagUID & obstacle )
{
   auto bodyIt = globalBodyStorage.find( globalBodySystemID );
   if( bodyIt != globalBodyStorage.end() )
   {
      mapBody< BoundaryHandling_T >( *bodyIt, block, blockStorage, boundaryHandlingPtr, obstacle );
   }

}


} // namespace pe_coupling
} // namespace walberla
