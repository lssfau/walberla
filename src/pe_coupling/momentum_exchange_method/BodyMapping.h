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

#include "core/debug/Debug.h"
#include "core/cell/Cell.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/FlagField.h"
#include "field/iterators/IteratorMacros.h"

#include "pe/rigidbody/BodyIterators.h"

#include "pe_coupling/mapping/BodyBBMapping.h"

namespace walberla {
namespace pe_coupling {

template< typename BoundaryHandling_T >
class BodyMapping
{
public:

   typedef typename BoundaryHandling_T::FlagField FlagField_T;
   typedef typename FlagField_T::iterator         FlagFieldIterator;
   typedef typename BoundaryHandling_T::flag_t    flag_t;
   typedef Field< pe::BodyID, 1 >                 BodyField_T;

   inline BodyMapping( const shared_ptr<StructuredBlockStorage> & blockStorage,
                       const BlockDataID & boundaryHandlingID,
                       const BlockDataID & bodyStorageID, const BlockDataID & bodyFieldID,
                       const FlagUID & obstacle, const FlagUID & formerObstacle );

   void operator()( IBlock * const block ) const;

protected:

   shared_ptr<StructuredBlockStorage> blockStorage_;

   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyStorageID_;
   const BlockDataID bodyFieldID_;

   const FlagUID obstacle_;
   const FlagUID formerObstacle_;

}; // class BodyMapping



template< typename BoundaryHandling_T >
inline BodyMapping< BoundaryHandling_T >::BodyMapping( const shared_ptr<StructuredBlockStorage> & blockStorage,
                                                       const BlockDataID & boundaryHandlingID,
                                                       const BlockDataID & bodyStorageID, const BlockDataID & bodyFieldID,
                                                       const FlagUID & obstacle, const FlagUID & formerObstacle ) :

   blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ), bodyStorageID_(bodyStorageID), bodyFieldID_( bodyFieldID ),
   obstacle_( obstacle ), formerObstacle_( formerObstacle )
{}



template< typename BoundaryHandling_T >
void BodyMapping< BoundaryHandling_T >::operator()( IBlock * const block ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );
   FlagField_T *        flagField        = boundaryHandling->getFlagField();
   BodyField_T *        bodyField        = block->getData< BodyField_T >( bodyFieldID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyField );

   WALBERLA_ASSERT_EQUAL( flagField->xyzSize(), bodyField->xyzSize() );

   WALBERLA_ASSERT( flagField->flagExists( obstacle_ ) );

   const flag_t       obstacle = flagField->getFlag( obstacle_ );
   const flag_t formerObstacle = flagField->flagExists( formerObstacle_ ) ? flagField->getFlag( formerObstacle_ ) :
                                                                            flagField->registerFlag( formerObstacle_ );

   for( auto bodyIt = pe::BodyIterator::begin(*block, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
   {
      if( bodyIt->isFixed() /*|| !isMOBody( *bodyIt )*/ )
         continue;

      // policy: every body manages only its own flags

      CellInterval cellBB = getCellBB( *bodyIt, *block, *blockStorage_, flagField->nrOfGhostLayers() );

      Vector3<real_t> startCellCenter = blockStorage_->getBlockLocalCellCenter( *block, cellBB.min() );
      const real_t dx = blockStorage_->dx( blockStorage_->getLevel(*block) );
      const real_t dy = blockStorage_->dy( blockStorage_->getLevel(*block) );
      const real_t dz = blockStorage_->dz( blockStorage_->getLevel(*block) );

      real_t cz = startCellCenter[2];
      for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
      {
         real_t cy = startCellCenter[1];
         for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
         {
            real_t cx = startCellCenter[0];
            for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
            {

               flag_t & cellFlagPtr = flagField->get(x,y,z);

               if( bodyIt->containsPoint(cx,cy,cz) )
               {
                  if( !isFlagSet( cellFlagPtr, obstacle ) )
                  {

                     if( isFlagSet( cellFlagPtr, formerObstacle ) )
                     {
                        boundaryHandling->setBoundary( obstacle, x, y, z );
                        removeFlag( cellFlagPtr, formerObstacle );
                     }
                     else
                     {
                        boundaryHandling->forceBoundary( obstacle, x, y, z );
                     }
                  }

                  (*bodyField)(x,y,z) = *bodyIt;
               }
               else
               {
                  if( isFlagSet( cellFlagPtr, obstacle ) && ((*bodyField)(x, y, z) == *bodyIt) )
                  {
                     boundaryHandling->removeBoundary( obstacle, x, y, z );
                     addFlag( cellFlagPtr, formerObstacle );
                     // entry at (*bodyField)(x,y,z) should still point to the previous body.
                     // If during initialization the overlap between neighboring blocks
                     // was chosen correctly/large enough, the body should still be on this block.
                     // The body information is needed in the PDF restoration step.
                  }
               }

               cx += dx;
            }
            cy += dy;
         }
         cz += dz;
      }
   }
}



////////////////////////////
// FREE MAPPING FUNCTIONS //
////////////////////////////

template< typename BoundaryHandling_T >
void mapMovingBody( const pe::BodyID body, IBlock & block, StructuredBlockStorage & blockStorage,
                    const BlockDataID & boundaryHandlingID, const BlockDataID & bodyFieldID, const FlagUID & obstacle )
{
   typedef Field< pe::BodyID, 1 > BodyField_T;

   WALBERLA_ASSERT_EQUAL( &block.getBlockStorage(), &(blockStorage.getBlockStorage()) );

   if( body->isFixed() /*|| !body->isFinite()*/ )
      return;

   BoundaryHandling_T * boundaryHandling = block.getData< BoundaryHandling_T >( boundaryHandlingID );
   auto *               flagField        = boundaryHandling->getFlagField();
   BodyField_T *        bodyField        = block.getData< BodyField_T >( bodyFieldID );

   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyField );
   WALBERLA_ASSERT_EQUAL( flagField->xyzSize(), bodyField->xyzSize() );
   WALBERLA_ASSERT( flagField->flagExists( obstacle ) );

   const auto obstacleFlag = flagField->getFlag( obstacle );

   CellInterval cellBB = getCellBB( body, block, blockStorage, flagField->nrOfGhostLayers() );

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
               (*bodyField)(x,y,z) = body;
            }

            cx += dx;
         }
         cy += dy;
      }
      cz += dz;
   }
}



template< typename BoundaryHandling_T >
void mapMovingBodies( StructuredBlockStorage & blockStorage, const BlockDataID & boundaryHandlingID, const BlockDataID & bodyStorageID,
                      const BlockDataID & bodyFieldID, const FlagUID & obstacle )
{
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
       for (auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt)
           mapMovingBody< BoundaryHandling_T >( *bodyIt, *blockIt, blockStorage, boundaryHandlingID, bodyFieldID, obstacle );
   }
}

template< typename BoundaryHandling_T >
void mapMovingGlobalBodies( StructuredBlockStorage & blockStorage, const BlockDataID & boundaryHandlingID,
                            pe::BodyStorage & globalBodyStorage, const BlockDataID & bodyFieldID, const FlagUID & obstacle )
{
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
      for( auto bodyIt = globalBodyStorage.begin(); bodyIt != globalBodyStorage.end(); ++bodyIt)
      {
         mapMovingBody< BoundaryHandling_T >( *bodyIt, *blockIt, blockStorage, boundaryHandlingID, bodyFieldID, obstacle );
      }
   }
}

template< typename BoundaryHandling_T >
void mapMovingGlobalBody( const id_t globalBodySystemID,
                          StructuredBlockStorage & blockStorage, const BlockDataID & boundaryHandlingID,
                          pe::BodyStorage & globalBodyStorage, const BlockDataID & bodyFieldID, const FlagUID & obstacle )
{
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
      auto bodyIt = globalBodyStorage.find( globalBodySystemID );
      if( bodyIt != globalBodyStorage.end() )
      {
         mapMovingBody< BoundaryHandling_T >( *bodyIt, *blockIt, blockStorage, boundaryHandlingID, bodyFieldID, obstacle );
      }
   }
}


} // namespace pe_coupling
} // namespace walberla
