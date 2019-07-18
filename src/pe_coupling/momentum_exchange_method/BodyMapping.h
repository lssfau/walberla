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

#include "core/ConcatIterator.h"
#include "core/debug/Debug.h"
#include "core/cell/Cell.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/FlagField.h"
#include "field/iterators/IteratorMacros.h"

#include "pe/rigidbody/BodyIterators.h"

#include "pe_coupling/mapping/BodyBBMapping.h"
#include "pe_coupling/utility/BodySelectorFunctions.h"

#include <functional>

namespace walberla {
namespace pe_coupling {

/*!\brief Maps the moving bodies into the simulation domain and updates the mapping
 *
 * Cells that are inside a body, will be marked with the 'obstacle' flag.
 * 'Inside' means that the cell center is contained in the body.
 * Thus, the body has to provide a valid containsPoint(x,y,z) function.
 *
 * Cells that in the last time step were inside the body but are now outside of it, i.e. the body has moved,
 * will be marked with the 'formerObstacle' flag.
 * The 'formerObstacle' flag is used in a second step by the PDFReconstruction class (see PDFReconstruction.h) to
 * re-initialize the missing PDFs. Afterwards, the 'formerObstacle' flag is removed and the 'fluid' flag is set.
 *
 * It is not required that the mapping has been initialized with one of the free functions from below.
 *
 * The 'mappingBodySelectorFct' can be used to decide which bodies should be mapped or not.
 */
template< typename LatticeModel_T, typename BoundaryHandling_T, typename Destroyer_T = NaNDestroyer<LatticeModel_T>, bool conserveMomentum = false >
class BodyMapping
{
public:

   typedef lbm::PdfField< LatticeModel_T >        PdfField_T;
   typedef typename BoundaryHandling_T::FlagField FlagField_T;
   typedef typename BoundaryHandling_T::flag_t    flag_t;
   typedef Field< pe::BodyID, 1 >                 BodyField_T;

   BodyMapping( const shared_ptr<StructuredBlockStorage> & blockStorage,
                const BlockDataID & pdfFieldID,
                const BlockDataID & boundaryHandlingID,
                const BlockDataID & bodyStorageID,
                const shared_ptr<pe::BodyStorage> & globalBodyStorage,
                const BlockDataID & bodyFieldID,
                const Destroyer_T & destroyer,
                const FlagUID & obstacle, const FlagUID & formerObstacle,
                const std::function<bool(pe::BodyID)> & mappingBodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), pdfFieldID_( pdfFieldID ), boundaryHandlingID_( boundaryHandlingID ),
     bodyStorageID_(bodyStorageID), globalBodyStorage_( globalBodyStorage ), bodyFieldID_( bodyFieldID ),
     destroyer_( destroyer ), obstacle_( obstacle ), formerObstacle_( formerObstacle ),
     mappingBodySelectorFct_( mappingBodySelectorFct )
   {}

   BodyMapping( const shared_ptr<StructuredBlockStorage> & blockStorage,
                const BlockDataID & pdfFieldID,
                const BlockDataID & boundaryHandlingID,
                const BlockDataID & bodyStorageID,
                const shared_ptr<pe::BodyStorage> & globalBodyStorage,
                const BlockDataID & bodyFieldID,
                const FlagUID & obstacle, const FlagUID & formerObstacle,
                const std::function<bool(pe::BodyID)> & mappingBodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), pdfFieldID_( pdfFieldID ), boundaryHandlingID_( boundaryHandlingID ),
     bodyStorageID_(bodyStorageID), globalBodyStorage_( globalBodyStorage ), bodyFieldID_( bodyFieldID ),
     destroyer_( Destroyer_T() ), obstacle_( obstacle ), formerObstacle_( formerObstacle ),
     mappingBodySelectorFct_( mappingBodySelectorFct )
   {}

   void operator()( IBlock * const block )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );

      PdfField_T *         pdfField         = block->getData< PdfField_T >( pdfFieldID_ );
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );
      FlagField_T *        flagField        = boundaryHandling->getFlagField();
      BodyField_T *        bodyField        = block->getData< BodyField_T >( bodyFieldID_ );

      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
      WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );
      WALBERLA_ASSERT_NOT_NULLPTR( bodyField );

      WALBERLA_ASSERT_EQUAL( flagField->xyzSize(), bodyField->xyzSize() );

      WALBERLA_ASSERT( flagField->flagExists( obstacle_ ) );

      const flag_t       obstacle = flagField->getFlag( obstacle_ );
      const flag_t formerObstacle = flagField->flagExists( formerObstacle_ ) ? flagField->getFlag( formerObstacle_ ) :
                                    flagField->registerFlag( formerObstacle_ );

      const real_t dx = blockStorage_->dx( blockStorage_->getLevel(*block) );
      const real_t dy = blockStorage_->dy( blockStorage_->getLevel(*block) );
      const real_t dz = blockStorage_->dz( blockStorage_->getLevel(*block) );

      for( auto bodyIt = pe::BodyIterator::begin(*block, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( mappingBodySelectorFct_(bodyIt.getBodyID()) )
            mapBodyAndUpdateMapping(bodyIt.getBodyID(), block, pdfField, boundaryHandling, flagField , bodyField, obstacle, formerObstacle, dx, dy, dz);
      }
      for( auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt)
      {
         if( mappingBodySelectorFct_(bodyIt.getBodyID()))
            mapBodyAndUpdateMapping(bodyIt.getBodyID(), block, pdfField, boundaryHandling, flagField , bodyField, obstacle, formerObstacle, dx, dy, dz);
      }
   }

private:

   void mapBodyAndUpdateMapping(pe::BodyID body, IBlock * const block,
                                PdfField_T * pdfField, BoundaryHandling_T * boundaryHandling, FlagField_T * flagField, BodyField_T * bodyField,
                                const flag_t & obstacle, const flag_t & formerObstacle,
                                real_t dx, real_t dy, real_t dz)
   {
      // policy: every body manages only its own flags

      CellInterval cellBB = getCellBB( body, *block, *blockStorage_, flagField->nrOfGhostLayers() );

      WALBERLA_ASSERT_LESS_EQUAL(body->getLinearVel().length(), real_t(1),
            "Velocity is above 1 (" << body->getLinearVel() << "), which violates the assumption made in the getCellBB() function. The coupling might thus not work properly. Body:\n" << *body);

      Vector3<real_t> startCellCenter = blockStorage_->getBlockLocalCellCenter( *block, cellBB.min() );

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

               if( body->containsPoint(cx,cy,cz) )
               {
                  // cell is inside body
                  if( !isFlagSet( cellFlagPtr, obstacle ) )
                  {
                     // cell is not yet an obstacle cell
                     if( isFlagSet( cellFlagPtr, formerObstacle ) )
                     {
                        // cell was marked as former obstacle, e.g. by another body that has moved away
                        boundaryHandling->setBoundary( obstacle, x, y, z );
                        removeFlag( cellFlagPtr, formerObstacle );
                     }
                     else
                     {
                        // set obstacle flag
                        boundaryHandling->forceBoundary( obstacle, x, y, z );

                        if( conserveMomentum && pdfField->isInInnerPart(Cell(x,y,z)) ) {
                           // this removes the fluid information (PDFs) from the simulation, together with the containing momentum
                           // to ensure momentum conservation, the momentum of this fluid cell has to be added to the particle
                           // see Aidun, C. K., Lu, Y., & Ding, E. J. (1998). Direct analysis of particulate suspensions with inertia using the discrete Boltzmann equation. Journal of Fluid Mechanics, 373, 287-311.

                           // force = momentum / dt, with dt = 1
                           Vector3<real_t> momentum = pdfField->getMomentumDensity(x, y, z);
                           body->addForceAtPos(momentum[0], momentum[1], momentum[2], cx, cy, cz);
                        }

                        // invalidate PDF values
                        destroyer_( x, y, z, block, pdfField );
                     }
                  }
                  // let pointer from body field point to this body
                  (*bodyField)(x,y,z) = body;

                  WALBERLA_ASSERT(isFlagSet( cellFlagPtr, obstacle ), "Flag mapping incorrect for body\n" << *body );
                  WALBERLA_ASSERT_EQUAL((*bodyField)(x,y,z), body, "Body field does not point to correct body\n" << *body << ".");
               }
               else
               {
                  // cell is outside body
                  if( isFlagSet( cellFlagPtr, obstacle ) && ((*bodyField)(x, y, z) == body) )
                  {
                     // cell was previously occupied by this body
                     boundaryHandling->removeBoundary( obstacle, x, y, z );
                     addFlag( cellFlagPtr, formerObstacle );
                     // entry at (*bodyField)(x,y,z) should still point to the previous body.
                     // If during initialization the overlap between neighboring blocks
                     // was chosen correctly/large enough, the body should still be on this block.
                     // The body information is needed in the PDF restoration step.
                     // There, the flag will be removed and replaced by a domain flag after restoration.
                  }
               }

               cx += dx;
            }
            cy += dy;
         }
         cz += dz;
      }
   }

   shared_ptr<StructuredBlockStorage> blockStorage_;

   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyStorageID_;
   shared_ptr<pe::BodyStorage> globalBodyStorage_;
   const BlockDataID bodyFieldID_;

   Destroyer_T destroyer_;

   const FlagUID obstacle_;
   const FlagUID formerObstacle_;

   std::function<bool(pe::BodyID)> mappingBodySelectorFct_;

}; // class BodyMapping


////////////////////////////
// FREE MAPPING FUNCTIONS //
////////////////////////////


// general mapping functions for a given single (moving) body on a given single block
template< typename BoundaryHandling_T >
void mapMovingBody( pe::BodyID body, IBlock & block, StructuredBlockStorage & blockStorage,
                    const BlockDataID & boundaryHandlingID, const BlockDataID & bodyFieldID, const FlagUID & obstacle )
{
   typedef Field< pe::BodyID, 1 > BodyField_T;

   WALBERLA_ASSERT_EQUAL( &block.getBlockStorage(), &(blockStorage.getBlockStorage()) );

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
               (*bodyField)(x,y,z) = body;
            }

            cx += dx;
         }
         cy += dy;
      }
      cz += dz;
   }
}

/*!\brief Mapping function to map all bodies - with certain properties - onto all blocks to the "moving obstacle" boundary condition, that requires the additional bodyField
 *
 * All bodies (from bodyStorage and globalBodyStorage) are iterated and mapped to all blocks.
 * Cells that are inside the bodies are set to 'obstacle' and a pointer to the body is stored in the bodyField.
 *
 * Whether or not a body is mapped depends on the return value of the 'mappingBodySelectorFct'.
 */
template< typename BoundaryHandling_T >
void mapMovingBodies( StructuredBlockStorage & blockStorage, const BlockDataID & boundaryHandlingID,
                      const BlockDataID & bodyStorageID, pe::BodyStorage & globalBodyStorage,
                      const BlockDataID & bodyFieldID, const FlagUID & obstacle,
                      const std::function<bool(pe::BodyID)> & mappingBodySelectorFct = selectAllBodies )
{
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
      for (auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt)
      {
         if( mappingBodySelectorFct(bodyIt.getBodyID()))
            mapMovingBody< BoundaryHandling_T >( bodyIt.getBodyID(), *blockIt, blockStorage, boundaryHandlingID, bodyFieldID, obstacle );
      }
      for( auto bodyIt = globalBodyStorage.begin(); bodyIt != globalBodyStorage.end(); ++bodyIt)
      {
         if( mappingBodySelectorFct(bodyIt.getBodyID()))
            mapMovingBody< BoundaryHandling_T >( bodyIt.getBodyID(), *blockIt, blockStorage, boundaryHandlingID, bodyFieldID, obstacle );
      }
   }
}


} // namespace pe_coupling
} // namespace walberla
