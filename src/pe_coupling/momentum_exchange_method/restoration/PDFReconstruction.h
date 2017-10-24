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
//! \file PDFRestoration.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================


#pragma once

#include "pe/rigidbody/BodyIterators.h"
#include "pe_coupling/mapping/BodyBBMapping.h"

#include "core/debug/Debug.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/Field.h"
#include "field/FlagField.h"
#include "field/iterators/IteratorMacros.h"
#include "Reconstructor.h"


namespace walberla {
namespace pe_coupling {


//**************************************************************************************************************************************
/*!
*   \brief Class to manage the reconstruction of PDFs that is needed when cells are becoming uncovered by moving obstacles.
*
*   Due to the explicit mapping of bodies into the fluid domain via flags, the PDFs of cells that turned from obstacle to fluid
*   are missing and must be reconstructed in order to continue with the simulation.
*   This class is to be used as a sweep in a LBM simulation with moving obstacles and calls for each cell that is tagged as
*   'formerObstacle' the specified reconstructor (see Reconstructor.h for the available variants).
*   After the successful reconstruction of all PDFs, the flags are updated to 'fluid'.
*   For small obstacle fractions, an optimized variant is available that only looks for 'formerObstacle' cells in the vicinity
*   of available bodies. It is activated via the 'optimizeForSmallObstacleFraction' argument in the constructor.
*
*/
//**************************************************************************************************************************************

template< typename LatticeModel_T, typename BoundaryHandling_T, typename Reconstructer_T >
class PDFReconstruction
{
public:

   typedef typename BoundaryHandling_T::FlagField FlagField_T;
   typedef typename FlagField_T::iterator         FlagFieldIterator;
   typedef typename BoundaryHandling_T::flag_t    flag_t;
   typedef Field< pe::BodyID, 1 >                 BodyField_T;

   inline PDFReconstruction( const shared_ptr<StructuredBlockStorage> & blockStorage,
                             const BlockDataID & boundaryHandlingID,
                             const BlockDataID & bodyStorageID, const BlockDataID & bodyFieldID,
                             const Reconstructer_T & reconstructor,
                             const FlagUID & formerObstacle, const FlagUID & fluid,
                             const bool optimizeForSmallObstacleFraction = false ) :
      blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ), bodyStorageID_(bodyStorageID), bodyFieldID_( bodyFieldID ),
      reconstructor_ ( reconstructor ), formerObstacle_( formerObstacle ), fluid_( fluid ),
      optimizeForSmallObstacleFraction_( optimizeForSmallObstacleFraction )
   {}

   void operator()( IBlock * const block );

protected:

   shared_ptr<StructuredBlockStorage> blockStorage_;

   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyStorageID_;
   const BlockDataID bodyFieldID_;

   Reconstructer_T reconstructor_;

   const FlagUID formerObstacle_;
   const FlagUID fluid_;

   const bool optimizeForSmallObstacleFraction_;

};


template< typename LatticeModel_T, typename BoundaryHandling_T, typename Reconstructer_T >
void PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructer_T >
::operator()( IBlock * const block )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );
   FlagField_T *        flagField        = boundaryHandling->getFlagField();
   BodyField_T *        bodyField        = block->getData< BodyField_T >( bodyFieldID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyField );

   WALBERLA_ASSERT_EQUAL( bodyField->xyzSize(), flagField->xyzSize() );

   WALBERLA_ASSERT( flagField->flagExists( formerObstacle_ ) );
   WALBERLA_ASSERT( flagField->flagExists( fluid_ ) );

   const flag_t formerObstacle = flagField->getFlag( formerObstacle_ );
   const flag_t fluid          = flagField->getFlag( fluid_ );

   // reconstruct all missing PDFs (only inside the domain)
   if( optimizeForSmallObstacleFraction_ )
   {
      const uint_t numberOfGhostLayersToInclude = uint_t(0);

      for( auto bodyIt = pe::BodyIterator::begin(*block, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( bodyIt->isFixed() )
            continue;

         CellInterval cellBB = getCellBB( *bodyIt, *block, blockStorage_, numberOfGhostLayersToInclude );

         for( auto cellIt = cellBB.begin(); cellIt != cellBB.end(); ++cellIt )
         {
            Cell cell( *cellIt );
            if( isFlagSet( flagField->get(cell), formerObstacle ) )
            {
               reconstructor_( cell.x(), cell.y(), cell.z(), block );
            }
         }
      }
   }
   else
   {
      auto xyzFieldSize = flagField->xyzSize();
      for( auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt )
      {
         Cell cell( *fieldIt );
         if( isFlagSet( flagField->get(cell), formerObstacle ) )
         {
            reconstructor_( cell.x(), cell.y(), cell.z(), block );
         }
      }
   }

   // update the flags from formerObstacle to fluid (inside domain & in ghost layers)
   if( optimizeForSmallObstacleFraction_ )
   {
      const uint_t numberOfGhostLayersToInclude = flagField->nrOfGhostLayers();

      for( auto bodyIt = pe::BodyIterator::begin(*block, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( bodyIt->isFixed() )
            continue;

         CellInterval cellBB = getCellBB( *bodyIt, *block, blockStorage_, numberOfGhostLayersToInclude );

         for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
         {
            for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
            {
               for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
               {
                  if( isFlagSet( flagField->get(x,y,z), formerObstacle ) )
                  {
                     boundaryHandling->setDomain( fluid, x, y, z );
                     removeFlag( flagField->get(x,y,z), formerObstacle );
                     (*bodyField)(x,y,z) = NULL;
                  }
               }
            }
         }
      }
   }
   else
   {
      CellInterval cells = flagField->xyzSizeWithGhostLayer();
      for( cell_idx_t z = cells.zMin(); z <= cells.zMax(); ++z )
      {
         for( cell_idx_t y = cells.yMin(); y <= cells.yMax(); ++y )
         {
            for( cell_idx_t x = cells.xMin(); x <= cells.xMax(); ++x )
            {
               if( isFlagSet( flagField->get(x,y,z), formerObstacle ) )
               {
                  boundaryHandling->setDomain( fluid, x, y, z );
                  removeFlag( flagField->get(x,y,z), formerObstacle );
                  (*bodyField)(x,y,z) = NULL;
               }
            }
         }
      }
   }
}


} // namespace pe_coupling
} // namespace walberla
