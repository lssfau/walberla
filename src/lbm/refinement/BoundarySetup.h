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
//! \file LinearExplosion.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/Block.h"
#include "blockforest/BlockNeighborhoodSection.h"

#include "boundary/Boundary.h"

#include "core/cell/CellInterval.h"
#include "core/debug/Debug.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagUID.h"

#include "stencil/D3Q27.h"



namespace walberla {
namespace lbm {
namespace refinement {


namespace internal {
inline shared_ptr< BoundaryConfiguration > defaultBoundaryConfiguration( const Cell &, const Vector3< real_t > & )
{
   return shared_ptr< BoundaryConfiguration >( new BoundaryConfiguration );
}
}


//**********************************************************************************************************************
/*!
*   \brief Function for consistently setting boundaries across block borders (even if there is a level transition)
*
*   Boundaries are set not just in the inside of the block that is passed to this function, but also in the block's
*   ghost layers.
*   The user must provide a callback function "isBoundary": For any point in the 3D simulation space, this function must
*   either return true or false. Returning true means at this point in space, a boundary of type "flag" must be set.
*   If the boundary requires a corresponding boundary configuration parameter, a second callback function "parameter"
*   must be provided: If a boundary is set at a certain cell, this cell is passed to the callback function together
*   with the point in the 3D simulation space that, prior to that, was passed to "isBoundary". The function then must
*   return a shared pointer to a BoundaryConfiguration object.
*
*   \param blocks              block storage data structure
*   \param block               block for which boundaries of type "flag" will be set
*   \param boundaryHandlingId  block data ID for the boundary handling object
*   \param flag                flag UID that corresponds to the boundary which will be set
*   \param isBoundary          callback function for determining if for a certain point in the 3D simulation space a
*                              boundary corresponding to "flag" must be set
*   \param parameter           callback function for retrieving boundary configuration parameters
*/
//**********************************************************************************************************************
template< typename BoundaryHandling_T, bool ForceBoundary >
void consistentlySetBoundary( StructuredBlockStorage & blocks, Block & block, const BlockDataID & boundaryHandlingId, const FlagUID & flag,
                              const std::function< bool ( const Vector3< real_t > & ) > & isBoundary,
                              const std::function< shared_ptr< BoundaryConfiguration > ( const Cell &, const Vector3< real_t > & ) > & parameter = internal::defaultBoundaryConfiguration );



/// For documentation see function "template< typename BoundaryHandling_T, bool ForceBoundary > void consistentlySetBoundary(...)"
template< typename BoundaryHandling_T >
void consistentlySetBoundary( StructuredBlockStorage & blocks, Block & block, const BlockDataID & boundaryHandlingId, const FlagUID & flag,
                              const std::function< bool ( const Vector3< real_t > & ) > & isBoundary,
                              const std::function< shared_ptr< BoundaryConfiguration > ( const Cell &, const Vector3< real_t > & ) > & parameter = internal::defaultBoundaryConfiguration )
{
   consistentlySetBoundary< BoundaryHandling_T, false >( blocks, block, boundaryHandlingId, flag, isBoundary, parameter );
}



/// For documentation see function "template< typename BoundaryHandling_T, bool ForceBoundary > void consistentlySetBoundary(...)"
template< typename BoundaryHandling_T >
void consistentlyForceBoundary( StructuredBlockStorage & blocks, Block & block, const BlockDataID & boundaryHandlingId, const FlagUID & flag,
                                const std::function< bool ( const Vector3< real_t > & ) > & isBoundary,
                                const std::function< shared_ptr< BoundaryConfiguration > ( const Cell &, const Vector3< real_t > & ) > & parameter = internal::defaultBoundaryConfiguration )
{
   consistentlySetBoundary< BoundaryHandling_T, true >( blocks, block, boundaryHandlingId, flag, isBoundary, parameter );
}



template< typename BoundaryHandling_T, bool ForceBoundary >
void consistentlySetBoundary( StructuredBlockStorage & blocks, Block & block, const BlockDataID & boundaryHandlingId, const FlagUID & flag,
                              const std::function< bool ( const Vector3< real_t > & ) > & isBoundary,
                              const std::function< shared_ptr< BoundaryConfiguration > ( const Cell &, const Vector3< real_t > & ) > & parameter )
{
   const uint_t level = blocks.getLevel(block);

   BoundaryHandling_T * boundaryHandling = block.template getData< BoundaryHandling_T >( boundaryHandlingId );

   const auto * const flagField = boundaryHandling->getFlagField();
   int ghostLayers = int_c( flagField->nrOfGhostLayers() );

   CellInterval cells = flagField->xyzSize();
   cells.expand( cell_idx_c(ghostLayers) );

   std::vector< CellInterval > coarseRegions;
   for( auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir )
   {
      const auto index = blockforest::getBlockNeighborhoodSectionIndex( dir.cx(), dir.cy(), dir.cz() );
      if( block.neighborhoodSectionHasLargerBlock( index ) )
      {
         CellInterval coarseRegion( cells );
         for( uint_t i = 0; i != 3; ++i )
         {
            const auto c = stencil::c[i][*dir];

            if( c == -1 ) coarseRegion.max()[i] = coarseRegion.min()[i] + cell_idx_c( 2 * ghostLayers - 1 );
            else if( c == 1 ) coarseRegion.min()[i] = coarseRegion.max()[i] - cell_idx_c( 2 * ghostLayers - 1 );
         }
         coarseRegions.push_back( coarseRegion );
      }
   }

   for( auto cell = cells.begin(); cell != cells.end(); ++cell )
   {
      bool inCoarseRegion( false );
      for( auto region = coarseRegions.begin(); region != coarseRegions.end() && !inCoarseRegion; ++region )
         inCoarseRegion = region->contains( *cell );

      if( !inCoarseRegion )
      {
         Vector3< real_t > center;
         blocks.getBlockLocalCellCenter( block, *cell, center );
         blocks.mapToPeriodicDomain( center );

         if( isBoundary(center) )
         {
            if( ForceBoundary )
            {
               auto p = parameter( *cell, center );
               boundaryHandling->forceBoundary( flag, cell->x(), cell->y(), cell->z(), *p );
            }
            else
            {
               auto p = parameter( *cell, center );
               boundaryHandling->setBoundary( flag, cell->x(), cell->y(), cell->z(), *p );
            }
         }
      }
      else
      {
         Cell globalCell( *cell );
         blocks.transformBlockLocalToGlobalCell( globalCell, block );

         Cell coarseCell( globalCell );
         for( uint_t i = 0; i < 3; ++i )
         {
            if( coarseCell[i] < cell_idx_t(0) )
            {
               coarseCell[i] = -( ( cell_idx_t(1) - coarseCell[i] ) >> 1 );
            }
            else
            {
               coarseCell[i] >>= 1;
            }
         }

         Vector3< real_t > coarseCenter;
         blocks.getCellCenter( coarseCenter, coarseCell, level - uint_t(1) );
         blocks.mapToPeriodicDomain( coarseCenter );

         if( isBoundary(coarseCenter) )
         {
            if( ForceBoundary )
            {
               auto p = parameter( *cell, coarseCenter );            
               boundaryHandling->forceBoundary( flag, cell->x(), cell->y(), cell->z(), *p );
            }
            else
            {
               auto p = parameter( *cell, coarseCenter );
               boundaryHandling->setBoundary( flag, cell->x(), cell->y(), cell->z(), *p );
            }
         }
      }
   }
}



} // namespace refinement
} // namespace lbm
} // namespace walberla
