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
//! \file BoundaryFromCellInterval.impl.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Implementations for BoundaryFromCellInterval
//
//======================================================================================================================

#include "core/logging/Logging.h"

namespace walberla {
namespace geometry {
namespace initializer {



template< typename BoundaryHandlerT >
BoundaryFromCellInterval<BoundaryHandlerT>::BoundaryFromCellInterval( StructuredBlockStorage & blocks, BlockDataID & boundaryHandlerID )
   : structuredBlockStorage_( blocks ),
     boundaryHandlerID_     ( boundaryHandlerID )
{}


template< typename BoundaryHandlerT >
void BoundaryFromCellInterval<BoundaryHandlerT>::init( const CellInterval & globalCellInterval,
                                                       BoundarySetter<BoundaryHandlerT> & boundarySetter)
{
   for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
   {
      CellInterval localCellInterval;
      structuredBlockStorage_.transformGlobalToBlockLocalCellInterval( localCellInterval, *blockIt, globalCellInterval );
      boundarySetter.configure( *blockIt, boundaryHandlerID_ );

      auto flagField = boundarySetter.getFlagField();
      localCellInterval.intersect(flagField->xyzSizeWithGhostLayer());
      boundarySetter.set( localCellInterval );
   }
}


template< typename BoundaryHandlerT >
void BoundaryFromCellInterval<BoundaryHandlerT>::init( const Config::BlockHandle & blockHandle )
{
   if( !blockHandle )
      return;

   CellInterval globalCellInterval;
   if ( blockHandle.isDefined("CellInterval") )
      globalCellInterval = blockHandle.getParameter<CellInterval>("CellInterval");
   else
   {
      const Vector3<cell_idx_t> min = blockHandle.getParameter< Vector3<cell_idx_t> >( "min");
      const Vector3<cell_idx_t> max = blockHandle.getParameter< Vector3<cell_idx_t> >( "max");

      globalCellInterval.min() = Cell( min[0], min[1], min[2] );
      globalCellInterval.max() = Cell( max[0], max[1], max[2] );
   }

   BoundarySetter<BoundaryHandlerT> boundarySetter;
   boundarySetter.setConfigBlock( blockHandle );

   init( globalCellInterval, boundarySetter );
}


template< typename BoundaryHandlerT >
void BoundaryFromCellInterval<BoundaryHandlerT>::init( const CellInterval & ci, const BoundaryUID & uid,
                                                       const shared_ptr<BoundaryConfiguration> & bcConfig )
{
   BoundarySetter<BoundaryHandlerT> boundarySetter;
   boundarySetter.setBoundaryConfig( uid,bcConfig );
   init( ci, boundarySetter );
}

template< typename BoundaryHandlerT >
void BoundaryFromCellInterval<BoundaryHandlerT>::init( const CellInterval & ci, const FlagUID & uid )
{
   BoundarySetter<BoundaryHandlerT> boundarySetter;
   boundarySetter.setFlagUID( uid );
   init( ci, boundarySetter );
}



} // namespace initializer
} // namespace geometry
} // namespace walberla
