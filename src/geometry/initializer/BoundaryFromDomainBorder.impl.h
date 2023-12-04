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
//! \file BoundaryFromDomainBorder.impl.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/StringUtility.h"

#include <algorithm>

namespace walberla {
namespace geometry {
namespace initializer {


//*******************************************************************************************************************
/*! Constructor for BoundaryFromDomainBorder
*
* \param blocks                 the structured block storage
* \param handlerBlockDataID     the id of the BoundaryHandler, which type is the template parameter Handling
* \param globalDomain           This parameter can be used if the simulation domain should be smaller than the
*                               number of total cells the StructuredBlockStorage has. Instead of using
*                               blockStorage.getDomainCellBB() here a custom cell interval can be given.
*                               If the empty CellInterval is given (default), the complete domain
*                               is used.
*/
//*******************************************************************************************************************
template<typename Handling>
BoundaryFromDomainBorder<Handling>::BoundaryFromDomainBorder( StructuredBlockStorage & blocks,
                                                              BlockDataID handlerBlockDataID,
                                                              CellInterval globalDomain )
   : blocks_ ( blocks), handlerBlockDataID_( handlerBlockDataID )
{
   if ( globalDomain.empty() )
      globalDomain_ = blocks.getDomainCellBB();
   else
      globalDomain_ = globalDomain;
}


template<typename Handling>
void BoundaryFromDomainBorder<Handling>::init( BlockStorage & /*blockstorage*/, const Config::BlockHandle & blockHandle )
{
   init( blockHandle );
}

template<typename Handling>
void BoundaryFromDomainBorder<Handling>::init( const Config::BlockHandle & blockHandle )
{
   BoundarySetter<Handling> boundarySetter;
   boundarySetter.setConfigBlock( blockHandle );

   const std::string directionStr = blockHandle.getParameter<std::string>( "direction" );
   cell_idx_t  wallDistance            = blockHandle.getParameter<cell_idx_t>( "walldistance", 0 );
   cell_idx_t  ghostLayersToInitialize = blockHandle.getParameter<cell_idx_t>( "ghostLayersToInitialize", std::numeric_limits<cell_idx_t>::max() );

   std::vector<std::string> directionStrings = string_split( directionStr, "," );

   bool atLeastOneBoundarySet = false;
   using stencil::D3Q7;
   for( auto dirIt = D3Q7::beginNoCenter(); dirIt != D3Q7::end(); ++dirIt )
   {
      const bool isAll = string_icompare( directionStr, "all" ) == 0;
      const bool isInDirectionStrings = std::find( directionStrings.begin(), directionStrings.end(),
                                             stencil::dirToString[*dirIt] ) != directionStrings.end();

      if( isAll || isInDirectionStrings )
      {
         init( *dirIt, wallDistance, ghostLayersToInitialize, boundarySetter );
         atLeastOneBoundarySet = true;
      }
   }

   if ( ! atLeastOneBoundarySet )
      WALBERLA_ABORT( "Invalid Direction " << directionStr << ". Allowed values: all, N,S,W,E,T,B ")
}


template<typename Handling>
void BoundaryFromDomainBorder<Handling>::init( stencil::Direction direction,
                                               cell_idx_t wallDistance,
                                               cell_idx_t ghostLayersToInitialize,
                                               BoundarySetter<Handling> & boundarySetter )
{
   for( auto blockIt = blocks_.begin(); blockIt != blocks_.end(); ++blockIt )
   {
      boundarySetter.configure( *blockIt, handlerBlockDataID_ );

      const cell_idx_t gls = std::min( ghostLayersToInitialize, cell_idx_c( boundarySetter.getFlagField()->nrOfGhostLayers() ) );

      CellInterval ci;

      const cell_idx_t wd = wallDistance;

      CellInterval dBB = blocks_.getDomainCellBB( blocks_.getLevel(*blockIt) );

      for( uint_t dim = 0; dim< 3; ++dim )
         switch ( stencil::c[dim][direction] )
         {
            case -1: ci.min()[dim] =          wd;        ci.max()[dim] =  wd ;                   break;
            case  0: ci.min()[dim] =       - gls;        ci.max()[dim] =  dBB.max()[dim] + gls;  break;
            case  1: ci.min()[dim] = dBB.max()[dim]-wd;  ci.max()[dim] =  dBB.max()[dim] - wd;   break;
         }


      auto blockCellBB = blocks_.getBlockCellBB( *blockIt );
      blockCellBB.expand( cell_idx_c( boundarySetter.getFlagField()->nrOfGhostLayers() )  );
      ci.intersect( blockCellBB );
      blocks_.transformGlobalToBlockLocalCellInterval( ci, *blockIt );
      boundarySetter.set( ci );
   }
}

template<typename Handling>
void BoundaryFromDomainBorder<Handling>::init( FlagUID flagUID, stencil::Direction direction,
                                               cell_idx_t wallDistance, cell_idx_t ghostLayersToInitialize )
{
   BoundarySetter<Handling> boundarySetter;
   boundarySetter.setFlagUID( flagUID );
   init( direction, wallDistance, ghostLayersToInitialize, boundarySetter );
}

template<typename Handling>
void BoundaryFromDomainBorder<Handling>::init( BoundaryUID boundaryUID, stencil::Direction direction,
                                               const shared_ptr<BoundaryConfiguration> & conf,
                                               cell_idx_t wallDistance, cell_idx_t ghostLayersToInitialize )
{
   BoundarySetter<Handling> boundarySetter;
   boundarySetter.setBoundaryConfig( boundaryUID, conf );
   init( direction, wallDistance, ghostLayersToInitialize, boundarySetter );
}


template<typename Handling>
void BoundaryFromDomainBorder<Handling>::initAllBorders( FlagUID flagUID, cell_idx_t wallDistance,
                                                         cell_idx_t ghostLayersToInitialize )
{
   using stencil::D3Q7;
   for( auto dirIt = D3Q7::beginNoCenter(); dirIt != D3Q7::end(); ++dirIt )
      init( flagUID, *dirIt, wallDistance, ghostLayersToInitialize );

}

template<typename Handling>
void BoundaryFromDomainBorder<Handling>::initAllBorders( BoundaryUID boundaryUID, const shared_ptr<BoundaryConfiguration> & conf,
                                                         cell_idx_t wallDistance, cell_idx_t ghostLayersToInitialize )
{
   using stencil::D3Q7;
   for( auto dirIt = D3Q7::beginNoCenter(); dirIt != D3Q7::end(); ++dirIt )
      init( boundaryUID, *dirIt, conf, wallDistance, ghostLayersToInitialize );
}



} // namespace initializer
} // namespace geometry
} // namespace walberla

