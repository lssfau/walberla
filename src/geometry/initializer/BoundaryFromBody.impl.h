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
//! \file BoundaryFromBody.impl.h
//! \ingroup geometry
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//! \brief Implementations for BoundaryFromBody
//
//======================================================================================================================


#include "geometry/bodies/AABBBody.h"
#include "geometry/bodies/BodyFromConfig.h"
#include "geometry/bodies/Cylinder.h"
#include "geometry/bodies/Ellipsoid.h"
#include "geometry/bodies/Sphere.h"
#include "geometry/bodies/Torus.h"

#include "core/Abort.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"


namespace walberla {
namespace geometry {
namespace initializer {



template< typename BoundaryHandlerT >
BoundaryFromBody<BoundaryHandlerT>::BoundaryFromBody( StructuredBlockStorage & blocks, BlockDataID & boundaryHandlerID )
   : structuredBlockStorage_( blocks ),
     boundaryHandlerID_( boundaryHandlerID )
{}


template< typename BoundaryHandlerT >
void BoundaryFromBody<BoundaryHandlerT>::init( const Config::BlockHandle & blockHandle )
{
   if( !blockHandle )
      return;

   BoundarySetter<BoundaryHandlerT> boundarySetter;
   std::set< std::string> ignoreBlocks;
   ignoreBlocks.insert( "first"  );
   ignoreBlocks.insert( "second" );
   boundarySetter.setConfigBlock( blockHandle, ignoreBlocks );
   
   init ( *bodyFromConfig ( blockHandle ), boundarySetter );
}



template< typename BoundaryHandlerT >
template< typename Body >
void BoundaryFromBody<BoundaryHandlerT>::init( const Body & body, BoundarySetter<BoundaryHandlerT>  & boundarySetter )
{
   for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
   {
      IBlock * block = &(*blockIt);

      const uint_t level = structuredBlockStorage_.getLevel(*block);
      const real_t dx = structuredBlockStorage_.dx() / real_c(1 << level);
      const real_t dy = structuredBlockStorage_.dy() / real_c(1 << level);
      const real_t dz = structuredBlockStorage_.dz() / real_c(1 << level);

      boundarySetter.configure( *block, boundaryHandlerID_ );
      auto ff = boundarySetter.getFlagField();
      auto gl = cell_idx_c( ff->nrOfGhostLayers() );

      // If Block (extended with ghost layers) does not intersect body - skip the complete block
      AABB blockBB = block->getAABB();
      blockBB.extend( math::Vector3< real_t >( dx * real_c( gl ), dy * real_c( gl ), dz * real_c( gl ) ) );
      if( fastOverlapCheck( body, blockBB ) == geometry::COMPLETELY_OUTSIDE )
         continue;

      AABB firstCellBB;
      structuredBlockStorage_.getBlockLocalCellAABB( *block, ff->beginWithGhostLayer().cell(), firstCellBB );
      Vector3<real_t> firstCellMidpoint;
      for( uint_t i = 0; i < 3; ++i )
         firstCellMidpoint[i] = firstCellBB.min(i) + real_t(0.5) * firstCellBB.size(i);

      Vector3<real_t> currentMidpoint;
      currentMidpoint[2] = firstCellMidpoint[2];
      for ( cell_idx_t z = -gl; z < cell_idx_c( ff->zSize() ) + gl; ++z, currentMidpoint[2] += dz )
      {
         currentMidpoint[1] = firstCellMidpoint[1];
         for ( cell_idx_t y = -gl; y < cell_idx_c( ff->ySize() ) + gl; ++y, currentMidpoint[1] += dy )
         {
            currentMidpoint[0] = firstCellMidpoint[0];
            for( cell_idx_t x = -gl; x < cell_idx_c( ff->xSize() ) + gl; ++x, currentMidpoint[0] += dx )
            {
               if ( !contains(body, currentMidpoint) )
                  continue;

               boundarySetter.set( x,y,z );
            }
         }
      }
   }
}


template< typename BoundaryHandlerT >
template< typename Body >
void BoundaryFromBody<BoundaryHandlerT>::init( const Body & body, const BoundaryUID & uid, const shared_ptr<BoundaryConfiguration> & bcConfig )
{
   BoundarySetter<BoundaryHandlerT> boundarySetter;
   boundarySetter.setBoundaryConfig( uid,bcConfig );
   init( body, boundarySetter );
}

template< typename BoundaryHandlerT >
template< typename Body >
void BoundaryFromBody<BoundaryHandlerT>::init( const Body & body, const FlagUID & uid )
{
   BoundarySetter<BoundaryHandlerT> boundarySetter;
   boundarySetter.setFlagUID( uid );
   init( body, boundarySetter );
}


} // namespace initializer
} // namespace geometry
} // namespace walberla
