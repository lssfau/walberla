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
//! \file VorticityBasedLevelDetermination.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockForest.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/BlockDataID.h"
#include "field/GhostLayerField.h"

#include <vector>

namespace walberla {
namespace lbm {
namespace refinement {


/*!\brief Level determination for refinement check based on (scaled) vorticity values
 *
 * If (scaled) vorticity magnitude is below lowerLimit in all cells of a block, that block could be coarsened.
 * If the (scaled) vorticity value is above the upperLimit for at least one cell, that block gets marked for refinement.
 * Else, the block remains on the current level.
 *
 * The scaling originates from neglecting the actual mesh size on the block to obtain different vorticity values for
 * different mesh sizes.
 */
template< typename Filter_T, bool Pseudo2D = false >
class VorticityBasedLevelDetermination // used as a 'BlockForest::RefreshMinTargetLevelDeterminationFunction'
{

public:

   typedef GhostLayerField< Vector3<real_t>, 1 >  VectorField_T;

   VorticityBasedLevelDetermination( const ConstBlockDataID & fieldID, const Filter_T & filter,
                                     const real_t upperLimit, const real_t lowerLimit, const uint_t maxLevel ) :
         fieldID_( fieldID ), filter_( filter ),
         upperLimit_( upperLimit * upperLimit ), lowerLimit_( lowerLimit * lowerLimit ), maxLevel_( maxLevel )
   {}

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > & blocksAlreadyMarkedForRefinement,
                    const BlockForest & forest );

private:

   ConstBlockDataID fieldID_;

   Filter_T filter_;

   real_t upperLimit_;
   real_t lowerLimit_;

   uint_t maxLevel_;

}; // class VorticityRefinement

template< typename Filter_T, bool Pseudo2D >
void VorticityBasedLevelDetermination< Filter_T, Pseudo2D >::operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                                                                         std::vector< const Block * > &, const BlockForest & )
{
   for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
   {
      const Block * const block = it->first;
      const VectorField_T * u = block->template getData< VectorField_T >( fieldID_ );

      if( u == nullptr )
      {
         it->second = uint_t(0);
         continue;
      }

      WALBERLA_ASSERT_GREATER_EQUAL(u->nrOfGhostLayers(), uint_t(1));

      CellInterval interval = u->xyzSize();
      Cell expand( cell_idx_c(-1), cell_idx_c(-1), Pseudo2D ? cell_idx_t(0) : cell_idx_c(-1) );
      interval.expand( expand );

      const cell_idx_t one( cell_idx_t(1) );
      const real_t half( real_c(0.5) );

      bool refine( false );
      bool coarsen( true );

      filter_( *block );

      const cell_idx_t xSize = cell_idx_c( interval.xSize() );
      const cell_idx_t ySize = cell_idx_c( interval.ySize() );
      const cell_idx_t zSize = cell_idx_c( interval.zSize() );

      for( cell_idx_t z = cell_idx_t(0); z < zSize; ++z ) {
         for (cell_idx_t y = cell_idx_t(0); y < ySize; ++y) {
            for (cell_idx_t x = cell_idx_t(0); x < xSize; ++x) {
               if( filter_(x,y,z) && filter_(x+one,y,z) && filter_(x-one,y,z) && filter_(x,y+one,z) && filter_(x,y-one,z) &&
                   ( Pseudo2D || (filter_(x,y,z+one) && filter_(x,y,z-one)) ) )
               {
                  const Vector3< real_t > xa = u->get(x+one,y,z);
                  const Vector3< real_t > xb = u->get(x-one,y,z);
                  const Vector3< real_t > ya = u->get(x,y+one,z);
                  const Vector3< real_t > yb = u->get(x,y-one,z);
                  const Vector3< real_t > za = Pseudo2D ? Vector3< real_t >(0) : u->get(x,y,z+one);
                  const Vector3< real_t > zb = Pseudo2D ? Vector3< real_t >(0) : u->get(x,y,z-one);

                  const real_t duzdy = half * ( ya[2] - yb[2] );
                  const real_t duydz = half * ( za[1] - zb[1] );
                  const real_t duxdz = half * ( za[0] - zb[0] );
                  const real_t duzdx = half * ( xa[2] - xb[2] );
                  const real_t duydx = half * ( xa[1] - xb[1] );
                  const real_t duxdy = half * ( ya[0] - yb[0] );

                  // since no dx was used in the gradient calculation, this value is actually curl * dx
                  // this is needed to have changing curl values on different grid levels
                  const Vector3< real_t > curl( duzdy - duydz, duxdz - duzdx, duydx - duxdy );
                  const auto curlSqr = curl.sqrLength();

                  if( curlSqr > lowerLimit_ )
                  {
                     coarsen = false;
                     if( curlSqr > upperLimit_ )
                        refine = true;
                  }
               }
            }
         }
      }

      if( refine && block->getLevel() < maxLevel_ )
      {
         WALBERLA_ASSERT( !coarsen );
         it->second = block->getLevel() + uint_t(1);
      }
      if( coarsen && block->getLevel() > uint_t(0) )
      {
         WALBERLA_ASSERT( !refine );
         it->second = block->getLevel() - uint_t(1);
      }
   }
}

} // namespace refinement
} // namespace lbm
} // namespace walberla
