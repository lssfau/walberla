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
//! \author Lukas Werner <lks.werner@fau.de>
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


/*!\brief Level determination for refinement check based on local curl
 *
 * If (scaled) vorticity magnitude is below lowerLimit in all cells of a block, that block could be coarsened.
 * If the (scaled) vorticity value is above the upperLimit for at least one cell, that block gets marked for refinement.
 * Else, the block remains on the current level.
 *
 * The scaling originates from neglecting the actual mesh size on the block to obtain different vorticity values for
 * different mesh sizes.
 *
 * Parametes upperLimit corresponds to sigma_c, coarsenFactor to c, lowerLimit to c*sigma_c, lengthScaleWeight to r.
 */
template< typename Filter_T >
class CurlBasedLevelDetermination // used as a 'BlockForest::RefreshMinTargetLevelDeterminationFunction'
{

public:

   typedef GhostLayerField< Vector3<real_t>, 1 >  VectorField_T;

   CurlBasedLevelDetermination(const ConstBlockDataID & fieldID, const StructuredBlockForest & structuredBlockForest,
         const Filter_T & filter, const uint_t maxLevel,
         const real_t upperLimit, const real_t lowerLimit, const real_t lengthScaleWeight = real_t(2)) :
         fieldID_(fieldID), structuredBlockForest_(structuredBlockForest), filter_(filter), maxLevel_(maxLevel),
         upperLimitSqr_(upperLimit*upperLimit), lowerLimitSqr_(lowerLimit*lowerLimit),
         lengthScaleWeight_(lengthScaleWeight)
   {
      WALBERLA_CHECK_FLOAT_UNEQUAL(lengthScaleWeight_, real_t(0)); // else std::pow(x, y/0) is calculated further below
   }

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > & blocksAlreadyMarkedForRefinement,
                    const BlockForest & forest );

private:

   ConstBlockDataID fieldID_;
   const StructuredBlockForest & structuredBlockForest_;

   Filter_T filter_;

   uint_t maxLevel_;

   real_t upperLimitSqr_;
   real_t lowerLimitSqr_;

   real_t lengthScaleWeight_;
};

template< typename Filter_T >
void CurlBasedLevelDetermination< Filter_T >::operator()( std::vector< std::pair< const Block *,
      uint_t > > & minTargetLevels, std::vector< const Block * > &, const BlockForest & forest) {

   for(auto & minTargetLevel : minTargetLevels) {
      const Block * const block = minTargetLevel.first;
      const VectorField_T * u = block->template getData< VectorField_T >(fieldID_);

      if(u == nullptr) {
         minTargetLevel.second = uint_t(0);
         continue;
      }

      WALBERLA_ASSERT_GREATER_EQUAL(u->nrOfGhostLayers(), uint_t(1));

      CellInterval interval = u->xyzSize();
      Cell expand(cell_idx_c(-1), cell_idx_c(-1), cell_idx_c(-1));
      interval.expand(expand);

      const auto one = cell_idx_t(1);

      const real_t dx = structuredBlockForest_.dx(forest.getLevel(*block));
      const real_t dy = structuredBlockForest_.dy(forest.getLevel(*block));
      const real_t dz = structuredBlockForest_.dz(forest.getLevel(*block));

      const auto halfInvDx = real_t(0.5) * real_t(1) / dx;
      const auto halfInvDy = real_t(0.5) * real_t(1) / dy;
      const auto halfInvDz = real_t(0.5) * real_t(1) / dz;

      bool refine = false;
      bool coarsen = true;

      filter_(*block);

      const cell_idx_t xSize = cell_idx_c(interval.xSize());
      const cell_idx_t ySize = cell_idx_c(interval.ySize());
      const cell_idx_t zSize = cell_idx_c(interval.zSize());

      const real_t lengthScale = std::cbrt(dx*dy*dz);
      const real_t weightedLengthScale = std::pow(lengthScale, (lengthScaleWeight_+1)/lengthScaleWeight_);
      const real_t weightedLengthScaleSqr = weightedLengthScale*weightedLengthScale;

      for (auto z = cell_idx_t(0); z < zSize; ++z) {
         for (auto y = cell_idx_t(0); y < ySize; ++y) {
            for (auto x = cell_idx_t(0); x < xSize; ++x) {
               if (filter_(x,y,z) && filter_(x+one,y,z) && filter_(x-one,y,z) && filter_(x,y+one,z) && filter_(x,y-one,z)
                  && filter_(x,y,z+one) && filter_(x,y,z-one)) {
                  const Vector3< real_t > xa = u->get(x+one,y,z);
                  const Vector3< real_t > xb = u->get(x-one,y,z);
                  const Vector3< real_t > ya = u->get(x,y+one,z);
                  const Vector3< real_t > yb = u->get(x,y-one,z);
                  const Vector3< real_t > za = u->get(x,y,z+one);
                  const Vector3< real_t > zb = u->get(x,y,z-one);

                  const real_t duzdy = halfInvDy * (ya[2] - yb[2]);
                  const real_t duydz = halfInvDz * (za[1] - zb[1]);
                  const real_t duxdz = halfInvDz * (za[0] - zb[0]);
                  const real_t duzdx = halfInvDx * (xa[2] - xb[2]);
                  const real_t duydx = halfInvDx * (xa[1] - xb[1]);
                  const real_t duxdy = halfInvDy * (ya[0] - yb[0]);

                  const Vector3< real_t > curl( duzdy - duydz, duxdz - duzdx, duydx - duxdy );
                  const auto curlSqr = curl.sqrLength();

                  const auto curlSensorSqr = curlSqr * weightedLengthScaleSqr;

                  if (curlSensorSqr > lowerLimitSqr_) {
                     // curl is not small enough to coarsen, i.e. stay at least on the current level
                     coarsen = false;
                     if (curlSensorSqr > upperLimitSqr_) {
                        // curl is big enough for refinement
                        refine = true;
                     }
                  }
               }
            }
         }
      }

      if (refine && block->getLevel() < maxLevel_) {
         WALBERLA_ASSERT(!coarsen);
         minTargetLevel.second = block->getLevel() + uint_t(1);
      }
      if (coarsen && block->getLevel() > uint_t(0)) {
         WALBERLA_ASSERT(!refine);
         minTargetLevel.second = block->getLevel() - uint_t(1);
      }
   }
}

} // namespace refinement
} // namespace lbm
} // namespace walberla
