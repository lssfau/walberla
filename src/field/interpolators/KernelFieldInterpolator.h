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
//! \file KernelFieldInterpolator.h
//! \ingroup field
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagField.h"

#include "stencil/D3Q27.h"

#include <numeric>

namespace walberla {
namespace field {

namespace kernelweights
{

// corresponds to the smoothed dirac delta function from Roma et al - An Adaptive Version of the Immersed Boundary Method
// f(r) != 0 for abs(r) < 1.5 -> requires three neighborhood cells
inline real_t smoothedDeltaFunction( const real_t & r )
{
   real_t rAbs = std::fabs(r);
   if( rAbs <= real_t(0.5) )
   {
      return (real_t(1) + std::sqrt( - real_t(3) * r * r + real_t(1) ) ) * (real_t(1) / real_t(3));
   } else if ( rAbs < real_t(1.5) )
   {
      return (real_t(5) - real_t(3) * rAbs - std::sqrt( - real_t(3) * ( real_t(1) - rAbs ) * ( real_t(1) - rAbs ) + real_t(1) ) ) * ( real_t(1) / real_t(6) );
   } else
   {
      return real_t(0);
   }

}

// X: Lagrangian position, x: Eulerian position (usually cell center), global coordinates
// dx, dy, dz: mesh spacing
inline real_t kernelWeightFunction( const real_t & X, const real_t & Y, const real_t & Z,
                                    const real_t & x, const real_t & y, const real_t & z,
                                    const real_t & dx = real_t(1), const real_t & dy = real_t(1), const real_t & dz = real_t(1) )
{
   return smoothedDeltaFunction( ( X - x ) / dx ) * smoothedDeltaFunction( ( Y - y ) / dy ) * smoothedDeltaFunction( ( Z - z ) / dz );
}

// X: Lagrangian position, x: Eulerian position (usually cell center), global coordinates
// dx, dy, dz: mesh spacing
inline real_t kernelWeightFunction( const Vector3<real_t> & X, const Vector3<real_t> & x,
                                    const real_t & dx = real_t(1), const real_t & dy = real_t(1), const real_t & dz = real_t(1) )
{
   return smoothedDeltaFunction( ( X[0] - x[0] ) / dx ) * smoothedDeltaFunction( ( X[1] - x[1] ) / dy ) * smoothedDeltaFunction( ( X[2] - x[2] ) / dz );
}

} // namespace kernelweights


/*! Interpolator for walberla::field::GhostLayerField with kernel strategy
 *
 * \ingroup field
 *
 * This interpolator uses a smoothed dirac kernel function to interpolate values.
 * The applied weights are given in the namespace field::kernelweights.
 * Needs the full neighborhood of the containing cell, i.e. 27 cells.
 * Never construct this interpolator directly, but use the functionality from the FieldInterpolatorCreator.h instead.
 */
template< typename Field_T, typename FlagField_T >
class KernelFieldInterpolator
{
public:

   static const uint_t F_SIZE = Field_T::F_SIZE;

   typedef Field_T                                      BaseField_T;
   typedef typename FlagField_T::flag_t                 flag_t;
   typedef KernelFieldInterpolator<Field_T,FlagField_T> OwnType;

   KernelFieldInterpolator( const weak_ptr<StructuredBlockStorage> & blockStorage, const IBlock & block,
                            const BaseField_T & baseField, const FlagField_T & flagField,
                            const flag_t & evaluationMask )
   : blockStorage_( blockStorage ), block_( block ), baseField_( baseField ), flagField_( flagField ), evaluationMask_( evaluationMask )
   {
      WALBERLA_ASSERT(baseField.nrOfGhostLayers() > uint_t(0), "field for kernel interpolation needs at least one ghost layer");
   }


   inline bool operator==( const OwnType & other ){ return baseField_ == other.baseField_; }

   template< typename ForwardIterator_T >
   inline void get( const Vector3<real_t> & position, ForwardIterator_T interpolationResultBegin )
   {
      get( position[0], position[1], position[2], interpolationResultBegin );
   }

   template< typename ForwardIterator_T >
   inline void get( const real_t & x, const real_t & y, const real_t & z, ForwardIterator_T interpolationResultBegin )
   {

      WALBERLA_ASSERT(block_.getAABB().contains(x,y,z),
                      "Interpolation position <" << x << ", " << y << ", " << z << "> is not contained inside the block of this interpolator with AABB " << block_.getAABB() << " !");

      WALBERLA_CHECK( !blockStorage_.expired() );
      auto blockStorage = blockStorage_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockStorage);

      Cell centerCell = blockStorage->getBlockLocalCell( block_, x, y, z );

      const real_t dx = blockStorage->dx( blockStorage->getLevel( block_ ) );
      const real_t dy = blockStorage->dy( blockStorage->getLevel( block_ ) );
      const real_t dz = blockStorage->dz( blockStorage->getLevel( block_ ) );

      const uint_t neighborhoodSize = uint_t(1);

      CellInterval cellNeighborhood( centerCell[0] - cell_idx_c(neighborhoodSize), centerCell[1] - cell_idx_c(neighborhoodSize), centerCell[2] - cell_idx_c(neighborhoodSize),
                                     centerCell[0] + cell_idx_c(neighborhoodSize), centerCell[1] + cell_idx_c(neighborhoodSize), centerCell[2] + cell_idx_c(neighborhoodSize) );

      const uint_t kernelSizeOneDirection = uint_t(2) * neighborhoodSize + uint_t(1);
      // store the calculated weighting factors of the available cells, i.e. cells flagged by the evaluation mask
      std::vector<real_t> weights( kernelSizeOneDirection * kernelSizeOneDirection * kernelSizeOneDirection, real_t(0) );

      // store the calculated weighting factors of the UNavailable cells, i.e. cells not included in the evaluation mask
      std::vector<real_t> weightsUnavailable( kernelSizeOneDirection * kernelSizeOneDirection * kernelSizeOneDirection, real_t(0) );

      cell_idx_t cellIdx0x = cellNeighborhood.xMin();
      cell_idx_t cellIdx0y = cellNeighborhood.yMin();
      cell_idx_t cellIdx0z = cellNeighborhood.zMin();
      uint_t nx = kernelSizeOneDirection;
      uint_t nxy = kernelSizeOneDirection * kernelSizeOneDirection;
      Vector3<real_t> cellCenter0 = blockStorage->getBlockLocalCellCenter( block_, Cell(cellIdx0x, cellIdx0y, cellIdx0z) ); // = cell in neighborhood with smallest x-, y-, and z-indices

      // calculate kernel weights of all cells in neighborhood
      for( uint_t k = uint_t(0); k < kernelSizeOneDirection; ++k)
      {
         for( uint_t j = uint_t(0); j < kernelSizeOneDirection; ++j)
         {
            for( uint_t i = uint_t(0); i < kernelSizeOneDirection; ++i)
            {
               Vector3<real_t> curCellCenter( cellCenter0[0] + real_c(i) * dx, cellCenter0[1] + real_c(j) * dy, cellCenter0[2] + real_c(k) * dz);
               if( flagField_.isPartOfMaskSet( cellIdx0x + cell_idx_c(i), cellIdx0y + cell_idx_c(j), cellIdx0z + cell_idx_c(k), evaluationMask_ ) )
               {
                  weights[k*nxy+j*nx+i] = kernelweights::kernelWeightFunction( x, y, z, curCellCenter[0], curCellCenter[1], curCellCenter[2], dx, dy, dz );
               } 
               else
               {
                  weightsUnavailable[k*nxy+j*nx+i] = kernelweights::kernelWeightFunction( x, y, z, curCellCenter[0], curCellCenter[1], curCellCenter[2], dx, dy, dz );
               }
            }
         }
      }

/*
      // zero entries in weights array (and non-zero entries in weightsUnavailable) indicate non-fluid nodes
      // weight correction to incorporate extrapolation for kernel interpolation, center element can not be redistributed
      // without this part, the kernel interpolator can not extrapolate values
      // this is based on experimental findings by the author and general correctness is not ensured, thus it is commented out
      for( auto stenIt = stencil::D3Q27::beginNoCenter(); stenIt != stencil::D3Q27::end(); ++stenIt )
      {
         int cx = stenIt.cx();
         int cy = stenIt.cy();
         int cz = stenIt.cz();
         int i = cx + 1;
         int j = cy + 1;
         int k = cz + 1;

         if( weightsUnavailable[k*nxy+j*nx+i] > real_t(0) )
         {
            // check if we can distribute the non-fluid weight to nearby fluid weights that lie in one line inside the neighborhood
            // it should generally not matter in which direction it is distributed
            // check x direction
            if( weights[k*nxy+j*nx+i-cx] > real_t(0) && weights[k*nxy+j*nx+i-2*cx] > real_t(0) )
            {
               weights[k*nxy+j*nx+i-cx] += real_t(2)*weightsUnavailable[k*nxy+j*nx+i];
               weights[k*nxy+j*nx+i-2*cx] -= weightsUnavailable[k*nxy+j*nx+i];
               weightsUnavailable[k*nxy+j*nx+i] = real_t(0);
               continue;
            }
            // check y direction
            if( weights[k*nxy+(j-cy)*nx+i] > real_t(0) && weights[k*nxy+(j-2*cy)*nx+i] > real_t(0) )
            {
               weights[k*nxy+(j-cy)*nx+i] += real_t(2)*weightsUnavailable[k*nxy+j*nx+i];
               weights[k*nxy+(j-2*cy)*nx+i] -= weightsUnavailable[k*nxy+j*nx+i];
               weightsUnavailable[k*nxy+j*nx+i] = real_t(0);
               continue;
            }
            // check z direction
            if( weights[(k-cz)*nxy+j*nx+i] > real_t(0) && weights[(k-2*cz)*nxy+j*nx+i] > real_t(0) )
            {
               weights[(k-cz)*nxy+j*nx+i] += real_t(2)*weightsUnavailable[k*nxy+j*nx+i];
               weights[(k-2*cz)*nxy+j*nx+i] -= weightsUnavailable[k*nxy+j*nx+i];
               weightsUnavailable[k*nxy+j*nx+i] = real_t(0);
               continue;
            }
         }
      }
*/
      // scale available weights by the total amount of unavailable weights such that afterwards sum over all weights is 1
      const real_t sumOfWeightsUnavailable = std::accumulate(weightsUnavailable.begin(), weightsUnavailable.end(), real_t(0) );
      const real_t sumOfWeights = std::accumulate(weights.begin(), weights.end(), real_t(0) );

      // check if at least one cell was available, to prevent division by 0
      if( sumOfWeights <= real_t(0) )
         return;

      const real_t scalingFactor = real_t(1) + sumOfWeightsUnavailable / sumOfWeights;

      for( auto weightIt = weights.begin(); weightIt != weights.end(); ++weightIt )
      {
         *weightIt *= scalingFactor;
      }

      // interpolate value to interpolation position using the previously calculated weights
      // values to not be included have a weight of 0
      for( uint_t k = uint_t(0); k < kernelSizeOneDirection; ++k)
      {
         for( uint_t j = uint_t(0); j < kernelSizeOneDirection; ++j)
         {
            for( uint_t i = uint_t(0); i < kernelSizeOneDirection; ++i)
            {
               if ( weights[k*nxy+j*nx+i] > real_t(0) )
               {
                  Cell curCell( cellIdx0x + cell_idx_c(i), cellIdx0y + cell_idx_c(j), cellIdx0z + cell_idx_c(k) );
                  addWeightedCellValue( interpolationResultBegin, curCell, weights[k*nxy+j*nx+i] );
               }
            }
         }
      }
   }

private:

   template< typename ForwardIterator_T >
   void addWeightedCellValue( ForwardIterator_T interpolationResultBegin, const Cell & curCell, const real_t & weighting )
   {
      for( uint_t f = uint_t(0); f < F_SIZE; ++f )
      {
         WALBERLA_ASSERT( !math::isnan(baseField_( curCell, f)), "NaN found in component " << f << " when interpolating from cell " << curCell );

         *interpolationResultBegin += weighting * baseField_( curCell, f);
         ++interpolationResultBegin;
      }
   }

   weak_ptr<StructuredBlockStorage> blockStorage_;
   const IBlock & block_;
   const BaseField_T & baseField_;
   const FlagField_T & flagField_;
   flag_t evaluationMask_;

};


} // namespace field
} // namespace walberla
