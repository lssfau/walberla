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
//! \file KernelDistributor.h
//! \ingroup field
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/interpolators/KernelFieldInterpolator.h"
#include "field/GhostLayerField.h"

#include <vector>

namespace walberla {
namespace field {

/*! Distributor for walberla::field::GhostLayerField with kernel strategy
 *
 * \ingroup field
 *
 * This distributor uses a smoothed dirac kernel function to distribute values to the field at the given position.
 * The applied weights are given in the namespace field::kernelweights.
 * Needs the full neighborhood of the containing cell, i.e. 27 cells.
 * Never construct this distributor directly, but use the functionality from DistributorCreators.h instead.
 */
template< typename Field_T, typename FlagField_T >
class KernelDistributor
{
public:

   static const uint_t F_SIZE = Field_T::F_SIZE;

   using BaseField_T = Field_T;
   using flag_t = typename FlagField_T::flag_t;
   using OwnType = KernelDistributor<Field_T, FlagField_T>;

   KernelDistributor( const weak_ptr<StructuredBlockStorage> & blockStorage, const IBlock & block,
                      BaseField_T & baseField, const FlagField_T & flagField,
                      const flag_t & evaluationMask )
   : blockStorage_( blockStorage ), block_( block ), baseField_( baseField ), flagField_( flagField ), evaluationMask_( evaluationMask )
   {
      WALBERLA_ASSERT(baseField.nrOfGhostLayers() > uint_t(0), "field for kernel distribution needs at least one ghost layer");
   }

   inline bool operator==( const OwnType & other ){ return baseField_ == other.baseField_; }

   template< typename ForwardIterator_T >
   inline void distribute( const Vector3<real_t> & position, ForwardIterator_T distributeValueBegin )
   {
      distribute( position[0], position[1], position[2], distributeValueBegin );
   }

   template< typename ForwardIterator_T >
   inline void distribute( const real_t & x, const real_t & y, const real_t & z, ForwardIterator_T distributeValueBegin )
   {
      WALBERLA_ASSERT(block_.getAABB().contains(x,y,z),
                      "Distribution position <" << x << ", " << y << ", " << z << "> is not contained inside the block of this distributor with AABB " << block_.getAABB() << " !");

      WALBERLA_CHECK( !blockStorage_.expired() );
      auto blockStorage = blockStorage_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockStorage);

      Cell centerCell = blockStorage->getBlockLocalCell( block_, x, y, z );

      const real_t dx = blockStorage->dx( blockStorage->getLevel( block_ ) );
      const real_t dy = blockStorage->dy( blockStorage->getLevel( block_ ) );
      const real_t dz = blockStorage->dz( blockStorage->getLevel( block_ ) );
      
      const uint_t neighborhoodSize = cell_idx_t(1);

      CellInterval cellNeighborhood( centerCell[0] - cell_idx_c(neighborhoodSize), centerCell[1] - cell_idx_c(neighborhoodSize), centerCell[2] - cell_idx_c(neighborhoodSize),
                                     centerCell[0] + cell_idx_c(neighborhoodSize), centerCell[1] + cell_idx_c(neighborhoodSize), centerCell[2] + cell_idx_c(neighborhoodSize) );

      const uint_t kernelSizeOneDirection = uint_t(2) * neighborhoodSize + uint_t(1);
      std::vector<real_t> weights( kernelSizeOneDirection * kernelSizeOneDirection * kernelSizeOneDirection, real_t(0) );

      uint_t counter = uint_t(0);
      real_t sumOfWeights = real_t(0);
      real_t sumOfWeightsUnavailable = real_t(0);

      // get distribution weights and count available cells in surrounding cells
      for( auto cellIt = cellNeighborhood.begin(); cellIt != cellNeighborhood.end(); ++cellIt )
      {
         Vector3<real_t> curCellCenter = blockStorage->getBlockLocalCellCenter( block_, *cellIt );
         if( flagField_.isPartOfMaskSet( *cellIt, evaluationMask_ ) )
         {
            weights[counter] = kernelweights::kernelWeightFunction( x, y, z, curCellCenter[0], curCellCenter[1], curCellCenter[2], dx, dy, dz );
            sumOfWeights += weights[counter];
         }
         else
         {
            weights[counter] = real_t(0);
            sumOfWeightsUnavailable += kernelweights::kernelWeightFunction( x, y, z, curCellCenter[0], curCellCenter[1], curCellCenter[2], dx, dy, dz );
         }
         ++counter;
      }

      // check if at least one cell was available, to prevent division by 0
      if( sumOfWeights <= real_t(0) )
         return;

      // scale domain weights if some non-domain cells are in neighborhood
      const real_t scalingFactor = real_t(1) + sumOfWeightsUnavailable / sumOfWeights;

      // distribute the values to the neighboring domain cells with the corresponding (scaled) weighting
      counter = uint_t(0);
      for( auto cellIt = cellNeighborhood.begin(); cellIt != cellNeighborhood.end(); ++cellIt )
      {
         if ( weights[counter] > real_t(0) )
         {
            addWeightedCellValue( distributeValueBegin, *cellIt, scalingFactor * weights[counter] );
         }
         ++counter;
      }
   }

private:

   template< typename ForwardIterator_T >
   void addWeightedCellValue( ForwardIterator_T distributeValueBegin, const Cell & curCell, const real_t & weighting )
   {
      for( uint_t f = uint_t(0); f < F_SIZE; ++f )
      {
         baseField_( curCell, f) += weighting * (*distributeValueBegin);
         ++distributeValueBegin;
      }
   }

   weak_ptr<StructuredBlockStorage> blockStorage_;
   const IBlock & block_;
   BaseField_T & baseField_;
   const FlagField_T & flagField_;
   flag_t evaluationMask_;

};


} // namespace field
} // namespace walberla
