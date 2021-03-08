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
//! \file NearestNeighborFieldInterpolator.h
//! \ingroup field
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagField.h"

namespace walberla {
namespace field {

/*! Interpolator for walberla::field::Field with nearest neighbor strategy
 *
 * \ingroup field
 *
 * This interpolator obtains the value from the nearest cell, flagged in the evaluation mask, of the interpolation position.
 * This is usually the cell that contains the interpolation position.
 * If this cell is a cell that is not within the evaluation mask, the direct neighborhood (8 cells) will be searched for an available cell.
 * Never construct this interpolator directly, but use the functionality from FieldInterpolatorCreators.h instead.
 */
template< typename Field_T, typename FlagField_T >
class NearestNeighborFieldInterpolator
{
public:

   static const uint_t F_SIZE = Field_T::F_SIZE;

   using BaseField_T = Field_T;
   using flag_t = typename FlagField_T::flag_t;
   using OwnType = NearestNeighborFieldInterpolator<Field_T, FlagField_T>;

   NearestNeighborFieldInterpolator( const weak_ptr<StructuredBlockStorage> & blockStorage, const IBlock & block,
                                     const BaseField_T & baseField, const FlagField_T & flagField,
                                     const flag_t & evaluationMask )
   : blockStorage_( blockStorage ), block_( block ), baseField_( baseField ), flagField_( flagField ), evaluationMask_( evaluationMask )
   {}


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

      Cell nearestCell = blockStorage->getBlockLocalCell( block_, x, y, z );

      if( flagField_.isPartOfMaskSet( nearestCell, evaluationMask_ ) )
      {
         for( uint_t f = uint_t(0); f < F_SIZE; ++f )
         {
            WALBERLA_ASSERT( !math::isnan( baseField_( nearestCell, f) ), "NaN found in component " << f << " when interpolating from cell " << nearestCell );
            *interpolationResultBegin = baseField_( nearestCell, f);
            ++interpolationResultBegin;
         }
      }
      else
      {
         // look for available cell in direct neighborhood

         CellInterval fieldXYZSize = baseField_.xyzSize();

         Vector3<real_t> nearestCellCenter = blockStorage->getBlockLocalCellCenter( block_, nearestCell );
         const cell_idx_t xNeighbor = cell_idx_c( floor( x - nearestCellCenter[0] ) );
         const cell_idx_t yNeighbor = cell_idx_c( floor( y - nearestCellCenter[1] ) );
         const cell_idx_t zNeighbor = cell_idx_c( floor( z - nearestCellCenter[2] ) );

         const cell_idx_t xMin = nearestCell.x() + xNeighbor;
         const cell_idx_t yMin = nearestCell.y() + yNeighbor;
         const cell_idx_t zMin = nearestCell.z() + zNeighbor;

         for( cell_idx_t zC = zMin; zC <= zMin + cell_idx_t(1); ++zC)
         {
            for( cell_idx_t yC = yMin; yC <= yMin + cell_idx_t(1); ++yC)
            {
               for( cell_idx_t xC = xMin; xC <= xMin + cell_idx_t(1); ++xC)
               {
                  Cell curCell(xC,yC,zC);
                  if( fieldXYZSize.contains(curCell) )
                  {
                     if (flagField_.isPartOfMaskSet(curCell, evaluationMask_))
                     {
                        for (uint_t f = uint_t(0); f < F_SIZE; ++f)
                        {
                           WALBERLA_ASSERT(!math::isnan(baseField_(curCell, f)),
                                           "NaN found in component " << f << " when interpolating from cell "
                                                                     << curCell);
                           *interpolationResultBegin = baseField_(curCell, f);
                           ++interpolationResultBegin;
                        }
                        return;
                     }
                  }
               }
            }
         }
      }
   }

private:

   weak_ptr<StructuredBlockStorage> blockStorage_;
   const IBlock & block_;
   const BaseField_T & baseField_;
   const FlagField_T & flagField_;
   flag_t evaluationMask_;
};

} // namespace field
} // namespace walberla
