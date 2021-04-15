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
//! \file TrilinearFieldInterpolator.h
//! \ingroup field
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagField.h"
#include "field/interpolators/NearestNeighborFieldInterpolator.h"

namespace walberla {
namespace field {

/*! Interpolator for walberla::field::GhostLayerField with trilinear strategy
 *
 * \ingroup field
 *
 * This interpolator uses trilinear interpolation to obtain the interpolated value from the interpolation position.
 * If at least one of the cells, that are required for trilinear interpolation, i.e. the 8 closest cells,
 * is not contained in the evaluation mask, the nearest-neighbor interpolation (see NearestNeighborFieldInterpolator.h)
 * will be used instead.
 * Never construct this interpolator directly, but use the functionality from FieldInterpolatorCreators.h instead.
 */
template< typename Field_T, typename FlagField_T >
class TrilinearFieldInterpolator
{
public:

   static const uint_t F_SIZE = Field_T::F_SIZE;

   using BaseField_T = Field_T;
   using flag_t = typename FlagField_T::flag_t;
   using OwnType = TrilinearFieldInterpolator<Field_T, FlagField_T>;

   TrilinearFieldInterpolator( const weak_ptr<StructuredBlockStorage> & blockStorage, const IBlock & block,
                               const BaseField_T & baseField, const FlagField_T & flagField,
                               const flag_t & evaluationMask )
   : blockStorage_( blockStorage ), block_( block ), baseField_( baseField ), flagField_( flagField ), evaluationMask_( evaluationMask ),
     nearestNeighborInterpolator_( blockStorage, block, baseField, flagField, evaluationMask )
   {
      WALBERLA_ASSERT(baseField.nrOfGhostLayers() > uint_t(0), "field for trilinear interpolation needs at least one ghost layer");
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

      const real_t dx = blockStorage->dx( blockStorage->getLevel( block_ ) );
      const real_t dy = blockStorage->dy( blockStorage->getLevel( block_ ) );
      const real_t dz = blockStorage->dz( blockStorage->getLevel( block_ ) );

      Cell containingCell = blockStorage->getBlockLocalCell( block_, x, y, z );
      Vector3<real_t> containingCellCenter = blockStorage->getBlockLocalCellCenter( block_, containingCell );

      const cell_idx_t xNeighbor1 = cell_idx_c( floor( x - containingCellCenter[0] ) );
      const cell_idx_t xNeighbor2 = xNeighbor1 + cell_idx_t(1);
      const cell_idx_t yNeighbor1 = cell_idx_c( floor( y - containingCellCenter[1] ) );
      const cell_idx_t yNeighbor2 = yNeighbor1 + cell_idx_t(1);
      const cell_idx_t zNeighbor1 = cell_idx_c( floor( z - containingCellCenter[2] ) );
      const cell_idx_t zNeighbor2 = zNeighbor1 + cell_idx_t(1);

      // define the 8 nearest cells required for the trilinear interpolation
      // the cell 'ccc' is the one with the smallest x-, y-, and z-indices
      Cell ccc( containingCell.x() + xNeighbor1, containingCell.y() + yNeighbor1, containingCell.z() + zNeighbor1 );
      Cell hcc( containingCell.x() + xNeighbor2, containingCell.y() + yNeighbor1, containingCell.z() + zNeighbor1 );
      Cell chc( containingCell.x() + xNeighbor1, containingCell.y() + yNeighbor2, containingCell.z() + zNeighbor1 );
      Cell hhc( containingCell.x() + xNeighbor2, containingCell.y() + yNeighbor2, containingCell.z() + zNeighbor1 );
      Cell cch( containingCell.x() + xNeighbor1, containingCell.y() + yNeighbor1, containingCell.z() + zNeighbor2 );
      Cell hch( containingCell.x() + xNeighbor2, containingCell.y() + yNeighbor1, containingCell.z() + zNeighbor2 );
      Cell chh( containingCell.x() + xNeighbor1, containingCell.y() + yNeighbor2, containingCell.z() + zNeighbor2 );
      Cell hhh( containingCell.x() + xNeighbor2, containingCell.y() + yNeighbor2, containingCell.z() + zNeighbor2 );

      // check if all neighboring cells are domain cells
      if( flagField_.isPartOfMaskSet( ccc, evaluationMask_ ) &&
          flagField_.isPartOfMaskSet( hcc, evaluationMask_ ) &&
          flagField_.isPartOfMaskSet( chc, evaluationMask_ ) &&
          flagField_.isPartOfMaskSet( hhc, evaluationMask_ ) &&
          flagField_.isPartOfMaskSet( cch, evaluationMask_ ) &&
          flagField_.isPartOfMaskSet( hch, evaluationMask_ ) &&
          flagField_.isPartOfMaskSet( chh, evaluationMask_ ) &&
          flagField_.isPartOfMaskSet( hhh, evaluationMask_ ) )
      {
         // trilinear interpolation can be applied

         const real_t inv_totalVolume = real_t(1) / ( dx * dy * dz );
         Vector3<real_t> cccCellCenter = blockStorage->getBlockLocalCellCenter( block_, ccc );

         // weighting = volume of opposing volume element / total volume
         real_t weighting(0.0);
         // x: (c)---dxc--(X)--------------dxh-------------(h)
         const real_t dxc = x - cccCellCenter[0];
         const real_t dxh = cccCellCenter[0] + dx - x;
         const real_t dyc = y - cccCellCenter[1];
         const real_t dyh = cccCellCenter[1] + dy - y;
         const real_t dzc = z - cccCellCenter[2];
         const real_t dzh = cccCellCenter[2] + dz - z;


         weighting = dxh * dyh * dzh * inv_totalVolume;
         addWeightedCellValue( interpolationResultBegin, ccc, weighting );

         weighting = dxc * dyh * dzh * inv_totalVolume;
         addWeightedCellValue( interpolationResultBegin, hcc, weighting );

         weighting = dxh * dyc * dzh * inv_totalVolume;
         addWeightedCellValue( interpolationResultBegin, chc, weighting );

         weighting = dxc * dyc * dzh * inv_totalVolume;
         addWeightedCellValue( interpolationResultBegin, hhc, weighting );

         weighting = dxh * dyh * dzc * inv_totalVolume;
         addWeightedCellValue( interpolationResultBegin, cch, weighting );

         weighting = dxc * dyh * dzc * inv_totalVolume;
         addWeightedCellValue( interpolationResultBegin, hch, weighting );

         weighting = dxh * dyc * dzc * inv_totalVolume;
         addWeightedCellValue( interpolationResultBegin, chh, weighting );

         weighting = dxc * dyc * dzc * inv_totalVolume;
         addWeightedCellValue( interpolationResultBegin, hhh, weighting );

      }
      else
      {
         // revert to nearest neighbor interpolation
         nearestNeighborInterpolator_.get( x, y, z, interpolationResultBegin );
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

   NearestNeighborFieldInterpolator<Field_T,FlagField_T> nearestNeighborInterpolator_;

};


} // namespace field
} // namespace walberla
