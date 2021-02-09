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
//! \file FieldInterpolatorCreators.h
//! \ingroup field
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/Set.h"

#include "blockforest/BlockDataHandling.h"

#include "domain_decomposition/StructuredBlockStorage.h"

namespace walberla {
namespace field {

//**********************************************************************************************************************
/*! FieldInterpolatorCreators
 *
 * \ingroup field
 *
 * Interpolator_T: A field interpolator that has a constructor
 *  ( const weak_ptr<StructuredBlockStorage> & blockStorage, const IBlock & block, const BaseField_T & baseField,
 *    const FlagField_T & flagField, const flag_t & evaluationMask )
 * and getter functions:
 *  template< typename ForwardIterator_T > inline void get( const Vector3<real_t> & position, ForwardIterator_T interpolationResultBegin )
 *  template< typename ForwardIterator_T > inline void get( const real_t & x, const real_t & y, const real_t & z, ForwardIterator_T interpolationResultBegin )
 *
 * See TrilinearFieldInterpolator for an example implementation.
 *
 * A field interpolator is aware of the flag field (FlagField_T) and uses only values that contain flags from a given mask.
 *
 * Field Interpolators can be used to sample the underlying field at certain positions.
 * E.g. the fluid velocity can be interpolated to a desired global position from a velocity field.
 *
 */
//**********************************************************************************************************************
template< typename Interpolator_T, typename FlagField_T >
class InterpolatorHandling : public blockforest::AlwaysInitializeBlockDataHandling< Interpolator_T >
{
public:

   InterpolatorHandling( const weak_ptr<StructuredBlockStorage> & blockStorage,
                         const ConstBlockDataID & interpolatedFieldID,
                         const ConstBlockDataID & flagFieldID,
                         const Set< FlagUID > & cellsToEvaluate ) :
   blockStorage_( blockStorage ), interpolatedFieldID_( interpolatedFieldID ), flagFieldID_( flagFieldID ), cellsToEvaluate_( cellsToEvaluate )
   {}

   Interpolator_T * initialize( IBlock * const block ) override
   {
      typedef typename Interpolator_T::BaseField_T InterpolatedField_T;
      typedef typename FlagField_T::flag_t flag_t;
      const InterpolatedField_T * interpolatedField = block->getData< InterpolatedField_T >( interpolatedFieldID_ );
      const FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );

      WALBERLA_ASSERT_NOT_NULLPTR( interpolatedField );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );

      const flag_t evaluationMask = flagField->getMask( cellsToEvaluate_ );

      return new Interpolator_T( blockStorage_, *block, *interpolatedField, *flagField, evaluationMask );
   }

private:

   weak_ptr<StructuredBlockStorage> blockStorage_;
   ConstBlockDataID interpolatedFieldID_;
   ConstBlockDataID flagFieldID_;
   Set< FlagUID > cellsToEvaluate_;
   
}; // class InterpolatorHandling



template< typename Interpolator_T, typename FlagField_T  >
inline BlockDataID addFieldInterpolator( const shared_ptr< StructuredBlockStorage > & blocks,
                                         const ConstBlockDataID & interpolatedFieldID,
                                         const ConstBlockDataID & flagFieldID,
                                         const Set< FlagUID > & cellsToEvaluate,
                                         const std::string & identifier = std::string(),
                                         const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                         const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< InterpolatorHandling< Interpolator_T, FlagField_T > >( blocks, interpolatedFieldID, flagFieldID, cellsToEvaluate ), identifier, requiredSelectors, incompatibleSelectors );
}


} // namespace field
} // namespace walberla
