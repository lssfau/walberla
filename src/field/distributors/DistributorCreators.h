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
//! \file DistributorCreators.h
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
/*! DistributorCreators
 *
 * \ingroup field
 *
 * Distributor_T: A distributor that has a constructor
 *  ( const weak_ptr<StructuredBlockStorage> & blockStorage, const IBlock & block, const BaseField_T & baseField,
 *    const FlagField_T & flagField, const flag_t & evaluationMask )
 * and distribution functions:
 *  template< typename ForwardIterator_T > inline void distribute( const Vector3<real_t> & position, ForwardIterator_T distributeValueBegin )
 *  template< typename ForwardIterator_T > inline void distribute( const real_t & x, const real_t & y, const real_t & z, ForwardIterator_T distributeValueBegin )
 *
 * See NearestNeighborDistributor for an example implementation.
 *
 * A distributor is aware of the flag field (FlagField_T) and distributes values only to cells flagged by a given mask.
 *
 * Distributors are used to spread a given value to the corresponding destination field.
 * E.g. if a certain force has to be applied at some specific position onto the fluid, a distributor can be used
 * to do so by distributing this force value (and conservation fo this force value is ensured) onto the force field.
 *
 */
//**********************************************************************************************************************
template< typename Distributor_T, typename FlagField_T >
class DistributorHandling : public blockforest::AlwaysInitializeBlockDataHandling< Distributor_T >
{
public:

   DistributorHandling( const weak_ptr<StructuredBlockStorage> & blockStorage,
                        const BlockDataID & distributionDestinationFieldID,
                        const ConstBlockDataID & flagFieldID,
                        const Set< FlagUID > & cellsToEvaluate ) :
   blockStorage_( blockStorage ), distributionDestinationFieldID_( distributionDestinationFieldID ), flagFieldID_( flagFieldID ), cellsToEvaluate_( cellsToEvaluate )
   {}

   Distributor_T * initialize( IBlock * const block ) override
   {
      typedef typename Distributor_T::BaseField_T DistributionDestinationField_T;
      typedef typename FlagField_T::flag_t flag_t;
      DistributionDestinationField_T * distributionDestinationField = block->getData< DistributionDestinationField_T >( distributionDestinationFieldID_ );
      const FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );

      WALBERLA_ASSERT_NOT_NULLPTR( distributionDestinationField );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );

      const flag_t evaluationMask = flagField->getMask( cellsToEvaluate_ );

      return new Distributor_T( blockStorage_, *block, *distributionDestinationField, *flagField, evaluationMask );
   }

private:

   weak_ptr<StructuredBlockStorage> blockStorage_;
   BlockDataID distributionDestinationFieldID_;
   ConstBlockDataID flagFieldID_;
   Set< FlagUID > cellsToEvaluate_;
   
}; // class DistributorHandling



template< typename Distributor_T, typename FlagField_T  >
inline BlockDataID addDistributor( const shared_ptr< StructuredBlockStorage > & blocks,
                                   const BlockDataID & distributionDestinationFieldID,
                                   const ConstBlockDataID & flagFieldID,
                                   const Set< FlagUID > & cellsToEvaluate,
                                   const std::string & identifier = std::string(),
                                   const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   return blocks->addBlockData( make_shared< DistributorHandling< Distributor_T, FlagField_T > >( blocks, distributionDestinationFieldID, flagFieldID, cellsToEvaluate ), identifier, requiredSelectors, incompatibleSelectors );
}


} // namespace field
} // namespace walberla
