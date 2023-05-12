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
//! \file EvaluationFilter.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/Set.h"
#include "domain_decomposition/IBlock.h"

#include "FlagUID.h"



namespace walberla {
namespace field {



//**********************************************************************************************************************
/*!
*   \section docEvaluationFilter Evaluation Filter
*
*   Evaluation filters always work on block local cell coordinates.
*
*   The concept for evaluation filters looks like as follows (class with two member functions):
*
*   1. void operator()( const IBlock & block )
*      -> called every time a new block is processed
*
*   2. bool operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
*      -> called every time a new cell is evaluated
*         might be called in parallel (must be threat-safe!) and must return true or false
*         (true selects the cell for further processing, false skips the cell in question)
*/
//**********************************************************************************************************************



class DefaultEvaluationFilter
{
public:
   void operator()( const IBlock & ) {}
   bool operator()( const cell_idx_t, const cell_idx_t, const cell_idx_t ) const { return true; }
};



template< typename FlagField_T >
class FlagFieldEvaluationFilter
{
public:

   FlagFieldEvaluationFilter( const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate ) :
      flagFieldId_( flagFieldId ), flagField_( nullptr ), cellsToEvaluate_( cellsToEvaluate ) {}

   void operator()( const IBlock & block )
   {
      flagField_ = block.template getData< const FlagField_T >( flagFieldId_ );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField_ )
      evaluationMask_ = flagField_->getMask( cellsToEvaluate_ );
   }

   bool operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( flagField_ )
      return flagField_->isPartOfMaskSet( x, y, z, evaluationMask_ );
   }

private:

   ConstBlockDataID flagFieldId_;
   FlagField_T const * flagField_;

   Set< FlagUID > cellsToEvaluate_;
   typename FlagField_T::flag_t evaluationMask_;
};



} // namespace field
} // namespace walberla
