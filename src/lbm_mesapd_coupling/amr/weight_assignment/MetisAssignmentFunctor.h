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
//! \file MetisAssignmentFunctor.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/PhantomBlockForest.h"
#include "blockforest/PhantomBlock.h"
#include "blockforest/loadbalancing/DynamicParMetis.h"

#include <functional>
#include <vector>

namespace walberla {
namespace lbm_mesapd_coupling {
namespace amr {

class MetisAssignmentFunctor
{
public:

   using WeightEvaluationFct = std::function<real_t(const PhantomBlock *)>;

   explicit MetisAssignmentFunctor( const WeightEvaluationFct& weightEvaluationFct)
         : weightEvaluationFctVector_(1, weightEvaluationFct) {}

   explicit MetisAssignmentFunctor( const std::vector< WeightEvaluationFct > & weightEvaluationFctVector )
         : weightEvaluationFctVector_(weightEvaluationFctVector) {}

   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & phantomBlockForest);

   uint_t getNumberOfConstrains() const { return weightEvaluationFctVector_.size(); }

   inline void   setBlockBaseWeight( const real_t blockBaseWeight ){ blockBaseWeight_ = blockBaseWeight; }
   inline real_t getBlockBaseWeight() const { return blockBaseWeight_; }

   inline void setWeightMultiplicator( const real_t weightMultiplicator ){ weightMultiplicator_ = weightMultiplicator; }

private:
   std::vector< WeightEvaluationFct > weightEvaluationFctVector_;
   real_t blockBaseWeight_ = real_t(1);
   real_t weightMultiplicator_ = real_t(1);
};



} // namespace amr
} // namespace lbm_mesapd_coupling
} // namespace walberla

