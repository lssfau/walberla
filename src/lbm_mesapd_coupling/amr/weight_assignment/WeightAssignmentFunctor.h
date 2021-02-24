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
//! \file WeightAssignmentFunctor.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/loadbalancing/PODPhantomData.h"
#include "blockforest/PhantomBlockForest.h"
#include "blockforest/PhantomBlock.h"

#include <functional>

namespace walberla {
namespace lbm_mesapd_coupling  {
namespace amr {



class WeightAssignmentFunctor
{
public:
   using PhantomBlockWeight = walberla::blockforest::PODPhantomWeight<double>;
   using WeightEvaluationFct = std::function<real_t(const PhantomBlock *)>;

   explicit WeightAssignmentFunctor( const WeightEvaluationFct & weightEvaluationFct ) :
         weightEvaluationFct_(weightEvaluationFct) {}

   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & )
   {
      for (auto &it : blockData) {
         const PhantomBlock * block = it.first;
         //only change of one level is supported!
         WALBERLA_ASSERT_LESS( std::abs(int_c(block->getLevel()) - int_c(block->getSourceLevel())), 2 );

         real_t blockWeight = std::max(weightEvaluationFct_(block), blockBaseWeight_);
         it.second = PhantomBlockWeight( double_c( blockWeight ) );
      }
   }

   inline void   setBlockBaseWeight( const real_t blockBaseWeight ) { blockBaseWeight_ = blockBaseWeight; }
   inline real_t getBlockBaseWeight() const { return blockBaseWeight_; }

private:
   WeightEvaluationFct weightEvaluationFct_;
   real_t blockBaseWeight_ = real_t(1);
};


} // namespace amr
} // namespace lbm_mesapd_coupling
} // namespace walberla
