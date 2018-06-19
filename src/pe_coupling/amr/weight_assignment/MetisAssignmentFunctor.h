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

#include "pe_coupling/amr/InfoCollection.h"

#include "blockforest/loadbalancing/DynamicParMetis.h"

#include <functional>
#include <vector>

namespace walberla {
namespace pe_coupling {
namespace amr {

class MetisAssignmentFunctor
{
public:

   MetisAssignmentFunctor( shared_ptr<InfoCollection>& ic,
                           const std::function<real_t(const BlockInfo&)> & weightEvaluationFct,
                           const uint_t ncon = 1)
      : ic_(ic), ncon_(ncon), weightEvaluationFct_(ncon, weightEvaluationFct) {}

   MetisAssignmentFunctor( shared_ptr<InfoCollection>& ic,
                           const std::vector< std::function<real_t(const BlockInfo&)> > & weightEvaluationFct )
      : ic_(ic), ncon_(weightEvaluationFct.size()), weightEvaluationFct_(weightEvaluationFct) {}

   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & phantomBlockForest);

   inline void setWeightEvaluationFct( const std::function<real_t(const BlockInfo &)> & weightEvaluationFct, const uint_t con = 0 ) { weightEvaluationFct_[con] = weightEvaluationFct; }
   inline void setWeightEvaluationFcts( const std::vector< std::function<real_t(const BlockInfo &)> > & weightEvaluationFct) { weightEvaluationFct_ = weightEvaluationFct; }

   uint_t getNcon() const { return ncon_; }

   inline void   setBlockBaseWeight( const real_t blockBaseWeight ){ blockBaseWeight_ = blockBaseWeight; }
   inline real_t getBlockBaseWeight() const { return blockBaseWeight_; }

private:

   shared_ptr<InfoCollection> ic_;
   uint_t ncon_;
   std::vector< std::function<real_t(const BlockInfo&)> > weightEvaluationFct_;
   real_t blockBaseWeight_ = real_t(1);
};

} // namespace amr
} // namespace pe_coupling
} // namespace walberla

