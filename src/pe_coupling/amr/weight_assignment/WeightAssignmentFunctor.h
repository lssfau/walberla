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

#include "pe_coupling/amr/BlockInfo.h"
#include "pe_coupling/amr/InfoCollection.h"

#include "blockforest/loadbalancing/PODPhantomData.h"

#include <functional>

namespace walberla {
namespace pe_coupling {
namespace amr {

class WeightAssignmentFunctor
{
public:
   using PhantomBlockWeight = walberla::blockforest::PODPhantomWeight<double>;

   WeightAssignmentFunctor( const shared_ptr<InfoCollection>& ic,
                            const std::function<real_t(const BlockInfo&)> & weightEvaluationFct ) :
         ic_(ic), weightEvaluationFct_(weightEvaluationFct) {}

   void operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & );

   inline void setWeightEvaluationFct( const std::function<real_t(const BlockInfo &)> & weightEvaluationFct ) { weightEvaluationFct_ = weightEvaluationFct;}

   inline void   setBlockBaseWeight( const real_t blockBaseWeight ){blockBaseWeight_ = blockBaseWeight;}
   inline real_t getBlockBaseWeight() const { return blockBaseWeight_;}

private:
   shared_ptr<InfoCollection> ic_;
   std::function<real_t(const BlockInfo&)> weightEvaluationFct_;
   real_t blockBaseWeight_ = real_t(1);
};

} // namespace amr
} // namespace pe_coupling
} // namespace walberla
