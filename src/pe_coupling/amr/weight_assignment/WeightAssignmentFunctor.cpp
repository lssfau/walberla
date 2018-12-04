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
//! \file WeightAssignmentFunctor.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "WeightAssignmentFunctor.h"

namespace walberla {
namespace pe_coupling {
namespace amr {

void WeightAssignmentFunctor::operator()( std::vector< std::pair< const PhantomBlock *, walberla::any > > & blockData, const PhantomBlockForest & )
{
   for (auto &it : blockData) {
      const PhantomBlock * block = it.first;
      //only change of one level is supported!
      WALBERLA_ASSERT_LESS( std::abs(int_c(block->getLevel()) - int_c(block->getSourceLevel())), 2 );

      BlockInfo blockInfo;
      pe_coupling::getBlockInfoFromInfoCollection(block, ic_, blockInfo);

      real_t blockWeight = std::max(weightEvaluationFct_(blockInfo), blockBaseWeight_);

      it.second = PhantomBlockWeight( double_c( blockWeight ) );

   }
}

} // namespace amr
} // namespace pe_coupling
} // namespace walberla
