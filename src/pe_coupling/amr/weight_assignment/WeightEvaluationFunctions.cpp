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
//! \file WeightEvaluationFunctions.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "pe_coupling/amr/weight_assignment/WeightEvaluationFunctions.h"

namespace walberla {
namespace pe_coupling {
namespace amr {


real_t defaultWeightEvaluationFunction(const BlockInfo& blockInfo)
{
   return real_c(blockInfo.numberOfCells);
}

} // namespace amr
} // namespace pe_coupling
} // namespace walberla
