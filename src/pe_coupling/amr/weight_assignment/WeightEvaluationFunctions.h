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
//! \file WeightEvaluationFunctions.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe_coupling/amr/BlockInfo.h"

namespace walberla {
namespace pe_coupling {
namespace amr {


/*
 * Examples of weight evaluation functions, useful for coupled LBM-PE simulations:
 *  - defaultWeightEvaluationFunction: weight is just the number of cells (thus constant on each block)
 */

real_t defaultWeightEvaluationFunction(const BlockInfo& blockInfo);

} // namespace amr
} // namespace pe_coupling
} // namespace walberla
