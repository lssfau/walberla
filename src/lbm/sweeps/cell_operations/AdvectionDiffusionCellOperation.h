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
//! \file AdvectionDiffusionCellOperation.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once



namespace walberla {
namespace lbm {

template< typename LM_AdvDiff, typename LM_Hydro, class Enable = void >
class AdvectionDiffusionCellOperation;

} // namespace lbm
} // namespace walberla

#include "lbm/srt/cell_operations/AdvectionDiffusionCellOperation.impl.h"
//#include "lbm/trt/cell_operations/AdvectionDiffusionCellOperation.impl.h"
//#include "lbm/mrt/cell_operations/AdvectionDiffusionCellOperation.impl.h"
