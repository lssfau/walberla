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
//! \file SplitSweep.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once



namespace walberla {
namespace lbm {

template< typename LatticeModel_T, typename FlagField_T, class Enable = void >
class SplitSweep
{
   static_assert( never_true<LatticeModel_T>::value, "For your current LB lattice model, there is yet no implementation for class 'lbm::SplitSweep'!" );
};

} // namespace lbm
} // namespace walberla

#include "lbm/srt/SplitSweep.impl.h"
#include "lbm/trt/SplitSweep.impl.h"
#include "lbm/mrt/SplitSweep.impl.h"
