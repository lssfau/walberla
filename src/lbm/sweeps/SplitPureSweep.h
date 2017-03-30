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
//! \file SplitPureSweep.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once



namespace walberla {
namespace lbm {

template< typename LatticeModel_T, class Enable = void >
class SplitPureSweep
{
   static_assert( never_true<LatticeModel_T>::value, "Instantiating 'lbm::SplitPureSweep' failed, possible reasons:\n"
                                                     " - For your current LB lattice model, there is yet no implementation for class 'lbm::SplitPureSweep'.\n"
                                                     " - You are providing more than just one template argument to class 'lbm::SplitPureSweep'.\n"
                                                     "   'lbm::SplitPureSweep' only needs the type of the lattice model - no further template arguments are required!" );
};

} // namespace lbm
} // namespace walberla

#include "lbm/srt/SplitPureSweep.impl.h"
#include "lbm/trt/SplitPureSweep.impl.h"
#include "lbm/mrt/SplitPureSweep.impl.h"
