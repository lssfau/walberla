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
//! \file TimingNode.cpp
//! \ingroup core
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "TimingNode.h"



namespace walberla {
namespace timing {

// Explicit instantiation
template struct TimingNode<WcPolicy>;
template struct TimingNode<CpuPolicy>;

} // namespace timing
} // namespace walberla
