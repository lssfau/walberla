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
//! \file BlockDataHandling.cpp
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockDataHandling.h"
#include "BlockStorage.h"


namespace walberla {
namespace domain_decomposition {
namespace internal {


BlockDataHandlingAdder::operator BlockDataID()
{
   return storage_.addBlockData( dataHandling_, identifier_ );
}


} // namespace internal
} // namespace domain_decomposition
} // namespace walberla
