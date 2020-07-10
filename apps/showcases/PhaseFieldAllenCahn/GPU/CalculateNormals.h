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
//! \file CalculateNormals.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/StructuredBlockForest.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"

#include "field/FlagField.h"

namespace walberla
{
void calculate_normals(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID normalsFieldID,
                       ConstBlockDataID flagFieldID, FlagUID fluidFlagUID, FlagUID boundaryFlagUID);
}