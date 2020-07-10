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
//! \file CalculateNormals.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "CalculateNormals.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"

#include "field/FlagField.h"

namespace walberla
{
using FlagField_T    = FlagField< uint8_t >;
using NormalsField_T = GhostLayerField< int8_t, 3 >;

void calculate_normals(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID normalsFieldID,
                       ConstBlockDataID flagFieldID, FlagUID domainFlagUID, FlagUID boundaryFlagUID)
{
   for (auto& block : *blocks)
   {
      CellInterval globalCellBB = blocks->getBlockCellBB(block);
      CellInterval blockLocalCellBB;
      blocks->transformGlobalToBlockLocalCellInterval(blockLocalCellBB, block, globalCellBB);

      auto* normalsField = block.getData< NormalsField_T >(normalsFieldID);
      auto* flagField    = block.getData< FlagField_T >(flagFieldID);
      auto boundaryFlag  = flagField->getFlag(boundaryFlagUID);
      auto domainFlag    = flagField->getFlag(domainFlagUID);

      // clang-format off
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(normalsField,

         if( x < blockLocalCellBB.xMax() ){
            if(flagField->get(x, y, z) == boundaryFlag && flagField->get(x + 1, y, z) == domainFlag)
               normalsField->get(x, y, z, 0) = 1;
         }

         if( x > blockLocalCellBB.xMin() ){
            if(flagField->get(x, y, z) == boundaryFlag && flagField->get(x - 1, y, z) == domainFlag)
               normalsField->get(x, y, z, 0) = - 1;
         }

         if( y < blockLocalCellBB.yMax() ){
            if(flagField->get(x, y, z) == boundaryFlag && flagField->get(x, y + 1, z) == domainFlag)
               normalsField->get(x, y, z, 1) = 1;
         }

         if( y > blockLocalCellBB.yMin() ){
            if(flagField->get(x, y, z) == boundaryFlag && flagField->get(x, y - 1, z) == domainFlag)
               normalsField->get(x, y, z, 1) = - 1;
         }

         if( z < blockLocalCellBB.zMax() ){
            if(flagField->get(x, y, z) == boundaryFlag && flagField->get(x, y, z + 1) == domainFlag)
               normalsField->get(x, y, z, 2) = 1;
         }

         if( z > blockLocalCellBB.zMin() ){
            if(flagField->get(x, y, z) == boundaryFlag && flagField->get(x, y, z - 1) == domainFlag)
               normalsField->get(x, y, z, 2) = - 1;
         }

      )
      // clang-format on
   }
}
} // namespace walberla
