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
//! \file InitializerFunctions.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "InitializerFunctions.h"

namespace walberla
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using FlagField_T  = FlagField< uint8_t >;

void initPhaseFieldDroplet(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID,
                           const real_t dropletRadius                         = real_c(10.0),
                           const Vector3< real_t > dropletMidPoint = Vector3< real_t >(0.0, 0.0, 0.0),
                           const real_t W)
{
   for (auto& block : *blocks)
   {
      auto phaseField = block.getData< ScalarField_T >(phaseFieldID);
      // clang-format off
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(phaseField, Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         real_t Ri = sqrt((real_c(globalCell[0]) - dropletMidPoint[0]) * (real_c(globalCell[0]) - dropletMidPoint[0]) +
                          (real_c(globalCell[1]) - dropletMidPoint[1]) * (real_c(globalCell[1]) - dropletMidPoint[1]) +
                          (real_c(globalCell[2]) - dropletMidPoint[2]) * (real_c(globalCell[2]) - dropletMidPoint[2]));
         phaseField->get(x, y, z)        = real_c(0.5) - real_c(0.5) * real_c(tanh(2.0 * (Ri - dropletRadius) / W));
      )
      // clang-format on
   }
}

void initMicroChannel(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID, BlockDataID temperatureFieldID,
                           const real_t Th, const real_t T0, const real_t Tc, const real_t W)
{
   auto halfX          = real_c(blocks->getDomainCellBB().xMax()) / real_c(2.0);
   auto halfY          = real_c((blocks->getDomainCellBB().yMax())) / real_c(2.0);
   for (auto& block : *blocks)
   {
      auto phaseField = block.getData< ScalarField_T >(phaseFieldID);
      // clang-format off
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(phaseField, Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         phaseField->get(x, y, z)     = real_c(0.5) + real_c(0.5) * real_c(tanh(   (real_c(globalCell[1]) - halfY) / (W / real_c(2.0))   ));
      )
      // clang-format on

      auto temperatureField = block.getData< ScalarField_T >(temperatureFieldID);
      // clang-format off
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(temperatureField, Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         if(globalCell[1] == -1)
         {
            auto X = (globalCell[0] < 0) ? 0.0 : globalCell[0];
            temperatureField->get(x, y, z)     = Th + T0 * cos( (math::pi / halfX) * (X - halfX) );
         }
         else
         {
            temperatureField->get(x, y, z)     = Tc;
         }
      )
      // clang-format on
   }
}

} // namespace walberla
