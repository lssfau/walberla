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
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "python_coupling/DictWrapper.h"

namespace walberla
{
using PhaseField_T = GhostLayerField< real_t, 1 >;
using FlagField_T  = FlagField< uint8_t >;

void initPhaseField_sphere(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID,
                           const real_t R                         = 10,
                           const Vector3< real_t > bubbleMidPoint = Vector3< real_t >(0.0, 0.0, 0.0),
                           const bool bubble = true, const real_t W = 5)
{
   for (auto& block : *blocks)
   {
      auto phaseField = block.getData< PhaseField_T >(phaseFieldID);
      // clang-format off
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(phaseField, Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         real_t Ri = sqrt((globalCell[0] - bubbleMidPoint[0]) * (globalCell[0] - bubbleMidPoint[0]) +
                          (globalCell[1] - bubbleMidPoint[1]) * (globalCell[1] - bubbleMidPoint[1]) +
                          (globalCell[2] - bubbleMidPoint[2]) * (globalCell[2] - bubbleMidPoint[2]));
         if (bubble) phaseField->get(x, y, z) = 0.5 + 0.5 * tanh(2.0 * (Ri - R) / W);
         else phaseField->get(x, y, z)        = 0.5 - 0.5 * tanh(2.0 * (Ri - R) / W);
      )
      // clang-format on
   }
}

void initPhaseField_RTI(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID,
                        const real_t W = 5)
{
   auto X              = blocks->getDomainCellBB().xMax();
   auto halfY          = (blocks->getDomainCellBB().yMax()) / 2.0;
   double perturbation = 0.05;

   for (auto& block : *blocks)
   {
      auto phaseField = block.getData< PhaseField_T >(phaseFieldID);
      // clang-format off
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(phaseField, Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         real_t tmp =
         perturbation * X * (cos((2.0 * math::pi * globalCell[0]) / X) + cos((2.0 * math::pi * globalCell[2]) / X));
         phaseField->get(x, y, z) = 0.5 + 0.5 * tanh(((globalCell[1] - halfY) - tmp) / (W / 2.0));
      )
      // clang-format on
   }
}
} // namespace walberla
