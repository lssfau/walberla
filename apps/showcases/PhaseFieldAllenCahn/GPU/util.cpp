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
//! \file util.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "core/cell/Cell.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include <queue>
namespace walberla
{
using flag_t          = walberla::uint8_t;
using PhaseField_T    = GhostLayerField< real_t, 1 >;
using VelocityField_T = GhostLayerField< real_t, 3 >;
using FlagField_T     = FlagField< flag_t >;

void calc_total_velocity(const shared_ptr< StructuredBlockStorage >& blocks, std::array< real_t, 5 >& total_velocity,
                         BlockDataID phaseFieldID, BlockDataID velocityFieldID, ConstBlockDataID flagFieldID,
                         FlagUID fluidFlagUID)
{
   for (auto& block : *blocks)
   {
      auto phaseField = block.getData< PhaseField_T >(phaseFieldID);
      auto velField   = block.getData< VelocityField_T >(velocityFieldID);
      auto flagField  = block.getData< FlagField_T >(flagFieldID);
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
         phaseField,

         auto fluidFlag = flagField->getFlag(fluidFlagUID); if (flagField->get(x, y, z) == fluidFlag) {
            if (phaseField->get(x, y, z) < real_c(0.5))
            {
               real_t invC = 1 - phaseField->get(x, y, z);
               real_t U    = velField->get(x, y, z, 0);
               real_t V    = velField->get(x, y, z, 1);
               real_t W    = velField->get(x, y, z, 2);
               total_velocity[0] += real_c(sqrt((U * invC) * (U * invC) + (V * invC) * (V * invC) + (W * invC) * (W * invC)));
               total_velocity[1] += U * invC;
               total_velocity[2] += V * invC;
               total_velocity[3] += W * invC;
               total_velocity[4] += invC;
            }
         })
   }
}

void flood_fill(PhaseField_T& phaseField, VelocityField_T& velocityField, CellInterval boundingBox, real_t& volume,
                uint_t& nrOfCells, std::array< real_t, 3 >& center_of_mass, std::array< real_t, 4 >& total_velocity)
{
   Cell startCell(boundingBox.xSize() / 2, boundingBox.ySize() / 2, boundingBox.zSize() / 2);
   Field< bool, 1 > visit(phaseField.xSize(), phaseField.ySize(), phaseField.zSize(), false, field::fzyx);
   using namespace stencil;

   volume            = real_c(0.0);
   nrOfCells         = uint_c(0);
   center_of_mass[0] = real_c( 0.0);
   center_of_mass[1] = real_c(0.0);
   center_of_mass[2] = real_c(0.0);

   while (phaseField.get(startCell) > 0.5 && startCell.x() > 0)
      --startCell.x();

   if (phaseField.get(startCell) > 0.5) WALBERLA_ABORT("startCell for flood fill was not suitable")

   std::queue< Cell > cellQueue;
   cellQueue.push(startCell);
   visit.get(startCell) = true;

   real_t invC = 1 - phaseField.get(startCell);
   real_t v_U  = velocityField.get(startCell, 0);
   real_t v_V  = velocityField.get(startCell, 1);
   real_t v_W  = velocityField.get(startCell, 2);

   nrOfCells++;
   volume += invC;

   total_velocity[0] += real_c(sqrt((v_U * invC) * (v_U * invC) + (v_V * invC) * (v_V * invC) + (v_W * invC) * (v_W * invC)));
   total_velocity[1] += v_U * invC;
   total_velocity[2] += v_V * invC;
   total_velocity[3] += v_W * invC;

   center_of_mass[0] += real_c(startCell.x() + boundingBox.xMin());
   center_of_mass[1] += real_c(startCell.y() + boundingBox.yMin());
   center_of_mass[2] += real_c(startCell.z() + boundingBox.xMin());

   const int DIRS[6] = { N, S, E, W, T, B };

   CellInterval sizeInterval = phaseField.xyzSize();
   while (!cellQueue.empty())
   {
      Cell& cell = cellQueue.front();
      cellQueue.pop();

      for (int i : DIRS)
      {
         Cell neighborCell(cell.x() + cx[i], cell.y() + cy[i], cell.z() + cz[i]);

         if (!sizeInterval.contains(neighborCell)) { continue; }

         if (phaseField.get(neighborCell) < 0.5 && !visit.get(neighborCell))
         {
            invC = 1 - phaseField.get(neighborCell);
            v_U  = velocityField.get(neighborCell, 0);
            v_V  = velocityField.get(neighborCell, 1);
            v_W  = velocityField.get(neighborCell, 2);

            nrOfCells++;
            volume += invC;

            total_velocity[0] +=
               real_c(sqrt((v_U * invC) * (v_U * invC) + (v_V * invC) * (v_V * invC) + (v_W * invC) * (v_W * invC)));
            total_velocity[1] += v_U * invC;
            total_velocity[2] += v_V * invC;
            total_velocity[3] += v_W * invC;

            center_of_mass[0] += real_c(neighborCell.x() + boundingBox.xMin());
            center_of_mass[1] += real_c(neighborCell.y() + boundingBox.yMin());
            center_of_mass[2] += real_c(neighborCell.z() + boundingBox.xMin());

            visit.get(neighborCell) = true;
            cellQueue.push(neighborCell);
         }
      }
   }
   center_of_mass[0] = center_of_mass[0] / real_t(nrOfCells);
   center_of_mass[1] = center_of_mass[1] / real_t(nrOfCells);
   center_of_mass[2] = center_of_mass[2] / real_t(nrOfCells);
}
} // namespace walberla
