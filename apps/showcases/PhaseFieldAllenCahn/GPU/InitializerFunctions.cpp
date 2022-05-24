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
                           const real_t R                         = real_c(10.0),
                           const Vector3< real_t > bubbleMidPoint = Vector3< real_t >(0.0, 0.0, 0.0),
                           const bool bubble = true, const real_t W = real_c(5.0))
{
   for (auto& block : *blocks)
   {
      auto phaseField = block.getData< PhaseField_T >(phaseFieldID);
      // clang-format off
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(phaseField, Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         real_t Ri = sqrt((real_c(real_c(globalCell[0])) - bubbleMidPoint[0]) * (real_c(real_c(globalCell[0])) - bubbleMidPoint[0]) +
                          (real_c(globalCell[1]) - bubbleMidPoint[1]) * (real_c(globalCell[1]) - bubbleMidPoint[1]) +
                          (real_c(globalCell[2]) - bubbleMidPoint[2]) * (real_c(globalCell[2]) - bubbleMidPoint[2]));
         if (bubble) phaseField->get(x, y, z) = real_c(0.5) + real_c(0.5) * real_c(tanh(real_c(2.0) * (Ri - R) / W));
         else phaseField->get(x, y, z)        = real_c(0.5) - real_c(0.5) * tanh(real_c(2.0) * (Ri - R) / W);
      )
      // clang-format on
   }
}

void init_Taylor_bubble(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID,
                        const real_t D = 5.0, const real_t H = 2.0, const real_t DT = 20.0, const real_t Donut_x0 = 40.0)
{
   auto Mx = real_c(blocks->getDomainCellBB().xMax()) / real_c(2.0);
   auto Mz = real_c(blocks->getDomainCellBB().zMax()) / real_c(2.0);

   for (auto& block : *blocks)
   {
      auto phaseField = block.getData< PhaseField_T >(phaseFieldID);
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
         phaseField, Cell globalCell; blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));

         real_t Ri = D * real_c(sqrt(pow(H, 2) - pow(DT - sqrt(pow(real_c(real_c(globalCell[0])) - Mx, 2) + pow(real_c(globalCell[2]) - Mz, 2)), 2)));

         real_t shifter           = atan2((real_c(globalCell[2]) - Mz), (real_c(globalCell[0]) - Mx));
         if (shifter < 0) shifter = shifter + 2 * math::pi;
         if ((real_c(globalCell[1]) < Donut_x0 + Ri * sin(shifter / 2.0)) && (real_c(globalCell[1]) > Donut_x0 - Ri)) {
            phaseField->get(x, y, z) = 0.0;
         } else { phaseField->get(x, y, z) = 1.0; })
   }
}

void init_bubble_field(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID, real_t R,
                       real_t W = 5.0)
{
   Vector3< real_t > bubbleMidPoint;

   auto X = real_c(blocks->getDomainCellBB().xMax());
   auto Y = real_c(blocks->getDomainCellBB().yMax());
   auto Z = real_c(blocks->getDomainCellBB().zMax());

   // 20 percent from the top are filled with the gas phase
   real_t gas_top = Y - Y / real_c(5.0);

   // Diameter of the bubble
   real_t D = R * real_c(2.0);

   // distance in between the bubbles
   auto dist = real_c(4.0);
   auto nx  = uint_c(floor(X / (D + dist * W)));
   auto nz  = uint_c(floor(Z / (D + dist * W)));

   // fluctuation of the bubble radii
   std::vector< std::vector< real_t > > fluctuation_radius(nx, std::vector< real_t >(nz, real_c(0.0)));
   std::vector< std::vector< real_t > > fluctuation_pos(nx, std::vector< real_t >(nz, real_c(0.0)));

   real_t max_fluctuation_radius = R / real_c(5.0);
   real_t max_fluctuation_pos    = (dist * W) / real_c(3.0);
   for (unsigned int i = 0; i < nx; ++i)
   {
      for (unsigned int j = 0; j < nz; ++j)
      {
         fluctuation_radius[i][j] = math::realRandom< real_t >(-max_fluctuation_radius, max_fluctuation_radius);
         fluctuation_pos[i][j]    = math::realRandom< real_t >(-max_fluctuation_pos, max_fluctuation_pos);
      }
   }

   for (auto& block : *blocks)
   {
      auto phaseField = block.getData< PhaseField_T >(phaseFieldID);
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
         phaseField, Cell globalCell; blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         for (unsigned int i = 0; i < nx; ++i) {
            for (unsigned int j = 0; j < nz; ++j)
            {
               bubbleMidPoint[0] = real_c(i + 1) * (D + (dist * W)) - (D + (dist * W)) / real_c(2.0) + fluctuation_pos[i][j];
               bubbleMidPoint[1] = R + W + 4;
               bubbleMidPoint[2] = real_c(j + 1) * (D + (dist * W)) - (D + (dist * W)) / real_c(2.0) + fluctuation_pos[i][j];

               real_t Ri = sqrt((real_c(globalCell[0]) - bubbleMidPoint[0]) * (real_c(globalCell[0]) - bubbleMidPoint[0]) +
                                (real_c(globalCell[1]) - bubbleMidPoint[1]) * (real_c(globalCell[1]) - bubbleMidPoint[1]) +
                                (real_c(globalCell[2]) - bubbleMidPoint[2]) * (real_c(globalCell[2]) - bubbleMidPoint[2]));
               if (real_c(globalCell[0]) >= real_c(i) * (D + dist * W) && real_c(globalCell[0]) <= real_c(i + 1) * (D + dist * W) &&
                   real_c(globalCell[2]) >= real_c(j) * (D + dist * W) && real_c(globalCell[2]) <= real_c(j + 1) * (D + dist * W))
                  phaseField->get(x, y, z) = real_c(0.5) + real_c(0.5) * real_c(tanh(2.0 * (Ri - (R - fluctuation_radius[i][j])) / W));

               if (real_c(globalCell[0]) > real_c(nx) * (D + dist * W)) phaseField->get(x, y, z) = 1.0;
               if (real_c(globalCell[2]) > real_c(nz) * (D + dist * W)) phaseField->get(x, y, z) = 1.0;
            }
         }

         if (real_c(globalCell[1]) > gas_top) {
            phaseField->get(x, y, z) = real_c(0.5) + real_c(0.5) * real_c(tanh(2.0 * (gas_top + 10 - real_c(globalCell[1])) / W));
         })
   }
}

void initPhaseField_RTI(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID,
                        const real_t W = 5.0, const bool pipe = true)
{
   auto X              = real_c(blocks->getDomainCellBB().xMax());
   auto Z              = real_c(blocks->getDomainCellBB().zMax());
   auto halfY          = real_c((blocks->getDomainCellBB().yMax())) / real_c(2.0);
   auto perturbation = real_c(0.05);

   if (pipe)
   {
      for (auto& block : *blocks)
      {
         auto phaseField = block.getData< PhaseField_T >(phaseFieldID);
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
            phaseField, Cell globalCell; blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
            real_t R     = sqrt((real_c(globalCell[0]) - X / 2) * (real_c(globalCell[0]) - X / 2) +
                                (real_c(globalCell[2]) - Z / 2) * (real_c(globalCell[2]) - Z / 2));
            if (R > X) R = X; real_t tmp = perturbation * X * real_c(cos((2.0 * math::pi * R) / X));
            phaseField->get(x, y, z)     = real_c(0.5) + real_c(0.5) * real_c(tanh(((real_c(globalCell[1]) - halfY) + tmp)) / (W / 2.0));)
      }
   }
   else
   {
      for (auto& block : *blocks)
      {
         auto phaseField = block.getData< PhaseField_T >(phaseFieldID);
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
            phaseField, Cell globalCell; blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
            real_t tmp = perturbation * X *
                         real_c(cos((real_c(2.0) * math::pi * real_c(globalCell[0]))) / X + cos((2.0 * math::pi * real_c(globalCell[2])) / X));
            phaseField->get(x, y, z) = real_c(0.5) + real_c(0.5) * real_c(tanh(((real_c(globalCell[1]) - halfY) - tmp)) / (W / 2.0));)
      }
   }
}

void initTubeWithCylinder(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID flagFieldID,
                          field::FlagUID boundaryFlagUID, real_t const R_in, real_t const eccentricity,
                          real_t const start_transition, real_t const length_transition,
                          bool const eccentricity_or_pipe_ratio)
{
   if (eccentricity_or_pipe_ratio)
   {
      auto Mx = real_c(blocks->getDomainCellBB().xMax()) / real_c(2.0);
      auto Mz = real_c(blocks->getDomainCellBB().zMax()) / real_c(2.0);

      auto R_outer = real_c(blocks->getDomainCellBB().xMax()) / real_c(2.0) + real_c(1.0);

      real_t const shift = eccentricity * Mx / real_c(2.0);

      for (auto& block : *blocks)
      {
         auto flagField    = block.template getData< FlagField_T >(flagFieldID);
         auto boundaryFlag = flagField->getOrRegisterFlag(boundaryFlagUID);
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flagField, Cell globalCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
            real_t R1;
            if (real_c(globalCell[1]) <= start_transition) {
               R1 = sqrt((real_c(globalCell[0]) - Mx) * (real_c(globalCell[0]) - Mx) + (real_c(globalCell[2]) - Mz) * (real_c(globalCell[2]) - Mz));
            } else if (real_c(globalCell[1]) > start_transition && real_c(globalCell[1]) < start_transition + length_transition) {
               real_t tmp       = math::pi * (real_c(globalCell[1]) - start_transition) / (length_transition);
               real_t shift_tmp = shift * real_c(0.5) * (real_c(1.0) - real_c(cos(tmp)));
               R1               = sqrt((real_c(globalCell[0]) - Mx - shift_tmp) * (real_c(globalCell[0]) - Mx - shift_tmp) +
                                       (real_c(globalCell[2]) - Mz) * (real_c(globalCell[2]) - Mz));
            } else {
               R1 = sqrt((real_c(globalCell[0]) - Mx - shift) * (real_c(globalCell[0]) - Mx - shift) +
                         (real_c(globalCell[2]) - Mz) * (real_c(globalCell[2]) - Mz));
            }

            real_t R2 = sqrt((real_c(globalCell[0]) - Mx) * (real_c(globalCell[0]) - Mx) + (real_c(globalCell[2]) - Mz) * (real_c(globalCell[2]) - Mz));
            if (R1 < R_in) addFlag(flagField->get(x, y, z), boundaryFlag);
            if (R2 > R_outer) addFlag(flagField->get(x, y, z), boundaryFlag);)
      }
   }
   else
   {
      auto Mx = real_c(blocks->getDomainCellBB().xMax()) / real_c(2.0);
      auto Mz = real_c(blocks->getDomainCellBB().zMax()) / real_c(2.0);

      auto R_outer = real_c(blocks->getDomainCellBB().xMax()) / real_c(2.0) + real_c(1.0);

      auto const shift = real_c(eccentricity * R_in);
      real_t R_tmp;

      for (auto& block : *blocks)
      {
         auto flagField    = block.template getData< FlagField_T >(flagFieldID);
         auto boundaryFlag = flagField->getOrRegisterFlag(boundaryFlagUID);
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
            flagField, Cell globalCell; blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
            if (real_c(globalCell[1]) <= start_transition) {
               R_tmp = R_in;
            } else if (real_c(globalCell[1]) > start_transition && real_c(globalCell[1]) < start_transition + length_transition) {
               real_t tmp       = math::pi * (real_c(globalCell[1]) - start_transition) / (length_transition);
               real_t shift_tmp = shift * real_c(0.5) * (real_c(1.0) - real_c(cos(tmp)));
               R_tmp = R_in + shift_tmp;
            } else {
               R_tmp = R_in + shift;
            }

            real_t R2 = real_c(sqrt((real_c(globalCell[0]) - Mx) * (real_c(globalCell[0]) - Mx) + (real_c(globalCell[2]) - Mz) * (real_c(globalCell[2]) - Mz)));
            if (R2 < R_tmp) addFlag(flagField->get(x, y, z), boundaryFlag);
            if (R2 > R_outer) addFlag(flagField->get(x, y, z), boundaryFlag);)
      }
   }
}
} // namespace walberla
