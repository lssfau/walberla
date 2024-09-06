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
//! \file InitializerFunctions.h
//! \author Ravi Ayyala Somayajula <ravi.k.ayyala@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/all.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"

#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

// #include "python_coupling/DictWrapper.h"
#include "GeneralInfoHeader.h"
#include "InitializerFunctions.h"
#pragma once

namespace walberla
{

void initConcentrationField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& ConcentrationFieldID,
                            const math::AABB& domainAABB, Vector3< uint_t > domainSize)
{
   const real_t radius = real_c(domainSize[1] / 10);
   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);

      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(ConcentrationField, {
         Cell globalCell;
         const auto cellAABB = blocks->getBlockLocalCellAABB(block, Cell(x, y, z));
         auto cellCenter     = cellAABB.center();
         const real_t posX   = x; // cellCenter[0];
         const real_t posY   = y; // cellCenter[1];
         const real_t posZ   = z; // cellCenter[2];

         real_t distance =
            real_c(sqrt(pow((domainAABB.center()[0] - posX), 2) + pow((domainAABB.center()[1] - posY), 2) +
                        pow((domainAABB.center()[2] - posZ), 2)));

         if (distance <= radius) { ConcentrationField->get(x, y, z) = real_t(1.05); }
         else { ConcentrationField->get(x, y, z) = real_t(1); }
      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

} // initConcentrationField

void initConcentrationFieldGaussian(const shared_ptr< StructuredBlockStorage >& blocks,
                                    BlockDataID& ConcentrationFieldID, const math::AABB& domainAABB,
                                    Vector3< uint_t > domainSize, const real_t sigma_0, const real_t sigma_D,
                                    const Vector3< real_t > uInflow, const Vector3< real_t > x_0)
{
   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);
      Block& b                = dynamic_cast< Block& >(block);
      uint_t level            = b.getLevel();
      CellInterval xyz        = ConcentrationField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);

         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         // WALBERLA_LOG_INFO("posx " << pos[0]);
         ConcentrationField->get(*cellIt) = std::exp(
            -(std::pow((pos[0] - x_0[0]), 2) + std::pow((pos[1] - x_0[1]), 2) /*+ std::pow((posZ-x_0[2]),2)*/) /
            (2 * sigma_0 * sigma_0));
         ConcentrationField->get(*cellIt) = std::max(ConcentrationField->get(*cellIt), 1e-15);
      }
   }
}

void initFluidField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& FluidFieldID,
                    const Vector3< real_t > uInflow)
{
   for (auto& block : *blocks)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Velocity init reached here");
      auto FluidVelocityField = block.getData< VelocityField_fluid_T >(FluidFieldID);

      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(FluidVelocityField, {
         Cell globalCell;
         const auto cellAABB                 = blocks->getBlockLocalCellAABB(block, Cell(x, y, z));
         FluidVelocityField->get(x, y, z, 0) = uInflow[0];
         FluidVelocityField->get(x, y, z, 1) = uInflow[1];
         // FluidVelocityField->get(x,y,z,2) = uInflow[2];
      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }
}

void analyticalSolGaussian(const shared_ptr< StructuredBlockStorage >& blocks,
                           BlockDataID& AnalyticalConcentrationFieldID, const math::AABB& domainAABB,
                           Vector3< uint_t > domainSize, const real_t sigma_0, const real_t diffusivity,
                           const Vector3< real_t > uInflow, const Vector3< real_t > x_0, const real_t time)
{
   real_t sigma_D2 = (2 * diffusivity * time);

   for (auto& block : *blocks)
   {
      auto AnalyticalConcentrationField = block.getData< DensityField_concentration_T >(AnalyticalConcentrationFieldID);
      Block& b                          = dynamic_cast< Block& >(block);
      uint_t level                      = b.getLevel();
      CellInterval xyz                  = AnalyticalConcentrationField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, *cellIt);

         Vector3< real_t > pos = blocks->getCellCenter(globalCell, level);
         // WALBERLA_LOG_INFO("posx " << pos[0]);
         real_t prefactor = (sigma_0 * sigma_0) / (sigma_0 * sigma_0 + sigma_D2);
         AnalyticalConcentrationField->get(*cellIt) =
            prefactor * std::exp(-(std::pow((pos[0] - x_0[0] - uInflow[0] * time), 2) +
                                   std::pow((pos[1] - x_0[1] - uInflow[1] * time),
                                            2) /*+ std::pow((posZ - x_0[2] - uInflow[2]*time),2)*/) /
                                 (2 * (sigma_0 * sigma_0 + sigma_D2)));
      }
   }
}

std::vector< real_t > computeErrorL2(const shared_ptr< StructuredBlockStorage >& blocks,
                                     BlockDataID& NumericalSolFieldID, BlockDataID& AnalyticalSolFieldID,
                                     BlockDataID& ErrorFieldID, const math::AABB& domainAABB)
{
   real_t Linf{ 0.0 };
   real_t L1{ 0.0 };
   real_t L2{ 0.0 };
   real_t analytical_squared{ 0.0 };
   uint_t cells{ 0 };
   std::vector< real_t > Errors{ 0, 0, 0 };

   for (auto block = blocks->begin(); block != blocks->end(); ++block)
   {
      auto numericalSolutionField  = block->getData< DensityField_concentration_T >(NumericalSolFieldID);
      auto analyticalSolutionField = block->getData< DensityField_concentration_T >(AnalyticalSolFieldID);
      auto ErrorField              = block->getData< DensityField_concentration_T >(ErrorFieldID);
      Block& b                     = dynamic_cast< Block& >(*block);
      uint_t level                 = b.getLevel();
      CellInterval xyz             = analyticalSolutionField->xyzSize();

      for (auto cellIt = xyz.begin(); cellIt != xyz.end(); ++cellIt)
      {
         real_t currErr = (numericalSolutionField->get(*cellIt) - analyticalSolutionField->get(*cellIt));

         // potentialErrorField ->get(x, y, z) = real_c(fabs(currErr));
         ErrorField->get(*cellIt) = fabs(currErr);
         L1 += std::abs(currErr);
         L2 += currErr * currErr;
         analytical_squared += analyticalSolutionField->get(*cellIt) * analyticalSolutionField->get(*cellIt);
         Linf = std::max(Linf, std::abs(currErr));
         cells += 1;
      }
   }
   mpi::allReduceInplace(L1, mpi::SUM);
   mpi::allReduceInplace(L2, mpi::SUM);
   mpi::allReduceInplace(analytical_squared, mpi::SUM);
   mpi::allReduceInplace(Linf, mpi::MAX);
   mpi::allReduceInplace(cells, mpi::SUM);
   Errors[0] = Linf;
   Errors[1] = L1 / cells;
   Errors[2] = std::sqrt(L2 / analytical_squared);
   return Errors;
}

} // namespace walberla