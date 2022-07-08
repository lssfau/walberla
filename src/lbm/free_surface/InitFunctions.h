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
//! \file InitFunctions.h
//! \ingroup free_surface
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Initialization functions.
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "domain_decomposition/BlockDataID.h"

#include <functional>

#include "FlagInfo.h"
#include "InterfaceFromFillLevel.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Initialize fill level with "value" in cells belonging to boundary and obstacles such that the bubble model does not
 * detect obstacles as gas cells.
 **********************************************************************************************************************/
template< typename BoundaryHandling_T, typename Stencil_T, typename ScalarField_T >
void initFillLevelsInBoundaries(const std::weak_ptr< StructuredBlockForest >& blockForestPtr,
                                const ConstBlockDataID& handlingID, const BlockDataID& fillFieldID,
                                real_t value = real_c(1))
{
   const auto blockForest = blockForestPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField           = blockIt->getData< ScalarField_T >(fillFieldID);
      const BoundaryHandling_T* const handling = blockIt->getData< const BoundaryHandling_T >(handlingID);

      // set fill level to "value" in every cell belonging to boundary
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(fillField, {
         if (handling->isBoundary(x, y, z)) { fillField->get(x, y, z) = value; }
      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }
}

/***********************************************************************************************************************
 * Clear and initialize flags in every cell according to the fill level.
 **********************************************************************************************************************/
template< typename BoundaryHandling_T, typename Stencil_T, typename FlagField_T, typename ScalarField_T >
void initFlagsFromFillLevels(const std::weak_ptr< StructuredBlockForest >& blockForestPtr,
                             const FlagInfo< FlagField_T >& flagInfo, const BlockDataID& handlingID,
                             const ConstBlockDataID& fillFieldID)
{
   const auto blockForest = blockForestPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(fillFieldID);
      BoundaryHandling_T* const handling   = blockIt->getData< BoundaryHandling_T >(handlingID);

      // clear all flags in the boundary handling
      handling->removeFlag(flagInfo.gasFlag);
      handling->removeFlag(flagInfo.liquidFlag);
      handling->removeFlag(flagInfo.interfaceFlag);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         // set flags only in non-boundary and non-obstacle cells
         if (!handling->isBoundary(fillFieldIt.x(), fillFieldIt.y(), fillFieldIt.z()))
         {
            if (*fillFieldIt <= real_c(0))
            {
               // set gas flag
               handling->forceFlag(flagInfo.gasFlag, fillFieldIt.x(), fillFieldIt.y(), fillFieldIt.z());
            }
            else
            {
               if (*fillFieldIt < real_c(1.0))
               {
                  // set interface flag
                  handling->forceFlag(flagInfo.interfaceFlag, fillFieldIt.x(), fillFieldIt.y(), fillFieldIt.z());
               }
               else
               {
                  // check if the cell is an interface cell due to direct neighboring gas cells
                  if (isInterfaceFromFillLevel< Stencil_T >(fillFieldIt))
                  {
                     // set interface flag
                     handling->forceFlag(flagInfo.interfaceFlag, fillFieldIt.x(), fillFieldIt.y(), fillFieldIt.z());
                  }
                  else
                  {
                     // set liquid flag
                     handling->forceFlag(flagInfo.liquidFlag, fillFieldIt.x(), fillFieldIt.y(), fillFieldIt.z());
                  }
               }
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }
}

/***********************************************************************************************************************
 * Initialize the hydrostatic pressure in the direction in which a force is acting in ALL cells (regardless of a cell's
 * flag). The velocity remains unchanged.
 *
 * The force vector must have only one component, i.e., the direction of the force can only be in x-, y- or z-axis.
 * The variable fluidHeight determines the height at which the density is equal to reference density (=1).
 **********************************************************************************************************************/
template< typename PdfField_T >
void initHydrostaticPressure(const std::weak_ptr< StructuredBlockForest >& blockForestPtr,
                             const BlockDataID& pdfFieldID, const Vector3< real_t >& force, real_t fluidHeight)
{
   // count number of non-zero components of the force vector
   uint_t numForceComponents = uint_c(0);
   if (!realIsEqual(force[0], real_c(0), real_c(1e-14))) { ++numForceComponents; }
   if (!realIsEqual(force[1], real_c(0), real_c(1e-14))) { ++numForceComponents; }
   if (!realIsEqual(force[2], real_c(0), real_c(1e-14))) { ++numForceComponents; }

   WALBERLA_CHECK_EQUAL(numForceComponents, uint_c(1),
                        "The current implementation of the hydrostatic pressure initialization does not support "
                        "forces that have none or multiple components, i. e., a force that points in a direction other "
                        "than the x-, y- or z-axis.");

   const auto blockForest = blockForestPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      PdfField_T* const pdfField = blockIt->getData< PdfField_T >(pdfFieldID);

      CellInterval local = pdfField->xyzSizeWithGhostLayer(); // block-, i.e., process-local cell interval

      for (auto cellIt = local.begin(); cellIt != local.end(); ++cellIt)
      {
         // transform the block-local coordinate to global coordinates
         Cell global;
         blockForest->transformBlockLocalToGlobalCell(global, *blockIt, *cellIt);

         // get the current global coordinate, i.e., height of the fluid in the relevant direction
         cell_idx_t coordinate = cell_idx_c(0);
         real_t forceComponent = real_c(0);
         if (!realIsEqual(force[0], real_c(0), real_c(1e-14)))
         {
            coordinate     = global.x();
            forceComponent = force[0];
         }
         else
         {
            if (!realIsEqual(force[1], real_c(0), real_c(1e-14)))
            {
               coordinate     = global.y();
               forceComponent = force[1];
            }
            else
            {
               if (!realIsEqual(force[2], real_c(0), real_c(1e-14)))
               {
                  coordinate     = global.z();
                  forceComponent = force[2];
               }
               else
               {
                  WALBERLA_ABORT(
                     "The current implementation of the hydrostatic pressure initialization does not support "
                     "forces that have none or multiple components, i. e., a force that points in a direction other "
                     "than the x-, y- or z-axis.")
               }
            }
         }

         // initialize the (hydrostatic) pressure, i.e., LBM density
         // Bernoulli: p = p0 + density * gravity * height
         // => LBM (density=1): rho = rho0 + gravity * height = 1 + 1/cs^2 * g * h = 1 + 3 * g * h
         // shift global cell by 0.5 since density is set for cell center
         const real_t rho =
            real_c(1) + real_c(3) * forceComponent * (real_c(coordinate) + real_c(0.5) - std::ceil(fluidHeight));

         const Vector3< real_t > velocity = pdfField->getVelocity(*cellIt);

         pdfField->setDensityAndVelocity(*cellIt, velocity, rho);
      }
   }
}

/***********************************************************************************************************************
 * Set density in non-liquid and non-interface cells to 1.
 **********************************************************************************************************************/
template< typename FlagField_T, typename PdfField_T >
void setDensityInNonFluidCellsToOne(const std::weak_ptr< StructuredBlockForest >& blockForestPtr,
                                    const FlagInfo< FlagField_T >& flagInfo, const ConstBlockDataID& flagFieldID,
                                    const BlockDataID& pdfFieldID)
{
   const auto blockForest = blockForestPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      PdfField_T* const pdfField         = blockIt->getData< PdfField_T >(pdfFieldID);
      const FlagField_T* const flagField = blockIt->getData< const FlagField_T >(flagFieldID);

      WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
         if (!flagInfo.isLiquid(*flagFieldIt) && !flagInfo.isInterface(*flagFieldIt))
         {
            // set density in gas cells to 1
            pdfField->setDensityAndVelocity(pdfFieldIt.cell(), Vector3< real_t >(real_c(0)), real_c(1));
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }
}
} // namespace free_surface
} // namespace walberla
