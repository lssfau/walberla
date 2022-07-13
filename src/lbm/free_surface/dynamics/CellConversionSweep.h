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
//! \file CellConversionSweep.h
//! \ingroup dynamics
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Convert cells to/from interface cells.
//
//======================================================================================================================

#pragma once

#include "field/FlagField.h"

#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/bubble_model/BubbleModel.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Convert cells to/from interface
 * - expects that a previous sweep has set the convertToLiquidFlag and convertToGasFlag flags
 * - notifies the bubble model that conversions occurred
 * - all cells that have been converted, have the "converted" flag set
 * - PDF field is required in order to initialize the velocity in new gas cells with information from surrounding cells
 *
 * A) Interface -> Liquid/Gas
 *     - always converts interface to liquid
 *     - converts interface to gas only if no newly created liquid cell is in neighborhood
 * B) Liquid/Gas -> Interface (to obtain a closed interface layer)
 *     - "old" liquid cells in the neighborhood of newly converted cells in A) are converted to interface
 *     - "old" gas cells in the neighborhood of newly converted cells in A) are converted to interface
 *     The term "old" refers to cells that were not converted themselves in the same time step.
 * C) GAS -> INTERFACE (due to inflow boundary condition)
 * D) LIQUID/GAS -> INTERFACE (due to wetting; only when using local triangulation for curvature computation)
 *
 * For gas cells that were converted to interface in B), the flag "convertedFromGasToInterface" is set to signal another
 * sweep (i.e. "PdfRefillingSweep.h") that the this cell's PDFs need to be reinitialized.
 *
 * For gas cells that were converted to interface in C), the cell's PDFs are reinitialized with equilibrium constructed
 * with the inflow velocity.
 *
 * More information can be found in the dissertation of N. Thuerey, 2007, section 4.3.
 * ********************************************************************************************************************/
template< typename LatticeModel_T, typename BoundaryHandling_T, typename ScalarField_T >
class CellConversionSweep
{
 public:
   using FlagField_T = typename BoundaryHandling_T::FlagField;
   using flag_t      = typename FlagField_T::flag_t;
   using PdfField_T  = lbm::PdfField< LatticeModel_T >;
   using Stencil_T   = typename LatticeModel_T::Stencil;

   CellConversionSweep(BlockDataID handlingID, BlockDataID pdfFieldID, const FlagInfo< FlagField_T >& flagInfo,
                       BubbleModelBase* bubbleModel)
      : handlingID_(handlingID), pdfFieldID_(pdfFieldID), bubbleModel_(bubbleModel), flagInfo_(flagInfo)
   {}

   void operator()(IBlock* const block)
   {
      BoundaryHandling_T* const handling = block->getData< BoundaryHandling_T >(handlingID_);
      PdfField_T* const pdfField         = block->getData< PdfField_T >(pdfFieldID_);

      FlagField_T* const flagField = handling->getFlagField();

      // A) INTERFACE -> LIQUID/GAS
      // convert interface cells that have filled/emptied to liquid/gas (cflagInfo_. dissertation of N. Thuerey, 2007,
      // section 4.3)
      // the conversion is also performed in the first ghost layer, since B requires an up-to-date first ghost layer;
      // explicitly avoid OpenMP, as bubble IDs are set here
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(flagField, uint_c(1), omp critical, {
         const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);

         // 1) convert interface cells to liquid
         if (isFlagSet(flagFieldPtr, flagInfo_.convertToLiquidFlag))
         {
            handling->removeFlag(flagInfo_.interfaceFlag, x, y, z);
            handling->setFlag(flagInfo_.liquidFlag, x, y, z);
            handling->setFlag(flagInfo_.convertedFlag, x, y, z);

            if (flagField->isInInnerPart(flagFieldPtr.cell()))
            {
               // register and detect splitting of bubbles
               bubbleModel_->reportInterfaceToLiquidConversion(block, flagFieldPtr.cell());
            }
         }

         // 2) convert interface cells to gas only if no newly converted liquid cell from 1) is in neighborhood to
         // ensure a closed interface
         if (isFlagSet(flagFieldPtr, flagInfo_.convertToGasFlag) &&
             !isFlagInNeighborhood< Stencil_T >(flagFieldPtr, flagInfo_.convertToLiquidFlag))
         {
            handling->removeFlag(flagInfo_.interfaceFlag, x, y, z);
            handling->setFlag(flagInfo_.gasFlag, x, y, z);
            handling->setFlag(flagInfo_.convertedFlag, x, y, z);
         }
      }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP

      // B) LIQUID/GAS -> INTERFACE
      // convert those liquid/gas cells to interface that are in the neighborhood of the newly created liquid/gas cells
      // from A; this maintains a closed interface layer;
      // explicitly avoid OpenMP, as bubble IDs are set here
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP(flagField, omp critical, {
         const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);

         // only consider "old" liquid cells, i.e., cells that have not been converted in this time step
         if (flagInfo_.isLiquid(flagFieldPtr) && !flagInfo_.hasConverted(flagFieldPtr))
         {
            const flag_t newGasFlagMask = flagInfo_.convertedFlag | flagInfo_.gasFlag; // flag newly converted gas cell

            // the state of ghost layer cells becomes relevant here
            for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d != LatticeModel_T::Stencil::end(); ++d)
               if (isMaskSet(flagFieldPtr.neighbor(*d), newGasFlagMask)) // newly converted gas cell is in neighborhood
               {
                  // convert the current cell to interface
                  handling->removeFlag(flagInfo_.liquidFlag, x, y, z);
                  handling->setFlag(flagInfo_.interfaceFlag, x, y, z);
                  handling->setFlag(flagInfo_.convertedFlag, x, y, z);

                  if (flagField->isInInnerPart(flagFieldPtr.cell()))
                  {
                     // register and detect merging of bubbles
                     bubbleModel_->reportLiquidToInterfaceConversion(block, flagFieldPtr.cell());
                  }

                  // current cell was already converted to interface, flags of other neighbors are not relevant
                  break;
               }
         }
         // only consider "old" gas cells, i.e., cells that have not been converted in this time step
         else
         {
            if (flagInfo_.isGas(flagFieldPtr) && !flagInfo_.hasConverted(flagFieldPtr))
            {
               const flag_t newLiquidFlagMask =
                  flagInfo_.convertedFlag | flagInfo_.liquidFlag; // flag of newly converted liquid cell

               // the state of ghost layer cells is relevant here
               for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d != LatticeModel_T::Stencil::end(); ++d)

                  // newly converted liquid cell is in neighborhood
                  if (isMaskSet(flagFieldPtr.neighbor(*d), newLiquidFlagMask))
                  {
                     // convert the current cell to interface
                     handling->removeFlag(flagInfo_.gasFlag, x, y, z);
                     handling->setFlag(flagInfo_.interfaceFlag, x, y, z);
                     handling->setFlag(flagInfo_.convertedFlag, x, y, z);
                     handling->setFlag(flagInfo_.convertFromGasToInterfaceFlag, x, y, z);

                     // current cell was already converted to interface, flags of other neighbors are not relevant
                     break;
                  }
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS_XYZ_OMP

      // C) GAS -> INTERFACE (due to inflow boundary condition)
      // convert gas cells to interface cells near inflow boundaries;
      // explicitly avoid OpenMP, such that cell conversions are performed sequentially
      convertedFromGasToInterfaceDueToInflow.clear();
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP(flagField, omp critical, {
         const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);

         if (flagInfo_.isConvertToInterfaceForInflow(flagFieldPtr) && !flagInfo_.hasConverted(flagFieldPtr))
         {
            // newly converted liquid cell is in neighborhood
            handling->removeFlag(flagInfo_.convertToInterfaceForInflowFlag, x, y, z);

            // convert the current cell to interface
            handling->removeFlag(flagInfo_.gasFlag, x, y, z);
            handling->setFlag(flagInfo_.interfaceFlag, x, y, z);
            handling->setFlag(flagInfo_.convertedFlag, x, y, z);
            convertedFromGasToInterfaceDueToInflow.insert(flagFieldPtr.cell());
         }
      }) // WALBERLA_FOR_ALL_CELLS_XYZ_OMP

      // D) LIQUID/GAS -> INTERFACE (due to wetting; only active when using local triangulation for curvature
      // computation)
      // convert liquid/gas to interface cells where the interface cell is required for a smooth
      // continuation of the wetting surface (see dissertation of S. Donath, 2011, section 6.3.5.3);
      // explicitly avoid OpenMP, as bubble IDs are set here
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP(flagField, omp critical, {
         const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);

         // only consider wetting and non-interface cells
         if (flagInfo_.isKeepInterfaceForWetting(flagFieldPtr) && !flagInfo_.isInterface(flagFieldPtr))
         {
            // convert liquid cell to interface
            if (isFlagSet(flagFieldPtr, flagInfo_.liquidFlag))
            {
               handling->removeFlag(flagInfo_.liquidFlag, x, y, z);
               handling->setFlag(flagInfo_.interfaceFlag, x, y, z);
               handling->setFlag(flagInfo_.convertedFlag, x, y, z);
               handling->removeFlag(flagInfo_.keepInterfaceForWettingFlag, x, y, z);
               if (flagField->isInInnerPart(flagFieldPtr.cell()))
               {
                  // register and detect merging of bubbles
                  bubbleModel_->reportLiquidToInterfaceConversion(block, flagFieldPtr.cell());
               }
            }
            else
            {
               // convert gas cell to interface
               if (isFlagSet(flagFieldPtr, flagInfo_.gasFlag))
               {
                  handling->removeFlag(flagInfo_.gasFlag, x, y, z);
                  handling->setFlag(flagInfo_.interfaceFlag, x, y, z);
                  handling->setFlag(flagInfo_.convertedFlag, x, y, z);
                  handling->removeFlag(flagInfo_.keepInterfaceForWettingFlag, x, y, z);
                  handling->setFlag(flagInfo_.convertFromGasToInterfaceFlag, x, y, z);
               }
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS_XYZ_OMP

      // initialize PDFs of interface cells that were created due to an inflow boundary; the PDFs are set to equilibrium
      // with density=1 and velocity of the inflow boundary
      initializeFromInflow(convertedFromGasToInterfaceDueToInflow, flagField, pdfField, handling);
   }

 protected:
   /********************************************************************************************************************
    * Initializes PDFs in cells that are converted to interface due to a neighboring inflow boundary.
    *
    * The PDFs of these cells are set to equilibrium values using density=1 and the average velocity of neighboring
    * inflow boundaries. An inflow cell is used for averaging only if the velocity actually flows towards the newly
    * created interface cell. In other words, the velocity direction is compared to the converted cell's direction with
    * respect to the inflow location.
    *
    * REMARK: The inflow boundary condition must implement function "getValue()" that returns the prescribed velocity
    * (see e.g. UBB).
    *******************************************************************************************************************/
   void initializeFromInflow(const std::set< Cell >& cells, FlagField_T* flagField, PdfField_T* pdfField,
                             BoundaryHandling_T* handling)
   {
      for (auto setIt = cells.begin(); setIt != cells.end(); ++setIt)
      {
         const Cell& cell = *setIt;

         Vector3< real_t > u(real_c(0.0));
         uint_t numNeighbors = uint_c(0);

         // get UBB inflow boundary
         auto ubbInflow = handling->template getBoundaryCondition<
            typename FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::UBB_Inflow_T >(
            handling->getBoundaryUID(
               FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::ubbInflowFlagID));

         for (auto i = LatticeModel_T::Stencil::beginNoCenter(); i != LatticeModel_T::Stencil::end(); ++i)
         {
            using namespace stencil;
            const Cell neighborCell(cell[0] + i.cx(), cell[1] + i.cy(), cell[2] + i.cz());

            const flag_t neighborFlag = flagField->get(neighborCell);

            // neighboring cell is inflow
            if (isPartOfMaskSet(neighborFlag, flagInfo_.inflowFlagMask))
            {
               // get direction towards cell containing inflow boundary
               const Vector3< int > dir = Vector3< int >(-i.cx(), -i.cy(), -i.cz());

               // get velocity from UBB boundary
               const Vector3< real_t > inflowVel =
                  ubbInflow.getValue(cell[0] + i.cx(), cell[1] + i.cy(), cell[2] + i.cz());

               // skip directions in which the corresponding velocity component is zero
               if (realIsEqual(inflowVel[0], real_c(0), real_c(1e-14)) && dir[0] != 0) { continue; }
               if (realIsEqual(inflowVel[1], real_c(0), real_c(1e-14)) && dir[1] != 0) { continue; }
               if (realIsEqual(inflowVel[2], real_c(0), real_c(1e-14)) && dir[2] != 0) { continue; }

               // skip directions in which the corresponding velocity component is in opposite direction
               if (inflowVel[0] > real_c(0) && dir[0] < 0) { continue; }
               if (inflowVel[1] > real_c(0) && dir[1] < 0) { continue; }
               if (inflowVel[2] > real_c(0) && dir[2] < 0) { continue; }

               // use inflow velocity to get average velocity
               u += inflowVel;
               numNeighbors++;
            }
         }
         if (numNeighbors > uint_c(0)) { u /= real_c(numNeighbors); } // else: velocity is zero

         pdfField->setDensityAndVelocity(cell, u, real_c(1)); // set density=1
      }
   }

   BlockDataID handlingID_;
   BlockDataID pdfFieldID_;
   BubbleModelBase* bubbleModel_;

   FlagInfo< FlagField_T > flagInfo_;

   std::set< Cell > convertedFromGasToInterfaceDueToInflow;
}; // class CellConversionSweep

} // namespace free_surface
} // namespace walberla
