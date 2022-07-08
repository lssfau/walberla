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
//! \file FindInterfaceCellConversion.h
//! \ingroup dynamics
//! \author Martin Bauer
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Find and mark interface cells for conversion to gas/liquid.
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/timing/TimingPool.h"

#include "lbm/free_surface/FlagInfo.h"

#include <type_traits>

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Finds interface cell conversions to gas/liquid and sets corresponding conversion flag that marks this cell for
 * conversion.
 *
 * This version uses the cached OR'ed flag neighborhood given in neighborFlags. See free_surface::getOredNeighborhood.
 *
 * This function decides by looking at the fill level which interface cells should be converted to gas or liquid cells.
 * It does not change the state directly, it only sets "conversion suggestions" by setting flags called
 * 'convertToGasFlag' and 'convertToLiquidFlag'.
 *
 * An interface is first classified as one of the following types (same classification as in free_surface::advectMass)
 * - _to gas_   : if cell has no liquid cell in neighborhood, this cell should become gas
 * - _to liquid_: if cell has no gas cell in neighborhood, this cell should become liquid
 * - _pure interface_: if not '_to gas_' and not '_to liquid_'
 *
 * This classification up to now depends only on the neighborhood flags, not on the fill level.
 *
 *  _Pure interface_ cells are marked for to-gas-conversion if their fill level is lower than
 *  "0 - cellConversionThreshold" and marked for to-liquid-conversion if the fill level is higher than
 *  "1 + cellConversionThreshold". The threshold is introduced to prevent oscillating cell conversions.
 *  The value of the offset is chosen heuristically (see dissertation of N. Thuerey, 2007: 1e-3, dissertation of
 *  T. Pohl, 2008: 1e-2).
 *
 * Additionally, interface cells without fluid neighbors are marked for conversion to gas if an:
 * - interface cell is almost empty (fill < cellConversionForceThreshold) AND
 * - interface cell has no neighboring interface cells
 * OR
 * - interface cell has neighboring cell with outflow boundary condition
 *
 * Similarly, interface cells without gas neighbors are marked for conversion to liquid if:
 * - interface cell is almost full (fill > 1.0-cellConversionForceThreshold) AND
 * - interface cell has no neighboring interface cells
 *
 * The value of cellConversionForceThreshold is chosen heuristically: 1e-1 (see dissertation of N. Thuerey, 2007,
 * section 4.4).
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename BoundaryHandling_T, typename ScalarField_T, typename FlagField_T,
          typename ScalarIt_T, typename FlagIt_T >
void findInterfaceCellConversion(const BoundaryHandling_T& handling, const ScalarIt_T& fillFieldIt,
                                 FlagIt_T& flagFieldIt, const typename FlagField_T::flag_t& neighborFlags,
                                 const FlagInfo< FlagField_T >& flagInfo, real_t cellConversionThreshold,
                                 real_t cellConversionForceThreshold)
{
   static_assert(std::is_same< typename FlagField_T::value_type, typename FlagIt_T::value_type >::value,
                 "The given flagFieldIt does not seem to be an iterator of the provided FlagField_T.");

   static_assert(std::is_floating_point< typename ScalarIt_T::value_type >::value,
                 "The given fillFieldIt has to be a floating point value.");

   cellConversionThreshold      = std::abs(cellConversionThreshold);
   cellConversionForceThreshold = std::abs(cellConversionForceThreshold);

   WALBERLA_ASSERT_LESS(cellConversionThreshold, cellConversionForceThreshold);

   // in the neighborhood of inflow boundaries, convert gas cells to interface cells depending on the direction of the
   // prescribed inflow velocity
   if (field::isFlagSet(neighborFlags, flagInfo.inflowFlagMask) && field::isFlagSet(flagFieldIt, flagInfo.gasFlag))
   {
      // get UBB inflow boundary
      auto ubbInflow = handling->template getBoundaryCondition<
         typename FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::UBB_Inflow_T >(
         handling->getBoundaryUID(
            FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::ubbInflowFlagID));

      for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d != LatticeModel_T::Stencil::end(); ++d)
      {
         if (field::isMaskSet(flagFieldIt.neighbor(*d), flagInfo.inflowFlagMask))
         {
            // get direction of cell containing inflow boundary
            const Vector3< int > dir = Vector3< int >(-d.cx(), -d.cy(), -d.cz());

            // get velocity from UBB inflow boundary
            const Vector3< real_t > inflowVel =
               ubbInflow.getValue(flagFieldIt.x() + d.cx(), flagFieldIt.y() + d.cy(), flagFieldIt.z() + d.cz());

            // skip directions in which the corresponding velocity component is zero
            if (realIsEqual(inflowVel[0], real_c(0), real_c(1e-14)) && dir[0] != 0) { continue; }
            if (realIsEqual(inflowVel[1], real_c(0), real_c(1e-14)) && dir[1] != 0) { continue; }
            if (realIsEqual(inflowVel[2], real_c(0), real_c(1e-14)) && dir[2] != 0) { continue; }

            // skip directions in which the corresponding velocity component is in opposite direction
            if (inflowVel[0] > real_c(0) && dir[0] < 0) { continue; }
            if (inflowVel[1] > real_c(0) && dir[1] < 0) { continue; }
            if (inflowVel[2] > real_c(0) && dir[2] < 0) { continue; }

            // set conversion flag to remaining cells
            field::addFlag(flagFieldIt, flagInfo.convertToInterfaceForInflowFlag);
         }
      }
      return;
   }

   // only interface cells are converted directly (except for cells near inflow boundaries, see above)
   if (!field::isFlagSet(flagFieldIt, flagInfo.interfaceFlag)) { return; }

   // interface cell is empty and should be converted to gas
   if (*fillFieldIt < -cellConversionThreshold)
   {
      if (field::isFlagSet(flagFieldIt, flagInfo.keepInterfaceForWettingFlag))
      {
         field::removeFlag(flagFieldIt, flagInfo.keepInterfaceForWettingFlag);
      }

      field::addFlag(flagFieldIt, flagInfo.convertToGasFlag);
      return;
   }

   // interface cell is full and should be converted to liquid
   if (*fillFieldIt > real_c(1.0) + cellConversionThreshold)
   {
      if (field::isFlagSet(flagFieldIt, flagInfo.keepInterfaceForWettingFlag))
      {
         field::removeFlag(flagFieldIt, flagInfo.keepInterfaceForWettingFlag);
      }

      field::addFlag(flagFieldIt, flagInfo.convertToLiquidFlag);
      return;
   }

   // interface cell has no liquid neighbor and should be converted to gas (see dissertation of N. Thuerey, 2007
   // section 4.4)
   if (!field::isFlagSet(neighborFlags, flagInfo.liquidFlag) &&
       !field::isFlagSet(flagFieldIt, flagInfo.keepInterfaceForWettingFlag) &&
       !field::isFlagSet(neighborFlags, flagInfo.inflowFlagMask))
   {
      // interface cell is almost empty
      if (*fillFieldIt < cellConversionForceThreshold && field::isFlagSet(neighborFlags, flagInfo.interfaceFlag))
      {
         // mass is not necessarily lost as it can be distributed to a neighboring interface cell
         field::addFlag(flagFieldIt, flagInfo.convertToGasFlag);
         return;
      }

      // interface cell has no other interface neighbors; conversion might lead to loss in mass (depending on the excess
      // mass distribution model)
      if (!field::isFlagSet(neighborFlags, flagInfo.interfaceFlag))
      {
         field::addFlag(flagFieldIt, flagInfo.convertToGasFlag);
         return;
      }
   }

   // interface cell has no gas neighbor and should be converted to liquid (see dissertation of N. Thuerey, 2007
   // section 4.4)
   if (!field::isFlagSet(neighborFlags, flagInfo.gasFlag) &&
       !field::isFlagSet(flagFieldIt, flagInfo.keepInterfaceForWettingFlag))
   {
      // interface cell is almost full
      if (*fillFieldIt > real_c(1.0) - cellConversionForceThreshold &&
          field::isFlagSet(neighborFlags, flagInfo.interfaceFlag))
      {
         // mass is not necessarily gained as it can be taken from a neighboring interface cell
         field::addFlag(flagFieldIt, flagInfo.convertToLiquidFlag);
         return;
      }

      // interface cell has no other interface neighbors; conversion might lead to gain in mass (depending on the excess
      // mass distribution model)
      if (!field::isFlagSet(neighborFlags, flagInfo.interfaceFlag))
      {
         field::addFlag(flagFieldIt, flagInfo.convertToLiquidFlag);
         return;
      }
   }
}

/***********************************************************************************************************************
 * Triggers findInterfaceCellConversion() for each cell of the given fields and recomputes the OR'ed flag neighborhood
 * info.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename BoundaryHandling_T, typename ScalarField_T, typename FlagField_T >
void findInterfaceCellConversions(const BoundaryHandling_T& handling, const ScalarField_T* fillField,
                                  FlagField_T* flagField, const FlagInfo< FlagField_T >& flagInfo,
                                  real_t cellConversionThreshold, real_t cellConversionForceThreshold)
{
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(flagField, uint_c(1), omp critical, {
      const typename FlagField_T::Ptr flagFieldPtr(*flagField, x, y, z);
      const typename ScalarField_T::ConstPtr fillFieldPtr(*fillField, x, y, z);

      if (field::isFlagSet(flagFieldPtr, flagInfo.interfaceFlag))
      {
         const typename FlagField_T::value_type neighborFlags =
            field::getOredNeighborhood< typename LatticeModel_T::Stencil >(flagFieldPtr);

         (findInterfaceCellConversion< LatticeModel_T, BoundaryHandling_T, ScalarField_T,
                                       FlagField_T >) (handling, fillFieldPtr, flagFieldPtr, neighborFlags, flagInfo,
                                                       cellConversionThreshold, cellConversionForceThreshold);
      }
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}

/***********************************************************************************************************************
 * Triggers findInterfaceCellConversion() for each cell of the given fields using the cached OR'ed flag neighborhood
 * given in neighField.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename BoundaryHandling_T, typename ScalarField_T, typename FlagField_T,
          typename NeighField_T >
void findInterfaceCellConversions(const BoundaryHandling_T& handling, const ScalarField_T* fillField,
                                  FlagField_T* flagField, const NeighField_T* neighField,
                                  const FlagInfo< FlagField_T >& flagInfo, real_t cellConversionThreshold,
                                  real_t cellConversionForceThreshold)
{
   WALBERLA_ASSERT_EQUAL_2(flagField->xyzSize(), fillField->xyzSize());
   WALBERLA_ASSERT_EQUAL_2(flagField->xyzSize(), neighField->xyzSize());

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(flagField, uint_c(1), omp critical, {
      const typename FlagField_T::Ptr flagFieldPtr(*flagField, x, y, z);
      const typename ScalarField_T::ConstPtr fillFieldPtr(*fillField, x, y, z);

      (findInterfaceCellConversion< LatticeModel_T, BoundaryHandling_T, ScalarField_T,
                                    FlagField_T >) (handling, fillFieldPtr, flagFieldPtr, neighField->get(x, y, z),
                                                    flagInfo, cellConversionThreshold, cellConversionForceThreshold);
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP
}

} // namespace free_surface
} // namespace walberla
