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
//! \file ExcessMassDistributionSweep.impl.h
//! \ingroup dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Distribute excess mass, i.e., mass that is undistributed after conversions from interface to liquid or gas.
//
//======================================================================================================================

#pragma once

#include "core/logging/Logging.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"
#include "field/GhostLayerField.h"

#include "lbm/field/PdfField.h"
#include "lbm/free_surface/FlagInfo.h"

#include "ExcessMassDistributionSweep.h"

namespace walberla
{
namespace free_surface
{
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExcessMassDistributionSweepInterfaceEvenly< LatticeModel_T, FlagField_T, ScalarField_T,
                                                 VectorField_T >::operator()(IBlock* const block)
{
   using Base_T = ExcessMassDistributionSweepBase_T;

   const FlagField_T* const flagField = block->getData< const FlagField_T >(Base_T::flagFieldID_);
   ScalarField_T* const fillField     = block->getData< ScalarField_T >(Base_T::fillFieldID_);
   const lbm::PdfField< LatticeModel_T >* const pdfField =
      block->getData< const lbm::PdfField< LatticeModel_T > >(Base_T::pdfFieldID_);

   // disable OpenMP to avoid mass distribution to neighboring cells before they have distributed their excess mass
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(fillField, uint_c(1), omp critical, {
      const Cell cell(x, y, z);

      if (flagField->isFlagSet(cell, Base_T::flagInfo_.convertedFlag))
      {
         // identify cells that were converted to gas/liquid in this time step
         const bool newGas = flagField->isMaskSet(cell, Base_T::flagInfo_.convertedFlag | Base_T::flagInfo_.gasFlag);
         const bool newLiquid =
            flagField->isMaskSet(cell, Base_T::flagInfo_.convertedFlag | Base_T::flagInfo_.liquidFlag);

         if (newGas || newLiquid)
         {
            // a cell can not be converted to both gas and liquid
            WALBERLA_ASSERT(!(newGas && newLiquid));

            // calculate excess fill level
            const real_t excessFill = newGas ? fillField->get(cell) : (fillField->get(cell) - real_c(1.0));

            distributeMassEvenly(fillField, flagField, pdfField, cell, excessFill);

            if (newGas) { fillField->get(cell) = real_c(0.0); }
            else { fillField->get(cell) = real_c(1.0); }
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExcessMassDistributionSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   getNumberOfLiquidAndInterfaceNeighbors(const FlagField_T* flagField, const Cell& cell, uint_t& liquidNeighbors,
                                          uint_t& interfaceNeighbors, uint_t& newInterfaceNeighbors)
{
   newInterfaceNeighbors = uint_c(0);
   interfaceNeighbors    = uint_c(0);
   liquidNeighbors       = uint_c(0);

   for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d != LatticeModel_T::Stencil::end(); ++d)
   {
      const Cell neighborCell = Cell(cell.x() + d.cx(), cell.y() + d.cy(), cell.z() + d.cz());
      auto neighborFlags      = flagField->get(neighborCell);

      if (isFlagSet(neighborFlags, flagInfo_.interfaceFlag))
      {
         ++interfaceNeighbors;
         if (isFlagSet(neighborFlags, flagInfo_.convertedFlag)) { ++newInterfaceNeighbors; }
      }
      else
      {
         if (isFlagSet(neighborFlags, flagInfo_.liquidFlag)) { ++liquidNeighbors; }
      }
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExcessMassDistributionSweepBase< LatticeModel_T, FlagField_T, ScalarField_T,
                                      VectorField_T >::getNumberOfInterfaceNeighbors(const FlagField_T* flagField,
                                                                                     const Cell& cell,
                                                                                     uint_t& newInterfaceNeighbors,
                                                                                     uint_t& interfaceNeighbors)
{
   interfaceNeighbors    = uint_c(0);
   newInterfaceNeighbors = uint_c(0);

   for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d != LatticeModel_T::Stencil::end(); ++d)
   {
      const Cell neighborCell = Cell(cell.x() + d.cx(), cell.y() + d.cy(), cell.z() + d.cz());
      auto neighborFlags      = flagField->get(neighborCell);

      if (isFlagSet(neighborFlags, flagInfo_.interfaceFlag))
      {
         ++interfaceNeighbors;
         if (isFlagSet(neighborFlags, flagInfo_.convertedFlag)) { ++newInterfaceNeighbors; }
      }
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
template< typename PdfField_T >
void ExcessMassDistributionSweepInterfaceEvenly< LatticeModel_T, FlagField_T, ScalarField_T,
                                                 VectorField_T >::distributeMassEvenly(ScalarField_T* fillField,
                                                                                       const FlagField_T* flagField,
                                                                                       const PdfField_T* pdfField,
                                                                                       const Cell& cell,
                                                                                       real_t excessFill)
{
   using Base_T = ExcessMassDistributionSweepBase_T;

   bool useEvenlyAll = Base_T::excessMassDistributionModel_.getModelType() ==
                       ExcessMassDistributionModel::ExcessMassModel::EvenlyAllInterface;
   bool useEvenlyNew = Base_T::excessMassDistributionModel_.getModelType() ==
                       ExcessMassDistributionModel::ExcessMassModel::EvenlyNewInterface;
   bool useEvenlyOld = Base_T::excessMassDistributionModel_.getModelType() ==
                       ExcessMassDistributionModel::ExcessMassModel::EvenlyOldInterface;

   // get number of interface neighbors
   uint_t newInterfaceNeighbors = uint_c(0);
   uint_t interfaceNeighbors    = uint_c(0);
   Base_T::getNumberOfInterfaceNeighbors(flagField, cell, newInterfaceNeighbors, interfaceNeighbors);
   const uint_t oldInterfaceNeighbors = interfaceNeighbors - newInterfaceNeighbors;

   if (interfaceNeighbors == uint_c(0))
   {
      WALBERLA_LOG_WARNING(
         "No interface cell is in the neighborhood to distribute excess mass to. Mass is lost/gained.");
      return;
   }

   // get density of the current cell
   const real_t density = pdfField->getDensity(cell);

   // compute mass to be distributed to neighboring cells
   real_t deltaMass = real_c(0);
   if ((useEvenlyOld && oldInterfaceNeighbors > uint_c(0)) || newInterfaceNeighbors == uint_c(0))
   {
      useEvenlyOld = true;
      useEvenlyAll = false;
      useEvenlyNew = false;

      deltaMass = excessFill / real_c(oldInterfaceNeighbors) * density;
   }
   else
   {
      if (useEvenlyNew || oldInterfaceNeighbors == uint_c(0))
      {
         useEvenlyOld = false;
         useEvenlyAll = false;
         useEvenlyNew = true;

         deltaMass = excessFill / real_c(newInterfaceNeighbors) * density;
      }
      else
      {
         useEvenlyOld = false;
         useEvenlyAll = true;
         useEvenlyNew = false;

         deltaMass = excessFill / real_c(interfaceNeighbors) * density;
      }
   }

   // distribute the excess mass
   for (auto pushDir = LatticeModel_T::Stencil::beginNoCenter(); pushDir != LatticeModel_T::Stencil::end(); ++pushDir)
   {
      const Cell neighborCell = Cell(cell.x() + pushDir.cx(), cell.y() + pushDir.cy(), cell.z() + pushDir.cz());

      // do not push mass in the direction of the second ghost layer
      // - inner domain: 0 to *Size()-1
      // - inner domain incl. first ghost layer: -1 to *Size()
      if (neighborCell.x() < cell_idx_c(-1) || neighborCell.y() < cell_idx_c(-1) || neighborCell.z() < cell_idx_c(-1) ||
          neighborCell.x() > cell_idx_c(fillField->xSize()) || neighborCell.y() > cell_idx_c(fillField->ySize()) ||
          neighborCell.z() > cell_idx_c(fillField->zSize()))
      {
         continue;
      }

      // only push mass to neighboring interface cells
      if (flagField->isFlagSet(neighborCell, Base_T::flagInfo_.interfaceFlag))
      {
         // get density of neighboring interface cell
         const real_t neighborDensity = pdfField->getDensity(neighborCell);

         if (flagField->isFlagSet(neighborCell, Base_T::flagInfo_.convertedFlag) && (useEvenlyAll || useEvenlyNew))
         {
            // push mass to newly converted interface cell
            fillField->get(neighborCell) += deltaMass / neighborDensity;
         }
         else
         {
            if (!flagField->isFlagSet(neighborCell, Base_T::flagInfo_.convertedFlag) && (useEvenlyOld || useEvenlyAll))
            {
               // push mass to old, i.e., non-newly converted interface cells
               fillField->getNeighbor(cell, *pushDir) += deltaMass / neighborDensity;
            }
         }
      }
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExcessMassDistributionSweepInterfaceWeighted< LatticeModel_T, FlagField_T, ScalarField_T,
                                                   VectorField_T >::operator()(IBlock* const block)
{
   using Base_T = ExcessMassDistributionSweepBase_T;

   const FlagField_T* const flagField = block->getData< const FlagField_T >(Base_T::flagFieldID_);
   ScalarField_T* const fillField     = block->getData< ScalarField_T >(Base_T::fillFieldID_);
   const lbm::PdfField< LatticeModel_T >* const pdfField =
      block->getData< const lbm::PdfField< LatticeModel_T > >(Base_T::pdfFieldID_);
   const VectorField_T* const normalField = block->getData< const VectorField_T >(normalFieldID_);

   // disable OpenMP to avoid mass distribution to neighboring cells before they have distributed their excess mass
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(fillField, uint_c(1), omp critical, {
      const Cell cell(x, y, z);

      if (flagField->isFlagSet(cell, Base_T::flagInfo_.convertedFlag))
      {
         // identify cells that were converted to gas/liquid in this time step
         const bool newGas = flagField->isMaskSet(cell, Base_T::flagInfo_.convertedFlag | Base_T::flagInfo_.gasFlag);
         const bool newLiquid =
            flagField->isMaskSet(cell, Base_T::flagInfo_.convertedFlag | Base_T::flagInfo_.liquidFlag);

         if (newGas || newLiquid)
         {
            // a cell can not be converted to both gas and liquid
            WALBERLA_ASSERT(!(newGas && newLiquid));

            // calculate excess fill level
            const real_t excessFill = newGas ? fillField->get(cell) : (fillField->get(cell) - real_c(1.0));

            distributeMassWeighted(fillField, flagField, pdfField, normalField, cell, newLiquid, excessFill);

            if (newGas) { fillField->get(cell) = real_c(0.0); }
            else { fillField->get(cell) = real_c(1.0); }
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
template< typename PdfField_T >
void ExcessMassDistributionSweepInterfaceWeighted< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   distributeMassWeighted(ScalarField_T* fillField, const FlagField_T* flagField, const PdfField_T* pdfField,
                          const VectorField_T* normalField, const Cell& cell, bool isNewLiquid, real_t excessFill)
{
   using Base_T = ExcessMassDistributionSweepBase_T;

   bool useWeightedAll = Base_T::excessMassDistributionModel_.getModelType() ==
                         ExcessMassDistributionModel::ExcessMassModel::WeightedAllInterface;
   bool useWeightedNew = Base_T::excessMassDistributionModel_.getModelType() ==
                         ExcessMassDistributionModel::ExcessMassModel::WeightedNewInterface;
   bool useWeightedOld = Base_T::excessMassDistributionModel_.getModelType() ==
                         ExcessMassDistributionModel::ExcessMassModel::WeightedOldInterface;

   // get number of interface neighbors
   uint_t newInterfaceNeighbors = uint_c(0);
   uint_t interfaceNeighbors    = uint_c(0);
   Base_T::getNumberOfInterfaceNeighbors(flagField, cell, newInterfaceNeighbors, interfaceNeighbors);
   const uint_t oldInterfaceNeighbors = interfaceNeighbors - newInterfaceNeighbors;

   if (interfaceNeighbors == uint_c(0))
   {
      WALBERLA_LOG_WARNING(
         "No interface cell is in the neighborhood to distribute excess mass to. Mass is lost/gained.");
      return;
   }

   // check applicability of the chosen model
   if ((useWeightedOld && oldInterfaceNeighbors > uint_c(0)) || newInterfaceNeighbors == uint_c(0))
   {
      useWeightedOld = true;
      useWeightedAll = false;
      useWeightedNew = false;
   }
   else
   {
      if (useWeightedNew || oldInterfaceNeighbors == uint_c(0))
      {
         useWeightedOld = false;
         useWeightedAll = false;
         useWeightedNew = true;
      }
      else
      {
         useWeightedOld = false;
         useWeightedAll = true;
         useWeightedNew = false;
      }
   }

   // get normal-direction-based weights of the excess mass
   std::vector< real_t > weights(LatticeModel_T::Stencil::Size, real_c(0));
   getExcessMassWeights(flagField, normalField, cell, isNewLiquid, useWeightedOld, useWeightedAll, useWeightedNew,
                        weights);

   // get the sum of all weights
   real_t weightSum = real_c(0);
   for (const auto& w : weights)
   {
      weightSum += w;
   }

   // if there are either no old or no newly converted interface cells in normal direction, distribute mass to whatever
   // interface cell is available in normal direction
   if (realIsEqual(weightSum, real_c(0), real_c(1e-14)) && (useWeightedOld || useWeightedNew))
   {
      useWeightedOld = false;
      useWeightedAll = true;
      useWeightedNew = false;

      // recompute mass weights since other type of interface cells are now considered also
      getExcessMassWeights(flagField, normalField, cell, isNewLiquid, useWeightedOld, useWeightedAll, useWeightedNew,
                           weights);

      // update sum of all weights
      for (const auto& w : weights)
      {
         weightSum += w;
      }

      // if no interface cell is available in normal direction, distribute mass evenly to all neighboring interface
      // cells
      if (realIsEqual(weightSum, real_c(0), real_c(1e-14)))
      {
         WALBERLA_LOG_WARNING_ON_ROOT(
            "Excess mass can not be distributed with a weighted approach since no interface cell is available in "
            "normal direction. Distributing excess mass evenly among all surrounding interface cells.");

         // manually set weights to 1 to get equal weight in any direction
         for (auto& w : weights)
         {
            w = real_c(1);
         }

         // weight sum is now the number of neighboring interface cells
         weightSum = real_c(interfaceNeighbors);
      }
   }

   WALBERLA_ASSERT_GREATER(
      weightSum, real_c(0),
      "Sum of all weights is zero in ExcessMassDistribution. This means that no neighboring interface cell is "
      "available for distributing the excess mass to. This error should have been caught earlier.");

   const real_t excessMass = excessFill * pdfField->getDensity(cell);

   // distribute the excess mass
   for (auto pushDir = LatticeModel_T::Stencil::beginNoCenter(); pushDir != LatticeModel_T::Stencil::end(); ++pushDir)
   {
      const Cell neighborCell = Cell(cell.x() + pushDir.cx(), cell.y() + pushDir.cy(), cell.z() + pushDir.cz());

      // do not push mass in the direction of the second ghost layer
      // - inner domain: 0 to *Size()-1
      // - inner domain incl. first ghost layer: -1 to *Size()
      if (neighborCell.x() < cell_idx_c(-1) || neighborCell.y() < cell_idx_c(-1) || neighborCell.z() < cell_idx_c(-1) ||
          neighborCell.x() > cell_idx_c(fillField->xSize()) || neighborCell.y() > cell_idx_c(fillField->ySize()) ||
          neighborCell.z() > cell_idx_c(fillField->zSize()))
      {
         continue;
      }

      // only push mass to neighboring interface cells
      if (flagField->isFlagSet(neighborCell, Base_T::flagInfo_.interfaceFlag))
      {
         // get density of neighboring interface cell
         const real_t neighborDensity = pdfField->getDensity(neighborCell);

         if (flagField->isFlagSet(neighborCell, Base_T::flagInfo_.convertedFlag) && (useWeightedAll || useWeightedNew))
         {
            // push mass to newly converted interface cell
            const real_t deltaMass = excessMass * weights[pushDir.toIdx()] / weightSum;
            fillField->get(neighborCell) += deltaMass / neighborDensity;
         }
         else
         {
            if (!flagField->isFlagSet(neighborCell, Base_T::flagInfo_.convertedFlag) &&
                (useWeightedOld || useWeightedAll))
            {
               // push mass to old, i.e., non-newly converted interface cells
               const real_t deltaMass = excessMass * weights[pushDir.toIdx()] / weightSum;
               fillField->getNeighbor(cell, *pushDir) += deltaMass / neighborDensity;
            }
         }
      }
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExcessMassDistributionSweepInterfaceWeighted< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   getExcessMassWeights(const FlagField_T* flagField, const VectorField_T* normalField, const Cell& cell,
                        bool isNewLiquid, bool useWeightedOld, bool useWeightedAll, bool useWeightedNew,
                        std::vector< real_t >& weights)
{
   using Base_T = ExcessMassDistributionSweepBase_T;

   // iterate all neighboring cells
   for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d != LatticeModel_T::Stencil::end(); ++d)
   {
      const Cell neighborCell = Cell(cell.x() + d.cx(), cell.y() + d.cy(), cell.z() + d.cz());
      auto neighborFlags      = flagField->get(neighborCell);

      if (isFlagSet(neighborFlags, Base_T::flagInfo_.interfaceFlag))
      {
         // compute dot product of normal direction and lattice direction to neighboring cell
         const real_t n_dot_ci =
            normalField->get(cell) * Vector3< real_t >(real_c(d.cx()), real_c(d.cy()), real_c(d.cz()));

         if (useWeightedAll || (useWeightedOld && !isFlagSet(neighborFlags, Base_T::flagInfo_.convertedFlag)))
         {
            computeWeightWithNormal(n_dot_ci, isNewLiquid, d, weights);
         }
         else
         {
            if (useWeightedNew && isFlagSet(neighborFlags, Base_T::flagInfo_.convertedFlag))
            {
               computeWeightWithNormal(n_dot_ci, isNewLiquid, d, weights);
            }
            else { weights[d.toIdx()] = real_c(0); }
         }
      }
      else
      {
         // no interface cell in this direction, weight is zero
         weights[d.toIdx()] = real_c(0);
      }
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExcessMassDistributionSweepInterfaceWeighted< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   computeWeightWithNormal(real_t n_dot_ci, bool isNewLiquid, typename LatticeModel_T::Stencil::iterator dir,
                           std::vector< real_t >& weights)
{
   // dissertation of N. Thuerey, 2007, equation (4.9)
   if (isNewLiquid)
   {
      if (n_dot_ci > real_c(0)) { weights[dir.toIdx()] = n_dot_ci; }
      else { weights[dir.toIdx()] = real_c(0); }
   }
   else // cell was converted from interface to gas
   {
      if (n_dot_ci < real_c(0)) { weights[dir.toIdx()] = -n_dot_ci; }
      else { weights[dir.toIdx()] = real_c(0); }
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExcessMassDistributionSweepInterfaceAndLiquid< LatticeModel_T, FlagField_T, ScalarField_T,
                                                    VectorField_T >::operator()(IBlock* const block)
{
   using Base_T = ExcessMassDistributionSweepBase_T;

   const FlagField_T* const flagField = block->getData< const FlagField_T >(Base_T::flagFieldID_);
   ScalarField_T* const fillField     = block->getData< ScalarField_T >(Base_T::fillFieldID_);
   const lbm::PdfField< LatticeModel_T >* const pdfField =
      block->getData< const lbm::PdfField< LatticeModel_T > >(Base_T::pdfFieldID_);

   ScalarField_T* const srcExcessMassField = block->getData< ScalarField_T >(excessMassFieldID_);
   ScalarField_T* const dstExcessMassField = excessMassFieldClone_.get(block);

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(dstExcessMassField, uint_c(1), {
      dstExcessMassField->get(x, y, z) = real_c(0);
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   // disable OpenMP to avoid mass distribution to neighboring cells before they have distributed their excess mass
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(fillField, uint_c(1), omp critical, {
      const Cell cell(x, y, z);

      if (flagField->isFlagSet(cell, Base_T::flagInfo_.convertedFlag))
      {
         // identify cells that were converted to gas/liquid in this time step
         const bool newGas = flagField->isMaskSet(cell, Base_T::flagInfo_.convertedFlag | Base_T::flagInfo_.gasFlag);
         const bool newLiquid =
            flagField->isMaskSet(cell, Base_T::flagInfo_.convertedFlag | Base_T::flagInfo_.liquidFlag);

         if (newGas || newLiquid)
         {
            // a cell can not be converted to both gas and liquid
            WALBERLA_ASSERT(!(newGas && newLiquid));

            // calculate excess fill level
            const real_t excessFill = newGas ? fillField->get(cell) : (fillField->get(cell) - real_c(1.0));

            // store excess mass such that it can be distributed below (no += here because cell was an interface cell
            // that can not have an excess mass stored in the field; any excess mass is added to the interface cell's
            // fill level)
            srcExcessMassField->get(cell) = excessFill * pdfField->getDensity(cell);

            if (newGas) { fillField->get(cell) = real_c(0.0); }
            else { fillField->get(cell) = real_c(1.0); }
         }
      }

      if (!realIsEqual(srcExcessMassField->get(cell), real_c(0), real_c(1e-14)))
      {
         distributeMassInterfaceAndLiquid(fillField, dstExcessMassField, flagField, pdfField, cell,
                                          srcExcessMassField->get(cell));
      }
   }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   srcExcessMassField->swapDataPointers(dstExcessMassField);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
template< typename PdfField_T >
void ExcessMassDistributionSweepInterfaceAndLiquid< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   distributeMassInterfaceAndLiquid(ScalarField_T* fillField, ScalarField_T* dstExcessMassField,
                                    const FlagField_T* flagField, const PdfField_T* pdfField, const Cell& cell,
                                    real_t excessMass)
{
   using Base_T = ExcessMassDistributionSweepBase_T;

   // get number of liquid and interface neighbors
   uint_t liquidNeighbors       = uint_c(0);
   uint_t interfaceNeighbors    = uint_c(0);
   uint_t newInterfaceNeighbors = uint_c(0);
   Base_T::getNumberOfLiquidAndInterfaceNeighbors(flagField, cell, liquidNeighbors, interfaceNeighbors,
                                                  newInterfaceNeighbors);
   const uint_t liquidAndInterfaceNeighbors = liquidNeighbors + interfaceNeighbors;

   if (liquidAndInterfaceNeighbors == uint_c(0))
   {
      WALBERLA_LOG_WARNING(
         "No liquid or interface cell is in the neighborhood to distribute excess mass to. Mass is lost/gained.");
      return;
   }

   // check if there are neighboring new interface cells
   const bool preferNewInterface = Base_T::excessMassDistributionModel_.getModelType() ==
                                      ExcessMassDistributionModel::ExcessMassModel::EvenlyNewInterfaceFallbackLiquid &&
                                   newInterfaceNeighbors > uint_c(0);

   // check if there are neighboring interface cells
   const bool preferInterface = (Base_T::excessMassDistributionModel_.getModelType() ==
                                    ExcessMassDistributionModel::ExcessMassModel::EvenlyAllInterfaceFallbackLiquid ||
                                 !preferNewInterface) &&
                                interfaceNeighbors > uint_c(0) &&
                                Base_T::excessMassDistributionModel_.getModelType() !=
                                   ExcessMassDistributionModel::ExcessMassModel::EvenlyAllInterfaceAndLiquid;

   // compute mass to be distributed to neighboring cells
   real_t deltaMass;
   if (preferNewInterface) { deltaMass = excessMass / real_c(newInterfaceNeighbors); }
   else
   {
      if (preferInterface) { deltaMass = excessMass / real_c(interfaceNeighbors); }
      else { deltaMass = excessMass / real_c(liquidAndInterfaceNeighbors); }
   }

   // distribute the excess mass
   for (auto pushDir = LatticeModel_T::Stencil::beginNoCenter(); pushDir != LatticeModel_T::Stencil::end(); ++pushDir)
   {
      const Cell neighborCell = Cell(cell.x() + pushDir.cx(), cell.y() + pushDir.cy(), cell.z() + pushDir.cz());

      // do not push mass to cells in the ghost layer (done by the process from which the ghost layer is synchronized)
      // - inner domain: 0 to *Size()-1
      // - inner domain incl. first ghost layer: -1 to *Size()
      if (neighborCell.x() <= cell_idx_c(-1) || neighborCell.y() <= cell_idx_c(-1) ||
          neighborCell.z() <= cell_idx_c(-1) || neighborCell.x() >= cell_idx_c(fillField->xSize()) ||
          neighborCell.y() >= cell_idx_c(fillField->ySize()) || neighborCell.z() >= cell_idx_c(fillField->zSize()))
      {
         continue;
      }

      // distribute excess mass to newly converted interface cell
      if (flagField->isMaskSet(neighborCell, Base_T::flagInfo_.convertedFlag | Base_T::flagInfo_.interfaceFlag))
      {
         // get density of neighboring interface cell
         const real_t neighborDensity = pdfField->getDensity(neighborCell);

         // add excess mass directly to fill level for newly converted neighboring interface cells
         fillField->get(neighborCell) += deltaMass / neighborDensity;
      }
      else
      {
         // distribute excess mass to old interface cell
         if (flagField->isFlagSet(neighborCell, Base_T::flagInfo_.interfaceFlag) && !preferNewInterface)
         {
            // get density of neighboring interface cell
            const real_t neighborDensity = pdfField->getDensity(neighborCell);

            // add excess mass directly to fill level for neighboring interface cells
            fillField->get(neighborCell) += deltaMass / neighborDensity;
         }
         else
         {
            // distribute excess mass to liquid cell
            if (flagField->isFlagSet(neighborCell, Base_T::flagInfo_.liquidFlag) && !preferInterface &&
                !preferNewInterface)
            {
               // add excess mass to excessMassField for neighboring liquid cells
               dstExcessMassField->get(neighborCell) += deltaMass;
            }
         }
      }
   }
}

} // namespace free_surface
} // namespace walberla
