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
//! \file PdfRefillingSweep.impl.h
//! \ingroup dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \author Michael Zikeli
//! \brief Sweeps for refilling cells (i.e. reinitializing PDFs) after the cell was converted from gas to interface.
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/cell/Cell.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"

#include "domain_decomposition/IBlock.h"

#include "field/GhostLayerField.h"

#include "lbm/field/Equilibrium.h"
#include "lbm/field/PdfField.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/dynamics/PdfRefillingModel.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"

#include "stencil/D3Q6.h"

#include <set>
#include <vector>

#include "PdfRefillingSweep.h"

namespace walberla
{
namespace free_surface
{
template< typename LatticeModel_T, typename FlagField_T >
real_t RefillingSweepBase< LatticeModel_T, FlagField_T >::getAverageDensityAndVelocity(
   const Cell& cell, const PdfField_T& pdfField, const FlagField_T& flagField, const FlagInfo< FlagField_T >& flagInfo,
   Vector3< real_t >& avgVelocity, std::vector< bool >& validStencilIndices)
{
   real_t rho = real_c(0.0);
   Vector3< real_t > u(real_c(0.0));
   uint_t numNeighbors = uint_c(0);

   // do not use data from ghost layer if optimized communication is used (see comment in PdfRefillingSweep.h at
   // ExtrapolationRefillingSweepBase)
   const CellInterval localDomain = useDataFromGhostLayers_ ? pdfField.xyzSizeWithGhostLayer() : pdfField.xyzSize();

   for (auto i = Stencil_T::beginNoCenter(); i != Stencil_T::end(); ++i)
   {
      const Cell neighborCell(cell[0] + i.cx(), cell[1] + i.cy(), cell[2] + i.cz());

      const flag_t neighborFlag        = flagField.get(neighborCell);
      const flag_t liquidInterfaceMask = flagInfo.interfaceFlag | flagInfo.liquidFlag;

      // only use neighboring cell if
      // - neighboring cell is part of the block-local domain
      // - neighboring cell is liquid or interface
      // - not newly converted from G->I
      const bool useNeighbor = isPartOfMaskSet(neighborFlag, liquidInterfaceMask) &&
                               !flagInfo.hasConvertedFromGasToInterface(flagField.get(neighborCell)) &&
                               localDomain.contains(neighborCell);

      // calculate the average of valid neighbor cells to calculate an average density and velocity.
      if (useNeighbor)
      {
         numNeighbors++;
         const typename PdfField_T::ConstPtr neighbor(pdfField, neighborCell[0], neighborCell[1], neighborCell[2]);
         Vector3< real_t > neighborU;
         real_t neighborRho;
         neighborRho = lbm::getDensityAndMomentumDensity(neighborU, pdfField.latticeModel(), neighbor);
         neighborU /= neighborRho;
         u += neighborU;
         rho += neighborRho;

         validStencilIndices[i.toIdx()] = true;
      }
   }

   // normalize the newly calculated velocity and density
   if (numNeighbors != uint_c(0))
   {
      u /= real_c(numNeighbors);
      rho /= real_c(numNeighbors);
   }
   else
   {
      u   = Vector3< real_t >(0.0);
      rho = real_c(1.0);
      WALBERLA_LOG_WARNING_ON_ROOT("There are no valid neighbors for the refilling of (block-local) cell: " << cell);
   }

   avgVelocity = u;
   return rho;
}

template< typename LatticeModel_T, typename FlagField_T >
std::vector< real_t > RefillingSweepBase< LatticeModel_T, FlagField_T >::getAveragePdfs(
   const Cell& cell, const PdfField_T& pdfField, const FlagField_T& flagField, const FlagInfo< FlagField_T >& flagInfo)
{
   uint_t numNeighbors = uint_c(0);

   // do not use data from ghost layer if optimized communication is used (see comment in PdfRefillingSweep.h at
   // ExtrapolationRefillingSweepBase)
   const CellInterval localDomain = useDataFromGhostLayers_ ? pdfField.xyzSizeWithGhostLayer() : pdfField.xyzSize();

   std::vector< real_t > pdfSum(Stencil_T::Size, real_c(0));

   for (auto i = Stencil_T::beginNoCenter(); i != Stencil_T::end(); ++i)
   {
      const Cell neighborCell(cell[0] + i.cx(), cell[1] + i.cy(), cell[2] + i.cz());

      const flag_t neighborFlag        = flagField.get(neighborCell);
      const flag_t liquidInterfaceMask = flagInfo.interfaceFlag | flagInfo.liquidFlag;

      // only use neighboring cell if
      // - neighboring cell is part of the block-local domain
      // - neighboring cell is liquid or interface
      // - not newly converted from G->I
      const bool useNeighbor = isPartOfMaskSet(neighborFlag, liquidInterfaceMask) &&
                               !flagInfo.hasConvertedFromGasToInterface(flagField.get(neighborCell)) &&
                               localDomain.contains(neighborCell);

      // calculate the average of valid neighbor cells to calculate an average of PDFs in each direction
      if (useNeighbor)
      {
         ++numNeighbors;
         for (auto pdfDir = Stencil_T::begin(); pdfDir != Stencil_T::end(); ++pdfDir)
         {
            pdfSum[pdfDir.toIdx()] += pdfField.get(neighborCell, *pdfDir);
         }
      }
   }

   // average the PDFs of all neighboring cells
   if (numNeighbors != uint_c(0))
   {
      for (auto& pdf : pdfSum)
      {
         pdf /= real_c(numNeighbors);
      }
   }
   else
   {
      // fall back to EquilibriumRefilling by setting PDFs according to equilibrium with velocity=0 and density=1
      lbm::Equilibrium< LatticeModel_T >::set(pdfSum, Vector3< real_t >(real_c(0)), real_c(1));

      WALBERLA_LOG_WARNING_ON_ROOT("There are no valid neighbors for the refilling of (block-local) cell: " << cell);
   }

   return pdfSum;
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
Vector3< cell_idx_t > ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   findCorrespondingLatticeDirection(const Vector3< real_t >& direction)
{
   if (direction == Vector3< real_t >(real_c(0))) { return direction; }

   stencil::Direction bestFittingDirection = stencil::C; // arbitrary default initialization
   real_t scalarProduct                    = real_c(0);

   for (auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
   {
      // compute inner product <dir,c_i>
      const real_t scalarProductTmp = direction[0] * stencil::cNorm[0][*dir] + direction[1] * stencil::cNorm[1][*dir] +
                                      direction[2] * stencil::cNorm[2][*dir];
      if (scalarProductTmp > scalarProduct)
      {
         // largest scalar product is the best fitting lattice direction, i.e., this direction has the smallest angle to
         // the given direction
         scalarProduct        = scalarProductTmp;
         bestFittingDirection = *dir;
      }
   }

   return Vector3< cell_idx_t >(stencil::cx[bestFittingDirection], stencil::cy[bestFittingDirection],
                                stencil::cz[bestFittingDirection]);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
Vector3< cell_idx_t >
   ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T,
                                    VectorField_T >::findExtrapolationDirection(const Cell& cell,
                                                                                const FlagField_T& flagField,
                                                                                const ScalarField_T& fillField)
{
   Vector3< real_t > normal = Vector3< real_t >(real_c(0));

   // get flag mask for obstacle boundaries
   const flag_t obstacleFlagMask = flagField.getMask(RefillingSweepBase_T::flagInfo_.getObstacleIDSet());

   // get flag mask for liquid, interface, and gas cells
   const flag_t liquidInterfaceGasMask = flagField.getMask(flagIDs::liquidInterfaceGasFlagIDs);

   const typename ScalarField_T::ConstPtr fillPtr(fillField, cell.x(), cell.y(), cell.z());
   const typename FlagField_T::ConstPtr flagPtr(flagField, cell.x(), cell.y(), cell.z());

   if (!isFlagInNeighborhood< Stencil_T >(flagPtr, obstacleFlagMask))
   {
      normal_computation::computeNormal<
         typename std::conditional< Stencil_T::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type >(
         normal, fillPtr, flagPtr, liquidInterfaceGasMask);
   }
   else
   {
      normal_computation::computeNormalNearSolidBoundary<
         typename std::conditional< Stencil_T::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type >(
         normal, fillPtr, flagPtr, liquidInterfaceGasMask, obstacleFlagMask);
   }

   // normalize normal (in contrast to the usual definition in FSLBM, it points from gas to fluid here)
   normal = normal.getNormalizedOrZero();

   // find the lattice direction that most closely resembles the normal direction
   // Note: the normal must point from gas to fluid because the extrapolation direction is also defined to point from
   // gas to fluid
   return findCorrespondingLatticeDirection(normal);
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
uint_t ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   getNumberOfExtrapolationCells(const Cell& cell, const FlagField_T& flagField, const PdfField_T& pdfField,
                                 const Vector3< cell_idx_t >& extrapolationDirection)
{
   // skip center cell
   if (extrapolationDirection == Vector3< cell_idx_t >(cell_idx_c(0))) { return uint_c(0); }

   // do not use data from ghost layer if optimized communication is used (see comment in PdfRefillingSweep.h at
   // ExtrapolationRefillingSweepBase)
   const CellInterval localDomain =
      (RefillingSweepBase_T::useDataFromGhostLayers_) ? pdfField.xyzSizeWithGhostLayer() : pdfField.xyzSize();

   // for extrapolation order n, n+1 applicable cells must be available in the desired direction
   const uint_t maxExtrapolationCells = extrapolationOrder_ + uint_c(1);

   for (uint_t numCells = uint_c(1); numCells <= maxExtrapolationCells; ++numCells)
   {
      // potential cell used for extrapolation
      const Cell checkCell = cell + Cell(cell_idx_c(numCells) * extrapolationDirection);

      const flag_t neighborFlag = flagField.get(checkCell);
      const flag_t liquidInterfaceMask =
         RefillingSweepBase_T::flagInfo_.interfaceFlag | RefillingSweepBase_T::flagInfo_.liquidFlag;

      // only use cell if the cell is
      // - inside domain
      // - liquid or interface
      // - not a cell that has also just been converted from gas to interface
      if (!localDomain.contains(checkCell) || !isPartOfMaskSet(neighborFlag, liquidInterfaceMask) ||
          RefillingSweepBase_T::flagInfo_.hasConvertedFromGasToInterface(flagField.get(checkCell)))
      {
         return numCells - uint_c(1);
      }
   }

   return maxExtrapolationCells;
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
std::vector< real_t > ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   getNonEquilibriumPdfsInCell(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField)
{
   std::vector< real_t > nonEquilibriumPartOfPdfs(Stencil_T::Size);

   Vector3< real_t > velocity;
   const real_t density = pdfField.getDensityAndVelocity(velocity, cell);

   // compute non-equilibrium of each PDF
   for (auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir)
   {
      nonEquilibriumPartOfPdfs[dir.toIdx()] =
         pdfField.get(cell, dir.toIdx()) - lbm::EquilibriumDistribution< LatticeModel_T >::get(*dir, velocity, density);
   }
   return nonEquilibriumPartOfPdfs;
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
std::vector< real_t >
   ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::getPdfsInCell(
      const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField)
{
   std::vector< real_t > Pdfs(Stencil_T::Size);

   for (auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir)
   {
      Pdfs[dir.toIdx()] = pdfField.get(cell, dir.toIdx());
   }
   return Pdfs;
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   applyQuadraticExtrapolation(
      const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField, const Vector3< cell_idx_t >& extrapolationDirection,
      bool includeThisCell,
      const std::function< std::vector< real_t >(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField) >&
         getPdfFunc)
{
   // store the PDFs of the cells participating in the extrapolation
   std::vector< real_t > pdfsXf(LatticeModel_T::Stencil::Size);   // cell + 1 * extrapolationDirection
   std::vector< real_t > pdfsXff(LatticeModel_T::Stencil::Size);  // cell + 2 * extrapolationDirection
   std::vector< real_t > pdfsXfff(LatticeModel_T::Stencil::Size); // cell + 3 * extrapolationDirection

   // determines if the PDFs of the current cell are also considered
   const real_t centerFactor = includeThisCell ? real_c(1) : real_c(0);

   // get PDFs
   pdfsXf   = getPdfFunc(cell + Cell(cell_idx_c(1) * extrapolationDirection), pdfField);
   pdfsXff  = getPdfFunc(cell + Cell(cell_idx_c(2) * extrapolationDirection), pdfField);
   pdfsXfff = getPdfFunc(cell + Cell(cell_idx_c(3) * extrapolationDirection), pdfField);

   // compute the resulting PDF values of "cell" according to second order extrapolation
   for (auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir)
   {
      pdfField.get(cell, dir.toIdx()) = centerFactor * pdfField.get(cell, dir.toIdx()) +
                                        real_c(3) * pdfsXf[dir.toIdx()] - real_c(3) * pdfsXff[dir.toIdx()] +
                                        real_c(1) * pdfsXfff[dir.toIdx()];
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   applyLinearExtrapolation(
      const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField, const Vector3< cell_idx_t >& extrapolationDirection,
      bool includeThisCell,
      const std::function< std::vector< real_t >(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField) >&
         getPdfFunc)
{
   // store the PDFs of the cells participating in the interpolation
   std::vector< real_t > pdfsXf(Stencil_T::Size);  // cell + 1 * extrapolationDirection
   std::vector< real_t > pdfsXff(Stencil_T::Size); // cell + 2 * extrapolationDirection

   // determines if the PDFs of the current cell are also considered
   const real_t centerFactor = includeThisCell ? real_c(1) : real_c(0);

   // get PDFs
   pdfsXf  = getPdfFunc(cell + Cell(cell_idx_c(1) * extrapolationDirection), pdfField);
   pdfsXff = getPdfFunc(cell + Cell(cell_idx_c(2) * extrapolationDirection), pdfField);

   // compute the resulting PDF values of "cell" according to first order extrapolation
   for (auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir)
   {
      pdfField.get(cell, dir.toIdx()) = centerFactor * pdfField.get(cell, dir.toIdx()) +
                                        real_c(2) * pdfsXf[dir.toIdx()] - real_c(1) * pdfsXff[dir.toIdx()];
   }
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExtrapolationRefillingSweepBase< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::
   applyConstantExtrapolation(
      const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField, const Vector3< cell_idx_t >& extrapolationDirection,
      bool includeThisCell,
      const std::function< std::vector< real_t >(const Cell& cell, lbm::PdfField< LatticeModel_T >& pdfField) >&
         getPdfFunc)
{
   // store the PDFs of the cells participating in the interpolation
   std::vector< real_t > pdfsXf(Stencil_T::Size); // cell + 1 * extrapolationDirection

   // determines if the PDFs of the current cell are also considered
   const real_t centerFactor = includeThisCell ? real_c(1) : real_c(0);

   // get PDFs
   pdfsXf = getPdfFunc(cell + Cell(cell_idx_c(1) * extrapolationDirection), pdfField);

   // compute the resulting PDF values of "cell" according to zeroth order extrapolation
   for (auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir)
   {
      pdfField.get(cell, dir.toIdx()) =
         centerFactor * pdfField.get(cell, dir.toIdx()) + real_c(1) * pdfsXf[dir.toIdx()];
   }
}

template< typename LatticeModel_T, typename FlagField_T >
void EquilibriumRefillingSweep< LatticeModel_T, FlagField_T >::operator()(IBlock* const block)
{
   PdfField_T* const pdfField         = block->getData< PdfField_T >(RefillingSweepBase_T::pdfFieldID_);
   const FlagField_T* const flagField = block->getData< const FlagField_T >(RefillingSweepBase_T::flagFieldID_);

   WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
      if (RefillingSweepBase_T::flagInfo_.hasConvertedFromGasToInterface(flagFieldIt))
      {
         const Cell cell = pdfFieldIt.cell();

         Vector3< real_t > avgVelocity;
         real_t avgDensity;
         avgDensity = RefillingSweepBase_T::getAverageDensityAndVelocity(cell, *pdfField, *flagField,
                                                                         RefillingSweepBase_T::flagInfo_, avgVelocity);

         pdfField->setDensityAndVelocity(cell, avgVelocity, avgDensity);
      }
   }) // WALBERLA_FOR_ALL_CELLS
}

template< typename LatticeModel_T, typename FlagField_T >
void AverageRefillingSweep< LatticeModel_T, FlagField_T >::operator()(IBlock* const block)
{
   PdfField_T* const pdfField         = block->getData< PdfField_T >(RefillingSweepBase_T::pdfFieldID_);
   const FlagField_T* const flagField = block->getData< const FlagField_T >(RefillingSweepBase_T::flagFieldID_);

   WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
      if (RefillingSweepBase_T::flagInfo_.hasConvertedFromGasToInterface(flagFieldIt))
      {
         const Cell cell = pdfFieldIt.cell();

         // compute average PDFs (in each direction) from all applicable neighboring cells
         const std::vector< real_t > pdfAverage =
            RefillingSweepBase_T::getAveragePdfs(cell, *pdfField, *flagField, RefillingSweepBase_T::flagInfo_);

         for (auto pdfDir = Stencil_T::begin(); pdfDir != Stencil_T::end(); ++pdfDir)
         {
            pdfField->get(cell, *pdfDir) = pdfAverage[pdfDir.toIdx()];
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void EquilibriumAndNonEquilibriumRefillingSweep< LatticeModel_T, FlagField_T, ScalarField_T,
                                                 VectorField_T >::operator()(IBlock* const block)
{
   PdfField_T* const pdfField =
      block->getData< PdfField_T >(ExtrapolationRefillingSweepBase_T::RefillingSweepBase_T::pdfFieldID_);
   const FlagField_T* const flagField =
      block->getData< const FlagField_T >(ExtrapolationRefillingSweepBase_T::RefillingSweepBase_T::flagFieldID_);
   const ScalarField_T* const fillField =
      block->getData< const ScalarField_T >(ExtrapolationRefillingSweepBase_T::fillFieldID_);

   // function to fetch relevant PDFs
   auto getPdfFunc = std::bind(&ExtrapolationRefillingSweepBase_T::getNonEquilibriumPdfsInCell, this,
                               std::placeholders::_1, std::placeholders::_2);

   WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
      if (RefillingSweepBase_T::flagInfo_.hasConvertedFromGasToInterface(flagFieldIt))
      {
         const Cell cell = pdfFieldIt.cell();

         // apply EquilibriumRefilling first
         Vector3< real_t > avgVelocity;
         real_t avgDensity;
         avgDensity = RefillingSweepBase_T::getAverageDensityAndVelocity(cell, *pdfField, *flagField,
                                                                         RefillingSweepBase_T::flagInfo_, avgVelocity);
         pdfField->setDensityAndVelocity(cell, avgVelocity, avgDensity);

         // find valid cells for extrapolation
         const Vector3< cell_idx_t > extrapolationDirection =
            ExtrapolationRefillingSweepBase_T::findExtrapolationDirection(cell, *flagField, *fillField);
         const uint_t numberOfCellsForExtrapolation = ExtrapolationRefillingSweepBase_T::getNumberOfExtrapolationCells(
            cell, *flagField, *pdfField, extrapolationDirection);

         // add non-equilibrium part of PDF (which might be obtained by extrapolation)
         if (numberOfCellsForExtrapolation >= uint_c(3))
         {
            ExtrapolationRefillingSweepBase_T::applyQuadraticExtrapolation(cell, *pdfField, extrapolationDirection,
                                                                           true, getPdfFunc);
         }
         else
         {
            if (numberOfCellsForExtrapolation >= uint_c(2))
            {
               ExtrapolationRefillingSweepBase_T::applyLinearExtrapolation(cell, *pdfField, extrapolationDirection,
                                                                           true, getPdfFunc);
            }
            else
            {
               if (numberOfCellsForExtrapolation >= uint_c(1))
               {
                  ExtrapolationRefillingSweepBase_T::applyConstantExtrapolation(cell, *pdfField, extrapolationDirection,
                                                                                true, getPdfFunc);
               }
               // else: do nothing here; this corresponds to using EquilibriumRefilling (done already at the beginning)
            }
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS
}

template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
void ExtrapolationRefillingSweep< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >::operator()(
   IBlock* const block)
{
   PdfField_T* const pdfField = block->getData< PdfField_T >(ExtrapolationRefillingSweepBase_T::pdfFieldID_);
   const FlagField_T* const flagField =
      block->getData< const FlagField_T >(ExtrapolationRefillingSweepBase_T::flagFieldID_);
   const ScalarField_T* const fillField =
      block->getData< const ScalarField_T >(ExtrapolationRefillingSweepBase_T::fillFieldID_);

   WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
      if (RefillingSweepBase_T::flagInfo_.hasConvertedFromGasToInterface(flagFieldIt))
      {
         const Cell cell = pdfFieldIt.cell();

         // find valid cells for extrapolation
         const Vector3< cell_idx_t > extrapolationDirection =
            ExtrapolationRefillingSweepBase_T::findExtrapolationDirection(cell, *flagField, *fillField);
         const uint_t numberOfCellsForExtrapolation = ExtrapolationRefillingSweepBase_T::getNumberOfExtrapolationCells(
            cell, *flagField, *pdfField, extrapolationDirection);

         // function to fetch relevant PDFs
         auto getPdfFunc = std::bind(&ExtrapolationRefillingSweepBase_T::getPdfsInCell, this, std::placeholders::_1,
                                     std::placeholders::_2);

         if (numberOfCellsForExtrapolation >= uint_c(3))
         {
            ExtrapolationRefillingSweepBase_T::applyQuadraticExtrapolation(cell, *pdfField, extrapolationDirection,
                                                                           false, getPdfFunc);
         }
         else
         {
            if (numberOfCellsForExtrapolation >= uint_c(2))
            {
               ExtrapolationRefillingSweepBase_T::applyLinearExtrapolation(cell, *pdfField, extrapolationDirection,
                                                                           false, getPdfFunc);
            }
            else
            {
               if (numberOfCellsForExtrapolation >= uint_c(1))
               {
                  ExtrapolationRefillingSweepBase_T::applyConstantExtrapolation(cell, *pdfField, extrapolationDirection,
                                                                                false, getPdfFunc);
               }
               else
               {
                  // if not enough cells for extrapolation are available, use EquilibriumRefilling
                  Vector3< real_t > avgVelocity;
                  real_t avgDensity;
                  avgDensity = RefillingSweepBase_T::getAverageDensityAndVelocity(
                     cell, *pdfField, *flagField, RefillingSweepBase_T::flagInfo_, avgVelocity);
                  pdfField->setDensityAndVelocity(cell, avgVelocity, avgDensity);
               }
            }
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS
}

template< typename LatticeModel_T, typename FlagField_T >
Vector3< real_t > GradsMomentsRefillingSweep< LatticeModel_T, FlagField_T >::getVelocityGradient(
   stencil::Direction direction, const Cell& cell, const PdfField_T* pdfField, const Vector3< real_t >& avgVelocity,
   const std::vector< bool >& validStencilIndices)
{
   stencil::Direction dir;
   stencil::Direction invDir;

   switch (direction)
   {
   case stencil::E:
   case stencil::W:
      dir    = stencil::E;
      invDir = stencil::W;
      break;
   case stencil::N:
   case stencil::S:
      dir    = stencil::N;
      invDir = stencil::S;
      break;
   case stencil::T:
   case stencil::B:
      dir    = stencil::T;
      invDir = stencil::B;
      break;
   default:
      WALBERLA_ABORT("Velocity gradient for GradsMomentsRefilling can not be computed in a direction other than in x-, "
                     "y-, or z-direction.");
   }

   Vector3< real_t > velocityGradient(real_c(0));

   // apply central finite differences if both neighboring cells are available
   if (validStencilIndices[Stencil_T::idx[dir]] && validStencilIndices[Stencil_T::idx[invDir]])
   {
      const Vector3< real_t > neighborVelocity1 = pdfField->getVelocity(cell + dir);
      const Vector3< real_t > neighborVelocity2 = pdfField->getVelocity(cell + invDir);
      velocityGradient[0] = real_c(0.5) * (neighborVelocity1[0] - neighborVelocity2[0]); // assuming dx = 1
      velocityGradient[1] = real_c(0.5) * (neighborVelocity1[1] - neighborVelocity2[1]); // assuming dx = 1
      velocityGradient[2] = real_c(0.5) * (neighborVelocity1[2] - neighborVelocity2[2]); // assuming dx = 1
   }
   else
   {
      // apply first order upwind scheme
      stencil::Direction upwindDirection = (avgVelocity[0] > real_c(0)) ? invDir : dir;

      stencil::Direction sourceDirection = stencil::C; // arbitrary default initialization

      if (validStencilIndices[Stencil_T::idx[upwindDirection]]) { sourceDirection = upwindDirection; }
      else
      {
         if (validStencilIndices[Stencil_T::idx[stencil::inverseDir[upwindDirection]]])
         {
            sourceDirection = stencil::inverseDir[upwindDirection];
         }
      }

      if (sourceDirection == dir)
      {
         auto neighborVelocity = pdfField->getVelocity(cell + sourceDirection);
         velocityGradient[0]   = neighborVelocity[0] - avgVelocity[0]; // assuming dx = 1
         velocityGradient[1]   = neighborVelocity[1] - avgVelocity[1]; // assuming dx = 1
         velocityGradient[2]   = neighborVelocity[2] - avgVelocity[2]; // assuming dx = 1
      }
      else
      {
         if (sourceDirection == invDir)
         {
            auto neighborVelocity = pdfField->getVelocity(cell + sourceDirection);
            velocityGradient[0]   = avgVelocity[0] - neighborVelocity[0]; // assuming dx = 1
            velocityGradient[1]   = avgVelocity[1] - neighborVelocity[1]; // assuming dx = 1
            velocityGradient[2]   = avgVelocity[2] - neighborVelocity[2]; // assuming dx = 1
         }
         // else: no stencil direction is valid, velocityGradient is zero
      }
   }

   return velocityGradient;
}

template< typename LatticeModel_T, typename FlagField_T >
void GradsMomentsRefillingSweep< LatticeModel_T, FlagField_T >::operator()(IBlock* const block)
{
   PdfField_T* pdfField               = block->getData< PdfField_T >(RefillingSweepBase_T::pdfFieldID_);
   const FlagField_T* const flagField = block->getData< const FlagField_T >(RefillingSweepBase_T::flagFieldID_);

   WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
      if (RefillingSweepBase_T::flagInfo_.hasConvertedFromGasToInterface(flagFieldIt))
      {
         const Cell cell = pdfFieldIt.cell();

         // get average density and velocity from valid neighboring cells and store the directions of valid neighbors
         std::vector< bool > validStencilIndices(Stencil_T::Size, false);
         Vector3< real_t > avgVelocity;
         real_t avgDensity;
         avgDensity = RefillingSweepBase_T::getAverageDensityAndVelocity(
            cell, *pdfField, *flagField, RefillingSweepBase_T::flagInfo_, avgVelocity, validStencilIndices);

         // get velocity gradients
         // - using a first order central finite differences (if two neighboring cells are available)
         // - using a first order upwind scheme (if only one neighboring cell is available)
         // - assuming a zero gradient if no valid neighboring cell is available
         // velocityGradient(u) =
         // | du1/dx1 du2/dx1 du3/dx1 |   | 0 1 2 |   | 0,0  0,1  0,2 |
         // | du1/dx2 du2/dx2 du3/dx2 | = | 3 4 5 | = | 1,0  1,1  1,2 |
         // | du1/dx3 du2/dx3 du3/dx3 |   | 6 7 8 |   | 2,0  2,1  2,2 |
         const Vector3< real_t > gradientX =
            getVelocityGradient(stencil::E, cell, pdfField, avgVelocity, validStencilIndices);
         const Vector3< real_t > gradientY =
            getVelocityGradient(stencil::N, cell, pdfField, avgVelocity, validStencilIndices);
         Vector3< real_t > gradientZ = Vector3< real_t >(real_c(0));
         if (Stencil_T::D == 3)
         {
            gradientZ = getVelocityGradient(stencil::T, cell, pdfField, avgVelocity, validStencilIndices);
         }

         Matrix3< real_t > velocityGradient(gradientX, gradientY, gradientZ);

         // compute non-equilibrium pressure tensor (equation (13) in Dorschner et al.); rho is not included in the
         // pre-factor here, but will be considered later
         Matrix3< real_t > pressureTensorNeq(real_c(0));

         // in equation (13) in Dorschner et al., 2*beta is used in the pre-factor; note that 2*beta=omega (relaxation
         // rate related to kinematic viscosity)
         const real_t preFac = -real_c(1) / (real_c(3) * relaxRate_);

         for (uint_t j = uint_c(0); j != uint_c(3); ++j)
         {
            for (uint_t i = uint_c(0); i != uint_c(3); ++i)
            {
               pressureTensorNeq(i, j) += preFac * (velocityGradient(i, j) + velocityGradient(j, i));
            }
         }

         // set PDFs according to equation (10) in Dorschner et al.; this is equivalent to setting the PDFs to
         // equilibrium f^{eq} and adding a contribution of the non-equilibrium pressure tensor P^{neq}
         for (auto q = Stencil_T::begin(); q != Stencil_T::end(); ++q)
         {
            const real_t velCi = lbm::internal::multiplyVelocityDirection(*q, avgVelocity);

            real_t contributionFromPneq = real_c(0);
            for (uint_t j = uint_c(0); j != uint_c(3); ++j)
            {
               for (uint_t i = uint_c(0); i != uint_c(3); ++i)
               {
                  // P^{neq}_{a,b} * c_{q,a} * c_{q,b}
                  contributionFromPneq +=
                     pressureTensorNeq(i, j) * real_c(stencil::c[i][*q]) * real_c(stencil::c[j][*q]);
               }
            }

            // - P^{neq}_{a,b} * cs^2 * delta_{a,b}
            contributionFromPneq -=
               (pressureTensorNeq(0, 0) + pressureTensorNeq(1, 1) + pressureTensorNeq(2, 2)) / real_c(3);

            // compute f^{eq} and add contribution from P^{neq}
            const real_t pdf = LatticeModel_T::w[q.toIdx()] * avgDensity *
                               (real_c(1) + real_c(3) * velCi - real_c(1.5) * avgVelocity.sqrLength() +
                                real_c(4.5) * velCi * velCi + real_c(4.5) * contributionFromPneq);
            pdfField->get(cell, *q) = pdf;
         }
      }
   }) // WALBERLA_FOR_ALL_CELLS
}

} // namespace free_surface
} // namespace walberla