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
//! \file StreamReconstructAdvectSweep.h
//! \ingroup free_surface
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Sweep for reconstruction of PDFs, streaming of PDFs (only in interface cells), advection of mass, update of
//!        bubble volumes and marking interface cells for conversion.
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "core/mpi/Reduce.h"
#include "core/timing/TimingPool.h"

#include "field/FieldClone.h"
#include "field/FlagField.h"

#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/bubble_model/BubbleModel.h"
#include "lbm/sweeps/StreamPull.h"

#include "PdfReconstructionModel.h"
#include "functionality/AdvectMass.h"
#include "functionality/FindInterfaceCellConversion.h"
#include "functionality/GetOredNeighborhood.h"
#include "functionality/ReconstructInterfaceCellABB.h"

namespace walberla
{
namespace free_surface
{
template< typename LatticeModel_T, typename BoundaryHandling_T, typename FlagField_T, typename FlagInfo_T,
          typename ScalarField_T, typename VectorField_T, bool useCodegen >
class StreamReconstructAdvectSweep
{
 public:
   using flag_t     = typename FlagInfo_T::flag_t;
   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   StreamReconstructAdvectSweep(real_t sigma, BlockDataID handlingID, BlockDataID fillFieldID, BlockDataID flagFieldID,
                                BlockDataID pdfField, ConstBlockDataID normalFieldID, ConstBlockDataID curvatureFieldID,
                                const FlagInfo_T& flagInfo, BubbleModelBase* bubbleModel,
                                const PdfReconstructionModel& pdfReconstructionModel, bool useSimpleMassExchange,
                                real_t cellConversionThreshold, real_t cellConversionForceThreshold)
      : sigma_(sigma), handlingID_(handlingID), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID),
        pdfFieldID_(pdfField), normalFieldID_(normalFieldID), curvatureFieldID_(curvatureFieldID), flagInfo_(flagInfo),
        bubbleModel_(bubbleModel), neighborhoodFlagFieldClone_(flagFieldID), fillFieldClone_(fillFieldID),
        pdfFieldClone_(pdfField), pdfReconstructionModel_(pdfReconstructionModel),
        useSimpleMassExchange_(useSimpleMassExchange), cellConversionThreshold_(cellConversionThreshold),
        cellConversionForceThreshold_(cellConversionForceThreshold)
   {}

   void operator()(IBlock* const block);

 protected:
   real_t sigma_; // surface tension

   BlockDataID handlingID_;
   BlockDataID fillFieldID_;
   BlockDataID flagFieldID_;
   BlockDataID pdfFieldID_;

   ConstBlockDataID normalFieldID_;
   ConstBlockDataID curvatureFieldID_;

   FlagInfo_T flagInfo_;
   bubble_model::BubbleModelBase* const bubbleModel_;

   // efficient clones of fields to provide temporary fields (for writing to)
   field::FieldClone< FlagField_T, true > neighborhoodFlagFieldClone_;
   field::FieldClone< ScalarField_T, true > fillFieldClone_;
   field::FieldClone< PdfField_T, true > pdfFieldClone_;

   PdfReconstructionModel pdfReconstructionModel_;
   bool useSimpleMassExchange_;
   real_t cellConversionThreshold_;
   real_t cellConversionForceThreshold_;
}; // class StreamReconstructAdvectSweep

template< typename LatticeModel_T, typename BoundaryHandling_T, typename FlagField_T, typename FlagInfo_T,
          typename ScalarField_T, typename VectorField_T, bool useCodegen >
void StreamReconstructAdvectSweep< LatticeModel_T, BoundaryHandling_T, FlagField_T, FlagInfo_T, ScalarField_T,
                                   VectorField_T, useCodegen >::operator()(IBlock* const block)
{
   // fetch data
   ScalarField_T* const fillSrcField = block->getData< ScalarField_T >(fillFieldID_);
   PdfField_T* const pdfSrcField     = block->getData< PdfField_T >(pdfFieldID_);

   const ScalarField_T* const curvatureField = block->getData< const ScalarField_T >(curvatureFieldID_);
   const VectorField_T* const normalField    = block->getData< const VectorField_T >(normalFieldID_);
   FlagField_T* const flagField              = block->getData< FlagField_T >(flagFieldID_);

   // temporary fields that act as destination fields
   PdfField_T* const pdfDstField            = pdfFieldClone_.get(block);
   FlagField_T* const neighborhoodFlagField = neighborhoodFlagFieldClone_.get(block);
   ScalarField_T* const fillDstField        = fillFieldClone_.get(block);

   // combine all neighbor flags using bitwise OR and write them to the neighborhood field
   // this is simply a pre-computation of often required values
   // IMPORTANT REMARK: the "OredNeighborhood" is also required for each cell in the first ghost layer; this requires
   // access to all first ghost layer cell's neighbors (i.e., to the second ghost layer)
   WALBERLA_CHECK_GREATER_EQUAL(flagField->nrOfGhostLayers(), uint_c(2));
   getOredNeighborhood< typename LatticeModel_T::Stencil >(flagField, neighborhoodFlagField);

   // explicitly avoid OpenMP due to bubble model update (reportFillLevelChange)
   WALBERLA_FOR_ALL_CELLS_OMP(
      pdfDstFieldIt, pdfDstField, pdfSrcFieldIt, pdfSrcField, fillSrcFieldIt, fillSrcField, fillDstFieldIt,
      fillDstField, flagFieldIt, flagField, neighborhoodFlagFieldIt, neighborhoodFlagField, normalFieldIt, normalField,
      curvatureFieldIt, curvatureField, omp critical, {
         if (flagInfo_.isInterface(flagFieldIt))
         {
            // get density (rhoGas) for interface PDF reconstruction
            const real_t bubbleDensity = bubbleModel_->getDensity(block, pdfSrcFieldIt.cell());
            const real_t rhoGas        = computeDeltaRhoLaplacePressure(sigma_, *curvatureFieldIt) + bubbleDensity;

            // reconstruct missing PDFs coming from gas neighbors according to the specified model (see dissertation of
            // N. Thuerey, 2007, section 4.2); reconstruction includes streaming of PDFs to interface cells (no LBM
            // stream required here)
            (reconstructInterfaceCellLegacy< LatticeModel_T, FlagField_T >) (flagField, pdfSrcFieldIt, flagFieldIt,
                                                                             normalFieldIt, flagInfo_, rhoGas,
                                                                             pdfDstFieldIt, pdfReconstructionModel_);

            // density before LBM stream (post-collision)
            const real_t oldRho = lbm::getDensity(pdfSrcField->latticeModel(), pdfSrcFieldIt);

            // density after LBM stream in reconstruction
            const real_t newRho = lbm::getDensity(pdfDstField->latticeModel(), pdfDstFieldIt);

            // compute mass advection using post-collision PDFs (explicitly not PDFs updated by stream above)
            const real_t deltaMass =
               (advectMass< LatticeModel_T, FlagField_T, typename ScalarField_T::iterator,
                            typename PdfField_T::iterator, typename FlagField_T::iterator,
                            typename FlagField_T::iterator, FlagInfo_T >) (flagField, fillSrcFieldIt, pdfSrcFieldIt,
                                                                           flagFieldIt, neighborhoodFlagFieldIt,
                                                                           flagInfo_, useSimpleMassExchange_);

            // update fill level after LBM stream and mass exchange
            *fillDstFieldIt        = (*fillSrcFieldIt * oldRho + deltaMass) / newRho;
            const real_t deltaFill = *fillDstFieldIt - *fillSrcFieldIt;

            // update the volume of bubbles
            bubbleModel_->reportFillLevelChange(block, fillSrcFieldIt.cell(), deltaFill);
         }
         else // treat non-interface cells
         {
            // manually adjust the fill level to avoid outdated fill levels being copied from fillSrcField
            if (flagInfo_.isGas(flagFieldIt)) { *fillDstFieldIt = real_c(0); }
            else
            {
               if (flagInfo_.isLiquid(flagFieldIt))
               {
                  const Cell cell = pdfSrcFieldIt.cell();
                  if constexpr (useCodegen)
                  {
                     auto lbmSweepGenerated = typename LatticeModel_T::Sweep(pdfFieldID_);
                     const CellInterval ci(cell, cell);
                     lbmSweepGenerated.streamInCellInterval(pdfSrcField, pdfDstField, ci);
                  }
                  else
                  {
                     lbm::StreamPull< LatticeModel_T >::execute(pdfSrcField, pdfDstField, cell[0], cell[1], cell[2]);
                  }

                  *fillDstFieldIt = real_c(1);
               }
               else // flag is e.g. obstacle or outflow
               {
                  *fillDstFieldIt = *fillSrcFieldIt;
               }
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS_XYZ_OMP

   pdfSrcField->swapDataPointers(pdfDstField);
   fillSrcField->swapDataPointers(fillDstField);

   BoundaryHandling_T* const handling = block->getData< BoundaryHandling_T >(handlingID_);

   // mark interface cell for conversion
   findInterfaceCellConversions< LatticeModel_T >(handling, fillSrcField, flagField, neighborhoodFlagField, flagInfo_,
                                                  cellConversionThreshold_, cellConversionForceThreshold_);
}

} // namespace free_surface
} // namespace walberla
