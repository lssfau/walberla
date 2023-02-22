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
//! \file AdvectSweep.h
//! \ingroup free_surface
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Sweep for mass advection and interface cell conversion marking (simplified StreamReconstructAdvectSweep).
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
#include "lbm/free_surface/dynamics/PdfReconstructionModel.h"
#include "lbm/free_surface/dynamics/functionality/AdvectMass.h"
#include "lbm/free_surface/dynamics/functionality/FindInterfaceCellConversion.h"
#include "lbm/free_surface/dynamics/functionality/GetOredNeighborhood.h"
#include "lbm/free_surface/dynamics/functionality/ReconstructInterfaceCellABB.h"
#include "lbm/sweeps/StreamPull.h"

namespace walberla
{
namespace free_surface
{
template< typename LatticeModel_T, typename BoundaryHandling_T, typename FlagField_T, typename FlagInfo_T,
          typename ScalarField_T, typename VectorField_T >
class AdvectSweep
{
 public:
   using flag_t     = typename FlagInfo_T::flag_t;
   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   AdvectSweep(BlockDataID handlingID, BlockDataID fillFieldID, BlockDataID flagFieldID, BlockDataID pdfField,
               const FlagInfo_T& flagInfo, const PdfReconstructionModel& pdfReconstructionModel,
               bool useSimpleMassExchange, real_t cellConversionThreshold, real_t cellConversionForceThreshold)
      : handlingID_(handlingID), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID), pdfFieldID_(pdfField),
        flagInfo_(flagInfo), neighborhoodFlagFieldClone_(flagFieldID), fillFieldClone_(fillFieldID),
        pdfFieldClone_(pdfField), pdfReconstructionModel_(pdfReconstructionModel),
        useSimpleMassExchange_(useSimpleMassExchange), cellConversionThreshold_(cellConversionThreshold),
        cellConversionForceThreshold_(cellConversionForceThreshold)
   {}

   void operator()(IBlock* const block);

 protected:
   BlockDataID handlingID_;
   BlockDataID fillFieldID_;
   BlockDataID flagFieldID_;
   BlockDataID pdfFieldID_;

   FlagInfo_T flagInfo_;

   // efficient clones of fields to provide temporary fields (for writing to)
   field::FieldClone< FlagField_T, true > neighborhoodFlagFieldClone_;
   field::FieldClone< ScalarField_T, true > fillFieldClone_;
   field::FieldClone< PdfField_T, true > pdfFieldClone_;

   PdfReconstructionModel pdfReconstructionModel_;
   bool useSimpleMassExchange_;
   real_t cellConversionThreshold_;
   real_t cellConversionForceThreshold_;
}; // class AdvectSweep

template< typename LatticeModel_T, typename BoundaryHandling_T, typename FlagField_T, typename FlagInfo_T,
          typename ScalarField_T, typename VectorField_T >
void AdvectSweep< LatticeModel_T, BoundaryHandling_T, FlagField_T, FlagInfo_T, ScalarField_T,
                  VectorField_T >::operator()(IBlock* const block)
{
   // fetch data
   ScalarField_T* const fillSrcField = block->getData< ScalarField_T >(fillFieldID_);
   PdfField_T* const pdfSrcField     = block->getData< PdfField_T >(pdfFieldID_);

   FlagField_T* const flagField = block->getData< FlagField_T >(flagFieldID_);

   // temporary fields that act as destination fields
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
      pdfSrcFieldIt, pdfSrcField, fillSrcFieldIt, fillSrcField, fillDstFieldIt, fillDstField, flagFieldIt, flagField,
      neighborhoodFlagFieldIt, neighborhoodFlagField, omp critical, {
         if (flagInfo_.isInterface(flagFieldIt))
         {
            const real_t rho = lbm::getDensity(pdfSrcField->latticeModel(), pdfSrcFieldIt);

            // compute mass advection using post-collision PDFs (explicitly not PDFs updated by stream above)
            const real_t deltaMass =
               (advectMass< LatticeModel_T, FlagField_T, typename ScalarField_T::iterator,
                            typename PdfField_T::iterator, typename FlagField_T::iterator,
                            typename FlagField_T::iterator, FlagInfo_T >) (flagField, fillSrcFieldIt, pdfSrcFieldIt,
                                                                           flagFieldIt, neighborhoodFlagFieldIt,
                                                                           flagInfo_, useSimpleMassExchange_);

            // update fill level after LBM stream and mass exchange
            *fillDstFieldIt = *fillSrcFieldIt + deltaMass / rho;
         }
         else // treat non-interface cells
         {
            // manually adjust the fill level to avoid outdated fill levels being copied from fillSrcField
            if (flagInfo_.isGas(flagFieldIt)) { *fillDstFieldIt = real_c(0); }
            else
            {
               if (flagInfo_.isLiquid(flagFieldIt)) { *fillDstFieldIt = real_c(1); }
               else // flag is e.g. obstacle or outflow
               {
                  *fillDstFieldIt = *fillSrcFieldIt;
               }
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS_XYZ_OMP

   fillSrcField->swapDataPointers(fillDstField);

   BoundaryHandling_T* const handling = block->getData< BoundaryHandling_T >(handlingID_);

   // mark interface cell for conversion
   findInterfaceCellConversions< LatticeModel_T >(handling, fillSrcField, flagField, neighborhoodFlagField, flagInfo_,
                                                  cellConversionThreshold_, cellConversionForceThreshold_);
}

} // namespace free_surface
} // namespace walberla
