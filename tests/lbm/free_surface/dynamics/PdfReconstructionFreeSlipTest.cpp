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
//! \file PdfReconstructionFreeSlipTest.cpp
//! \ingroup lbm/free_surface/dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test PDF reconstruction due to a free-slip cell near the free surface boundary.
//!
//! Initialize a 3x3 grid and test reconstruction of PDFs due to a free slip boundary condition.
//!     [L][L][L]       with L: liquid cell; I: interface cell; G: gas cell; f: free-slip cell
//!     [L][I][G]            F: free-slip cell of interest (for the following explanation)
//!     [f][f][F]
//! The PDF that streams from the free-slip cell in the lower right corner (F) into the interface cell (I) must be
//! reconstructed. This is because PDFs in the free-slip cells are specularly reflected. Therefore, this free-slip
//! cell's PDF in direction (-1,1) is the same as the the gas cell's PDF in direction (-1,-1). However, since PDFs in
//! gas cells are not available, this PDF must be reconstructed.
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/FieldClone.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/dynamics/PdfReconstructionModel.h"
#include "lbm/free_surface/dynamics/functionality/ReconstructInterfaceCellABB.h"
#include "lbm/lattice_model/D2Q9.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace free_surface
{
namespace PdfReconstructionFreeSlipTest
{
using LatticeModel_T = lbm::D2Q9< lbm::collision_model::SRT, true, lbm::force_model::None, 2 >;
using Stencil        = typename LatticeModel_T::Stencil;

using Communication_T = blockforest::SimpleCommunication< LatticeModel_T::CommunicationStencil >;

using PdfField_T    = lbm::PdfField< LatticeModel_T >;
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using flag_t                        = uint32_t;
using FlagField_T                   = FlagField< flag_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

void runSimulation()
{
   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(3), uint_c(3), uint_c(1));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          false, false, false);                                 // periodicity

   // create (dummy) lattice model
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::SRT(real_c(1.8)));

   // add pdf source and destination fields
   const BlockDataID pdfSrcFieldID =
      lbm::addPdfFieldToStorage(blockForest, "PDF source field", latticeModel, uint_c(0), field::fzyx);
   const BlockDataID pdfDstFieldID =
      lbm::addPdfFieldToStorage(blockForest, "PDF destination field", latticeModel, uint_c(0), field::fzyx);

   // add (dummy) normal field
   const BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normals", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add fill level field (MUST be initialized with 1, i.e., fluid everywhere for this test; otherwise the fluid
   // flag is not detected below by initFlagsFromFillLevel())
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(1), field::fzyx, uint_c(1));

   // central interface cell, in which the reconstruction and evaluation will be performed
   const Cell centralCell = Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0));

   // construct 3x3 grid with flags according to the description at the top of this file
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField   = blockIt->getData< ScalarField_T >(fillFieldID);
      VectorField_T* const normalField = blockIt->getData< VectorField_T >(normalFieldID);

      WALBERLA_FOR_ALL_CELLS(
         fillFieldIt, fillField, normalFieldIt, normalField,

         // initialize gas cell
         if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0))) { *fillFieldIt = real_c(0); }

         // initialize interface cells
         if (fillFieldIt.cell() == centralCell) { *fillFieldIt = real_c(0.5); }

         // initialize fluid cells
         if (fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0))) {
            *fillFieldIt = real_c(1);
         }) // WALBERLA_FOR_ALL_CELLS
   }

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfSrcFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();
   freeSurfaceBoundaryHandling->setFreeSlipAtBorder(stencil::S, cell_idx_c(0));

   // initial communication
   Communication_T(blockForest, pdfSrcFieldID, fillFieldID, flagFieldID)();

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   // perform reconstruction
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const PdfField_T* const pdfSrcField    = blockIt->getData< const PdfField_T >(pdfSrcFieldID);
      PdfField_T* const pdfDstField          = blockIt->getData< PdfField_T >(pdfDstFieldID);
      const VectorField_T* const normalField = blockIt->getData< const VectorField_T >(normalFieldID);
      const FlagField_T* const flagField     = blockIt->getData< const FlagField_T >(flagFieldID);

      WALBERLA_FOR_ALL_CELLS(
         pdfSrcFieldIt, pdfSrcField, pdfDstFieldIt, pdfDstField, flagFieldIt, flagField, normalFieldIt, normalField,
         if (pdfSrcFieldIt.cell() == centralCell) {
            // reconstruct with rhoGas!=1 such that reconstructed values differ from those in pdfSrcField
            reconstructInterfaceCellLegacy< LatticeModel_T >(flagField, pdfSrcFieldIt, flagFieldIt, normalFieldIt,
                                                             flagInfo, real_c(2), pdfDstFieldIt,
                                                             PdfReconstructionModel("OnlyMissing"));
         }) // WALBERLA_FOR_ALL_CELLS
   }

   // evaluate if the correct cells were reconstructed
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const PdfField_T* const pdfSrcField = blockIt->getData< const PdfField_T >(pdfSrcFieldID);
      const PdfField_T* const pdfDstField = blockIt->getData< const PdfField_T >(pdfDstFieldID);

      WALBERLA_FOR_ALL_CELLS(
         pdfSrcFieldIt, pdfSrcField, pdfDstFieldIt, pdfDstField, if (pdfSrcFieldIt.cell() == centralCell) {
            // the boundary handling is not executed so the only change must be the PDF coming from the free-slip
            // boundary cell that reflects the PDF from the gas cell
            WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]);   // C, (0,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]);   // N, (0,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
            WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]); // W, (-1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]);   // E, (1,0,0)
            WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]); // NW, (-1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]);   // NE, (1,1,0)

            // this is the PDF that must be reconstructed due to having
            WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]); // SW, (-1,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
         })
   }

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   runSimulation();

   return EXIT_SUCCESS;
}
} // namespace PdfReconstructionFreeSlipTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::PdfReconstructionFreeSlipTest::main(argc, argv); }