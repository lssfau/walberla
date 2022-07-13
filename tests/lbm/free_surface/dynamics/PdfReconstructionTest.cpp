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
//! \file PdfReconstructionTest.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test PDF reconstruction at free surface boundary.
//!
//! Initialize 3x3 grid similar to figure 3 in publication of Koerner et al., 2005 and test reconstruction of PDFs at
//! the free surface boundary with respect to the models specified in PdfReconstructionModel.h.
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
namespace PdfReconstructionTest
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

void runSimulation(const PdfReconstructionModel& pdfReconstructionModel)
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

   // add normal field
   const BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add fill level field (MUST be initialized with 1, i.e., fluid everywhere for this test; otherwise the fluid
   // flag is not detected below by initFlagsFromFillLevel())
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(1), field::fzyx, uint_c(1));

   // central interface cell, in which the reconstruction and evaluation will be performed
   const Cell centralCell = Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0));

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField   = blockIt->getData< ScalarField_T >(fillFieldID);
      VectorField_T* const normalField = blockIt->getData< VectorField_T >(normalFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, normalFieldIt, normalField, {
         // initialize gas cells as in figure 3 from Koerner et al., 2005
         if (fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
         {
            *fillFieldIt = real_c(0);
         }

         // initialize interface cells as in figure 3 from Koerner et al., 2005
         if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)) ||
             fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
         {
            *fillFieldIt = real_c(0.5);
         }

         // initialize fluid cell as in figure 3 from Koerner et al., 2005
         if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0))) { *fillFieldIt = real_c(1); }

         // initialize interface normal as in figure 3 from Koerner et al., 2005 (values estimated)
         // IMPORTANT REMARK: In this waLBerla's free surface implementation, the normal is defined to point from fluid
         // to gas, whereas in Koerner et al., the normal is defined to point from gas to fluid. Therefore, we
         // initialize the normal in opposite direction than in figure 3.
         if (normalFieldIt.cell() == centralCell)
         {
            *normalFieldIt = Vector3< real_t >(real_c(-0.83867), real_c(-0.54464), real_c(0));
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfSrcFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

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
         pdfSrcFieldIt, pdfSrcField, pdfDstFieldIt, pdfDstField, flagFieldIt, flagField, normalFieldIt, normalField, {
            if (pdfSrcFieldIt.cell() == centralCell)
            {
               // reconstruct with rhoGas!=1 such that reconstructed values differ from those in pdfSrcField
               reconstructInterfaceCellLegacy< LatticeModel_T >(flagField, pdfSrcFieldIt, flagFieldIt, normalFieldIt,
                                                                flagInfo, real_c(2), pdfDstFieldIt,
                                                                pdfReconstructionModel);
            }
         }) // WALBERLA_FOR_ALL_CELLS
   }

   const PdfReconstructionModel::ReconstructionModel reconstructionModel = pdfReconstructionModel.getModelType();
   const uint_t minReconstruct                               = pdfReconstructionModel.getNumMinReconstruct();
   const PdfReconstructionModel::FallbackModel fallbackModel = pdfReconstructionModel.getFallbackModel();

   // evaluate if the correct cells were reconstructed:
   // 1. equality/ inequality between pdfSrc and pdfDst was verified manually, i.e., by hand
   // 2. comparison with expected (reconstructed) values:
   //    - results only valid for rhoGas=2 (as specified above)
   //    - results obtained with a version that is believed to be correct
   //    - was included for detection of changes in PDF reconstruction boundary condition
   //    - makes 1. actually obsolete (1. was kept for as this is what the test was actually intended for, i.e., find
   //       out whether the correct PDFs are reconstructed)
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const PdfField_T* const pdfSrcField = blockIt->getData< const PdfField_T >(pdfSrcFieldID);
      const PdfField_T* const pdfDstField = blockIt->getData< const PdfField_T >(pdfDstFieldID);

      WALBERLA_FOR_ALL_CELLS(pdfSrcFieldIt, pdfSrcField, pdfDstFieldIt, pdfDstField, {
         if (pdfSrcFieldIt.cell() == centralCell)
         {
            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::NormalBasedKeepCenter ||
                (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                 minReconstruct > uint_c(3) &&
                 fallbackModel == PdfReconstructionModel::FallbackModel::NormalBasedKeepCenter))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(0.444444), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.111111), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.111111), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0277778), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]);   // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]);   // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]);   // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::NormalBasedReconstructCenter)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(1.333333), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.111111), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.111111), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0277778), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]);   // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]);   // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::All)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(1.333333), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.333333), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.333333), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0833333), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0833333), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            // in the setup here, 3 PDFs will always be reconstructed as they are coming from the gas-side and are
            // therefore missing
            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissing ||
                (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                 (minReconstruct == uint_c(0) || minReconstruct == uint_c(1) || minReconstruct == uint_c(2) ||
                  minReconstruct == uint_c(3)) &&
                 (fallbackModel == PdfReconstructionModel::FallbackModel::Smallest ||
                  fallbackModel == PdfReconstructionModel::FallbackModel::Largest ||
                  fallbackModel == PdfReconstructionModel::FallbackModel::NormalBasedKeepCenter)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(0.444444), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.111111), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.111111), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0277778), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0277778), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]);   // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]);   // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]);   // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]);   // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(4) && fallbackModel == PdfReconstructionModel::FallbackModel::Smallest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(1.333333), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.111111), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.111111), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0277778), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0277778), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]);   // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]);   // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]);   // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(4) && fallbackModel == PdfReconstructionModel::FallbackModel::Largest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(0.444444), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.111111), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.111111), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0833333), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0277778), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]);   // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]);   // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]);   // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(5) && fallbackModel == PdfReconstructionModel::FallbackModel::Smallest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(1.333333), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.111111), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.111111), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0277778), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]);   // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]);   // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(5) && fallbackModel == PdfReconstructionModel::FallbackModel::Largest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(0.444444), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.111111), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.333333), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0833333), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0277778), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]);   // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]);   // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(6) && fallbackModel == PdfReconstructionModel::FallbackModel::Smallest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(1.333333), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.111111), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.111111), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0833333), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0277778), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]);   // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]);   // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]);   // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(6) && fallbackModel == PdfReconstructionModel::FallbackModel::Largest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(0.444444), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.333333), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.333333), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0833333), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0277778), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]);   // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]);   // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(7) && fallbackModel == PdfReconstructionModel::FallbackModel::Smallest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(1.333333), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.333333), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.111111), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0833333), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0277778), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]);   // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]);   // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(7) && fallbackModel == PdfReconstructionModel::FallbackModel::Largest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(0.444444), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.333333), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.333333), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0277778), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0833333), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]);   // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]);   // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(8) && fallbackModel == PdfReconstructionModel::FallbackModel::Smallest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(1.333333), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.333333), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.333333), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0833333), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0277778), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]);   // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(8) && fallbackModel == PdfReconstructionModel::FallbackModel::Largest)
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(0.444444), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.333333), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.333333), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0833333), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0833333), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_EQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]);   // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }

            if (reconstructionModel == PdfReconstructionModel::ReconstructionModel::OnlyMissingMin &&
                minReconstruct == uint_c(9) &&
                (fallbackModel == PdfReconstructionModel::FallbackModel::Smallest ||
                 fallbackModel == PdfReconstructionModel::FallbackModel::Largest))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[0], real_c(1.333333), real_c(1e-6));  // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[1], real_c(0.333333), real_c(1e-6));  // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[2], real_c(0.333333), real_c(1e-6));  // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[3], real_c(0.333333), real_c(1e-6));  // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[4], real_c(0.333333), real_c(1e-6));  // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[5], real_c(0.0833333), real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[6], real_c(0.0833333), real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[7], real_c(0.0833333), real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfDstFieldIt[8], real_c(0.0833333), real_c(1e-6)); // SE, (1,-1,0)

               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[0], pdfDstFieldIt[0]); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[1], pdfDstFieldIt[1]); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[2], pdfDstFieldIt[2]); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[3], pdfDstFieldIt[3]); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[4], pdfDstFieldIt[4]); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[5], pdfDstFieldIt[5]); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[6], pdfDstFieldIt[6]); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[7], pdfDstFieldIt[7]); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_UNEQUAL(pdfSrcFieldIt[8], pdfDstFieldIt[8]); // SE, (1,-1,0)
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   PdfReconstructionModel model = PdfReconstructionModel("NormalBasedKeepCenter");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = PdfReconstructionModel("NormalBasedReconstructCenter");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = PdfReconstructionModel("All");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = PdfReconstructionModel("OnlyMissing");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   for (uint_t i = uint_c(0); i != uint_c(10); ++i)
   {
      model = PdfReconstructionModel("OnlyMissingMin-" + std::to_string(i) + "-smallest");
      WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
      runSimulation(model);

      model = PdfReconstructionModel("OnlyMissingMin-" + std::to_string(i) + "-largest");
      WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
      runSimulation(model);

      model = PdfReconstructionModel("OnlyMissingMin-" + std::to_string(i) + "-normalBasedKeepCenter");
      WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
      runSimulation(model);
   }

   return EXIT_SUCCESS;
}
} // namespace PdfReconstructionTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::PdfReconstructionTest::main(argc, argv); }