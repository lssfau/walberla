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
//! \file StreamInCellIntervalCodegenTest.cpp
//! \ingroup lbm/codegen
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test generated streamInCellInterval() by comparing with regular generated stream().
//
//! Initialize a periodic 2x1x1 PDF field with arbitrary values and stream the PDFs with both functions.
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "timeloop/SweepTimeloop.h"

#include "GeneratedLatticeModel.h"

namespace walberla
{
namespace StreamInCellIntervalCodegenTest
{
using LatticeModel_T         = lbm::GeneratedLatticeModel;
using Stencil_T              = LatticeModel_T::Stencil;
using PdfField_T             = lbm::PdfField< LatticeModel_T >;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   // define the domain size (2x1x1)
   const Vector3< uint_t > numBlocks(uint_t(1));
   const Vector3< uint_t > domainSize(uint_t(2), uint_t(1), uint_t(1));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // ease access to both cells
   const Cell cell0 = Cell(cell_idx_t(0), cell_idx_t(0), cell_idx_t(0));
   const Cell cell1 = Cell(cell_idx_t(1), cell_idx_t(0), cell_idx_t(0));

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, true, true);                                    // periodicity

   // relaxation rate
   real_t omega = real_c(1.8);

   // create lattice model
   LatticeModel_T latticeModel = LatticeModel_T(omega);

   // add PDF field (source for both functions, destination for regular stream)
   const BlockDataID pdfFieldID =
      lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, uint_c(1), field::fzyx);

   // add PDF field as destination for streamInInterval()
   const BlockDataID pdfFieldStreamIntervalID = lbm::addPdfFieldToStorage(
      blockForest, "PDF destination field for streamInInterval", latticeModel, uint_c(1), field::fzyx);

   // initialize PDF field with arbitrary values (should not be the same in both cells to improve validity of test)
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      PdfField_T* pdfField = blockIt->getData< PdfField_T >(pdfFieldID);

      pdfField->setDensityAndVelocity(cell0, Vector3< real_t >(real_c(0.1), real_c(0.1), real_c(0.1)), real_c(1.1));
      pdfField->setDensityAndVelocity(cell1, Vector3< real_t >(real_c(-0.1), real_c(-0.1), real_c(-0.1)), real_c(0.9));
   }

   // create communication for PDF fields
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication(blockForest);
   communication.addPackInfo(make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >(pdfFieldID));

   // communicate PDF fields to fill ghost layers
   communication();

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      PdfField_T* pdfField               = blockIt->getData< PdfField_T >(pdfFieldID);
      PdfField_T* pdfFieldStreamInterval = blockIt->getData< PdfField_T >(pdfFieldStreamIntervalID);

      // use streamInCellInterval() (does not change content of pdfField)
      auto lbmSweepGenerated = typename LatticeModel_T::Sweep(pdfFieldID);
      const CellInterval ci  = CellInterval(cell0, cell1);
      lbmSweepGenerated.streamInCellInterval(pdfField, pdfFieldStreamInterval, ci);

      // use regular stream(); overwrites content of pdfField and MUST therefore be called after streamInCellInterval()
      lbmSweepGenerated.stream(blockIt.get());
   }

   // check equality of streamInCellInterval() and stream()
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      PdfField_T* pdfField               = blockIt->getData< PdfField_T >(pdfFieldID);
      PdfField_T* pdfFieldStreamInterval = blockIt->getData< PdfField_T >(pdfFieldStreamIntervalID);

      WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, pdfFieldStreamIntervalIt, pdfFieldStreamInterval, {
         // check equality of each PDF
         for (uint_t i = uint_t(0); i != pdfField->F_SIZE; ++i)
         {
            WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[i], pdfFieldStreamIntervalIt[i]);
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   return EXIT_SUCCESS;
}
} // namespace StreamInCellIntervalCodegenTest
} // namespace walberla

int main(int argc, char** argv) { return walberla::StreamInCellIntervalCodegenTest::main(argc, argv); }
