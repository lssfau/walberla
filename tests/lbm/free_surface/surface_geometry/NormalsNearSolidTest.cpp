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
//! \file NormalsNearSolidTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test interface normal computation as influenced by obstacle cells, i.e., test narrower Parker-Youngs scheme.
//
//! The setup is similar to Figure 6.11 in the dissertation of S. Donath, 2011.
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/AddToStorage.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace free_surface
{
namespace NormalsNearSolidTest
{
// define types
using LatticeModel_T = lbm::D3Q27< lbm::collision_model::SRT, true >;
using Stencil_T      = LatticeModel_T::Stencil;
using flag_t         = uint32_t;

// define fields
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< flag_t >;

using Communication_T               = blockforest::SimpleCommunication< LatticeModel_T::CommunicationStencil >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(11), uint_c(3), uint_c(1));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          false, false, true);                                  // periodicity

   // create lattice model with omega=1
   LatticeModel_T latticeModel(real_c(1.0));

   // add fields
   BlockDataID pdfFieldID =
      lbm::addPdfFieldToStorage< LatticeModel_T >(blockForest, "PDF field", latticeModel, field::fzyx);
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(1.0), field::fzyx, uint_c(1));
   BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normals", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();

   // initialize domain as in Figure 6.11 in dissertation of S. Donath, 2011
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(fillFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         if (fillFieldIt.y() == cell_idx_c(0)) { *fillFieldIt = real_c(1); }
         if (fillFieldIt.y() == cell_idx_c(1)) { *fillFieldIt = real_c(fillFieldIt.x()) / real_c(25); }
         if (fillFieldIt.y() == cell_idx_c(2)) { *fillFieldIt = real_c(0); }
      }) // WALBERLA_FOR_ALL_CELLS
   }
   // set solid obstacle cells
   freeSurfaceBoundaryHandling->setNoSlipAtBorder(stencil::E, cell_idx_c(0));
   freeSurfaceBoundaryHandling->setNoSlipAtBorder(stencil::W, cell_idx_c(0));
   freeSurfaceBoundaryHandling->setNoSlipInCell(Cell(cell_idx_c(9), cell_idx_c(2), cell_idx_c(0)));

   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // initial communication (to initialize ghost layer in periodic z-direction)
   Communication_T(blockForest, pdfFieldID, fillFieldID, flagFieldID)();

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   // add sweep for computing interface normals
   NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalsSweep(
      normalFieldID, fillFieldID, flagFieldID, flagIDs::interfaceFlagID, flagIDs::liquidInterfaceGasFlagIDs,
      FreeSurfaceBoundaryHandling_T::noSlipFlagID, false, false, true, false);
   timeloop.add() << Sweep(normalsSweep, "Normals sweep");

   // perform a single time step
   timeloop.singleStep();

   // check correctness of computed interface normals; reference values have been obtained with a version of the code
   // that is assumed to be correct; results are also qualitatively verified with Figure 6.11 in dissertation of S.
   // Donath, 2011
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const FlagField_T* const flagField     = blockIt->getData< const FlagField_T >(flagFieldID);
      const VectorField_T* const normalField = blockIt->getData< const VectorField_T >(normalFieldID);
      const flag_t interfaceFlag             = flagField->getFlag(flagIDs::interfaceFlagID);

      WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, normalFieldIt, normalField, {
         if (isFlagSet(flagFieldIt, interfaceFlag))
         {
            // regular Parker-Youngs normal computation
            if (flagFieldIt.x() >= cell_idx_c(2) && flagFieldIt.x() <= cell_idx_c(7))
            {
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[0], real_c(-0.0399680383488715887), real_c(1e-6));
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[1], real_c(0.999200958721789267), real_c(1e-6));
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[2], real_c(0), real_c(1e-6));
            }

            // modified, i.e., narrower Parker-Youngs normal computation near solid boundaries
            if (flagFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[0], real_c(-0.0199960011996001309), real_c(1e-6));
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[1], real_c(0.999800059980007094), real_c(1e-6));
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[2], real_c(0), real_c(1e-6));
            }
            if (flagFieldIt.cell() == Cell(cell_idx_c(8), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[0], real_c(-0.129339184067768065), real_c(1e-6));
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[1], real_c(0.991600411186221775), real_c(1e-6));
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[2], real_c(-2.99157e-17), real_c(1e-6));
            }
            if (flagFieldIt.cell() == Cell(cell_idx_c(9), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[0], real_c(-0.0461047666084008420), real_c(1e-6));
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[1], real_c(0.998936609848685486), real_c(1e-6));
               WALBERLA_CHECK_FLOAT_EQUAL((*normalFieldIt)[2], real_c(-8.5311e-17), real_c(1e-6));
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   return EXIT_SUCCESS;
}
} // namespace NormalsNearSolidTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::NormalsNearSolidTest::main(argc, argv); }