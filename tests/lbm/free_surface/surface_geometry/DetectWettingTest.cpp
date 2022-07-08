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
//! \file DetectWettingTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test the DetectWettingSweep using a single interface cell surrounded by gas cells and obstacle cells.
//!
//! 3x4x3 domain where the first and fourth layer in y-direction are solid cells. The cell at (1,1,1) is an interface
//! cell with given normal and fill level. All remaining cells are gas cells that might be marked for conversion to
//! wetting cells by DetectWettingSweep.
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/surface_geometry/DetectWettingSweep.h"
#include "lbm/lattice_model/D3Q19.h"

#include "timeloop/SweepTimeloop.h"

#include <algorithm>
#include <vector>

namespace walberla
{
namespace free_surface
{
namespace DetectWettingTest
{
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::SRT >;
using PdfField_T     = lbm::PdfField< LatticeModel_T >;
using Stencil_T      = LatticeModel_T::Stencil;

using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using flag_t                        = uint32_t;
using FlagField_T                   = FlagField< flag_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

std::vector< Cell > detectWettingCells(const Vector3< real_t >& normal, const real_t fillLevel)
{
   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(3), uint_c(4), uint_c(3));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          false, false, false);                                 // periodicity

   // create lattice model (dummy, not relevant for this test)
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::SRT(real_c(1)));

   // add fields
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(0.0), field::fzyx, uint_c(1));
   const BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normal field", Vector3< real_t >(real_c(0.0)), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const BlockDataID boundaryHandlingID                               = freeSurfaceBoundaryHandling->getHandlingID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // initialize domain
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField   = blockIt->getData< ScalarField_T >(fillFieldID);
      VectorField_T* const normalField = blockIt->getData< VectorField_T >(normalFieldID);
      FreeSurfaceBoundaryHandling_T::BoundaryHandling_T* const boundaryHandling =
         blockIt->getData< FreeSurfaceBoundaryHandling_T::BoundaryHandling_T >(boundaryHandlingID);

      WALBERLA_FOR_ALL_CELLS_XYZ(fillField, {
         // set boundary cells at bottom and top of domain
         if (y == 0 || y == 3) { boundaryHandling->setBoundary(FreeSurfaceBoundaryHandling_T::noSlipFlagID, x, y, z); }

         // set interface cell in the center of the domain
         if (x == 1 && y == 1 && z == 1)
         {
            boundaryHandling->setFlag(flagInfo.interfaceFlag, x, y, z);
            fillField->get(x, y, z)   = fillLevel;
            normalField->get(x, y, z) = normal.getNormalized();
         }

         // set remaining domain to gas
         if ((y == 1 || y == 2) && !(x == 1 && y == 1 && z == 1))
         {
            boundaryHandling->setFlag(flagInfo.gasFlag, x, y, z);
         }
      }) // WALBERLA_FOR_ALL_CELLS_XYZ
   }

   // create timeloop
   SweepTimeloop timeloop(blockForest, 1);

   // add DetectWettingSweep
   DetectWettingSweep< Stencil_T, FreeSurfaceBoundaryHandling_T::BoundaryHandling_T, FlagField_T, ScalarField_T,
                       VectorField_T >
      detWetSweep(boundaryHandlingID, flagInfo, normalFieldID, fillFieldID);
   timeloop.add() << Sweep(detWetSweep, "Detect wetting sweep");

   timeloop.singleStep();

   std::vector< Cell > markedCells;

   // get cells that were marked by DetectWettingSweep
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const FlagField_T* const flagField = blockIt->getData< const FlagField_T >(flagFieldID);

      WALBERLA_FOR_ALL_CELLS_OMP(
         flagFieldIt, flagField, omp critical, if (flagInfo.isKeepInterfaceForWetting(flagFieldIt)) {
            markedCells.emplace_back(flagFieldIt.cell());
         }) // WALBERLA_FOR_ALL_CELLS_OMP
   }

   return markedCells;
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   Vector3< real_t > normal;
   real_t fillLevel;
   std::vector< Cell > expectedWettingCells;
   std::vector< Cell > computedWettingCells;
   bool vectorsEqual;

   // test various different normals and fill levels; expected results have been determined using ParaView
   normal               = Vector3< real_t >(real_c(1), real_c(1), real_c(1));
   fillLevel            = real_c(0.01);
   expectedWettingCells = std::vector< Cell >{ Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)) };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(-1), real_c(-1), real_c(-1));
   fillLevel            = real_c(0.01);
   expectedWettingCells = std::vector< Cell >{
      Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)), Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)),
      Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(1)), Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(2)),
      Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(2)), Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(2))
   };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(1), real_c(1), real_c(0));
   fillLevel            = real_c(0.01);
   expectedWettingCells = std::vector< Cell >{ Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)),
                                               Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(2)),
                                               Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(2)) };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(0), real_c(1), real_c(1));
   fillLevel            = real_c(0.01);
   expectedWettingCells = std::vector< Cell >{ Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)),
                                               Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)) };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(1), real_c(0), real_c(1));
   fillLevel            = real_c(0.01);
   expectedWettingCells = std::vector< Cell >{ Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)),
                                               Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)),
                                               Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(1)),
                                               Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)) };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(-1), real_c(0), real_c(-1));
   fillLevel            = real_c(0.01);
   expectedWettingCells = std::vector< Cell >{ Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)),
                                               Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)),
                                               Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(1)),
                                               Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(2)),
                                               Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(2)) };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(-1), real_c(0), real_c(-1));
   fillLevel            = real_c(0.5);
   expectedWettingCells = std::vector< Cell >{ Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),
                                               Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)),
                                               Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)),
                                               Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(1)),
                                               Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)) };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(1), real_c(1), real_c(1));
   fillLevel            = real_c(0.9);
   expectedWettingCells = std::vector< Cell >{
      Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)), Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)),
      Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(1)), Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(2)),
      Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(2)), Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(2))
   };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(1), real_c(1), real_c(0));
   fillLevel            = real_c(0.9);
   expectedWettingCells = std::vector< Cell >{
      Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)),
      Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)), Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(2)),
      Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(2)), Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(2))
   };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   normal               = Vector3< real_t >(real_c(0), real_c(1), real_c(0));
   fillLevel            = real_c(0.9);
   expectedWettingCells = std::vector< Cell >{
      Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)),
      Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)), Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(2)),
      Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(2)), Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(2))
   };
   WALBERLA_LOG_INFO("Testing wetting cells with normal=" << normal << " and fill level=" << fillLevel);
   computedWettingCells = detectWettingCells(normal, fillLevel);
   vectorsEqual         = std::is_permutation(computedWettingCells.begin(), computedWettingCells.end(),
                                              expectedWettingCells.begin(), expectedWettingCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   return EXIT_SUCCESS;
}
} // namespace DetectWettingTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::DetectWettingTest::main(argc, argv); }