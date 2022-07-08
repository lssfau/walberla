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
//! \file ObstacleFillLevelTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test the ObstacleFillLevelSweep with an obstacle cell that is surrounded by different other cells.
//!
//! 3x3x1 domain, with obstacle cell at (1,1,0) and given obstacle normal.
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/surface_geometry/ObstacleFillLevelSweep.h"
#include "lbm/lattice_model/D3Q19.h"

#include "timeloop/SweepTimeloop.h"

#include <unordered_map>
#include <unordered_set>

namespace walberla
{
namespace free_surface
{
namespace ObstacleFillLevelTest
{
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::SRT >;
using PdfField_T     = lbm::PdfField< LatticeModel_T >;
using Stencil_T      = LatticeModel_T::Stencil;

using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using flag_t                        = uint32_t;
using FlagField_T                   = FlagField< flag_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

real_t computeObstacleFillLevel(const std::unordered_map< Cell, real_t >& interfaceCells,
                                const std::unordered_set< Cell >& liquidCells,
                                const std::unordered_set< Cell >& gasCells,
                                const std::unordered_set< Cell >& solidCells,
                                const Vector3< real_t >& centerObstacleNormal)
{
   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(3), uint_c(3), uint_c(3));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, true, true);                                    // periodicity

   // create lattice model (dummy, not relevant for this test)
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::SRT(real_c(1)));

   // add fields
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(0.0), field::fzyx, uint_c(1));
   BlockDataID obstaclefillFieldID = field::addToStorage< ScalarField_T >(blockForest, "Obstacle fill level field",
                                                                          real_c(0.0), field::fzyx, uint_c(1));
   const BlockDataID obstacleNormalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Obstacle normals", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const BlockDataID boundaryHandlingID                               = freeSurfaceBoundaryHandling->getHandlingID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // initialize domain
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField           = blockIt->getData< ScalarField_T >(fillFieldID);
      VectorField_T* const obstacleNormalField = blockIt->getData< VectorField_T >(obstacleNormalFieldID);
      FreeSurfaceBoundaryHandling_T::BoundaryHandling_T* const boundaryHandling =
         blockIt->getData< FreeSurfaceBoundaryHandling_T::BoundaryHandling_T >(boundaryHandlingID);

      WALBERLA_FOR_ALL_CELLS_XYZ(fillField, {
         // set no slip at domain boundaries in z-direction to reduce problem to 2D
         if (z == cell_idx_c(0) || z == cell_idx_c(2))
         {
            boundaryHandling->setBoundary(FreeSurfaceBoundaryHandling_T::noSlipFlagID, x, y, z);

            // obstacle normal must be normalized to not trigger the assertion in ObstacleFillLevelSweep
            obstacleNormalField->get(x, y, z) = Vector3< real_t >(real_c(1), real_c(1), real_c(1)).getNormalized();
            continue;
         }

         // obstacle cell (to be evaluated) in the domain center
         if (x == cell_idx_c(1) && y == cell_idx_c(1) && z == cell_idx_c(1))
         {
            boundaryHandling->setBoundary(FreeSurfaceBoundaryHandling_T::noSlipFlagID, x, y, z);
            obstacleNormalField->get(x, y, z) = centerObstacleNormal.getNormalized();
            fillField->get(x, y, z)           = real_c(-5); // dummy value to check that the value did not change
            continue;
         }

         if (interfaceCells.find(Cell(x, y, z)) != interfaceCells.end())
         {
            boundaryHandling->setFlag(flagIDs::interfaceFlagID, x, y, z);
            fillField->get(x, y, z) = interfaceCells.find(Cell(x, y, z))->second;
            continue;
         }

         if (liquidCells.find(Cell(x, y, z)) != liquidCells.end())
         {
            boundaryHandling->setFlag(flagIDs::liquidFlagID, x, y, z);
            fillField->get(x, y, z) = real_c(1);
            continue;
         }

         if (gasCells.find(Cell(x, y, z)) != gasCells.end())
         {
            boundaryHandling->setFlag(flagIDs::gasFlagID, x, y, z);
            fillField->get(x, y, z) = real_c(0);
            continue;
         }

         if (solidCells.find(Cell(x, y, z)) != solidCells.end())
         {
            boundaryHandling->setFlag(FreeSurfaceBoundaryHandling_T::noSlipFlagID, x, y, z);

            // obstacle normal must be normalized to not trigger the assertion in ObstacleFillLevelSweep
            obstacleNormalField->get(x, y, z) = Vector3< real_t >(real_c(1), real_c(1), real_c(1)).getNormalized();
            continue;
         }
      }) // WALBERLA_FOR_ALL_CELLS_XYZ
   }

   // create timeloop
   SweepTimeloop timeloop(blockForest, 1);

   // add ObstacleFillLevelSweep
   ObstacleFillLevelSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > obstFillSweep(
      obstaclefillFieldID, fillFieldID, flagFieldID, obstacleNormalFieldID, flagIDs::liquidInterfaceGasFlagIDs,
      flagInfo.getObstacleIDSet());
   timeloop.add() << Sweep(obstFillSweep, "Obstacle fill level sweep");

   timeloop.singleStep();

   real_t obstacleFillLevel = real_c(0);

   // get fill level of (central) solid cell
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const ScalarField_T* const obstaclefillField = blockIt->getData< const ScalarField_T >(obstaclefillFieldID);

      obstacleFillLevel = obstaclefillField->get(cell_idx_c(1), cell_idx_c(1), cell_idx_c(1));
   }

   return obstacleFillLevel;
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   Vector3< real_t > centerObstacleNormal;
   std::unordered_map< Cell, real_t > interfaceCells;
   std::unordered_set< Cell > liquidCells;
   std::unordered_set< Cell > gasCells;
   std::unordered_set< Cell > solidCells;
   real_t expectedFillLevel;
   real_t computedFillLevel;

   // IMPORTANT REMARK:
   // the fill level of the obstacle center cell at (1, 1, 1) is going to be computed; therefore, Cell(1, 1, 1) must not
   // be set here

   // test case 1: interface neighbors at the top, no liquid or gas neighbors, remaining cells are solid, obstacle
   // normal points upwards
   WALBERLA_LOG_RESULT("Performing test case 1")
   centerObstacleNormal = Vector3< real_t >(real_c(0), real_c(1), real_c(0));
   solidCells.emplace(Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(1)));
   solidCells.emplace(Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(1)));
   solidCells.emplace(Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(1)));
   solidCells.emplace(Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)));
   solidCells.emplace(Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)));
   interfaceCells.emplace(Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   expectedFillLevel = real_c(0.5); // verified by hand
   computedFillLevel =
      computeObstacleFillLevel(interfaceCells, liquidCells, gasCells, solidCells, centerObstacleNormal);
   WALBERLA_CHECK_FLOAT_EQUAL(expectedFillLevel, computedFillLevel)

   interfaceCells.clear();
   liquidCells.clear();
   gasCells.clear();
   solidCells.clear();

   // test case 2: interface neighbors at the top, left cell is liquid, right cell is gas, bottom cells are solid,
   // obstacle normal points upwards
   WALBERLA_LOG_RESULT("Performing test case 2")
   centerObstacleNormal = Vector3< real_t >(real_c(0), real_c(1), real_c(0));
   solidCells.emplace(Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(1)));
   solidCells.emplace(Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(1)));
   solidCells.emplace(Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(1)));
   interfaceCells.emplace(Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(1)), real_c(0.75));
   interfaceCells.emplace(Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   liquidCells.emplace(Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)));
   gasCells.emplace(Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)));
   expectedFillLevel = real_c(0.5732233); // verified by hand
   computedFillLevel =
      computeObstacleFillLevel(interfaceCells, liquidCells, gasCells, solidCells, centerObstacleNormal);
   WALBERLA_CHECK_FLOAT_EQUAL(expectedFillLevel, computedFillLevel)

   interfaceCells.clear();
   liquidCells.clear();
   gasCells.clear();
   solidCells.clear();

   // test case 3: same as testcase 2 but obstacle normal points to the upper left corner
   WALBERLA_LOG_RESULT("Performing test case 3")
   centerObstacleNormal = Vector3< real_t >(real_c(-1), real_c(1), real_c(0));
   solidCells.emplace(Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(1)));
   solidCells.emplace(Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(1)));
   solidCells.emplace(Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(1)));
   interfaceCells.emplace(Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(1)), real_c(0.75));
   interfaceCells.emplace(Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   liquidCells.emplace(Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)));
   gasCells.emplace(Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)));
   expectedFillLevel = real_c(0.58009431); // verified by hand
   computedFillLevel =
      computeObstacleFillLevel(interfaceCells, liquidCells, gasCells, solidCells, centerObstacleNormal);
   WALBERLA_CHECK_FLOAT_EQUAL(expectedFillLevel, computedFillLevel)

   interfaceCells.clear();
   liquidCells.clear();
   gasCells.clear();
   solidCells.clear();

   // test case 4: only interface cells in neighborhood (all with same fill level)
   WALBERLA_LOG_RESULT("Performing test case 4")
   centerObstacleNormal = Vector3< real_t >(real_c(-1), real_c(1), real_c(0));
   interfaceCells.emplace(Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   interfaceCells.emplace(Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(1)), real_c(0.5));
   expectedFillLevel = real_c(0.5); // dummy value as initialized in code above
   computedFillLevel =
      computeObstacleFillLevel(interfaceCells, liquidCells, gasCells, solidCells, centerObstacleNormal);
   WALBERLA_CHECK_FLOAT_EQUAL(expectedFillLevel, computedFillLevel)

   return EXIT_SUCCESS;
}
} // namespace ObstacleFillLevelTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::ObstacleFillLevelTest::main(argc, argv); }