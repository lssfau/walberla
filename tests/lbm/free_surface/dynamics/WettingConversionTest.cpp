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
//! \file WettingConversionTest.cpp
//! \ingroup lbm/free_surface/dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test cell conversion initiated by wetting.
//!
//! Initialize drop as a cylinder section near a solid wall. Run a free surface LBM simulation with local triangulation
//! and wetting, and evaluate the converted interface cells for correctness.
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/lattice_model/D3Q19.h"

#include <algorithm>
#include <vector>
namespace walberla
{
namespace free_surface
{
namespace WettingConversionTest
{
// define types
using flag_t      = uint32_t;
using FlagField_T = FlagField< flag_t >;

using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

template< typename LatticeModel_T >
std::vector< Cell > runSimulation(uint_t timesteps, real_t contactAngle)
{
   using Stencil_T                     = typename LatticeModel_T::Stencil;
   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;
   using Communication_T = blockforest::SimpleCommunication< typename LatticeModel_T::CommunicationStencil >;

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(14), uint_c(7), uint_c(1));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          false, false, true);                                  // periodicity

   real_t relaxRate = real_c(1.8);

   // create lattice model
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::SRT(relaxRate));

   // add pdf field
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);

   // add fill level field (initialized with 0, i.e., gas everywhere)
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(0.0), field::fzyx, uint_c(1));

   // add dummy force field
   BlockDataID forceDensityFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Force density field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // add liquid drop (as cylinder section)
   Vector3< real_t > midpoint1(real_c(14) * real_c(0.5), real_c(1), real_c(0));
   Vector3< real_t > midpoint2(real_c(14) * real_c(0.5), real_c(1), real_c(5));
   geometry::Cylinder cylinder(midpoint1, midpoint2, real_c(14) * real_c(0.2));
   bubble_model::addBodyToFillLevelField< geometry::Cylinder >(*blockForest, fillFieldID, cylinder, false);

   // initialize bottom (in y-direction) of domain as no-slip boundary
   freeSurfaceBoundaryHandling->setNoSlipAtBorder(stencil::S, cell_idx_c(0));
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // initial communication
   Communication_T(blockForest, pdfFieldID, fillFieldID, flagFieldID, forceDensityFieldID)();

   // add (dummy) bubble model
   const bool disableSplits = true; // necessary if a gas bubble could split
   auto bubbleModel         = std::make_shared< bubble_model::BubbleModel< Stencil_T > >(blockForest, disableSplits);
   bubbleModel->initFromFillLevelField(fillFieldID);
   bubbleModel->setAtmosphere(Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)), real_c(1));

   // create timeloop
   SweepTimeloop timeloop(blockForest, timesteps);

   // add surface geometry handler
   SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > geometryHandler(
      blockForest, freeSurfaceBoundaryHandling, fillFieldID, "LocalTriangulation", true, true, contactAngle);
   geometryHandler.addSweeps(timeloop);

   const ConstBlockDataID curvatureFieldID = geometryHandler.getConstCurvatureFieldID();
   const ConstBlockDataID normalFieldID    = geometryHandler.getConstNormalFieldID();

   // add surface dynamics handler
   SurfaceDynamicsHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > dynamicsHandler(
      blockForest, pdfFieldID, flagFieldID, fillFieldID, forceDensityFieldID, normalFieldID, curvatureFieldID,
      freeSurfaceBoundaryHandling, bubbleModel, "NormalBasedKeepCenter", "EquilibriumRefilling", "EvenlyNewInterface",
      relaxRate, Vector3< real_t >(real_c(0)), real_c(1e-2), false, real_c(1e-3), real_c(1e-1));
   dynamicsHandler.addSweeps(timeloop);

   timeloop.run();

   // get interface cells after performing simulation
   std::vector< Cell > interfaceCells;
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const FlagField_T* const flagField = blockIt->getData< const FlagField_T >(flagFieldID);

      WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, omp critical, {
         if (flagInfo.isInterface(flagFieldIt)) { interfaceCells.emplace_back(flagFieldIt.cell()); }
      }) // WALBERLA_FOR_ALL_CELLS_OMP
   }

   return interfaceCells;
}

template< typename LatticeModel_T >
void testWettingConversion()
{
   uint_t timesteps;
   real_t contactAngle;
   std::vector< Cell > expectedInterfaceCells;
   std::vector< Cell > computedInterfaceCells;
   bool vectorsEqual;

   // test different contact angles; expected results have been determined using a version of the code that was assumed
   // to be correct
   timesteps              = uint_c(200);
   contactAngle           = real_c(120);
   expectedInterfaceCells = std::vector< Cell >{
      Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(9), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(9), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(4), cell_idx_c(3), cell_idx_c(0)), Cell(cell_idx_c(5), cell_idx_c(3), cell_idx_c(0)),
      Cell(cell_idx_c(6), cell_idx_c(3), cell_idx_c(0)), Cell(cell_idx_c(7), cell_idx_c(3), cell_idx_c(0)),
      Cell(cell_idx_c(8), cell_idx_c(3), cell_idx_c(0)), Cell(cell_idx_c(9), cell_idx_c(3), cell_idx_c(0))
   };
   WALBERLA_LOG_INFO("Testing interface conversion with wetting cells, contact angle=" << contactAngle << " degrees");
   computedInterfaceCells = runSimulation< LatticeModel_T >(timesteps, contactAngle);
   vectorsEqual           = std::is_permutation(computedInterfaceCells.begin(), computedInterfaceCells.end(),
                                                expectedInterfaceCells.begin(), expectedInterfaceCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }
   MPIManager::instance()->resetMPI();

   timesteps              = uint_c(200);
   contactAngle           = real_c(80);
   expectedInterfaceCells = std::vector< Cell >{
      Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(9), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(10), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(8), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(9), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(5), cell_idx_c(3), cell_idx_c(0)), Cell(cell_idx_c(6), cell_idx_c(3), cell_idx_c(0)),
      Cell(cell_idx_c(7), cell_idx_c(3), cell_idx_c(0)), Cell(cell_idx_c(8), cell_idx_c(3), cell_idx_c(0))
   };
   WALBERLA_LOG_INFO("Testing interface conversion with wetting cells, contact angle=" << contactAngle << " degrees");
   computedInterfaceCells = runSimulation< LatticeModel_T >(timesteps, contactAngle);
   vectorsEqual           = std::is_permutation(computedInterfaceCells.begin(), computedInterfaceCells.end(),
                                                expectedInterfaceCells.begin(), expectedInterfaceCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }
   MPIManager::instance()->resetMPI();

   timesteps              = uint_c(200);
   contactAngle           = real_c(45);
   expectedInterfaceCells = std::vector< Cell >{
      Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)),  Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(10), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(11), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)),  Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)),  Cell(cell_idx_c(6), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(7), cell_idx_c(2), cell_idx_c(0)),  Cell(cell_idx_c(8), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(9), cell_idx_c(2), cell_idx_c(0)),  Cell(cell_idx_c(10), cell_idx_c(2), cell_idx_c(0))
   };
   WALBERLA_LOG_INFO("Testing interface conversion with wetting cells, contact angle=" << contactAngle << " degrees");
   computedInterfaceCells = runSimulation< LatticeModel_T >(timesteps, contactAngle);
   vectorsEqual           = std::is_permutation(computedInterfaceCells.begin(), computedInterfaceCells.end(),
                                                expectedInterfaceCells.begin(), expectedInterfaceCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }
   MPIManager::instance()->resetMPI();

   timesteps              = uint_c(500);
   contactAngle           = real_c(1);
   expectedInterfaceCells = std::vector< Cell >{
      Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),  Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)),  Cell(cell_idx_c(10), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(11), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(12), cell_idx_c(1), cell_idx_c(0)),
      Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)),  Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)),  Cell(cell_idx_c(6), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(7), cell_idx_c(2), cell_idx_c(0)),  Cell(cell_idx_c(8), cell_idx_c(2), cell_idx_c(0)),
      Cell(cell_idx_c(9), cell_idx_c(2), cell_idx_c(0)),  Cell(cell_idx_c(10), cell_idx_c(2), cell_idx_c(0))
   };
   WALBERLA_LOG_INFO("Testing interface conversion with wetting cells, contact angle=" << contactAngle << " degrees");
   computedInterfaceCells = runSimulation< LatticeModel_T >(timesteps, contactAngle);
   vectorsEqual           = std::is_permutation(computedInterfaceCells.begin(), computedInterfaceCells.end(),
                                                expectedInterfaceCells.begin(), expectedInterfaceCells.end());
   if (!vectorsEqual) { WALBERLA_ABORT("Wrong cells converted."); }

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   WALBERLA_LOG_INFO("Testing with D3Q19 stencil.");
   testWettingConversion< lbm::D3Q19< lbm::collision_model::SRT, true, lbm::force_model::None, 2 > >();

   WALBERLA_LOG_INFO("Testing with D3Q27 stencil.");
   testWettingConversion< lbm::D3Q27< lbm::collision_model::SRT, true, lbm::force_model::None, 2 > >();

   return EXIT_SUCCESS;
}

} // namespace WettingConversionTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::WettingConversionTest::main(argc, argv); }