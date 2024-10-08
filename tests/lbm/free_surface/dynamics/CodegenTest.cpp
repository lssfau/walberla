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
//! \file CodegenTest.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test equivalence of generated LBM kernels in the free surface implementation.
//!
//! Simulates 100 time steps of a moving drop with diameter 2 in a periodic 4x4x4 domain. The drop moves due to a
//! constant body force.
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/InitFunctions.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/lattice_model/D3Q19.h"

#include <type_traits>

// include files generated by lbmpy
#include "GeneratedLatticeModel_FreeSurface.h"

namespace walberla
{
namespace free_surface
{
namespace CodegenTest
{

using ScalarField_T          = GhostLayerField< real_t, 1 >;
using VectorField_T          = GhostLayerField< Vector3< real_t >, 1 >;
using VectorFieldFlattened_T = GhostLayerField< real_t, 3 >;

template< bool useCodegen >
void runSimulation()
{
   using CollisionModel_T = lbm::collision_model::SRT;
   using ForceModel_T     = lbm::force_model::GuoField< VectorField_T >;
   using LatticeModel_T   = typename std::conditional< useCodegen, lbm::GeneratedLatticeModel_FreeSurface,
                                                     lbm::D3Q19< CollisionModel_T, true, ForceModel_T, 2 > >::type;
   using PdfField_T       = lbm::PdfField< LatticeModel_T >;

   using Communication_T = blockforest::SimpleCommunication< typename LatticeModel_T::CommunicationStencil >;

   using flag_t                        = uint32_t;
   using FlagField_T                   = FlagField< flag_t >;
   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(4), uint_c(4), uint_c(4));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, true, true);                                    // periodicity

   // physics parameters
   const real_t dropDiameter            = real_c(2);
   const real_t relaxationRate          = real_c(1.8);
   const real_t surfaceTension          = real_c(1e-5);
   const bool enableWetting             = false;
   const real_t contactAngle            = real_c(0);
   const Vector3< real_t > acceleration = Vector3< real_t >(real_c(1e-5), real_c(0), real_c(0));

   // model parameters
   const std::string pdfReconstructionModel      = "NormalBasedKeepCenter";
   const std::string pdfRefillingModel           = "EquilibriumRefilling";
   const std::string excessMassDistributionModel = "EvenlyAllInterface";
   const std::string curvatureModel              = "FiniteDifferenceMethod";
   const bool useSimpleMassExchange              = false;
   const real_t cellConversionThreshold          = real_c(1e-2);
   const real_t cellConversionForceThreshold     = real_c(1e-1);

   // add force field
   std::shared_ptr< BlockDataID > forceDensityFieldIDPtr = nullptr;

   std::shared_ptr< LatticeModel_T > latticeModel;

   // create lattice model
   if constexpr (useCodegen)
   {
      // add force field of type 'GhostLayerField< real_t, 3 >' as required by pystencils
      forceDensityFieldIDPtr = std::make_shared< BlockDataID >(field::addToStorage< VectorFieldFlattened_T >(
         blockForest, "Force density field (codegen)", real_c(0), field::fzyx, uint_c(1)));
      latticeModel           = std::make_shared< LatticeModel_T >(*forceDensityFieldIDPtr, relaxationRate);
   }
   else
   {
      forceDensityFieldIDPtr = std::make_shared< BlockDataID >(field::addToStorage< VectorField_T >(
         blockForest, "Force density field", Vector3< real_t >(0), field::fzyx, uint_c(1)));
      latticeModel =
         std::make_shared< LatticeModel_T >(CollisionModel_T(relaxationRate), ForceModel_T(*forceDensityFieldIDPtr));
   }

   WALBERLA_ASSERT_NOT_NULLPTR(forceDensityFieldIDPtr);

   // add various fields
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", *latticeModel, field::fzyx);
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(0.0), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // add drop to fill level field
   const geometry::Sphere sphereDrop(Vector3< real_t >(real_c(0.5) * real_c(domainSize[0]),
                                                       real_c(0.5) * real_c(domainSize[1]),
                                                       real_c(0.5) * real_c(domainSize[2])),
                                     real_c(dropDiameter) * real_c(0.5));
   bubble_model::addBodyToFillLevelField< geometry::Sphere >(*blockForest, fillFieldID, sphereDrop, false);

   // initialize flag field from fill level field
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   if constexpr (useCodegen)
   {
      initForceDensityFieldCodegen< PdfField_T, FlagField_T, VectorFieldFlattened_T, ScalarField_T >(
         blockForest, *forceDensityFieldIDPtr, fillFieldID, pdfFieldID, flagFieldID, flagInfo, acceleration);
   }
   else
   {
      // initialize force density field
      initForceDensityField< PdfField_T, FlagField_T, VectorField_T, ScalarField_T >(
         blockForest, *forceDensityFieldIDPtr, fillFieldID, pdfFieldID, flagFieldID, flagInfo, acceleration);
   }

   // initial communication
   Communication_T(blockForest, pdfFieldID, fillFieldID, flagFieldID, *forceDensityFieldIDPtr)();

   // add bubble model
   const std::shared_ptr< bubble_model::BubbleModelConstantPressure > bubbleModel =
      std::make_shared< bubble_model::BubbleModelConstantPressure >(real_c(1));

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(100));

   // Laplace pressure = 2 * surface tension * curvature; curvature computation is not necessary with 0 surface tension
   bool computeCurvature = false;
   if (!realIsEqual(surfaceTension, real_c(0), real_c(1e-14))) { computeCurvature = true; }

   // add surface geometry handler
   const SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > geometryHandler(
      blockForest, freeSurfaceBoundaryHandling, fillFieldID, curvatureModel, computeCurvature, enableWetting,
      contactAngle);

   const ConstBlockDataID curvatureFieldID = geometryHandler.getConstCurvatureFieldID();
   const ConstBlockDataID normalFieldID    = geometryHandler.getConstNormalFieldID();

   geometryHandler.addSweeps(timeloop);

   // add boundary handling for standard boundaries and free surface boundaries
   const SurfaceDynamicsHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T, useCodegen,
                                 VectorFieldFlattened_T >
      dynamicsHandler(blockForest, pdfFieldID, flagFieldID, fillFieldID, *forceDensityFieldIDPtr, normalFieldID,
                      curvatureFieldID, freeSurfaceBoundaryHandling, bubbleModel, pdfReconstructionModel,
                      pdfRefillingModel, excessMassDistributionModel, relaxationRate, acceleration, surfaceTension,
                      useSimpleMassExchange, cellConversionThreshold, cellConversionForceThreshold);

   dynamicsHandler.addSweeps(timeloop);

   timeloop.run();

   // check fill level (must be identical in lattice model from waLBerla and from lbmpy)
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const ScalarField_T* const fillField = blockIt->template getData< const ScalarField_T >(fillFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(1)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.512081), real_c(1e-4));
         }
         if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(1)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.529967), real_c(1e-4));
         }
         if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(1)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.512081), real_c(1e-4));
         }
         if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(1)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.529967), real_c(1e-4));
         }
         if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(2)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.512081), real_c(1e-4));
         }
         if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(2)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.529967), real_c(1e-4));
         }
         if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(2)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.512081), real_c(1e-4));
         }
         if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(2)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.529967), real_c(1e-4));
         }
      }); // WALBERLA_FOR_ALL_CELLS
   }
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment walberlaEnv(argc, argv);

   WALBERLA_LOG_INFO_ON_ROOT("Testing with lattice model from waLBerla.");
   runSimulation< false >();

   WALBERLA_LOG_INFO_ON_ROOT("Testing with lattice model generated by lbmpy.");
   runSimulation< true >();

   return EXIT_SUCCESS;
}

} // namespace CodegenTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::CodegenTest::main(argc, argv); }
