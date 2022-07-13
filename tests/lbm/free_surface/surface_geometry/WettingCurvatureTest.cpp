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
//! \file WettingCurvatureTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Initialize a drop as a cylinder section near a solid wall and evaluate the resulting wetting curvature.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"

#include "field/AddToStorage.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/bodies/Cylinder.h"
#include "geometry/initializer/OverlapFieldFromBody.h"

#include "lbm/free_surface/surface_geometry/ContactAngle.h"
#include "lbm/free_surface/surface_geometry/CurvatureSweep.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"
#include "lbm/free_surface/surface_geometry/ObstacleFillLevelSweep.h"
#include "lbm/free_surface/surface_geometry/ObstacleNormalSweep.h"
#include "lbm/free_surface/surface_geometry/SmoothingSweep.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/Initialization.h"
#include "vtk/VTKOutput.h"

namespace walberla
{
namespace free_surface
{
namespace WettingCurvatureTest
{
// define types
using Stencil_T = stencil::D3Q27;

using flag_t = uint32_t;

const FlagUID liquidFlagID("liquid");
const FlagUID interfaceFlagID("interface");
const FlagUID gasFlagID("gas");
const FlagUID solidID("solid");
const Set< FlagUID > liquidInterfaceGasFlagIDs =
   setUnion< FlagUID >(setUnion< FlagUID >(liquidFlagID, interfaceFlagID), gasFlagID);

// define fields
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< flag_t >;

real_t computeCurvature(real_t contactAngle, bool useTriangulation)
{
   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(6), uint_c(5), uint_c(5));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          false, false, false);                                 // periodicity

   // add fields
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(0.0), field::fzyx, uint_c(1));
   BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normals", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID curvatureFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Curvature", real_c(0), field::fzyx, uint_c(1));
   BlockDataID obstacleNormalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Obstacle normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(
      blockForest, "Flags", uint_c(2)); // flag field must have two ghost layers as in regular FSLBM application code
   BlockDataID smoothFillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Smooth fill levels", real_c(1.0), field::fzyx, uint_c(1));

   // add liquid drop
   Vector3< real_t > midpoint1(real_c(5), real_c(1), real_c(0));
   Vector3< real_t > midpoint2(real_c(5), real_c(1), real_c(5));
   geometry::Cylinder cylinder(midpoint1, midpoint2, real_c(5) * real_c(0.5));
   geometry::initializer::OverlapFieldFromBody(*blockForest, fillFieldID).init(cylinder, true);

   // set flags
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      FlagField_T* const flagField   = blockIt->getData< FlagField_T >(flagFieldID);
      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(fillFieldID);

      const auto gas       = flagField->registerFlag(gasFlagID);
      const auto liquid    = flagField->registerFlag(liquidFlagID);
      const auto interface = flagField->registerFlag(interfaceFlagID);
      const auto solid     = flagField->registerFlag(solidID);

      // initialize whole flag field as gas
      WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, { *flagFieldIt = gas; }) // WALBERLA_FOR_ALL_CELLS

      // initialize flag field according to fill level
      WALBERLA_FOR_ALL_CELLS_XYZ(
         flagField,
         if (y == 0) {
            flagField->get(x, y, z) = solid;
            fillField->get(x, y, z) = real_c(0);
         }

         if (fillField->get(x, y, z) <= real_c(0) && flagField->get(x, y, z) != solid) {
            flagField->get(x, y, z) = gas;
         }

         if (fillField->get(x, y, z) >= real_c(1)) { flagField->get(x, y, z) = liquid; }

         // not using else here to avoid overwriting solid flags
         if (fillField->get(x, y, z) < real_c(1) && fillField->get(x, y, z) > real_c(0)) {
            flagField->get(x, y, z) = interface;
         }) // WALBERLA_FOR_ALL_CELLS_XYZ
   }

   ContactAngle contactAngleObj = ContactAngle(contactAngle);

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   if (useTriangulation) // use local triangulation for curvature computation
   {
      // add sweep for computing obstacle normals in interface cells near obstacle cells
      ObstacleNormalSweep< Stencil_T, FlagField_T, VectorField_T > obstNormalsSweep(
         obstacleNormalFieldID, flagFieldID, interfaceFlagID, liquidInterfaceGasFlagIDs, solidID, true, false, false);
      timeloop.add() << Sweep(obstNormalsSweep, "Obstacle normals sweep");

      // add sweep for computing interface normals
      NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalsSweep(
         normalFieldID, fillFieldID, flagFieldID, interfaceFlagID, liquidInterfaceGasFlagIDs, solidID, false, false,
         true, false);
      timeloop.add() << Sweep(normalsSweep, "Normal sweep");

      // add sweep for computing curvature (including wetting)
      CurvatureSweepLocalTriangulation< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > curvatureSweep(
         blockForest, curvatureFieldID, normalFieldID, fillFieldID, flagFieldID, obstacleNormalFieldID, interfaceFlagID,
         solidID, true, contactAngleObj);
      timeloop.add() << Sweep(curvatureSweep, "Curvature sweep");
   }
   else // use finite differences for curvature computation
   {
      // add sweep for computing obstacle normals in obstacle cells
      ObstacleNormalSweep< Stencil_T, FlagField_T, VectorField_T > obstNormalsSweep(
         obstacleNormalFieldID, flagFieldID, interfaceFlagID, liquidInterfaceGasFlagIDs, solidID, false, true, true);
      timeloop.add() << Sweep(obstNormalsSweep, "Obstacle normals sweep");

      // add sweep for reflecting fill level into obstacle cells
      ObstacleFillLevelSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > obstacleFillLevelSweep(
         smoothFillFieldID, fillFieldID, flagFieldID, obstacleNormalFieldID, liquidInterfaceGasFlagIDs, solidID);
      timeloop.add() << Sweep(obstacleFillLevelSweep, "Obstacle fill level sweep");

      // add sweep for smoothing the fill level field (uses fill level values from obstacle cells set by
      // ObstacelFillLevelSweep)
      SmoothingSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > smoothingSweep(
         smoothFillFieldID, fillFieldID, flagFieldID, liquidInterfaceGasFlagIDs, solidID, true);
      timeloop.add() << Sweep(smoothingSweep, "Smoothing sweep");

      // add sweep for computing interface normals
      NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalsSweep(
         normalFieldID, smoothFillFieldID, flagFieldID, interfaceFlagID, liquidInterfaceGasFlagIDs, solidID, true, true,
         false, true);
      timeloop.add() << Sweep(normalsSweep, "Normal sweep");

      // add sweep for computing curvature (including wetting)
      CurvatureSweepFiniteDifferences< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > curvatureSweep(
         curvatureFieldID, normalFieldID, obstacleNormalFieldID, flagFieldID, interfaceFlagID,
         liquidInterfaceGasFlagIDs, solidID, true, contactAngleObj);
      timeloop.add() << Sweep(curvatureSweep, "Curvature sweep");
   }

   // run one time step
   timeloop.singleStep();

   // get the curvature in cell (2, 1, 2)
   real_t curvature = real_c(0);
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const ScalarField_T* const curvatureField = blockIt->getData< const ScalarField_T >(curvatureFieldID);
      curvature                                 = curvatureField->get(2, 1, 2);
   }

   //   auto vtkFlagField =
   //      field::createVTKOutput< FlagField_T >(flagFieldID, *blockForest, "flag_field", uint_c(1), uint_c(0));
   //   vtkFlagField();
   //   auto vtkFillField =
   //      field::createVTKOutput< ScalarField_T >(fillFieldID, *blockForest, "fill_field", uint_c(1), uint_c(0));
   //   vtkFillField();
   //   auto vtkCurvatureField =
   //      field::createVTKOutput< ScalarField_T >(curvatureFieldID, *blockForest, "curvature_field", uint_c(1),
   //      uint_c(0));
   //   vtkCurvatureField();
   //   auto vtkNormalField =
   //      field::createVTKOutput< VectorField_T >(normalFieldID, *blockForest, "normal_field", uint_c(1), uint_c(0));
   //   vtkNormalField();
   //   auto vtkObstNormalField = field::createVTKOutput< VectorField_T >(obstacleNormalFieldID, *blockForest,
   //                                                                     "obst_normal_field", uint_c(1), uint_c(0));
   //   vtkObstNormalField();
   //   auto vtkSmoothFillField = field::createVTKOutput< ScalarField_T >(smoothFillFieldID, *blockForest,
   //                                                                     "smooth_fill_field", uint_c(1), uint_c(0));
   //   vtkSmoothFillField();

   return curvature;
}

// the following values have been obtained with a version of the local triangulation curvature computation algorithm
// that is assumed to be correct (the angles computed in computeArtificalWallPoint() have been manually verified with
// ParaView for some of the values)
std::vector< std::pair< real_t, real_t > > expectedSolutionTriangulation()
{
   // pair contains: [0] contact angle, [1] expected curvature
   std::vector< std::pair< real_t, real_t > > testcases;

   testcases.emplace_back(real_c(0), real_c(-0.208893));
   testcases.emplace_back(real_c(1), real_c(-0.208832));
   testcases.emplace_back(real_c(10), real_c(-0.202813));
   testcases.emplace_back(real_c(30), real_c(-0.15704));
   testcases.emplace_back(real_c(45), real_c(-0.0994945));
   testcases.emplace_back(real_c(65), real_c(-0.00311236));
   testcases.emplace_back(real_c(65.5), real_c(-0.000512801));

   // IMPORTANT REMARK REGARDING THE CHANGE IN THE SIGN OF THE CURVATURE:
   // A change in the sign of the curvature would only be expected when the specified contact angle is equivalent to the
   // angle that is already present. However, the algorithm ensures that the constructed virtual wall point is valid
   // during curvature computation (such that the computed triangle by local triangulation is not degenerated).
   // Therefore, the algorithm internally changes the target contact angle (see lines 10 to 17 in algorithm 6.2 in
   // dissertation of S. Donath). The specified contact angle will then slowly be approached in subsequent time steps.

   testcases.emplace_back(real_c(65.8), real_c(0.00105015));
   testcases.emplace_back(real_c(66), real_c(0.00209344));
   testcases.emplace_back(real_c(66.5), real_c(0.00470618));
   testcases.emplace_back(real_c(67), real_c(0.00732526));
   testcases.emplace_back(real_c(70), real_c(0.0231628));
   testcases.emplace_back(real_c(75), real_c(0.0499505));
   testcases.emplace_back(real_c(77), real_c(0.0607714));
   testcases.emplace_back(real_c(78), real_c(0.0661992));

   // this contact angle is already reached by the initial setup
   // => the algorithm should take the route with if(realIsEqual(...)) in CurvatureSweepTR
   // => this route has been found to cause problems when using single precision (see comment in CurvatureSweepTR())
#ifdef WALBERLA_DOUBLE_ACCURACY
   testcases.emplace_back(real_c(78.40817854), real_c(0.0684177));
#endif

   testcases.emplace_back(real_c(80), real_c(0.0770832));
   testcases.emplace_back(real_c(90), real_c(0.131739));
   testcases.emplace_back(real_c(120), real_c(0.287743));
   testcases.emplace_back(real_c(160), real_c(0.431406));
   testcases.emplace_back(real_c(180), real_c(0.312837));

   return testcases;
}

// the following values have been obtained with a version of the finite difference curvature computation algorithm
// that is assumed to be correct
std::vector< std::pair< real_t, real_t > > expectedSolutionFiniteDifferences()
{
   std::vector< std::pair< real_t, real_t > > testcases;

   testcases.emplace_back(real_c(0), real_c(-0.0454018));
   testcases.emplace_back(real_c(1), real_c(-0.0474911));
   testcases.emplace_back(real_c(10), real_c(-0.0641205));
   testcases.emplace_back(real_c(30), real_c(-0.0810346));
   testcases.emplace_back(real_c(45), real_c(-0.0622399));
   testcases.emplace_back(real_c(65), real_c(0.0397875));
   testcases.emplace_back(real_c(65.5), real_c(0.0436974));
   testcases.emplace_back(real_c(65.8), real_c(0.046073));
   testcases.emplace_back(real_c(66), real_c(0.0476689));
   testcases.emplace_back(real_c(66.5), real_c(0.0517007));
   testcases.emplace_back(real_c(67), real_c(0.0557916));
   testcases.emplace_back(real_c(70), real_c(0.0814937));
   testcases.emplace_back(real_c(75), real_c(0.127884));
   testcases.emplace_back(real_c(77), real_c(0.147265));
   testcases.emplace_back(real_c(78), real_c(0.157051));
   testcases.emplace_back(real_c(78.40817854), real_c(0.161057));
   testcases.emplace_back(real_c(80), real_c(0.176711));
   testcases.emplace_back(real_c(90), real_c(0.271672));
   testcases.emplace_back(real_c(120), real_c(0.448388));
   testcases.emplace_back(real_c(160), real_c(0.487812));
   testcases.emplace_back(real_c(180), real_c(0.467033));

   return testcases;
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   // test with local triangulation curvature computation
   std::vector< std::pair< real_t, real_t > > testcases = expectedSolutionTriangulation();
   for (const auto& i : testcases)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Testing contact angle=" << i.first
                                                         << " with local triangulation curvature computation");
      const real_t curvature = computeCurvature(i.first, true);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(curvature, i.second, real_c(1e-5));
   }

   // test with finite difference curvature computation
   testcases = expectedSolutionFiniteDifferences();
   for (const auto& i : testcases)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Testing contact angle=" << i.first << " with finite difference curvature computation");
      const real_t curvature = computeCurvature(i.first, false);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(curvature, i.second, real_c(1e-6));
   }

   return EXIT_SUCCESS;
}
} // namespace WettingCurvatureTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::WettingCurvatureTest::main(argc, argv); }