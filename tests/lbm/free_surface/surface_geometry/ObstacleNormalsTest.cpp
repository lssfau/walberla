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
//! \file ObstacleNormalsTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test if mean obstacle normal computation is correct in setup with inclined plane and different inclination
//! angles.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"

#include "field/AddToStorage.h"

#include "lbm/free_surface/surface_geometry/ObstacleNormalSweep.h"

#include "stencil/D3Q27.h"
#include "stencil/D3Q7.h"

#include "timeloop/SweepTimeloop.h"

#include <cmath>

namespace walberla
{
namespace free_surface
{
namespace ObstacleNormalsTest
{
// define types
using Stencil_T = stencil::D3Q27;
using flag_t    = uint8_t;

const FlagUID Fluid_Flag("fluid");
const FlagUID FluidNearSolid_Flag("fluid near solid");
const FlagUID Solid_Flag("no slip");
const Set< FlagUID > All_Fluid_Flags = setUnion< FlagUID >(Fluid_Flag, FluidNearSolid_Flag);

// define fields
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< flag_t >;

void test(real_t degreeInclinationAngle)
{
   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(20), uint_c(3), uint_c(20));
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
   BlockDataID obstacleNormalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Obstacle normals", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(
      blockForest, "flags", uint_c(2)); // flag field must have two ghost layers as in regular FSLBM application code

   // initialize flag field as inclination, the lower part of the domain will be solid
   const real_t tanAngle = real_c(std::tan(degreeInclinationAngle * math::pi / real_c(180)));
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      FlagField_T* const flagField = blockIt->getData< FlagField_T >(flagFieldID);

      const auto fluidFlag          = flagField->registerFlag(Fluid_Flag);
      const auto fluidNearSolidFlag = flagField->registerFlag(FluidNearSolid_Flag);
      const auto solidFlag          = flagField->registerFlag(Solid_Flag);

      // create inclination in flag field
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
         flagField, uint_c(2),
         // all cells below the specified inclination angle are marked as solid
         if (real_c(x) * tanAngle <= real_c(z)) { flagField->get(x, y, z) = solidFlag; } else {
            flagField->get(x, y, z) = fluidFlag;
         }); // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

      // mark fluid cells that have neighboring solid cells ( only in these cells, the obstacle normal will be computed)
      WALBERLA_FOR_ALL_CELLS(
         flagFieldIt, flagField,
         if (isFlagSet(flagFieldIt, solidFlag) && isFlagInNeighborhood< stencil::D3Q7 >(flagFieldIt, fluidFlag)) {
            *flagFieldIt = fluidNearSolidFlag;
         }); // WALBERLA_FOR_ALL_CELLS
   }

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   // add obstacle normals sweep
   ObstacleNormalSweep< Stencil_T, FlagField_T, VectorField_T > obstNormalsSweep(
      obstacleNormalFieldID, flagFieldID, FluidNearSolid_Flag, All_Fluid_Flags, Solid_Flag, true, false, false);
   timeloop.add() << Sweep(obstNormalsSweep, "Obstacle normals sweep");

   // perform a single time step
   timeloop.singleStep();

   // compute analytical normal
   const real_t sinAngle = real_c(std::sin(degreeInclinationAngle * math::pi / real_c(180)));
   const real_t cosAngle = real_c(std::cos(degreeInclinationAngle * math::pi / real_c(180)));
   const Vector3< real_t > analyticalNormal(-sinAngle, real_c(0), cosAngle);

   // compare analytical normal with the average of all computed obstacle normals
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const VectorField_T* const obstacleNormalField = blockIt->getData< const VectorField_T >(obstacleNormalFieldID);
      const FlagField_T* const flagField             = blockIt->getData< const FlagField_T >(flagFieldID);

      const auto fluidNearSolidFlag = flagField->getFlag(FluidNearSolid_Flag);

      real_t averageObstNormalX = real_c(0);
      real_t averageObstNormalY = real_c(0);
      real_t averageObstNormalZ = real_c(0);
      uint_t cellCount          = uint_c(0);
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP(obstacleNormalField,
            omp parallel for schedule(static) reduction(+:averageObstNormalX) reduction(+:averageObstNormalY)
                                              reduction(+:averageObstNormalZ) reduction(+:cellCount), {
         // skip cells at the domain boundary; obstacle normal in these cells deviates further from analytical solution
         // since neighborhood for obstacle computation is too small
         if (x == cell_idx_c(0) || y == cell_idx_c(0) || z == cell_idx_c(0) ||
             x == cell_idx_c(domainSize[0] - uint_c(1)) || y == cell_idx_c(domainSize[1] - uint_c(1)) ||
             z == cell_idx_c(domainSize[2] - uint_c(1)))
         {
            continue;
         }

         if (!isFlagSet(flagField->get(x, y, z), fluidNearSolidFlag))
         {
            // obstacle normal must be zero in all cells that are not of type fluidNearSolid
            WALBERLA_CHECK_FLOAT_EQUAL(obstacleNormalField->get(x, y, z), Vector3< real_t >(real_c(0)));
         }
         else
         {
            ++cellCount;
            averageObstNormalX += obstacleNormalField->get(x, y, z)[0];
            averageObstNormalY += obstacleNormalField->get(x, y, z)[1];
            averageObstNormalZ += obstacleNormalField->get(x, y, z)[2];
         }}); // WALBERLA_FOR_ALL_CELLS_XYZ_OMP

      // compute average obstacle normal
      Vector3< real_t > averageObstNormal(averageObstNormalX, averageObstNormalY, averageObstNormalZ);
      averageObstNormal /= real_c(cellCount);

      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(averageObstNormal, analyticalNormal, real_c(0.1));
   }
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);
   // test with various inclination angles (only up to 87 degree, as otherwise in this setup no more fluid cells remain)
   for (uint_t angle = uint_c(1); angle <= uint_c(86); angle += uint_c(1))
   {
      WALBERLA_LOG_INFO_ON_ROOT("Testing inclination angle " << real_c(angle) << " degree.")
      test(real_c(angle));
   }

   return EXIT_SUCCESS;
}
} // namespace ObstacleNormalsTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::ObstacleNormalsTest::main(argc, argv); }