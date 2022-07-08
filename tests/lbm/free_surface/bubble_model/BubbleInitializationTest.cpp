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
//! \file BubbleInitializationTest.cpp
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test bubble initialization from fill level by evaluating initial bubble volumes (within and across blocks).
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "field/AddToStorage.h"
#include "field/Printers.h"

#include "geometry/bodies/Sphere.h"

#include "lbm/free_surface/bubble_model/BubbleModel.h"
#include "lbm/free_surface/bubble_model/Geometry.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

namespace walberla
{
namespace free_surface
{
namespace BubbleInitializationTest
{
// define field
using ScalarField_T = GhostLayerField< real_t, 1 >;

// class derived from BubbleModel to access its protected members for testing purposes
template< typename Stencil_T >
class BubbleModelTest : public bubble_model::BubbleModel< Stencil_T >
{
 public:
   BubbleModelTest(const std::shared_ptr< StructuredBlockForest >& blockForest)
      : bubble_model::BubbleModel< Stencil_T >(blockForest, true)
   {}

   real_t getBubbleInitVolume(bubble_model::BubbleID id)
   {
      WALBERLA_ASSERT_GREATER(bubble_model::BubbleModel< Stencil_T >::getBubbles().size(), id);
      return bubble_model::BubbleModel< Stencil_T >::getBubbles()[id].getInitVolume();
   }
}; // class BubbleModelTest

template< typename Stencil_T >
void testBubbleInitialization()
{
   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(2), uint_c(1), uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(56), uint_c(16), uint_c(1));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, false, false,                                   // periodicity
                                          true);                                                // global info

   // add fill level field
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(1.0), field::fzyx, uint_c(1));

   const real_t xQuarter = real_c(domainSize[0]) / real_c(4);
   const real_t yHalf    = real_c(domainSize[1]) / real_c(2);

   // add sphere in the left half (in x-direction) of the domain
   bubble_model::addBodyToFillLevelField(
      *blockForest, fillFieldID,
      geometry::Sphere(Vector3< real_t >(real_c(0.75) * xQuarter, yHalf, real_c(0)), xQuarter * real_c(0.25)), true);

   // add sphere in the center of the domain (across blockForest)
   bubble_model::addBodyToFillLevelField(
      *blockForest, fillFieldID,
      geometry::Sphere(Vector3< real_t >(real_c(domainSize[0]) * real_c(0.5), yHalf, real_c(0)), yHalf), true);

   // add sphere in the right half (in x-direction) of the domain
   bubble_model::addBodyToFillLevelField(
      *blockForest, fillFieldID,
      geometry::Sphere(Vector3< real_t >(real_c(3.25) * xQuarter, yHalf, real_c(0)), xQuarter * real_c(0.25)), true);

   // create bubble model
   BubbleModelTest< Stencil_T > bubbleModel(blockForest);
   bubbleModel.initFromFillLevelField(fillFieldID);

   // test correctness of bubble volumes (bubble IDs were determined empirically for this test)
   // left bubble
   WALBERLA_CHECK_LESS(
      std::abs(bubbleModel.getBubbleInitVolume(2) - (xQuarter * real_c(0.25) * xQuarter * real_c(0.25) * math::pi)),
      real_c(1.07));

   // center bubble
   WALBERLA_CHECK_LESS(std::abs(bubbleModel.getBubbleInitVolume(0) - (yHalf * yHalf * math::pi)), real_c(1.12));

   // right bubble
   WALBERLA_CHECK_LESS(
      std::abs(bubbleModel.getBubbleInitVolume(1) - (xQuarter * real_c(0.25) * xQuarter * real_c(0.25) * math::pi)),
      real_c(1.07));

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   auto mpiManager = MPIManager::instance();

   WALBERLA_CHECK_EQUAL(mpiManager->numProcesses(), 2);

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D2Q9 stencil.");
   testBubbleInitialization< stencil::D2Q9 >();

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D3Q19 stencil.");
   testBubbleInitialization< stencil::D3Q19 >();

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D3Q27 stencil.");
   testBubbleInitialization< stencil::D3Q27 >();

   return EXIT_SUCCESS;
}

} // namespace BubbleInitializationTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::BubbleInitializationTest::main(argc, argv); }