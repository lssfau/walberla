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
//! \file MovingSpheresTest.cpp
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//!
//! \brief Test volume conservation of moving bubbles, bubble merging and bubble splitting.
//!
//! A spherical bubble is moved towards a static spherical bubble in the center of the domain. The bubbles merge and
//! split again, as the movement of the spherical volume is continued after merging.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include "geometry/bodies/Sphere.h"

#include "lbm/free_surface/bubble_model/BubbleModel.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "BubbleBodyMover.h"

namespace walberla
{
namespace free_surface
{
namespace MovingSpheresTest
{
// class derived from BubbleModel to access its protected members for testing purposes
template< typename Stencil_T >
class BubbleModelTest : public bubble_model::BubbleModel< Stencil_T >
{
 public:
   BubbleModelTest(const std::shared_ptr< StructuredBlockForest >& blockForest)
      : bubble_model::BubbleModel< Stencil_T >(blockForest, true)
   {}

   static void testMovingSpheres();
}; // class BubbleModelTest

template< typename Stencil_T >
void BubbleModelTest< Stencil_T >::testMovingSpheres()
{
   auto mpiManager = MPIManager::instance();

   WALBERLA_CHECK_EQUAL(mpiManager->numProcesses(), 2);

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(2), uint_c(1), uint_c(1));
   Vector3< uint_t > domainSize(uint_c(60), uint_c(16), uint_c(16));

   if (Stencil_T::D == uint_c(2)) { domainSize[2] = uint_c(1); }

   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, false, false);                                  // periodicity

   // create bubble model
   std::shared_ptr< BubbleModelTest > bubbleModel = std::make_shared< BubbleModelTest >(blockForest);

   // create bubble body mover for moving bubbles
   bubble_model::BubbleBodyMover< geometry::Sphere, Stencil_T > bubbleSphereMover(blockForest, bubbleModel);

   Vector3< real_t > domainCenter(real_c(domainSize[0]) * real_c(0.5), real_c(domainSize[1]) * real_c(0.5),
                                  real_c(domainSize[2]) * real_c(0.5));

   // create a static spherical bubble in the center of the domain
   geometry::Sphere sphereCenter(domainCenter, real_c(5));
   auto doNotMoveSphere = [](geometry::Sphere&) {};
   bubbleSphereMover.addBody(sphereCenter, doNotMoveSphere);

   // create a moving spherical bubble in the left (in x-direction) half of the domain
   geometry::Sphere sphereLeft(Vector3< real_t >(real_c(8.0), domainCenter[1], domainCenter[2]), real_c(5));
   auto moveSphere = [](geometry::Sphere& sphere) {
      // midpoint of the sphere is shifted by 1 cell in positive x-direction at each call
      sphere.setMidpoint(sphere.midpoint() + Vector3< real_t >(real_c(1), real_c(0), real_c(0)));
   };
   bubbleSphereMover.addBody(sphereLeft, moveSphere);

   // initialize the just added bubbles
   bubbleSphereMover.initAddedBodies();

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(40));
   timeloop.addFuncBeforeTimeStep(bubbleSphereMover, "Move bubbles");
   timeloop.addFuncAfterTimeStep(std::bind(&bubble_model::BubbleModel< Stencil_T >::update, bubbleModel));

   real_t singleBubbleVolume = real_c(523.346);

   if (Stencil_T::D == uint_c(2)) { singleBubbleVolume = real_c(78.1987); }

   uint_t timestep = uint_c(0);

   // at time step 8, there should still be two bubbles
   for (; timestep < uint_c(8); ++timestep)
   {
      timeloop.singleStep();
   }
   WALBERLA_CHECK_EQUAL(bubbleModel->getBubbles().size(), 2);
   WALBERLA_CHECK_LESS(std::abs(bubbleModel->getBubbles()[0].getInitVolume() - singleBubbleVolume), real_c(0.5));
   WALBERLA_CHECK_LESS(std::abs(bubbleModel->getBubbles()[1].getInitVolume() - singleBubbleVolume), real_c(0.5));

   // at time step 12 the bubbles should have merged
   for (; timestep < uint_c(12); ++timestep)
   {
      timeloop.singleStep();
   }
   WALBERLA_CHECK_EQUAL(bubbleModel->getBubbles().size(), 1);
   WALBERLA_CHECK_LESS(std::abs(bubbleModel->getBubbles()[0].getInitVolume() - 2 * singleBubbleVolume), real_c(0.5));

   // at time step 35 the bubbles should have split again
   for (; timestep < uint_c(35); ++timestep)
   {
      timeloop.singleStep();
   }
   WALBERLA_CHECK_EQUAL(bubbleModel->getBubbles().size(), 2);
   WALBERLA_CHECK_LESS(std::abs(bubbleModel->getBubbles()[0].getInitVolume() - singleBubbleVolume), real_c(0.5));
   WALBERLA_CHECK_LESS(std::abs(bubbleModel->getBubbles()[1].getInitVolume() - singleBubbleVolume), real_c(0.5));

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D2Q9 stencil.")
   BubbleModelTest< stencil::D2Q9 >::testMovingSpheres();

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D3Q19 stencil.")
   BubbleModelTest< stencil::D3Q19 >::testMovingSpheres();

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D3Q27 stencil.")
   BubbleModelTest< stencil::D3Q27 >::testMovingSpheres();

   return EXIT_SUCCESS;
}
} // namespace MovingSpheresTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::MovingSpheresTest::main(argc, argv); }