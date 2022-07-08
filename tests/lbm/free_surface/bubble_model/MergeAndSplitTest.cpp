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
//! \file MergeAndSplitTest.cpp
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test bubble merging and splitting in a complex multi-process scenario.
//!
//! Initialize the fill levels, flags and bubble model from image MergeAndSplitTestUnconnected.png. This sets up a
//! complex scenario of 12 bubbles that are located on 10 processes. Bubble merging is tested by loading image
//! MergeAndSplitTestConnected.png. By loading the initial image in a second time step, bubble splitting is tested.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include "geometry/initializer/ScalarFieldFromGrayScaleImage.h"
#include "geometry/structured/GrayScaleImage.h"

#include "lbm/free_surface/bubble_model/BubbleModel.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <vector>

#include "BubbleModelTester.h"

namespace walberla
{
namespace free_surface
{
namespace MergeAndSplitTest
{
// class derived from BubbleModel to access its protected members for testing purposes
template< typename Stencil_T >
class BubbleModelTest : public bubble_model::BubbleModel< Stencil_T >
{
 public:
   BubbleModelTest(const std::shared_ptr< StructuredBlockForest >& blockStorage)
      : bubble_model::BubbleModel< Stencil_T >(blockStorage, true)
   {}

   static void testComplexMerge();
}; // class BubbleModelTest

// initialize and update a source and destination field for fill levels and flags from images; in an alternating manner,
// update the fields such that the bubble model must either handle bubble merging or splitting on every second call
template< typename Stencil_T >
class ImageMover : public bubble_model::BubbleModelTester< Stencil_T >
{
 public:
   ImageMover(const std::shared_ptr< StructuredBlockStorage >& blockStorage,
              const std::shared_ptr< bubble_model::BubbleModel< Stencil_T > >& bubbleModel)
      : bubble_model::BubbleModelTester< Stencil_T >(blockStorage, bubbleModel),
        imgInitializerSrc_(*blockStorage, bubble_model::BubbleModelTester< Stencil_T >::srcFillLevelFieldID_),
        imgInitializerDst_(*blockStorage, bubble_model::BubbleModelTester< Stencil_T >::dstFillLevelFieldID_),
        calls_(uint_c(0))
   {
      // load image of test scenario with bubbles that are not connected
      geometry::GrayScaleImage img("MergeAndSplitTestUnconnected.png");

      // initialize fill level field from image
      imgInitializerSrc_.init(img, 2, false);
      imgInitializerDst_.init(img, 2, false);

      // initialize flag field
      bubble_model::BubbleModelTester< Stencil_T >::srcFlagsFromSrcFills();
      bubble_model::BubbleModelTester< Stencil_T >::dstFlagsFromDstFills();

      // initialize bubble model
      bubble_model::BubbleModelTester< Stencil_T >::bubbleModel_->initFromFillLevelField(
         bubble_model::BubbleModelTester< Stencil_T >::dstFillLevelFieldID_);
   }

   // in an alternating manner, update fill levels such that the bubbles must either merge or split on every second call
   void updateDestinationFields() override
   {
      // load image of test scenario with bubbles that are not connected
      geometry::GrayScaleImage unconnected("MergeAndSplitTestUnconnected.png");

      // load image of test scenario with bubbles that are connected
      geometry::GrayScaleImage connected("MergeAndSplitTestConnected.png");

      if (uint_c(calls_) % uint_c(2) == uint_c(0))
      {
         // bubbles must merge (destination field is updated to be connected)
         imgInitializerSrc_.init(unconnected, 2, false);
         imgInitializerDst_.init(connected, 2, false);
      }
      else
      {
         // bubbles must split (destination field is updated to be unconnected)
         imgInitializerSrc_.init(connected, 2, false);
         imgInitializerDst_.init(unconnected, 2, false);
      }

      // update flag field
      bubble_model::BubbleModelTester< Stencil_T >::srcFlagsFromSrcFills();
      bubble_model::BubbleModelTester< Stencil_T >::dstFlagsFromDstFills();

      ++calls_;
   }

 private:
   geometry::initializer::ScalarFieldFromGrayScaleImage imgInitializerSrc_;
   geometry::initializer::ScalarFieldFromGrayScaleImage imgInitializerDst_;
   uint_t calls_;
}; // class ImageMover

template< typename Stencil_T >
void BubbleModelTest< Stencil_T >::testComplexMerge()
{
   auto mpiManager = MPIManager::instance();

   WALBERLA_CHECK_EQUAL(mpiManager->numProcesses(), 10);

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(2), uint_c(5), uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(100), uint_c(200), uint_c(1));
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

   // create imageMover that initializes the fill level field, and bubble model; updates the fill level field in an
   // alternating manner, such that bubbles must either merge or split on every second call
   ImageMover< Stencil_T > imageMover(blockForest, bubbleModel);

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(10));

   // add imageMover to timeloop
   timeloop.addFuncBeforeTimeStep(imageMover, "UpdateDomainFromImage");

   // update bubble model in timeloop
   timeloop.addFuncAfterTimeStep(std::bind(&bubble_model::BubbleModel< Stencil_T >::update, bubbleModel));

   // ensure correctness of initialization (number of bubbles)
   WALBERLA_CHECK_EQUAL(bubbleModel->getBubbles().size(), 12);

   // compute total volume of all bubbles
   real_t volumeBefore = real_c(0);
   for (auto b = bubbleModel->getBubbles().begin(); b != bubbleModel->getBubbles().end(); ++b)
   {
      volumeBefore += b->getCurrentVolume();
   }

   // ensure correctness of initialization (total volume of all bubbles)
   WALBERLA_CHECK_LESS(std::abs(volumeBefore - real_c(8905.51)), 0.1);

   // merge bubbles
   timeloop.singleStep();

   // there must be a single bubble in the system
   WALBERLA_CHECK_EQUAL(bubbleModel->getBubbles().size(), 1);

   // the total volume of this single bubble must be slightly larger than the initial volume of all bubbles
   real_t volumeAfter = bubbleModel->getBubbles()[0].getCurrentVolume();
   WALBERLA_CHECK_LESS(std::abs(volumeAfter - real_c(8918.41)), 0.1);

   // split bubbles
   timeloop.singleStep();

   // the number of bubbles must be as before merging
   WALBERLA_CHECK_EQUAL(bubbleModel->getBubbles().size(), 12);

   // compute total volume of all bubbles
   real_t volumeAfterSplit = real_c(0);
   for (auto b = bubbleModel->getBubbles().begin(); b != bubbleModel->getBubbles().end(); ++b)
   {
      volumeAfterSplit += b->getCurrentVolume();
   }

   // the total volume of all bubbles must be as before merging
   WALBERLA_CHECK_LESS(std::abs(volumeAfterSplit - real_c(8905.51)), 0.1);

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D2Q9 stencil.");
   BubbleModelTest< stencil::D2Q9 >::testComplexMerge();

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D3Q19 stencil.");
   BubbleModelTest< stencil::D3Q19 >::testComplexMerge();

   WALBERLA_LOG_INFO_ON_ROOT("Testing with D3Q27 stencil.");
   BubbleModelTest< stencil::D3Q27 >::testComplexMerge();

   return EXIT_SUCCESS;
}
} // namespace MergeAndSplitTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::MergeAndSplitTest::main(argc, argv); }