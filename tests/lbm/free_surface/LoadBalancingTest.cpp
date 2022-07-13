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
//! \file LoadBalancingTest.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test FSLBM load balancing in simplistic setup with 3x3x1 blocks on a 12x12x1 domain.
//
//! An initially less optimal weight distribution should be increased after performing load balancing.
//======================================================================================================================

#include "lbm/free_surface/LoadBalancing.h"

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/BlockStateDetectorSweep.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/lattice_model/D2Q9.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace free_surface
{
namespace LoadBalancingTest
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using LatticeModel_T = lbm::D2Q9< lbm::collision_model::SRT >;
using Stencil_T      = LatticeModel_T::Stencil;

using Communication_T = blockforest::SimpleCommunication< Stencil_T >;

using FlagField_T                   = FlagField< uint32_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment walberlaEnv(argc, argv);

   // define the domain size
   const Vector3< uint_t > domainSize(uint_c(12), uint_c(12), uint_c(1));
   const Vector3< uint_t > cellsPerBlock(uint_c(3), uint_c(3), uint_c(1));
   const Vector3< uint_t > periodicity(uint_c(0), uint_c(0), uint_c(0));

   Vector3< uint_t > numBlocks;
   numBlocks[0] = uint_c(std::ceil(domainSize[0] / cellsPerBlock[0]));
   numBlocks[1] = uint_c(std::ceil(domainSize[1] / cellsPerBlock[1]));
   numBlocks[2] = uint_c(std::ceil(domainSize[2] / cellsPerBlock[2]));

   uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());

   WALBERLA_CHECK_EQUAL(numProcesses, uint_c(4), "This test must be executed with four MPI processes.")

   WALBERLA_CHECK_LESS_EQUAL(numProcesses, numBlocks[0] * numBlocks[1] * numBlocks[2],
                             "The number of MPI processes is greater than the number of blocks as defined by "
                             "\"domainSize/cellsPerBlock\". This would result in unused MPI processes. Either decrease "
                             "the number of MPI processes or increase \"cellsPerBlock\".")

   // create non-uniform block forest (non-uniformity required for load balancing)
   const std::shared_ptr< StructuredBlockForest > blockForest =
      createNonUniformBlockForest(domainSize, cellsPerBlock, numBlocks, periodicity);

   // create (dummy) lattice model with dummy PDF field
   LatticeModel_T latticeModel  = LatticeModel_T(lbm::collision_model::SRT(real_c(1.0)));
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);

   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(0.0), field::fzyx, uint_c(1));

   // initialize fill level field
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(fillFieldID);
      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         // cell in block-local coordinates
         const Cell localCell = fillFieldIt.cell();

         // get cell in global coordinates
         Cell globalCell = fillFieldIt.cell();
         blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

         // set liquid cells
         if (globalCell[1] < cell_idx_c(5)) { *fillFieldIt = real_c(1); }

         // set interface cells
         if (globalCell[1] == cell_idx_c(5)) { *fillFieldIt = real_c(0.5); }

         // set gas cells
         if (globalCell[1] > cell_idx_c(5)) { *fillFieldIt = real_c(0); }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // create boundary handling for initializing flag field
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();
   freeSurfaceBoundaryHandling->setNoSlipAtAllBorders(cell_idx_c(-1));
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // communication after initialization
   Communication_T communication(blockForest, flagFieldID, fillFieldID);
   communication();

   Communication_T pdfCommunication(blockForest, pdfFieldID);
   pdfCommunication();

   // create bubble model
   std::shared_ptr< bubble_model::BubbleModelBase > bubbleModel =
      std::make_shared< bubble_model::BubbleModelConstantPressure >(real_c(1));

   // detect block states (detection performed during construction)
   BlockStateDetectorSweep< FlagField_T > blockStateDetector(blockForest, flagInfo, flagFieldID);

   // the initialization as chosen above results in the following block states
   // |G|G|G|   with G: onlyGasAndBoundary
   // |F|F|F|        F: fullFreeSurface
   // |F|F|F|        L: onlyLBM
   // |L|L|L|
   //
   // Note that the blocks in row 3 have also state F, although they only consist of gas cells. This is because there is
   // an interface cell in the ghost layer that is synchronized from the blocks of row 2. The BlockStateDetectorSweep
   // also checks the ghost layer, as e.g. during cell conversion, a block having an interface cell in its ghost layer
   // must also perform a corresponding conversion on its inner cells.

   uint_t blockWeightFullFreeSurface    = uint_c(50);
   uint_t blockWeightOnlyLBM            = uint_c(10);
   uint_t blockWeightOnlyGasAndBoundary = uint_c(5);

   // evaluate process loads
   ProcessLoadEvaluator< FlagField_T > loadEvaluator(blockForest, blockWeightFullFreeSurface, blockWeightOnlyLBM,
                                                     blockWeightOnlyGasAndBoundary, uint_c(1));

   // evaluate weight distribution on processes
   std::vector< real_t > weightSum = loadEvaluator.computeWeightSumPerProcess();

   // initial weight distribution before load balancing
   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_LOG_DEVEL("Checking initial weight distribution")
      WALBERLA_CHECK_FLOAT_EQUAL(weightSum[0], real_c(40), real_c(1e-14));
      WALBERLA_CHECK_FLOAT_EQUAL(weightSum[1], real_c(200), real_c(1e-14));
      WALBERLA_CHECK_FLOAT_EQUAL(weightSum[2], real_c(200), real_c(1e-14));
      WALBERLA_CHECK_FLOAT_EQUAL(weightSum[3], real_c(20), real_c(1e-14));
   }

   // perform load balancing
   LoadBalancer< FlagField_T, Stencil_T, Stencil_T > loadBalancer(
      blockForest, communication, pdfCommunication, bubbleModel, blockWeightFullFreeSurface, blockWeightOnlyLBM,
      blockWeightOnlyGasAndBoundary, uint_c(1), false);
   loadBalancer();

   // evaluate weight distribution on processes
   weightSum = loadEvaluator.computeWeightSumPerProcess();

   // check weight distribution after load balancing
   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_LOG_DEVEL("Checking weight distribution after load balancing")
      WALBERLA_CHECK_FLOAT_EQUAL(weightSum[0], real_c(140), real_c(1e-14));
      WALBERLA_CHECK_FLOAT_EQUAL(weightSum[1], real_c(100), real_c(1e-14));
      WALBERLA_CHECK_FLOAT_EQUAL(weightSum[2], real_c(100), real_c(1e-14));
      WALBERLA_CHECK_FLOAT_EQUAL(weightSum[3], real_c(120), real_c(1e-14));
   }

   return EXIT_SUCCESS;
}
} // namespace LoadBalancingTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::LoadBalancingTest::main(argc, argv); }
