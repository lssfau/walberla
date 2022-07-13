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
//! \file LoadBalancing.h
//! \ingroup free_surface
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Free surface-specific functionality for load balancing.
//
//======================================================================================================================

#include "blockforest/BlockForest.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/loadbalancing/DynamicCurve.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "core/logging/Logging.h"
#include "core/math/AABB.h"
#include "core/math/DistributedSample.h"
#include "core/math/Vector3.h"
#include "core/mpi/MPIManager.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/free_surface/BlockStateDetectorSweep.h"
#include "lbm/free_surface/bubble_model/BubbleModel.h"

#include <algorithm>
#include <numeric>

namespace walberla
{
namespace free_surface
{

template< typename FlagField_T >
class ProcessLoadEvaluator;

/***********************************************************************************************************************
 * Create non-uniform block forest to be used for load balancing.
 **********************************************************************************************************************/
std::shared_ptr< StructuredBlockForest > createNonUniformBlockForest(const Vector3< uint_t >& domainSize,
                                                                     const Vector3< uint_t >& cellsPerBlock,
                                                                     const Vector3< uint_t >& numBlocks,
                                                                     const Vector3< bool >& periodicity)
{
   WALBERLA_CHECK_EQUAL(domainSize[0], cellsPerBlock[0] * numBlocks[0],
                        "The domain size is not divisible by the specified \"cellsPerBlock\" in x-direction.");
   WALBERLA_CHECK_EQUAL(domainSize[1], cellsPerBlock[1] * numBlocks[1],
                        "The domain size is not divisible by the specified \"cellsPerBlock\" in y-direction.");
   WALBERLA_CHECK_EQUAL(domainSize[2], cellsPerBlock[2] * numBlocks[2],
                        "The domain size is not divisible by the specified \"cellsPerBlock\" in z-direction.");

   // create SetupBlockForest for allowing load balancing
   SetupBlockForest setupBlockForest;

   AABB domainAABB(real_c(0), real_c(0), real_c(0), real_c(domainSize[0]), real_c(domainSize[1]),
                   real_c(domainSize[2]));

   setupBlockForest.init(domainAABB, numBlocks[0], numBlocks[1], numBlocks[2], periodicity[0], periodicity[1],
                         periodicity[2]);

   // compute initial process distribution
   setupBlockForest.balanceLoad(blockforest::StaticLevelwiseCurveBalance(true),
                                uint_c(MPIManager::instance()->numProcesses()));

   WALBERLA_LOG_INFO_ON_ROOT(setupBlockForest);

   // define MPI communicator
   if (!MPIManager::instance()->rankValid()) { MPIManager::instance()->useWorldComm(); }

   // create BlockForest (will be encapsulated in StructuredBlockForest)
   const std::shared_ptr< BlockForest > blockForest =
      std::make_shared< BlockForest >(uint_c(MPIManager::instance()->rank()), setupBlockForest, false);

   // create StructuredBlockForest
   std::shared_ptr< StructuredBlockForest > structuredBlockForest =
      std::make_shared< StructuredBlockForest >(blockForest, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
   structuredBlockForest->createCellBoundingBoxes();

   return structuredBlockForest;
}

/***********************************************************************************************************************
 * Example for a load balancing implementation.
 *
 * IMPORTANT REMARK: The following implementation should be considered a demonstrator based on best-practices. That is,
 * it was not thoroughly benchmarked and does not guarantee a performance-optimal load balancing.
 **********************************************************************************************************************/
template< typename FlagField_T, typename CommunicationStencil_T, typename LatticeModelStencil_T >
class LoadBalancer
{
 public:
   LoadBalancer(const std::shared_ptr< StructuredBlockForest >& blockForestPtr,
                const blockforest::SimpleCommunication< CommunicationStencil_T >& communication,
                const blockforest::SimpleCommunication< LatticeModelStencil_T >& pdfCommunication,
                const std::shared_ptr< bubble_model::BubbleModelBase >& bubbleModel, uint_t blockWeightFullFreeSurface,
                uint_t blockWeightOnlyLBM, uint_t blockWeightOnlyGasAndBoundary, uint_t frequency,
                bool printStatistics = false)
      : blockForest_(blockForestPtr), communication_(communication), pdfCommunication_(pdfCommunication),
        bubbleModel_(bubbleModel), blockWeightFullFreeSurface_(blockWeightFullFreeSurface),
        blockWeightOnlyLBM_(blockWeightOnlyLBM), blockWeightOnlyGasAndBoundary_(blockWeightOnlyGasAndBoundary),
        frequency_(frequency), printStatistics_(printStatistics), executionCounter_(uint_c(0)),
        evaluator_(ProcessLoadEvaluator< FlagField_T >(blockForest_, blockWeightFullFreeSurface_, blockWeightOnlyLBM_,
                                                       blockWeightOnlyGasAndBoundary_, uint_c(1)))
   {
      BlockForest& blockForest = blockForest_->getBlockForest();

      // refinement is not implemented in FSLBM such that this can be set to false
      blockForest.recalculateBlockLevelsInRefresh(false);

      // rebalancing of blocks must be forced here, as it would normally be done when refinement levels change
      blockForest.alwaysRebalanceInRefresh(true);

      // depth of levels is not changing and must therefore not be communicated
      blockForest.allowRefreshChangingDepth(false);

      // refinement is not implemented in FSLBM such that this can be set to false
      blockForest.allowMultipleRefreshCycles(false);

      // leave algorithm when load balancing has been performed
      blockForest.checkForEarlyOutAfterLoadBalancing(true);

      // use PhantomWeight as defined below
      blockForest.setRefreshPhantomBlockDataAssignmentFunction(blockWeightAssignment);
      blockForest.setRefreshPhantomBlockDataPackFunction(phantomWeightsPack);
      blockForest.setRefreshPhantomBlockDataUnpackFunction(phantomWeightsUnpack);

      // assign load balancing function
      blockForest.setRefreshPhantomBlockMigrationPreparationFunction(
         blockforest::DynamicCurveBalance< PhantomWeight >(true, true));
   }

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      // only balance load in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {
         BlockForest& blockForest = blockForest_->getBlockForest();

         // balance load by updating the blockForest
         const uint_t modificationStamp = blockForest.getModificationStamp();
         blockForest.refresh();

         const uint_t newModificationStamp = blockForest.getModificationStamp();

         if (newModificationStamp != modificationStamp)
         {
            // communicate all fields
            communication_();
            pdfCommunication_();
            bubbleModel_->update();

            if (printStatistics_) { evaluator_(); }

            if (blockForest.getNumberOfBlocks() == uint_c(0))
            {
               WALBERLA_ABORT(
                  "Load balancing lead to a situation where there is a process with no blocks. This is "
                  "not supported yet. This can be avoided by either using smaller blocks or, equivalently, more blocks "
                  "per process.");
            }
         }
      }
      ++executionCounter_;
   }

 private:
   std::shared_ptr< StructuredBlockForest > blockForest_;
   blockforest::SimpleCommunication< CommunicationStencil_T > communication_;
   blockforest::SimpleCommunication< LatticeModelStencil_T > pdfCommunication_;
   std::shared_ptr< bubble_model::BubbleModelBase > bubbleModel_;

   uint_t blockWeightFullFreeSurface_;    // empirical choice, not thoroughly benchmarked: 50
   uint_t blockWeightOnlyLBM_;            // empirical choice, not thoroughly benchmarked: 10
   uint_t blockWeightOnlyGasAndBoundary_; // empirical choice, not thoroughly benchmarked: 5

   uint_t frequency_;
   bool printStatistics_;

   uint_t executionCounter_;

   ProcessLoadEvaluator< FlagField_T > evaluator_;

   class PhantomWeight // used as a 'PhantomBlockForest::PhantomBlockDataAssignmentFunction'
   {
    public:
      using weight_t = uint_t;
      PhantomWeight(const weight_t _weight) : weight_(_weight) {}
      weight_t weight() const { return weight_; }

    private:
      weight_t weight_;
   }; // class PhantomWeight

   std::function< void(mpi::SendBuffer& buffer, const PhantomBlock& block) > phantomWeightsPack =
      [](mpi::SendBuffer& buffer, const PhantomBlock& block) { buffer << block.getData< PhantomWeight >().weight(); };

   std::function< void(mpi::RecvBuffer& buffer, const PhantomBlock&, walberla::any& data) > phantomWeightsUnpack =
      [](mpi::RecvBuffer& buffer, const PhantomBlock&, walberla::any& data) {
         typename PhantomWeight::weight_t w;
         buffer >> w;
         data = PhantomWeight(w);
      };

   std::function< void(std::vector< std::pair< const PhantomBlock*, walberla::any > >& blockData,
                       const PhantomBlockForest&) >
      blockWeightAssignment =
         [this](std::vector< std::pair< const PhantomBlock*, walberla::any > >& blockData, const PhantomBlockForest&) {
            for (auto it = blockData.begin(); it != blockData.end(); ++it)
            {
               if (it->first->getState().contains(BlockStateDetectorSweep< FlagField_T >::fullFreeSurface))
               {
                  it->second = PhantomWeight(blockWeightFullFreeSurface_);
               }
               else
               {
                  if (it->first->getState().contains(BlockStateDetectorSweep< FlagField_T >::onlyLBM))
                  {
                     it->second = PhantomWeight(blockWeightOnlyLBM_);
                  }
                  else
                  {
                     if (it->first->getState().contains(BlockStateDetectorSweep< FlagField_T >::onlyGasAndBoundary))
                     {
                        it->second = PhantomWeight(blockWeightOnlyGasAndBoundary_);
                     }
                     else { WALBERLA_ABORT("Unknown block state"); }
                  }
               }
            }
         };

}; // class LoadBalancer

/***********************************************************************************************************************
 * Evaluates and prints statistics about the current load distribution situation:
 * - Average weight per process
 * - Maximum weight per process
 * - Minimum weight per process
 **********************************************************************************************************************/
template< typename FlagField_T >
class ProcessLoadEvaluator
{
 public:
   ProcessLoadEvaluator(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                        uint_t blockWeightFullFreeSurface, uint_t blockWeightOnlyLBM,
                        uint_t blockWeightOnlyGasAndBoundary, uint_t frequency)
      : blockForest_(blockForest), blockWeightFullFreeSurface_(blockWeightFullFreeSurface),
        blockWeightOnlyLBM_(blockWeightOnlyLBM), blockWeightOnlyGasAndBoundary_(blockWeightOnlyGasAndBoundary),
        frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % frequency_ != uint_c(0) && executionCounter_ != uint_c(1)) { return; }

      std::vector< real_t > weightSum = computeWeightSumPerProcess();

      print(weightSum);
   }

   std::vector< real_t > computeWeightSumPerProcess()
   {
      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      std::vector< real_t > weightSum(uint_c(MPIManager::instance()->numProcesses()), real_c(0));

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         if (blockForest->blockExistsLocally(blockIt->getId()))
         {
            if (blockIt->getState().contains(BlockStateDetectorSweep< FlagField_T >::fullFreeSurface))
            {
               weightSum[blockForest->getProcessRank(blockIt->getId())] += real_c(blockWeightFullFreeSurface_);
            }
            else
            {
               if (blockIt->getState().contains(BlockStateDetectorSweep< FlagField_T >::onlyLBM))
               {
                  weightSum[blockForest->getProcessRank(blockIt->getId())] += real_c(blockWeightOnlyLBM_);
               }
               else
               {
                  if (blockIt->getState().contains(BlockStateDetectorSweep< FlagField_T >::onlyGasAndBoundary))
                  {
                     weightSum[blockForest->getProcessRank(blockIt->getId())] += real_c(blockWeightOnlyGasAndBoundary_);
                  }
               }
            }
         }
      }

      mpi::reduceInplace< real_t >(weightSum, mpi::SUM, 0);

      return weightSum;
   }

   void print(const std::vector< real_t >& weightSum)
   {
      WALBERLA_ROOT_SECTION()
      {
         const std::vector< real_t >::const_iterator max = std::max_element(weightSum.cbegin(), weightSum.end());
         const std::vector< real_t >::const_iterator min = std::min_element(weightSum.cbegin(), weightSum.end());
         const real_t sum = std::accumulate(weightSum.cbegin(), weightSum.end(), real_c(0));
         const real_t avg = sum / real_c(MPIManager::instance()->numProcesses());

         WALBERLA_LOG_INFO("Load balancing:");
         WALBERLA_LOG_INFO("\t Average weight per process " << avg);
         WALBERLA_LOG_INFO("\t Maximum weight per process " << *max);
         WALBERLA_LOG_INFO("\t Minimum weight per process " << *min);
      }
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   uint_t blockWeightFullFreeSurface_;
   uint_t blockWeightOnlyLBM_;
   uint_t blockWeightOnlyGasAndBoundary_;

   uint_t frequency_;
   uint_t executionCounter_;
}; // class ProcessLoadEvaluator

} // namespace free_surface
} // namespace walberla
