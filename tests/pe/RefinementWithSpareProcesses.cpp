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
//! \file RefinementWithSpareProcesses.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/Initialization.h"
#include "blockforest/loadbalancing/DynamicCurve.h"
#include "blockforest/loadbalancing/DynamicDiffusive.h"
#include "blockforest/loadbalancing/DynamicParMetis.h"
#include "blockforest/loadbalancing/PODPhantomData.h"

#include "pe/basic.h"
#include "pe/amr/InfoCollection.h"
#include "pe/amr/level_determination/MinMaxLevelDetermination.h"
#include "pe/amr/weight_assignment/MetisAssignmentFunctor.h"
#include "pe/amr/weight_assignment/WeightAssignmentFunctor.h"
#include "pe/ccd/SimpleCCDDataHandling.h"

#include "core/debug/TestSubsystem.h"
#include "core/grid_generator/SCIterator.h"

#include <tuple>

#include <algorithm>
#include <limits>
#include <vector>

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Sphere> BodyTuple ;

int main( int /*argc*/, char ** /*argv*/, const std::string& LBAlgorithm )
{
   using namespace walberla::pe;

   //      logging::Logging::instance()->setStreamLogLevel( logging::Logging::DETAIL );
   //   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
   //   logging::Logging::instance()->includeLoggingToFile("SyncLog");

   shared_ptr<BodyStorage> globalStorage = make_shared<BodyStorage>();

   // create forest
   auto blockforest = blockforest::createBlockForest( math::AABB(0,0,0,10,10,10),
                                                      Vector3<uint_t>(2,2,2),
                                                      Vector3<bool>(false, false, false) );

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = blockforest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   blockforest->addBlockData(ccd::createHashGridsDataHandling( globalStorage, storageID ), "CCD");
   blockforest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   //***** SETUP LOADBALACING & REFINEMENT
   blockforest->recalculateBlockLevelsInRefresh( true );
   blockforest->alwaysRebalanceInRefresh( true );
   blockforest->reevaluateMinTargetLevelsAfterForcedRefinement( false );
   blockforest->allowRefreshChangingDepth( true );

   blockforest->allowMultipleRefreshCycles( false );
   blockforest->checkForEarlyOutInRefresh( true );
   blockforest->checkForLateOutInRefresh( true );

   auto ic = make_shared<InfoCollection>();

   blockforest->setRefreshMinTargetLevelDeterminationFunction( amr::MinMaxLevelDetermination(ic, 50, 100) );

   if (LBAlgorithm == "Morton")
   {
      blockforest->setRefreshPhantomBlockDataAssignmentFunction( amr::WeightAssignmentFunctor( ic, real_t(1) ) );
      blockforest->setRefreshPhantomBlockDataPackFunction( amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      blockforest->setRefreshPhantomBlockDataUnpackFunction( amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

      auto prepFunc = blockforest::DynamicCurveBalance< amr::WeightAssignmentFunctor::PhantomBlockWeight >( false, true, false );
      blockforest->setRefreshPhantomBlockMigrationPreparationFunction( prepFunc );
   } else if (LBAlgorithm == "Hilbert")
   {
      blockforest->setRefreshPhantomBlockDataAssignmentFunction( amr::WeightAssignmentFunctor( ic, real_t(1) ) );
      blockforest->setRefreshPhantomBlockDataPackFunction( amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      blockforest->setRefreshPhantomBlockDataUnpackFunction( amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

      auto prepFunc = blockforest::DynamicCurveBalance< amr::WeightAssignmentFunctor::PhantomBlockWeight >( true, true, false );
      blockforest->setRefreshPhantomBlockMigrationPreparationFunction( prepFunc );
   } else if (LBAlgorithm == "Metis")
   {
      auto assFunc = amr::MetisAssignmentFunctor( ic, real_t(1) );
      blockforest->setRefreshPhantomBlockDataAssignmentFunction( assFunc );
      blockforest->setRefreshPhantomBlockDataPackFunction( amr::MetisAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      blockforest->setRefreshPhantomBlockDataUnpackFunction( amr::MetisAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

      auto alg     = blockforest::DynamicParMetis::stringToAlgorithm(    "PART_GEOM_KWAY" );
      auto vWeight = blockforest::DynamicParMetis::stringToWeightsToUse( "VERTEX_WEIGHTS" );
      auto eWeight = blockforest::DynamicParMetis::stringToEdgeSource(   "EDGES_FROM_EDGE_WEIGHTS" );

      auto prepFunc = blockforest::DynamicParMetis( alg, vWeight, eWeight );
      prepFunc.setipc2redist(real_t(100000.0));
      blockforest->setRefreshPhantomBlockMigrationPreparationFunction( prepFunc );
   } else if (LBAlgorithm == "Diffusive")
   {
      blockforest->setRefreshPhantomBlockDataAssignmentFunction( amr::WeightAssignmentFunctor( ic, real_t(1) ) );
      blockforest->setRefreshPhantomBlockDataPackFunction( amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      blockforest->setRefreshPhantomBlockDataUnpackFunction( amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      auto prepFunc = blockforest::DynamicDiffusionBalance< amr::WeightAssignmentFunctor::PhantomBlockWeight >( 20, 12, false );
      prepFunc.adaptInflowWithGlobalInformation( true );
      prepFunc.adaptOutflowWithGlobalInformation( true );
      blockforest->setRefreshPhantomBlockMigrationPreparationFunction(prepFunc);
   } else
   {
      WALBERLA_ABORT("Unknown LBAlgorithm: " << LBAlgorithm);
   }

   for (auto& blk : *blockforest)
   {
      for (auto pt : grid_generator::SCGrid(blk.getAABB(), Vector3<real_t>(real_t(0.5), real_t(0.5), real_t(0.5)), real_t(1.0)) )
      {
         createSphere( *globalStorage, *blockforest, storageID, 0, pt, real_t(0.3) );
      }
   }

   WALBERLA_CHECK_GREATER_EQUAL(blockforest->size(), 0);
   WALBERLA_CHECK_LESS(blockforest->size(), 2);

   createWithNeighborhoodLocalShadow( *blockforest, storageID, *ic );
   blockforest->refresh();

   WALBERLA_CHECK_GREATER(blockforest->size(), 0);
   WALBERLA_CHECK_LESS(blockforest->size(), 5);

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   WALBERLA_LOG_DEVEL_ON_ROOT("*** TESTING MORTON ***");
   walberla::main( argc, argv, "Morton" );
   WALBERLA_LOG_DEVEL_ON_ROOT("*** TESTING HILBERT ***");
   walberla::main( argc, argv, "Hilbert" );
   WALBERLA_LOG_DEVEL_ON_ROOT("*** TESTING DIFFUSIVE ***");
   walberla::main( argc, argv, "Diffusive" );
#ifdef WALBERLA_BUILD_WITH_PARMETIS
   WALBERLA_LOG_DEVEL_ON_ROOT("*** TESTING PARMETIS ***");
   walberla::main( argc, argv, "Metis" );
#endif

   return EXIT_SUCCESS;
}
