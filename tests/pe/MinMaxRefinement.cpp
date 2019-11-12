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
//! \file MinMaxRefinement.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "blockforest/loadbalancing/InfoCollection.h"
#include <blockforest/loadbalancing/PODPhantomData.h>
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "timeloop/SweepTimeloop.h"
#include "vtk/VTKOutput.h"


#include "pe/basic.h"
#include "pe/amr/level_determination/MinMaxLevelDetermination.h"
#include "pe/amr/weight_assignment/WeightAssignmentFunctor.h"
#include "pe/ccd/SimpleCCDDataHandling.h"
#include "pe/synchronization/SyncNextNeighbors.h"
#include "pe/synchronization/ClearSynchronization.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"

#include <tuple>

#include <algorithm>
#include <limits>
#include <vector>

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Sphere, Plane> BodyTuple ;

int main( int argc, char ** argv )
{
   using namespace walberla::pe;

   debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   //      logging::Logging::instance()->setStreamLogLevel( logging::Logging::DETAIL );
   //   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
   //   logging::Logging::instance()->includeLoggingToFile("SyncLog");

   shared_ptr<BodyStorage> globalStorage = make_shared<BodyStorage>();

   // create forest
   shared_ptr< blockforest::StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            math::AABB(0,0,0,4,4,4),
            1,1,1,                                                 // number of blocks in x,y,z direction
            1,1,1,                                                  // how many cells per block (x,y,z)
            0,                                                      // max blocks per process
            false, false,                                           // include metis / force metis
            false, false, false );                                    // full periodicity

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");
   WALBERLA_UNUSED(fcdID);

   auto & blockforest = forest->getBlockForest();

   //***** SETUP LOADBALACING & REFINEMENT
   blockforest.recalculateBlockLevelsInRefresh( true );
   blockforest.alwaysRebalanceInRefresh( true );
   blockforest.reevaluateMinTargetLevelsAfterForcedRefinement( false );
   blockforest.allowRefreshChangingDepth( true );

   blockforest.allowMultipleRefreshCycles( false );
   blockforest.checkForEarlyOutInRefresh( true );
   blockforest.checkForLateOutInRefresh( true );

   auto infoCollection = make_shared<blockforest::InfoCollection>();

   amr::MinMaxLevelDetermination levelDetermination(infoCollection, 2, 5);
   blockforest.setRefreshMinTargetLevelDeterminationFunction( levelDetermination );

   blockforest.setRefreshPhantomBlockDataAssignmentFunction( amr::WeightAssignmentFunctor( infoCollection ) );
   blockforest.setRefreshPhantomBlockDataPackFunction( amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
   blockforest.setRefreshPhantomBlockDataUnpackFunction( amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

   blockforest.setRefreshPhantomBlockMigrationPreparationFunction(
            blockforest::DynamicCurveBalance< amr::WeightAssignmentFunctor::PhantomBlockWeight >( false, true, false ) );

   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(1,1,1), 1);
   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(1,1,3), 1);
   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(1,3,1), 1);
   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(1,3,3), 1);
   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(3,1,1), 1);
   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(3,1,3), 1);
   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(3,3,1), 1);
   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(3,3,3), 1);

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT( "Refinement 1" );
   createWithNeighborhoodLocalShadow(blockforest, storageID, *infoCollection);
   clearSynchronization( blockforest, storageID);
   forest->refresh();
   syncNextNeighbors<BodyTuple>(blockforest, storageID);

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
      ccd->reloadBodies();
   }

   WALBERLA_CHECK_EQUAL( blockforest.size(), 1);

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT( "Refinement 2" );
   blockforest.setRefreshMinTargetLevelDeterminationFunction( amr::MinMaxLevelDetermination(infoCollection, 9, 20) );
   createWithNeighborhoodLocalShadow(blockforest, storageID, *infoCollection);
   clearSynchronization( blockforest, storageID);
   forest->refresh();
   syncNextNeighbors<BodyTuple>(blockforest, storageID);

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
      ccd->reloadBodies();
   }

   WALBERLA_CHECK_EQUAL( blockforest.size(), mpi::MPIManager::instance()->worldRank() == 0 ? 1 : 0);
   WALBERLA_LOG_DEVEL( infoCollection->size() );

   for (unsigned int i = 0; i < 30; ++i)
   {
      createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(real_t(2.1), real_t(2.1), real_t(2.1)), 1);
   }

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT( "Refinement 3" );
   blockforest.setRefreshMinTargetLevelDeterminationFunction( amr::MinMaxLevelDetermination(infoCollection, 2, 3) );
   createWithNeighborhoodLocalShadow(blockforest, storageID, *infoCollection);
   clearSynchronization( blockforest, storageID);
   forest->refresh();
   syncNextNeighbors<BodyTuple>(blockforest, storageID);

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
      ccd->reloadBodies();
   }

   WALBERLA_LOG_DEVEL( infoCollection->size() );

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT( "Refinement 4" );
   createWithNeighborhoodLocalShadow(blockforest, storageID, *infoCollection);
   clearSynchronization( blockforest, storageID);
   forest->refresh();
   syncNextNeighbors<BodyTuple>(blockforest, storageID);

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
      ccd->reloadBodies();
   }

   WALBERLA_LOG_DEVEL( infoCollection->size() );

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT( "Refinement 5" );
   WALBERLA_LOG_DEVEL( "SIZE: " << blockforest.size() );
   createWithNeighborhoodLocalShadow(blockforest, storageID, *infoCollection);
   clearSynchronization( blockforest, storageID);
   forest->refresh();
   syncNextNeighbors<BodyTuple>(blockforest, storageID);

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
      ccd->reloadBodies();
   }

   WALBERLA_LOG_DEVEL( infoCollection->size() );

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
