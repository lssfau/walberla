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
//! \file ParallelEquivalence.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "timeloop/SweepTimeloop.h"
#include "vtk/VTKOutput.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "pe/basic.h"
#include "pe/ccd/SimpleCCDDataHandling.h"
#include "pe/synchronization/SyncNextNeighbors.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"

#include <tuple>

#include <algorithm>
#include <vector>

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Sphere, Plane> BodyTuple ;

class SimpleLB
{
public:
   SimpleLB( const weak_ptr< blockforest::StructuredBlockForest > & forest ) :
      forest_( forest )
   {}

   void operator()()
   {
      WALBERLA_LOG_INFO_ON_ROOT("call refresh");
      auto forest = forest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( forest );
      forest->getBlockForest().refresh();
   }

   bool operator()( std::vector< std::pair< const blockforest::PhantomBlock *, uint_t > > & targetProcess,
                    std::set< uint_t > & processesToRecvFrom,
                    const blockforest::PhantomBlockForest & /*phantomForest*/,
                    const uint_t /*iteration*/ ) const
   {
      WALBERLA_LOG_INFO_ON_ROOT("setting block ranks");
      auto rank = mpi::MPIManager::instance()->worldRank();
      auto numProcesses = mpi::MPIManager::instance()->numProcesses();
      uint_t recv = uint_c((rank < 1) ? numProcesses - 1  : rank - 1);
      uint_t send = uint_c(((rank + 1) >= numProcesses) ? 0  : rank + 1);
      WALBERLA_LOG_INFO("recv: " << recv << "\n send: " << send);
      processesToRecvFrom.insert( recv );
      for( auto it = targetProcess.begin(); it != targetProcess.end(); ++it)
      {
         it->second = send;
      }
      return false;
   }

private:

   weak_ptr< blockforest::StructuredBlockForest > forest_;
};

class RestoreRelations
{
public:
   RestoreRelations(const BlockDataID& storageID)
      : storageID_(storageID)
   {}

   void operator()( BlockForest & forest, const blockforest::PhantomBlockForest & /*phantomForest*/ )
   {
      for (auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt)
      {
         Block * block = dynamic_cast< Block * >( &(*blockIt) );
         for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID_); bodyIt != LocalBodyIterator::end(); ++bodyIt)
         {
            bodyIt->MPITrait.setOwner( Owner( int_c(block->getTargetProcess()[0]), block->getId().getID()));
         }
      }

   }
private:
   BlockDataID storageID_;
};



int main( int argc, char ** argv )
{
   using namespace walberla::pe;

   debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   logging::Logging::instance()->setStreamLogLevel( logging::Logging::DETAIL );
//   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
//   logging::Logging::instance()->includeLoggingToFile("SyncLog");

   shared_ptr<BodyStorage> globalStorage = make_shared<BodyStorage>();

   // create forest
   shared_ptr< blockforest::StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            math::AABB(0,0,0,20,10,10),
            2,1,1,                                                 // number of blocks in x,y,z direction
            1,1,1,                                                  // how many cells per block (x,y,z)
            0,                                                      // max blocks per process
            false, false,                                           // include metis / force metis
            false, false, false );                                    // full periodicity

//      vtk::writeDomainDecomposition( forest );

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   //cr::DEM    cr(globalStorage, forest->getBlockStorage(), storageID, ccdID, fcdID, NULL );
   cr::HCSITS cr(globalStorage, forest->getBlockStoragePointer(), storageID, ccdID, fcdID, nullptr );

   auto vtkOutput   = make_shared<DefaultBodyVTKOutput>(storageID, forest->getBlockStorage()) ;
   auto vtkWriter   = vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, "vtk", "simulation_step", false, false);

   auto & blockforest = forest->getBlockForest();

   blockforest.recalculateBlockLevelsInRefresh( true );
   blockforest.alwaysRebalanceInRefresh( true );
   blockforest.reevaluateMinTargetLevelsAfterForcedRefinement( false );
   blockforest.allowRefreshChangingDepth( false );

   blockforest.allowMultipleRefreshCycles( false );
   blockforest.checkForEarlyOutInRefresh( true );
   blockforest.checkForLateOutInRefresh( true );
   blockforest.addRefreshCallbackFunctionBeforeBlockDataIsPacked( RestoreRelations(storageID) );

   SimpleLB simpleLB( forest );

   //blockforest.setRefreshMinTargetLevelDeterminationFunction( regrid );

   blockforest.setRefreshPhantomBlockMigrationPreparationFunction( simpleLB );

   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(5,5,5), 1);
   createSphere(*globalStorage.get(), forest->getBlockStorage(), storageID, 0, Vec3(15,6,6), 1);

   SweepTimeloop timeloop( forest->getBlockStorage(), 1 );
   timeloop.addFuncBeforeTimeStep( simpleLB, "refreshFunctorName" );

   for (int i = 0; i < 1; ++i)
   {
      WALBERLA_LOG_INFO_ON_ROOT("timestep: " << i);
      timeloop.singleStep();
   }

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt)
      {
         WALBERLA_LOG_DEVEL(*bodyIt);
      }
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}