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
//! \file Callback.cpp
//! \brief checks callbacks of BodyStorage
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"
#include <pe/utility/DestroyBody.h>

#include "blockforest/Initialization.h"
#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Sphere> BodyTypeTuple ;

enum class State{ LOCALIZED0, SHADOW, MIGRATED, LOCALIZED1, REMOVED};
State state;

void addCallbackLocal(BodyID bd)
{
   bool shouldBeCalled = ((state == State::LOCALIZED0) && (mpi::MPIManager::instance()->worldRank() == 0)) ||
                         ((state == State::MIGRATED) && (mpi::MPIManager::instance()->worldRank() == 1));
   WALBERLA_CHECK(shouldBeCalled);
   WALBERLA_LOG_DEVEL("Add local body: " << bd->getSystemID() );
}
void removeCallbackLocal(BodyID bd)
{
   bool shouldBeCalled = ((state == State::MIGRATED) && (mpi::MPIManager::instance()->worldRank() == 0)) ||
                         ((state == State::REMOVED) && (mpi::MPIManager::instance()->worldRank() == 1));
   WALBERLA_CHECK(shouldBeCalled);
   WALBERLA_LOG_DEVEL("Remove local body: " << bd->getSystemID() );
}

void addCallbackShadow(BodyID bd)
{
   bool shouldBeCalled = ((state == State::SHADOW) && (mpi::MPIManager::instance()->worldRank() == 1)) ||
                         ((state == State::MIGRATED) && (mpi::MPIManager::instance()->worldRank() == 0));
   WALBERLA_CHECK(shouldBeCalled);
   WALBERLA_LOG_DEVEL("Add shadow body: " << bd->getSystemID() );
}
void removeCallbackShadow(BodyID bd)
{
   bool shouldBeCalled = ((state == State::MIGRATED) && (mpi::MPIManager::instance()->worldRank() == 1)) ||
                         ((state == State::LOCALIZED1) && (mpi::MPIManager::instance()->worldRank() == 0));
   WALBERLA_CHECK(shouldBeCalled);
   WALBERLA_LOG_DEVEL("Remove shadow body: " << bd->getSystemID() );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   auto forest = blockforest::createBlockForest( AABB(0,0,0,20,20,20), // simulation domain
                                                 Vector3<uint_t>(2,1,1), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );

   SetBodyTypeIDs<BodyTypeTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTypeTuple>(), "Storage");

   for (auto& currentBlock : *forest)
   {
      Storage * storage = currentBlock.getData< Storage >( storageID );
      BodyStorage& localStorage = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      localStorage.registerAddCallback( "Test", addCallbackLocal );
      localStorage.registerRemoveCallback( "Test", removeCallbackLocal );

      shadowStorage.registerAddCallback( "Test", addCallbackShadow );
      shadowStorage.registerRemoveCallback( "Test", removeCallbackShadow );
   }

   state = State::LOCALIZED0;
   BodyID sp = pe::createSphere(
                  *globalBodyStorage,
                  *forest,
                  storageID,
                  99,
                  Vec3(8,5,5),
                  real_c(1.0));
   if (sp != nullptr)
   {
      WALBERLA_CHECK_EQUAL( sp->getSystemID(), 1 );
   }
   syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT("=== localized on block 1 ===");
   WALBERLA_MPI_BARRIER();

   state = State::SHADOW;
   sp = getBody(*globalBodyStorage, *forest, storageID, 1);
   if (sp != nullptr)
   {
      sp->setPosition( Vec3(real_t(9.5), 5, 5) );
   }
   syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT("=== shadow on block 2 ===");
   WALBERLA_MPI_BARRIER();

   state = State::MIGRATED;
   sp = getBody(*globalBodyStorage, *forest, storageID, 1);
   if (sp != nullptr)
   {
      sp->setPosition( Vec3(real_t(10.5), 5, 5) );
   }
   syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT("=== migrated to block 2 ===");
   WALBERLA_MPI_BARRIER();

   state = State::LOCALIZED1;
   sp = getBody(*globalBodyStorage, *forest, storageID, 1);
   if (sp != nullptr)
   {
      sp->setPosition( Vec3(real_t(11.5), 5, 5) );
   }
   syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT("=== localized on block 2 ===");
   WALBERLA_MPI_BARRIER();

   state = State::REMOVED;
   destroyBodyBySID(*globalBodyStorage, *forest, storageID, 1);
   syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT("=== removed ===");
   WALBERLA_MPI_BARRIER();

   for (auto& currentBlock : *forest)
   {
      Storage * storage = currentBlock.getData< Storage >( storageID );
      BodyStorage& localStorage = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      localStorage.clearAddCallbacks();
      localStorage.clearRemoveCallbacks();

      shadowStorage.clearAddCallbacks();
      shadowStorage.clearRemoveCallbacks();
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
