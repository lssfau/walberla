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
//! \file ForceSync.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/basic.h"
#include "pe/ccd/SimpleCCDDataHandling.h"
#include "pe/synchronization/SyncNextNeighbors.h"

#include "pe/synchronization/SyncForces.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"

#include <tuple>

#include <algorithm>
#include <vector>

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Sphere, Plane> BodyTuple ;

int main( int argc, char ** argv )
{
   debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

//   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
//   logging::Logging::instance()->includeLoggingToFile("SyncLog");

   auto forest = blockforest::createUniformBlockGrid(
               uint_c( 2), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
               uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
               real_c(10),                         // dx: length of one cell in physical coordinates
               false,                                  // one block per process
               false, false, false );              // no periodicity

   shared_ptr<BodyStorage> globalStorage = make_shared<BodyStorage>();

//   WALBERLA_LOG_DEVEL("process: " << mpi::MPIManager::instance()->rank() << "\nnumber of blocks: " << forest->size());

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID_           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID1_              = forest->addBlockData(ccd::createSimpleCCDDataHandling( globalStorage, storageID_ ), "CCD1");
   auto ccdID2_              = forest->addBlockData(ccd::createHashGridsDataHandling( globalStorage, storageID_ ), "CCD2");

   walberla::id_t counter = 0;

   auto id1 = pe::createSphere( *globalStorage, forest->getBlockStorage(), storageID_,
                     ++counter, Vec3(8, 3, 3), real_c(1.6))->getSystemID();

   auto id2 = pe::createSphere( *globalStorage, forest->getBlockStorage(), storageID_,
                     ++counter, Vec3(11, 3, 3), real_c(1.6))->getSystemID();


   syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID_);
   syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID_);
   syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID_);

   for (auto it = forest->begin(); it != forest->end(); ++it){
      IBlock & currentBlock = *it;

      ccd::ICCD* ccd1 = currentBlock.getData< ccd::ICCD >( ccdID1_ );
      WALBERLA_LOG_DEVEL(ccd1->generatePossibleContacts().size());

      ccd::ICCD* ccd2 = currentBlock.getData< ccd::ICCD >( ccdID2_ );
      WALBERLA_LOG_DEVEL(ccd2->generatePossibleContacts().size());

      Storage * storage = currentBlock.getData< Storage >( storageID_ );
      BodyStorage& localStorage = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      for (auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt)
      {
         BodyID b = bodyIt.getBodyID();
         b->addForce( Vec3(1,0,0) );
      }
      for (auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); ++bodyIt)
      {
         BodyID b = bodyIt.getBodyID();
         b->addForce( Vec3(0,1,0) );
      }
   }

   reduceForces( forest->getBlockStorage(), storageID_);

   SphereID sp1 = static_cast<SphereID> (getBody(*globalStorage, forest->getBlockStorage(), storageID_, id1));
   SphereID sp2 = static_cast<SphereID> (getBody(*globalStorage, forest->getBlockStorage(), storageID_, id2));

   if ( (sp1 != nullptr) && (!sp1->isRemote() ))
      WALBERLA_ASSERT_FLOAT_EQUAL(sp1->getForce(), Vec3(1,0,0));
   if ( (sp2 != nullptr) && (!sp2->isRemote() ))
      WALBERLA_ASSERT_FLOAT_EQUAL(sp2->getForce(), Vec3(1,1,0));

//   for (auto it = forest->begin(); it != forest->end(); ++it){
//      IBlock & currentBlock = *it;

//      Storage * storage = currentBlock.getData< Storage >( storageID_ );
//      BodyStorage& localStorage = (*storage)[0];
//      BodyStorage& shadowStorage = (*storage)[1];

//      for (auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt)
//      {
//         BodyID b = bodyIt.getBodyID();
//         WALBERLA_LOG_DEVEL("LOCAL\n" << b << "\nForce: " << b->getForce());
//      }
//      for (auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); ++bodyIt)
//      {
//         BodyID b = bodyIt.getBodyID();
//         WALBERLA_LOG_DEVEL("SHADOW\n" << b << "\nForce: " << b->getForce());
//      }
//   }

   forest.reset();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}