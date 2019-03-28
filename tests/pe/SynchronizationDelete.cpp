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
//! \file SynchronizationDelete.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/communication/ParseMessage.h"
#include "pe/basic.h"
#include "pe/synchronization/SyncNextNeighbors.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"

#include <tuple>

namespace walberla {
using namespace walberla::pe;

using BodyTuple = std::tuple<Sphere> ;

void integrate(StructuredBlockForest& forest, BlockDataID storageID, const real_t dt)
{
   for (auto it = forest.begin(); it != forest.end(); ++it)
   {
      blockforest::Block& block = *(dynamic_cast<blockforest::Block*>(&(*it)));
      Storage& storage  = *(block.getData<Storage>(storageID));
      for (auto bd = storage[0].begin(); bd != storage[0].end(); ++bd)
      {
         bd->setPosition(bd->getPosition() + bd->getLinearVel() * dt);
      }
   }
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

//   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
//   logging::Logging::instance()->includeLoggingToFile("SyncLog");

   bool syncShadowOwners = false;
   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--syncShadowOwners" ) == 0 ) syncShadowOwners = true;
   }
   if (syncShadowOwners)
   {
      WALBERLA_LOG_DEVEL("running with syncShadowOwners");
   } else
   {
      WALBERLA_LOG_DEVEL("running with syncNextNeighbour");
   }

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage> ();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 3), uint_c( 3), uint_c( 3), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c(10),                         // dx: length of one cell in physical coordinates
            0,                                  // max blocks per process
            false, false,                       // include metis / force metis
            false, false, false );                 // full periodicity

   WALBERLA_LOG_DETAIL(forest->size() << " blocks located on this process!");

   for (auto it = forest->begin(); it != forest->end(); ++it)
   {
      blockforest::Block* block = dynamic_cast<blockforest::Block*>(&(*it));
      WALBERLA_UNUSED(block);
      WALBERLA_LOG_DETAIL("Block: " << Owner(int_c(block->getProcess()), block->getId().getID()) << "\n AABB: " << block->getAABB());
   }

   SetBodyTypeIDs<BodyTuple>::execute();

   BlockDataID storageID           = forest->addBlockData( createStorageDataHandling<BodyTuple>(), "BodyStorage");

   for (real_t x = real_c(-1); x < real_c(1.5); x += real_c(1.0))
      for (real_t y = real_c(-1); y < real_c(1.5); y += real_c(1.0))
         for (real_t z = real_c(-1); z < real_c(1.5); z += real_c(1.0))
         {
            if ((fabs(x) > 0.1) || (fabs(y) > 0.1) || (fabs(z) > 0.1))
            {
               SphereID sp = createSphere( *globalBodyStorage, forest->getBlockStorage(), storageID, 0, Vec3(15,15,15), 3);
               if (sp != nullptr) sp->setLinearVel(Vec3(x, y, z));
            }
         }

   for (auto it = forest->begin(); it != forest->end(); ++it)
   {
      blockforest::Block& block = *(dynamic_cast<blockforest::Block*>(&(*it)));
      Storage& storage  = *(block.getData<Storage>(storageID));
      if (block.getAABB().contains(Vec3(15,15,15)))
      {
         WALBERLA_CHECK_EQUAL( storage[0].size(), 26);
      } else
      {
         WALBERLA_CHECK_EQUAL( storage[0].size(), 0);
      }
   }

   std::function<void(void)> syncCall;
   if (!syncShadowOwners)
   {
      syncCall = std::bind( pe::syncNextNeighbors<BodyTuple>, std::ref(forest->getBlockForest()), storageID, static_cast<WcTimingTree*>(nullptr), real_c(0.0), false );
   } else
   {
      syncCall = std::bind( pe::syncShadowOwners<BodyTuple>, std::ref(forest->getBlockForest()), storageID, static_cast<WcTimingTree*>(nullptr), real_c(0.0), false );
   }

   for (int i = 0; i < 500; ++i){
      syncCall();
      integrate(*forest, storageID, real_c(0.1));
   }

   for (auto it = forest->begin(); it != forest->end(); ++it)
   {
      blockforest::Block& block = *(dynamic_cast<blockforest::Block*>(&(*it)));
      Storage& storage  = *(block.getData<Storage>(storageID));
      for (auto st = storage.begin(); st != storage.end(); ++st)
      {
         WALBERLA_CHECK_EQUAL(st->size(), 0);
      }
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}