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
//! \file Synchronization.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/basic.h"
#include "pe/communication/ParseMessage.h"
#include "pe/synchronization/SyncNextNeighbors.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"

#include <tuple>

namespace walberla {
using namespace walberla::pe;

using BodyTuple = std::tuple<Sphere> ;

void checkSphere(StructuredBlockForest& forest, BlockDataID storageID, walberla::id_t sid, Sphere& ref, const Vec3& newPos)
{
   for (auto it = forest.begin(); it != forest.end(); ++it)
   {
      blockforest::Block& block = *(dynamic_cast<blockforest::Block*>(&(*it)));
      Storage& storage  = *(block.getData<Storage>(storageID));
      BodyStorage& shadowStorage = storage[1];

      if (block.getAABB().contains( ref.getPosition() ))
      {
         WALBERLA_CHECK_EQUAL( storage[StorageType::LOCAL].size(), 1, "pos: " <<  ref.getPosition() << "\nradius: " << ref.getRadius() <<  "\ndomain: " << block.getAABB() );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 0, "pos: " << ref.getPosition() << "\nradius: " << ref.getRadius() <<  "\ndomain: " << block.getAABB() );
         SphereID bd = static_cast<SphereID> (storage[StorageType::LOCAL].find( sid ).getBodyID());
         WALBERLA_CHECK_NOT_NULLPTR(bd);
         checkVitalParameters(bd, &ref);
         WALBERLA_LOG_DEVEL("#shadows: " << bd->MPITrait.sizeShadowOwners() << " #block states set: " << bd->MPITrait.getBlockStateSize() << "\nowner domain: " << block.getAABB() << "\nowner: " << bd->MPITrait.getOwner());
         bd->setPosition( newPos );
      } else if (forest.periodicIntersect(block.getAABB(), ref.getAABB()) )
      {
         WALBERLA_CHECK_EQUAL( storage[StorageType::LOCAL].size(), 0, "pos: " << ref.getPosition() << "\nradius: " << ref.getRadius() <<  "\ndomain: " << block.getAABB() );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 1, "pos: " << ref.getPosition() << "\nradius: " << ref.getRadius() <<  "\ndomain: " << block.getAABB() );
         SphereID bd = static_cast<SphereID> (shadowStorage.find( sid ).getBodyID());
         WALBERLA_CHECK_NOT_NULLPTR(bd);
         auto backupPos =ref.getPosition();
         auto correctedPos = ref.getPosition();
         pe::communication::correctBodyPosition(forest.getDomain(), block.getAABB().center(), correctedPos);
         ref.setPosition(correctedPos);
         checkVitalParameters(bd, &ref);
         ref.setPosition(backupPos);
      } else
      {
         WALBERLA_CHECK_EQUAL( storage[StorageType::LOCAL].size(), 0, "pos: " << ref.getPosition() << "\nradius: " << ref.getRadius() <<  "\ndomain: " << block.getAABB() );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 0, "pos: " << ref.getPosition() << "\nradius: " << ref.getRadius() <<  "\ndomain: " << block.getAABB() );
      }
   }
   WALBERLA_LOG_PROGRESS("checked pos: " << ref.getPosition() << " | new pos: " << newPos);
   auto temp = newPos;
   forest.mapToPeriodicDomain(temp);
   ref.setPosition(temp);
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

//   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
//   logging::Logging::instance()->includeLoggingToFile("SyncLog");

   BodyStorage globalStorage;

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 3), uint_c( 3), uint_c( 3), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c(10),                         // dx: length of one cell in physical coordinates
            0,                                  // max blocks per process
            false, false,                       // include metis / force metis
            true, true, true );                 // full periodicity

   WALBERLA_LOG_DETAIL(forest->size() << " blocks located on this process!");

   for (auto it = forest->begin(); it != forest->end(); ++it)
   {
      blockforest::Block* block = dynamic_cast<blockforest::Block*>(&(*it));
      WALBERLA_UNUSED(block);
      WALBERLA_LOG_DETAIL("Block: " << Owner(int_c(block->getProcess()), block->getId().getID()) << "\n AABB: " << block->getAABB());
   }

   int rank = MPIManager::instance()->rank();

   SetBodyTypeIDs<BodyTuple>::execute();

   BlockDataID storageID           = forest->addBlockData( createStorageDataHandling<BodyTuple>(), "BodyStorage");

   MaterialID iron = Material::find("iron");
   walberla::id_t sid = 123;
   Sphere refSphere(1, 0, Vec3(15, 15, 15), Quat(), 3, iron, false, true, false);
   refSphere.setLinearVel(4, 5, 6);
   refSphere.setAngularVel( 1, 2, 3);
   Vec3 gpos = Vec3(15, 15, 15);

   SphereID sphere = createSphere( globalStorage, forest->getBlockStorage(), storageID, 0, gpos, 3);
   int sphereRank = -1;

   if (sphere != nullptr)
   {
      sphere->setLinearVel(4, 5, 6);
      sphere->setAngularVel( 1, 2, 3);
      sid = sphere->getSystemID();
      sphereRank = rank;
   }

   auto rankVec = mpi::allGather(sphereRank);
   for (auto it = rankVec.begin(); it != rankVec.end(); ++it)
   {
      if (*it != -1)
      {
         sphereRank = *it;
         break;
      }
   }
   mpi::broadcastObject(sid, sphereRank);
   WALBERLA_LOG_DETAIL("sphere with sid " << sid << " is located on rank " << sphereRank);

   WALBERLA_LOG_PROGRESS("*********************** [1 1 1] TEST ***********************");
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(19,19,19));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(21,21,21));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(25,25,25));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(29,29,29));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(31,31,31));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 5, 5, 5));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 9, 9, 9));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(11,11,11));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));


   WALBERLA_LOG_PROGRESS("*********************** [-1 1 1] TEST ***********************");
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(11,19,19));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 9,21,21));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 5,25,25));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 1,29,29));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(-1,31,31));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(25,05, 5));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(21, 9, 9));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(19,11,11));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));


   WALBERLA_LOG_PROGRESS("*********************** [-1 -1 1] TEST ***********************");
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(11,11,19));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 9, 9,21));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 5, 5,25));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 1, 1,29));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(-1,-1,31));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(25,25, 5));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(21,21, 9));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(19,19,11));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));


   WALBERLA_LOG_PROGRESS("*********************** [0 1 1] TEST ***********************");
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,19,19));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,21,21));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,25,25));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,29,29));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,31,31));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15, 5, 5));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15, 9, 9));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,11,11));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));


   WALBERLA_LOG_PROGRESS("*********************** [0 0 1] TEST ***********************");
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,19));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,21));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,25));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,29));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,31));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15, 5));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15, 9));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,11));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));

   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);
   syncNextNeighbors<BodyTuple>(forest->getBlockForest(), storageID);

   //*****************************************************************************************

   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   WALBERLA_LOG_PROGRESS("*********************** [1 1 1] TEST ***********************");
//   WALBERLA_LOG_DEVEL("BIG SYNC 0")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(19,19,19));

//   WALBERLA_LOG_DEVEL("BIG SYNC 1")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(21,21,21));

//   WALBERLA_LOG_DEVEL("BIG SYNC 2")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(25,25,25));

//   WALBERLA_LOG_DEVEL("BIG SYNC 3")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(29,29,29));

//   WALBERLA_LOG_DEVEL("BIG SYNC 4")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(31,31,31));

//   WALBERLA_LOG_DEVEL("BIG SYNC 5")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 5, 5, 5));

//   WALBERLA_LOG_DEVEL("BIG SYNC 6")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 9, 9, 9));

//   WALBERLA_LOG_DEVEL("BIG SYNC 7")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(11,11,11));

//   WALBERLA_LOG_DEVEL("BIG SYNC 8")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));


   WALBERLA_LOG_PROGRESS("*********************** [-1 1 1] TEST ***********************");
//   WALBERLA_LOG_DEVEL("BIG SYNC 9")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(11,19,19));

//   WALBERLA_LOG_DEVEL("BIG SYNC 10")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 9,21,21));

//   WALBERLA_LOG_DEVEL("BIG SYNC 11")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 5,25,25));

//   WALBERLA_LOG_DEVEL("BIG SYNC 12")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 1,29,29));

//   WALBERLA_LOG_DEVEL("BIG SYNC 13")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(-1,31,31));

//   WALBERLA_LOG_DEVEL("BIG SYNC 14")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(25,05, 5));

//   WALBERLA_LOG_DEVEL("BIG SYNC 15")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(21, 9, 9));

//   WALBERLA_LOG_DEVEL("BIG SYNC 16")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(19,11,11));

//   WALBERLA_LOG_DEVEL("BIG SYNC 17")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));


   WALBERLA_LOG_PROGRESS("*********************** [-1 -1 1] TEST ***********************");
//   WALBERLA_LOG_DEVEL("BIG SYNC 18")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(11,11,19));

//   WALBERLA_LOG_DEVEL("BIG SYNC 19")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 9, 9,21));

//   WALBERLA_LOG_DEVEL("BIG SYNC 20")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 5, 5,25));

//   WALBERLA_LOG_DEVEL("BIG SYNC 21")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3( 1, 1,29));

//   WALBERLA_LOG_DEVEL("BIG SYNC 22")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(-1,-1,31));

//   WALBERLA_LOG_DEVEL("BIG SYNC 23")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(25,25, 5));

//   WALBERLA_LOG_DEVEL("BIG SYNC 24")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(21,21, 9));

//   WALBERLA_LOG_DEVEL("BIG SYNC 25")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(19,19,11));

//   WALBERLA_LOG_DEVEL("BIG SYNC 26")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));


   WALBERLA_LOG_PROGRESS("*********************** [0 1 1] TEST ***********************");
//   WALBERLA_LOG_DEVEL("BIG SYNC 27")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,19,19));

//   WALBERLA_LOG_DEVEL("BIG SYNC 28")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,21,21));

//   WALBERLA_LOG_DEVEL("BIG SYNC 29")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,25,25));

//   WALBERLA_LOG_DEVEL("BIG SYNC 30")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,29,29));

//   WALBERLA_LOG_DEVEL("BIG SYNC 31")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,31,31));

//   WALBERLA_LOG_DEVEL("BIG SYNC 32")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15, 5, 5));

//   WALBERLA_LOG_DEVEL("BIG SYNC 33")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15, 9, 9));

//   WALBERLA_LOG_DEVEL("BIG SYNC 34")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,11,11));

//   WALBERLA_LOG_DEVEL("BIG SYNC 35")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));


   WALBERLA_LOG_PROGRESS("*********************** [0 0 1] TEST ***********************");
//   WALBERLA_LOG_DEVEL("BIG SYNC 36")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,19));

//   WALBERLA_LOG_DEVEL("BIG SYNC 37")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,21));

//   WALBERLA_LOG_DEVEL("BIG SYNC 38")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,25));

//   WALBERLA_LOG_DEVEL("BIG SYNC 39")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,29));

//   WALBERLA_LOG_DEVEL("BIG SYNC 40")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,31));

//   WALBERLA_LOG_DEVEL("BIG SYNC 41")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15, 5));

//   WALBERLA_LOG_DEVEL("BIG SYNC 42")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15, 9));

//   WALBERLA_LOG_DEVEL("BIG SYNC 43")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,11));

//   WALBERLA_LOG_DEVEL("BIG SYNC 44")
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   checkSphere(*forest, storageID, sid, refSphere, Vec3(15,15,15));

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}