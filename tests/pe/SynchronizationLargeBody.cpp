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
//! \file SynchronizationLargeBody.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/communication/ParseMessage.h"
#include "pe/basic.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"

#include <tuple>

namespace walberla {
using namespace walberla::pe;

using BodyTuple = std::tuple<Sphere> ;

// checkSphere without dx
void checkSphere(StructuredBlockForest& forest, BlockDataID storageID, walberla::id_t sid, SphereID ref, const Vec3& newPos)
{
   for (auto it = forest.begin(); it != forest.end(); ++it)
   {
      blockforest::Block& block = *(dynamic_cast<blockforest::Block*>(&(*it)));
      Storage& storage  = *(block.getData<Storage>(storageID));
      BodyStorage& shadowStorage = storage[1];

      if (block.getAABB().contains( ref->getPosition() ))
      {
         WALBERLA_CHECK_EQUAL( storage[StorageType::LOCAL].size(), 1 );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 0 );
         SphereID bd = static_cast<SphereID> (storage[0].find( sid ).getBodyID());
         WALBERLA_CHECK_NOT_NULLPTR(bd);
         checkVitalParameters(bd, ref);
         bd->setPosition( newPos );
      } else if (forest.periodicIntersect(block.getAABB(), ref->getAABB()) )
      {
         WALBERLA_CHECK_EQUAL( storage[0].size(), 0 );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 1 );
         SphereID bd = static_cast<SphereID> (shadowStorage.find( sid ).getBodyID());
         WALBERLA_CHECK_NOT_NULLPTR(bd);
         auto backupPos =ref->getPosition();
         auto correctedPos = ref->getPosition();
         pe::communication::correctBodyPosition(forest.getDomain(), block.getAABB().center(), correctedPos);
         ref->setPosition(correctedPos);
         checkVitalParameters(bd, ref);
         ref->setPosition(backupPos);
      } else
      {
         WALBERLA_CHECK_EQUAL( storage[0].size(), 0 );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 0 );
      }
   }
   //   WALBERLA_LOG_PROGRESS("checked pos: " << ref->getPosition() << " | new pos: " << newPos);
   auto temp = newPos;
   forest.mapToPeriodicDomain(temp);
   ref->setPosition(temp);
}

// checkSphere with dx
void checkSphere(StructuredBlockForest& forest, BlockDataID storageID, walberla::id_t sid, SphereID ref, const Vec3& newPos, const real_t dx)
{
   for (auto it = forest.begin(); it != forest.end(); ++it)
   {
      blockforest::Block& block = *(dynamic_cast<blockforest::Block*>(&(*it)));
      Storage& storage  = *(block.getData<Storage>(storageID));
      BodyStorage& shadowStorage = storage[1];

      if (block.getAABB().contains( ref->getPosition() ))
      {
         WALBERLA_CHECK_EQUAL( storage[StorageType::LOCAL].size(), 1 );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 0 );
         SphereID bd = static_cast<SphereID> (storage[0].find( sid ).getBodyID());
         WALBERLA_CHECK_NOT_NULLPTR(bd);
         checkVitalParameters(bd, ref);
         bd->setPosition( newPos );
      } else if (forest.periodicIntersect(block.getAABB(), ref->getAABB(), dx) )
      {
         WALBERLA_CHECK_EQUAL( storage[0].size(), 0 );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 1, "ref_sphere: " << ref << "\n" << block.getAABB() );
         SphereID bd = static_cast<SphereID> (shadowStorage.find( sid ).getBodyID());
         WALBERLA_CHECK_NOT_NULLPTR(bd);
         auto backupPos =ref->getPosition();
         auto correctedPos = ref->getPosition();
         pe::communication::correctBodyPosition(forest.getDomain(), block.getAABB().center(), correctedPos);
         ref->setPosition(correctedPos);
         checkVitalParameters(bd, ref);
         ref->setPosition(backupPos);
      } else
      {
         WALBERLA_CHECK_EQUAL( storage[0].size(), 0 );
         WALBERLA_CHECK_EQUAL( shadowStorage.size(), 0 );
      }
   }
   //   WALBERLA_LOG_PROGRESS("checked pos: " << ref->getPosition() << " | new pos: " << newPos);
   auto temp = newPos;
   forest.mapToPeriodicDomain(temp);
   ref->setPosition(temp);
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   //   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
   //   logging::Logging::instance()->includeLoggingToFile("SyncLog");
   //   logging::Logging::instance()->setStreamLogLevel( logging::Logging::PROGRESS );

   shared_ptr<BodyStorage> globalStorage = make_shared<BodyStorage> ();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 7), uint_c( 7), uint_c( 7), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c( 1),                         // dx: length of one cell in physical coordinates
            0,                                  // max blocks per process
            false, false,                       // include metis / force metis
            true, true, true,                   // full periodicity
            true);                              // keep global block information

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
   Vec3 gpos = Vec3(3.5, 3.5, 3.5);
   const real_t r = real_c(1.6);
   Sphere refSphere(1, 0, gpos, Quat(), r, iron, false, true, false);
   refSphere.setLinearVel(4, 5, 6);
   refSphere.setAngularVel( 1, 2, 3);


   SphereID sphere = createSphere( *globalStorage, forest->getBlockStorage(), storageID, 0, gpos, r);
   walberla::id_t sphereID = 789456123;
   if (sphere != nullptr) sphereID = sphere->getSystemID();
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

   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);

   for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
   {
      WALBERLA_LOG_PROGRESS("*********************** [" << dir.cx() << " " << dir.cy() << " "<< dir.cz() << "] TEST ***********************");
      Vec3 delta = Vec3( real_c(dir.cx()), real_c(dir.cy()), real_c(dir.cz()) ) / real_c(3.0);
      for (int i = 0; i < 21; ++i)
      {
         syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID);
         Vec3 pos = refSphere.getPosition() + delta;
         if (!forest->getDomain().contains( pos, real_c(0.5) ))
            forest->mapToPeriodicDomain(pos);
         checkSphere(*forest, storageID, sid, &refSphere, pos);
      }
   }
   WALBERLA_LOG_PROGRESS("TEST WITHOUT DX ... finished");

   //test with dx
   real_t dx = real_c(0.5);
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID, nullptr, dx);
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID, nullptr, dx);
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID, nullptr, dx);
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID, nullptr, dx);

   for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
   {
      WALBERLA_LOG_PROGRESS("*********************** [" << dir.cx() << " " << dir.cy() << " "<< dir.cz() << "] TEST ***********************");
      Vec3 delta = Vec3( real_c(dir.cx()), real_c(dir.cy()), real_c(dir.cz()) ) / real_c(3.0);
      for (int i = 0; i < 21; ++i)
      {
         syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID, nullptr, dx);
         Vec3 pos = refSphere.getPosition() + delta;
         if (!forest->getDomain().contains( pos, real_c(0.5) ))
            forest->mapToPeriodicDomain(pos);
         checkSphere(*forest, storageID, sid, &refSphere, pos, dx);
      }
   }
   syncShadowOwners<BodyTuple>(forest->getBlockForest(), storageID, nullptr, dx);
   WALBERLA_LOG_PROGRESS("TEST WITH DX ... finished");

   sphere = static_cast<SphereID> (getBody(*globalStorage, forest->getBlockStorage(), storageID, sphereID));

   if (sphere != nullptr)
   {
      //      WALBERLA_LOG_DEVEL("pos: " << sphere->getPosition());
      //      WALBERLA_LOG_DEVEL("aabb: " << sphere->getAABB());
      //      WALBERLA_LOG_DEVEL("shadows: " << sphere->MPITrait.sizeShadowOwners());
      //      for (auto shadowIt = sphere->MPITrait.beginShadowOwners(); shadowIt != sphere->MPITrait.endShadowOwners(); ++shadowIt)
      //      {
      //         WALBERLA_LOG_DEVEL(forest->getBlock(shadowIt->blockID_)->getAABB());
      //      }
      //      WALBERLA_LOG_DEVEL(math::AABB(Vec3(1,1,1), Vec3(2,2,2)).intersects(math::AABB(Vec3(1.9,1.9,1.9), Vec3(3,3,3))));
      //      WALBERLA_LOG_DEVEL(forest->periodicIntersect(math::AABB(Vec3(1,1,1), Vec3(2,2,2)), math::AABB(Vec3(1.9,1.9,1.9), Vec3(3,3,3))));
      if (!sphere->isRemote())
         WALBERLA_ASSERT_EQUAL(sphere->MPITrait.sizeShadowOwners(), 124);
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}