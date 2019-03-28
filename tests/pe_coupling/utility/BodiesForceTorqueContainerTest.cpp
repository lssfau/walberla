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
//! \file BodiesForceTorqueContainerTest.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"

#include "pe/basic.h"
#include "pe/utility/DestroyBody.h"

#include <pe_coupling/utility/all.h>

namespace force_torque_container_test
{

///////////
// USING //
///////////

using namespace walberla;

using BodyTypeTuple = std::tuple<pe::Sphere> ;


/*!\brief Test cases for the force torque container provided by the coupling module
 *
 * Spheres at different positions are created and force(s) and torque(s) are applied.
 * Then, they are store in the container and reset on the body.
 * Then, in some cases, the sphere is moved to cross block borders, and the stored forces/torques are reapplied by the container again.
 * The obtained force/torque values are compared against the originally set (and thus expected) values.
 */
//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   // uncomment to have logging
   //logging::Logging::instance()->setLogLevel(logging::Logging::LogLevel::DETAIL);

   const real_t dx     = real_t(1);
   const real_t radius = real_t(5);

   ///////////////////////////
   // DATA STRUCTURES SETUP //
   ///////////////////////////

   Vector3<uint_t> blocksPerDirection(uint_t(3), uint_t(1), uint_t(1));
   Vector3<uint_t> cellsPerBlock(uint_t(20), uint_t(20), uint_t(20));
   Vector3<bool> periodicity(true, false, false);

   // create fully periodic domain with refined blocks
   auto blocks = blockforest::createUniformBlockGrid( blocksPerDirection[0], blocksPerDirection[1], blocksPerDirection[2],
                                                      cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
                                                      dx,
                                                      0, false, false,
                                                      periodicity[0], periodicity[1], periodicity[2],
                                                      false );


   // pe body storage
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "Storage");
   auto sphereMaterialID = pe::createMaterial( "sphereMat", real_t(1) , real_t(0.3), real_t(0.2), real_t(0.2), real_t(0.24), real_t(200), real_t(200), real_t(0), real_t(0) );

   // pe coupling
   const real_t overlap = real_t( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );

   // sphere positions for test scenarios
   Vector3<real_t> positionInsideBlock(real_t(10), real_t(10), real_t(10));
   Vector3<real_t> positionAtBlockBorder(real_t(19.5), real_t(10), real_t(10));
   Vector3<real_t> positionAtBlockBorderUpdated(real_t(20.5), real_t(10), real_t(10));

   Vector3<real_t> positionAtBlockBorder2(real_t(20) + radius + overlap - real_t(0.5), real_t(10), real_t(10));
   Vector3<real_t> positionAtBlockBorderUpdated2(real_t(20) + radius + overlap + real_t(0.5), real_t(10), real_t(10));


   Vector3<real_t> testForce(real_t(2), real_t(1), real_t(0));
   Vector3<real_t> torqueOffset = Vector3<real_t>(real_t(1), real_t(0), real_t(0));

   pe_coupling::ForceTorqueOnBodiesResetter resetter(blocks, bodyStorageID);
   shared_ptr<pe_coupling::BodiesForceTorqueContainer> container1 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
   shared_ptr<pe_coupling::BodiesForceTorqueContainer> container2 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
   pe_coupling::BodyContainerSwapper containerSwapper(container1, container2);

   //////////////////
   // Inside block //
   //////////////////
   {
      std::string testIdentifier("Test: sphere inside block");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");
      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       positionInsideBlock, radius, sphereMaterialID, false, true, false);

      syncCall();

      uint_t applicationCounter( 0 );

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            ++applicationCounter;
            auto pos = bodyIt->getPosition();
            bodyIt->addForceAtPos(testForce, pos+torqueOffset);
         }
      }
      mpi::allReduceInplace(applicationCounter, mpi::SUM);

      container1->store();

      resetter();

      containerSwapper();

      container2->setOnBodies();

      Vector3<real_t> expectedForce = applicationCounter * testForce;
      Vector3<real_t> expectedTorque = applicationCounter * ( torqueOffset % testForce );

      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting force: " << expectedForce);
      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting torque: " << expectedTorque);

      Vector3<real_t> actingForce(real_t(0));
      Vector3<real_t> actingTorque(real_t(0));
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            actingForce += bodyIt->getForce();
            actingTorque += bodyIt->getTorque();

         }
      }

      mpi::allReduceInplace(actingForce[0], mpi::SUM);
      mpi::allReduceInplace(actingForce[1], mpi::SUM);
      mpi::allReduceInplace(actingForce[2], mpi::SUM);

      mpi::allReduceInplace(actingTorque[0], mpi::SUM);
      mpi::allReduceInplace(actingTorque[1], mpi::SUM);
      mpi::allReduceInplace(actingTorque[2], mpi::SUM);

      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[0], expectedForce[0], "Mismatch in force0");
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[1], expectedForce[1], "Mismatch in force1");
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[2], expectedForce[2], "Mismatch in force2");

         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[0], expectedTorque[0], "Mismatch in torque0");
         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[1], expectedTorque[1], "Mismatch in torque1");
         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[2], expectedTorque[2], "Mismatch in torque2");
      }

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

      container1->clear();
      container2->clear();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended");

   }

   /////////////////////
   // At block border //
   /////////////////////
   {
      std::string testIdentifier("Test: sphere at block border");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");
      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       positionAtBlockBorder, radius, sphereMaterialID, false, true, false);

      syncCall();

      uint_t applicationCounter( 0 );

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            ++applicationCounter;
            auto pos = bodyIt->getPosition();
            bodyIt->addForceAtPos(testForce, pos+torqueOffset);
         }
      }
      mpi::allReduceInplace(applicationCounter, mpi::SUM);

      container1->store();

      resetter();

      containerSwapper();

      container2->setOnBodies();

      Vector3<real_t> expectedForce = applicationCounter * testForce;
      Vector3<real_t> expectedTorque = applicationCounter * ( torqueOffset % testForce );

      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting force: " << expectedForce);
      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting torque: " << expectedTorque);

      Vector3<real_t> actingForce(real_t(0));
      Vector3<real_t> actingTorque(real_t(0));
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            actingForce += bodyIt->getForce();
            actingTorque += bodyIt->getTorque();
         }
      }

      mpi::allReduceInplace(actingForce[0], mpi::SUM);
      mpi::allReduceInplace(actingForce[1], mpi::SUM);
      mpi::allReduceInplace(actingForce[2], mpi::SUM);

      mpi::allReduceInplace(actingTorque[0], mpi::SUM);
      mpi::allReduceInplace(actingTorque[1], mpi::SUM);
      mpi::allReduceInplace(actingTorque[2], mpi::SUM);

      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[0], expectedForce[0], "Mismatch in force0");
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[1], expectedForce[1], "Mismatch in force1");
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[2], expectedForce[2], "Mismatch in force2");

         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[0], expectedTorque[0], "Mismatch in torque0");
         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[1], expectedTorque[1], "Mismatch in torque1");
         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[2], expectedTorque[2], "Mismatch in torque2");
      }

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

      container1->clear();
      container2->clear();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended");

   }

   ///////////////////////////////////////////
   // At block border with updated position //
   ///////////////////////////////////////////
   {
      std::string testIdentifier("Test: sphere at block border with updated position");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");
      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       positionAtBlockBorder, radius, sphereMaterialID, false, true, false);

      syncCall();

      uint_t applicationCounter( 0 );

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            ++applicationCounter;
            auto pos = bodyIt->getPosition();
            bodyIt->addForceAtPos(testForce, pos+torqueOffset);
         }
      }
      mpi::allReduceInplace(applicationCounter, mpi::SUM);

      container1->store();

      resetter();

      containerSwapper();

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) {
         for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
            bodyIt->setPosition(positionAtBlockBorderUpdated);
         }
      }
      syncCall();

      container2->setOnBodies();

      Vector3<real_t> expectedForce = applicationCounter * testForce;
      Vector3<real_t> expectedTorque = applicationCounter * ( torqueOffset % testForce );

      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting force: " << expectedForce);
      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting torque: " << expectedTorque);

      Vector3<real_t> actingForce(real_t(0));
      Vector3<real_t> actingTorque(real_t(0));
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            actingForce += bodyIt->getForce();
            actingTorque += bodyIt->getTorque();
         }
      }

      mpi::allReduceInplace(actingForce[0], mpi::SUM);
      mpi::allReduceInplace(actingForce[1], mpi::SUM);
      mpi::allReduceInplace(actingForce[2], mpi::SUM);

      mpi::allReduceInplace(actingTorque[0], mpi::SUM);
      mpi::allReduceInplace(actingTorque[1], mpi::SUM);
      mpi::allReduceInplace(actingTorque[2], mpi::SUM);

      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[0], expectedForce[0], "Mismatch in force0");
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[1], expectedForce[1], "Mismatch in force1");
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[2], expectedForce[2], "Mismatch in force2");

         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[0], expectedTorque[0], "Mismatch in torque0");
         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[1], expectedTorque[1], "Mismatch in torque1");
         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[2], expectedTorque[2], "Mismatch in torque2");
      }

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

      container1->clear();
      container2->clear();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended");

   }

   /////////////////////////////////////////////
   // At block border with updated position 2 //
   /////////////////////////////////////////////
   {
      std::string testIdentifier("Test: sphere at block border with updated position 2");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");
      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       positionAtBlockBorder2, radius, sphereMaterialID, false, true, false);

      syncCall();

      uint_t applicationCounter( 0 );

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            ++applicationCounter;
            auto pos = bodyIt->getPosition();
            bodyIt->addForceAtPos(testForce, pos+torqueOffset);
         }
      }
      mpi::allReduceInplace(applicationCounter, mpi::SUM);

      container1->store();

      resetter();

      containerSwapper();

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) {
         for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
            bodyIt->setPosition(positionAtBlockBorderUpdated2);
         }
      }
      syncCall();

      container2->setOnBodies();

      --applicationCounter; // sphere has vanished from one block

      // in this case, the complete force can no longer be applied onto the body since the part from the block
      // from which the body vanished is lost

      // HOWEVER: this case will never appear in a coupled simulation since NO FORCE will be acting a body
      // that is about to vanish from a block since it will not be mapped to the block from which it will vanish

      // If you are expecting a different behavior, you need to change the container mechanism by first all-reducing the
      // forces/torques to all processes/blocks that know this body such that every process/block has
      // the full information and then the process/block that owns this body sets the complete force

      Vector3<real_t> expectedForce = applicationCounter * testForce;
      Vector3<real_t> expectedTorque = applicationCounter * ( torqueOffset % testForce );

      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting force: " << expectedForce);
      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting torque: " << expectedTorque);

      Vector3<real_t> actingForce(real_t(0));
      Vector3<real_t> actingTorque(real_t(0));
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            actingForce += bodyIt->getForce();
            actingTorque += bodyIt->getTorque();
         }
      }

      mpi::allReduceInplace(actingForce[0], mpi::SUM);
      mpi::allReduceInplace(actingForce[1], mpi::SUM);
      mpi::allReduceInplace(actingForce[2], mpi::SUM);

      mpi::allReduceInplace(actingTorque[0], mpi::SUM);
      mpi::allReduceInplace(actingTorque[1], mpi::SUM);
      mpi::allReduceInplace(actingTorque[2], mpi::SUM);

      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[0], expectedForce[0], "Mismatch in force0");
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[1], expectedForce[1], "Mismatch in force1");
         WALBERLA_CHECK_FLOAT_EQUAL(actingForce[2], expectedForce[2], "Mismatch in force2");

         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[0], expectedTorque[0], "Mismatch in torque0");
         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[1], expectedTorque[1], "Mismatch in torque1");
         WALBERLA_CHECK_FLOAT_EQUAL(actingTorque[2], expectedTorque[2], "Mismatch in torque2");
      }

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

      container1->clear();
      container2->clear();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended");

   }

   return 0;

}

} //namespace force_torque_container_test

int main( int argc, char **argv ){
   force_torque_container_test::main(argc, argv);
}
