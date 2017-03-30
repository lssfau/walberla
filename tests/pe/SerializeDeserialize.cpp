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
//! \ingroup pe_tests
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
#include "core/grid_generator/SCIterator.h"

#include <boost/tuple/tuple.hpp>

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::blockforest;

typedef boost::tuple<Sphere> BodyTuple ;

void createDump()
{
   using namespace walberla::grid_generator;

   BodyStorage globalStorage;

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 2), uint_c( 2), uint_c( 2), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c(10),                         // dx: length of one cell in physical coordinates
            0,                                  // max blocks per process
            false, false,                       // include metis / force metis
            false, false, false );                 // full periodicity

   BlockDataID storageID           = forest->addBlockData( createStorageDataHandling<BodyTuple>() );

   for (auto it = SCIterator(forest->getDomain(), Vec3(-1,-1,-1), 3); it != SCIterator(); ++it)
   {
      createSphere( globalStorage, forest->getBlockStorage(), storageID, 0, *it, 1);
   }

   WALBERLA_LOG_DEVEL_ON_ROOT("dumping body storage");
   forest->saveBlockData("BodyStorageDump.dump", storageID);

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      BodyStorage& localStorage = (*(blockIt->getData< Storage >( storageID )))[0];
      WALBERLA_LOG_DEVEL("DUMPING BLOCK");
      WALBERLA_LOG_DEVEL("aabb: " << blockIt->getAABB());
      WALBERLA_LOG_DEVEL("#bodies: " << localStorage.size());
   }
}

void checkDump()
{
   using namespace walberla::grid_generator;

//   BodyStorage globalStorage;

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 2), uint_c( 2), uint_c( 2), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c(10),                         // dx: length of one cell in physical coordinates
            0,                                  // max blocks per process
            false, false,                       // include metis / force metis
            false, false, false );                 // full periodicity



   BlockDataID storageID = forest->loadBlockData("BodyStorageDump.dump", createStorageDataHandling<BodyTuple>());

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      BodyStorage& localStorage = (*(blockIt->getData< Storage >( storageID )))[0];
      WALBERLA_LOG_DEVEL("CHECKING BLOCK");
      WALBERLA_LOG_DEVEL("aabb: " << blockIt->getAABB());
      WALBERLA_LOG_DEVEL("#bodies: " << localStorage.size());
      auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID);
      for (auto it = SCIterator(forest->getDomain(), Vec3(-1,-1,-1), 3); it != SCIterator(); ++it)
      {
         if (blockIt->getAABB().contains(*it))
         {
            WALBERLA_CHECK_FLOAT_EQUAL( bodyIt->getPosition(), *it);
            ++bodyIt;
         }
      }
   }
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   SetBodyTypeIDs<BodyTuple>::execute();

   createDump();
   checkDump();

   return EXIT_SUCCESS;
}
