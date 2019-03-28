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
//! \file SerializeDeserialize.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/Initialization.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "pe/basic.h"
#include "pe/communication/ParseMessage.h"
#include "pe/synchronization/SyncNextNeighbors.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"
#include "core/grid_generator/SCIterator.h"

#include <tuple>

namespace walberla {
using namespace walberla::pe;

using BodyTuple = std::tuple<Sphere> ;

void createDump()
{
   using namespace walberla::grid_generator;

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   //! [Dump Blockforest]
   auto forest = blockforest::createBlockForest( math::AABB(0,0,0,60,60,60),
                                                 Vector3<uint_t>(2,2,2),                   // number of blocks
                                                 Vector3<bool>(false, false, false));      // periodicity
   forest->saveToFile("SerializeDeserialize.sbf");
   //! [Dump Blockforest]

   //! [Init Storage]
   auto storageID       = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   //! [Init Storage]
   auto ccdID           = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");

   for (auto it = SCIterator(forest->getDomain(), Vec3(-1,-1,-1), 3); it != SCIterator(); ++it)
   {
      createSphere( *globalBodyStorage, *forest, storageID, 0, *it, 1);
   }

   WALBERLA_LOG_DEVEL_ON_ROOT("dumping body storage");
   //! [Save Simulation Data]
   forest->saveBlockData("SerializeDeserialize.dump", storageID);
   //! [Save Simulation Data]

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      BodyStorage& localStorage = (*(blockIt->getData< Storage >( storageID )))[0];
      WALBERLA_LOG_DEVEL("DUMPING BLOCK (" << blockIt->getId() << ") " << blockIt->getAABB() );
      WALBERLA_LOG_DEVEL("#bodies: " << localStorage.size());

      ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
      WALBERLA_CHECK_EQUAL( ccd->getObservedBodyCount(), 1000 );
   }
}

void checkDump()
{
   using namespace walberla::grid_generator;

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   //! [Load Blockforest]
   auto forest = make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), "SerializeDeserialize.sbf", true, false );
   //! [Load Blockforest]

   //! [Load Storage]
   auto storageID    = forest->loadBlockData("SerializeDeserialize.dump", createStorageDataHandling<BodyTuple>(), "Storage");
   //! [Load Storage]
   auto ccdID        = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");

   //! [Reload CCD]
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
      ccd->reloadBodies();
   }
   //! [Reload CCD]

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
      WALBERLA_CHECK_EQUAL( ccd->getObservedBodyCount(), 1000 );

      BodyStorage& localStorage = (*(blockIt->getData< Storage >( storageID )))[0];
      WALBERLA_LOG_DEVEL("CHECKING BLOCK (" << blockIt->getId() << ") " << blockIt->getAABB() );
      WALBERLA_LOG_DEVEL("#bodies: " << localStorage.size());
      auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID);
      for (auto it = SCIterator(forest->getDomain(), Vec3(-1,-1,-1), 3); it != SCIterator(); ++it)
      {
         if (blockIt->getAABB().contains(*it))
         {
            WALBERLA_CHECK_FLOAT_EQUAL( bodyIt->getPosition(), *it, blockIt->getAABB());
            ++bodyIt;
         }
      }
   }
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   WALBERLA_MPI_SECTION()
   {
      walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   }

   SetBodyTypeIDs<BodyTuple>::execute();

   WALBERLA_LOG_DEVEL_ON_ROOT("*** DUMPING ***");
   createDump();
   WALBERLA_MPI_SECTION()
   {
      WALBERLA_MPI_BARRIER();
   }
   WALBERLA_LOG_DEVEL_ON_ROOT("*** CHECKING ***");
   checkDump();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
