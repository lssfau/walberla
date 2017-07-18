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
//! \file SaveLoadTest.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "core/math/IntegerFactorization.h"
#include "domain_decomposition/all.h"

#include "core/debug/TestSubsystem.h"

using namespace walberla;
using namespace walberla::blockforest;

int main( int argc, char ** argv )
{
   std::vector< uint64_t > dump;
   std::vector< uint64_t > check;

   walberla::debug::enterTestMode();

   WALBERLA_MPI_SECTION()
   {
      walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   }

   WALBERLA_LOG_DEVEL_ON_ROOT("*** DUMPING ***");

   dump.clear();

   auto proc = math::getFactors3D(uint_c( MPIManager::instance()->numProcesses() ));

   auto forestDump = createUniformBlockGrid( math::AABB(0,0,0,60,60,60), // domain
                                             2,2,2,                      // number of blocks
                                             1,1,1,                      // cells
                                             proc[0],proc[1],proc[2]);                     // on block per prozess
   forestDump->getBlockForest().saveToFile("SerializeDeserialize.sbf");

   for (auto blockIt = forestDump->begin(); blockIt != forestDump->end(); ++blockIt)
   {
      WALBERLA_LOG_DEVEL("DUMPING BLOCK (" << blockIt->getId() << ") " << blockIt->getAABB() );
      dump.push_back( blockIt->getId().getID() );
   }

   WALBERLA_MPI_SECTION()
   {
      WALBERLA_MPI_BARRIER();
   }

   WALBERLA_LOG_DEVEL_ON_ROOT("*** CHECKING ***");

   check.clear();

   auto forestCheck = shared_ptr< BlockForest >( new BlockForest( uint_c( MPIManager::instance()->rank() ), "SerializeDeserialize.sbf" ) );

   for (auto blockIt = forestCheck->begin(); blockIt != forestCheck->end(); ++blockIt)
   {
      WALBERLA_LOG_DEVEL("CHECKING BLOCK (" << blockIt->getId() << ") " << blockIt->getAABB() );
      check.push_back( blockIt->getId().getID() );
   }

   WALBERLA_CHECK_EQUAL(forestDump->getBlockIdBytes(), forestCheck->getBlockIdBytes());
   WALBERLA_CHECK_EQUAL(forestDump->getDepth(), forestCheck->getDepth());
   WALBERLA_CHECK_EQUAL(forestDump->getDomain(), forestCheck->getDomain());
   WALBERLA_CHECK_EQUAL(forestDump->getNumberOfBlocks(), forestCheck->getNumberOfBlocks());

   std::sort(dump.begin(), dump.end());
   std::sort(check.begin(), check.end());

   WALBERLA_CHECK_EQUAL( dump.size(), check.size() );
   for (size_t i = 0; i < dump.size(); ++i)
   {
      WALBERLA_CHECK_EQUAL(dump[i], check[i]);
   }

   return EXIT_SUCCESS;
}
