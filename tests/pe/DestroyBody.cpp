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
//! \file Destroy.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"
#include "pe/cr/PlainIntegrator.h"
#include "pe/ccd/SimpleCCDDataHandling.h"
#include "pe/utility/DestroyBody.h"

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "core/timing/TimingPool.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"

namespace walberla {
using namespace walberla::pe;

using BodyTuple = std::tuple<Sphere> ;

int main( int argc, char** argv )
{
    walberla::debug::enterTestMode();

    walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

    shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

    // create blocks
    shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
             uint_c( 2), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
             uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
             real_c(10),                         // dx: length of one cell in physical coordinates
             0,                                  // max blocks per process
             false, false,                       // include metis / force metis
             false, false, false );                 // full periodicity

    SetBodyTypeIDs<BodyTuple>::execute();

    auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
    Storage* firstStorage = nullptr, *secondStorage = nullptr;
    for (auto it = forest->begin(); it != forest->end(); ++it)
    {
       IBlock & currentBlock = *it;
       if (currentBlock.getAABB().contains( 5,5,5) )
          firstStorage = currentBlock.getData< Storage >( storageID );
       else
          secondStorage = currentBlock.getData< Storage >( storageID );
    }

    WALBERLA_CHECK_UNEQUAL(firstStorage, secondStorage);

    pe::createSphere(*globalBodyStorage, forest->getBlockStorage(), storageID, 0, Vec3(5,5,5), 2);
    auto sp = pe::createSphere(*globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(9,5,5), 2);

    WALBERLA_CHECK_EQUAL( firstStorage->at(0).size(), 2 );
    WALBERLA_CHECK_EQUAL( firstStorage->at(1).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(0).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(1).size(), 0 );

    syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID);

    WALBERLA_CHECK_EQUAL( firstStorage->at(0).size(), 2 );
    WALBERLA_CHECK_EQUAL( firstStorage->at(1).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(0).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(1).size(), 1 );

    destroyBodyBySID(*globalBodyStorage, forest->getBlockStorage(), storageID, sp->getSystemID());

    WALBERLA_CHECK_EQUAL( firstStorage->at(0).size(), 1 );
    WALBERLA_CHECK_EQUAL( firstStorage->at(1).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(0).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(1).size(), 0 );

    pe::createSphere(*globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(5,5,5), 2);
    pe::createSphere(*globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(9,5,5), 2);
    pe::createSphere(*globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(11,5,5), 2);
    pe::createSphere(*globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(15,5,5), 2);

    WALBERLA_CHECK_EQUAL( firstStorage->at(0).size(), 3 );
    WALBERLA_CHECK_EQUAL( firstStorage->at(1).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(0).size(), 2 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(1).size(), 0 );

    syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID);

    WALBERLA_CHECK_EQUAL( firstStorage->at(0).size(), 3 );
    WALBERLA_CHECK_EQUAL( firstStorage->at(1).size(), 1 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(0).size(), 2 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(1).size(), 1 );

    destroyBodyByUID(*globalBodyStorage, forest->getBlockStorage(), storageID, 1);

    WALBERLA_CHECK_EQUAL( firstStorage->at(0).size(), 1 );
    WALBERLA_CHECK_EQUAL( firstStorage->at(1).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(0).size(), 0 );
    WALBERLA_CHECK_EQUAL( secondStorage->at(1).size(), 0 );

    return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}