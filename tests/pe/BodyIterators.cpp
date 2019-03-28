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
//! \file BodyIterators.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"
#include "core/mpi/Reduce.h"

#include "pe/basic.h"

#include <cstdlib>
#include <iostream>


namespace walberla {
using namespace walberla::pe;

using BodyTuple = std::tuple<Sphere> ;

int main( int argc, char **argv )
{
    debug::enterTestMode();

    mpi::Environment env( argc, argv );

    if( MPIManager::instance()->numProcesses() != 2 )
    {
        std::cerr << "number of processes must be equal to 2!" << std::endl;
        return EXIT_FAILURE;
    }

    shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage> ();

    auto blocks = blockforest::createUniformBlockGrid(
                math::AABB( real_c(0), real_c(0), real_c(0), real_c(200), real_c(100), real_c(100) ),
                uint_c( 2 ), uint_c( 1 ), uint_c( 1 ),    // number of blocks in x/y/z direction
                uint_t( 10 ), uint_t( 10 ), uint_t( 10 ), // number of cells per block in x/y/z direction (not important for this test!)
                uint_c( 2 ), uint_c( 1 ), uint_c( 1 ),    // number of processes in x/y/z direction
                false, false, false );                    // NOT periodic

    SetBodyTypeIDs<BodyTuple>::execute();

    auto storageID           = blocks->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");

    // empty check
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
    {
        WALBERLA_CHECK_EQUAL( LocalBodyIterator::begin(*blockIt, storageID), LocalBodyIterator::end() );
        WALBERLA_CHECK_EQUAL( ShadowBodyIterator::begin(*blockIt, storageID), ShadowBodyIterator::end() );
        WALBERLA_CHECK_EQUAL( BodyIterator::begin(*blockIt, storageID), BodyIterator::end() );
    }

    uint_t sphereCount = 0;

    if (pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), storageID, 0, Vec3( real_t(25), real_t(25), real_t(50) ), 1) != nullptr)
        ++sphereCount;
    if (pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), storageID, 0, Vec3( real_t(99), real_t(25), real_t(50) ), 2) != nullptr)
        ++sphereCount;
    if (pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), storageID, 0, Vec3( real_t(101), real_t(25), real_t(50) ), 2) != nullptr)
        ++sphereCount;
    if (pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), storageID, 0, Vec3( real_t(125), real_t(25), real_t(50) ), 1) != nullptr)
        ++sphereCount;

    syncShadowOwners<BodyTuple>( blocks->getBlockForest(), storageID);

    uint_t counter = 0;

    // local bodies
    counter = uint_t(0);
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
        for( auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt )
        {
            ++counter;
            WALBERLA_CHECK( !(bodyIt->isRemote()) );
        }
    WALBERLA_CHECK_EQUAL( counter, sphereCount );
    counter = uint_t(0);
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
        for( auto bodyIt = LocalBodyIterator::begin<Sphere>(*blockIt, storageID); bodyIt != LocalBodyIterator::end<Sphere>(); ++bodyIt )
        {
            ++counter;
            WALBERLA_CHECK( !(bodyIt->isRemote()) );
        }
    WALBERLA_CHECK_EQUAL( counter, sphereCount );

    mpi::allReduceInplace( counter, mpi::SUM );
    WALBERLA_CHECK_EQUAL( counter, mpi::allReduce(sphereCount, mpi::SUM) );


    // shadow bodies
    counter = uint_t(0);
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
        for( auto bodyIt = ShadowBodyIterator::begin(*blockIt, storageID); bodyIt != ShadowBodyIterator::end(); ++bodyIt )
        {
            ++counter;
            WALBERLA_CHECK( bodyIt->isRemote() );
        }
    WALBERLA_CHECK_EQUAL( counter, 1 );
    counter = uint_t(0);
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
        for( auto bodyIt = ShadowBodyIterator::begin<Sphere>(*blockIt, storageID); bodyIt != ShadowBodyIterator::end<Sphere>(); ++bodyIt )
        {
            ++counter;
            WALBERLA_CHECK( bodyIt->isRemote() );
        }
    WALBERLA_CHECK_EQUAL( counter, 1 );

    mpi::allReduceInplace( counter, mpi::SUM );
    WALBERLA_CHECK_EQUAL( counter, 2 );

    // local&shadow bodies
    counter = uint_t(0);
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
        for( auto bodyIt = BodyIterator::begin(*blockIt, storageID); bodyIt != BodyIterator::end(); ++bodyIt )
        {
            ++counter;
        }
    WALBERLA_CHECK_EQUAL( counter, 3 );
    counter = uint_t(0);
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
        for( auto bodyIt = BodyIterator::begin<Sphere>(*blockIt, storageID); bodyIt != BodyIterator::end<Sphere>(); ++bodyIt )
        {
            ++counter;
        }
    WALBERLA_CHECK_EQUAL( counter, 3 );

    mpi::allReduceInplace( counter, mpi::SUM );
    WALBERLA_CHECK_EQUAL( counter, 6 );

    return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}