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
//! \file HashGrids.cpp
//! \brief checks equality of hash grids and simple ccd
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"
#include "pe/cr/PlainIntegrator.h"
#include "pe/ccd/SimpleCCDDataHandling.h"

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "core/timing/TimingPool.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Sphere, Plane> BodyTuple ;

int main( int argc, char** argv )
{
    if (argc != 2) WALBERLA_ABORT("Number of particles expected as first argument!");
    uint_t numParticles = uint_c(std::atoi(argv[1]));

    walberla::debug::enterTestMode();

    walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

    WALBERLA_LOG_DEVEL("running with " << numParticles << " particles (" << argv[1] << ")");

    MaterialID iron = Material::find("iron");

    shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

    // create blocks
    shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
             uint_c( 3), uint_c( 3), uint_c( 5), // number of blocks in x,y,z direction
             uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
             real_c(10),                         // dx: length of one cell in physical coordinates
             0,                                  // max blocks per process
             false, false,                       // include metis / force metis
             true, true, true );                 // full periodicity

    SetBodyTypeIDs<BodyTuple>::execute();

    auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
    auto sccdID              = forest->addBlockData(ccd::createSimpleCCDDataHandling( globalBodyStorage, storageID ), "SCCD");
    auto hccdID              = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "HCCD");
    auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");
    cr::PlainIntegrator cr(globalBodyStorage, forest->getBlockStoragePointer(), storageID, nullptr);

    pe::createPlane( *globalBodyStorage, 0, Vec3(0, +1, 0), Vec3(5, 0,5), iron);
    pe::createPlane( *globalBodyStorage, 0, Vec3(0, -1, 0), Vec3(5,30,5), iron);

    pe::createPlane( *globalBodyStorage, 0, Vec3(+1, 0, 0), Vec3( 0,5,5), iron);
    pe::createPlane( *globalBodyStorage, 0, Vec3(-1, 0, 0), Vec3(30,5,5), iron);

    pe::createPlane( *globalBodyStorage, 0, Vec3( 0, 0,+1), Vec3(5,5, 0), iron);
    pe::createPlane( *globalBodyStorage, 0, Vec3( 0, 0,-1), Vec3(5,5,30), iron);

    pe::createSphere( *globalBodyStorage, forest->getBlockStorage(), storageID, 999999999,
                     Vec3(1,1,1), real_c(1.1));


    syncShadowOwners<BodyTuple>( forest->getBlockForest(), storageID);
    syncShadowOwners<BodyTuple>( forest->getBlockForest(), storageID);

    for (auto it = forest->begin(); it != forest->end(); ++it)
    {
       IBlock & currentBlock = *it;

       ccd::ICCD* sccd = currentBlock.getData< ccd::ICCD >( sccdID );
       ccd::ICCD* hccd = currentBlock.getData< ccd::ICCD >( hccdID );
       fcd::IFCD* fcd  = currentBlock.getData< fcd::IFCD >( fcdID );
       Contacts cont1 = fcd->generateContacts( sccd->generatePossibleContacts() );
       auto tmp1 = cont1.size();
       Contacts cont2 = fcd->generateContacts( hccd->generatePossibleContacts() );
       auto tmp2 = cont2.size();
       WALBERLA_LOG_DETAIL_ON_ROOT("tracked particles: " << sccd->getObservedBodyCount() << "/" << hccd->getObservedBodyCount());
       WALBERLA_CHECK_EQUAL(tmp1, tmp2);
       WALBERLA_LOG_DETAIL_ON_ROOT("contacts on root: " << cont1.size());
    }

    math::seedRandomGenerator(1337);

    const real_t dv = 3;
    for (uint_t i = 0; i < numParticles; ++i)
    {
       SphereID sp = pe::createSphere(*globalBodyStorage, forest->getBlockStorage(), storageID, i,
                        Vec3(math::realRandom<real_t>(real_c(0), real_c(30)), math::realRandom<real_t>(real_c(0), real_c(30)), math::realRandom<real_t>(real_c(0), real_c(30))), real_c(0.4));
      if (sp != nullptr) sp->setLinearVel(Vec3(math::realRandom<real_t>(-dv, dv), math::realRandom<real_t>(-dv, dv), math::realRandom<real_t>(-dv, dv)));
    }

    pe::createSphere(*globalBodyStorage, forest->getBlockStorage(), storageID, 999999999,
                     Vec3(15,15,15), 7,
                     iron, true, false, true);

    syncShadowOwners<BodyTuple>( forest->getBlockForest(), storageID);
    for (int step=0; step < 100; ++step)
    {
       cr( real_c(0.1) );
       syncShadowOwners<BodyTuple>( forest->getBlockForest(), storageID);

       for (auto it = forest->begin(); it != forest->end(); ++it){
          IBlock & currentBlock = *it;

          ccd::ICCD* sccd = currentBlock.getData< ccd::ICCD >( sccdID );
          ccd::ICCD* hccd = currentBlock.getData< ccd::ICCD >( hccdID );
          fcd::IFCD* fcd  = currentBlock.getData< fcd::IFCD >( fcdID );
          Contacts cont1 = fcd->generateContacts( sccd->generatePossibleContacts() );
          auto tmp1 = cont1.size();
          Contacts cont2 = fcd->generateContacts( hccd->generatePossibleContacts() );
          auto tmp2 = cont2.size();
          WALBERLA_LOG_DETAIL_ON_ROOT("tracked particles: " << sccd->getObservedBodyCount() << "/" << hccd->getObservedBodyCount());
          WALBERLA_CHECK_EQUAL(tmp1, tmp2);
          WALBERLA_LOG_DETAIL_ON_ROOT("contacts on root: " << cont1.size());

          // check for correct ordering of bodies within contacts
          for (size_t i = 0; i < cont1.size(); ++i)
          {
             WALBERLA_CHECK_LESS(cont1[i].getBody1()->getSystemID(), cont1[i].getBody2()->getSystemID());
             WALBERLA_CHECK_LESS(cont2[i].getBody1()->getSystemID(), cont2[i].getBody2()->getSystemID());
          }
       }
    }

    forest.reset();

    return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}