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
//! \file SimpleCCD.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/Materials.h"
#include "pe/Types.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/ccd/SimpleCCD.h"
#include "pe/ccd/HashGrids.h"
#include "pe/fcd/SimpleFCD.h"

#include "core/DataTypes.h"
#include "core/UniqueID.h"
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

    SetBodyTypeIDs<BodyTuple>::execute();

    MaterialID iron = Material::find("iron");

    BodyStorage globalStorage;
    Storage storage;
    BodyStorage& bs = storage[0];
    ccd::SimpleCCD sccd(globalStorage, storage);
    ccd::HashGrids hg(globalStorage, storage[0], storage[1]);

    fcd::SimpleFCD<BodyTuple> s_fcd;
    fcd::SimpleFCD<BodyTuple> hg_fcd;

    math::seedRandomGenerator(1337);

    for (uint_t i = 0; i < 100; ++i)
      storage[0].add( std::make_unique<Sphere>(UniqueID<Sphere>::createGlobal(), 0, Vec3( math::realRandom(real_c(0), real_c(10)), math::realRandom(real_c(0), real_c(10)), math::realRandom(real_c(0), real_c(10))), Vec3(0,0,0), Quat(), real_t(1), iron, false, false, false) );

    sccd.generatePossibleContacts();

    WALBERLA_CHECK_EQUAL( sccd.getPossibleContacts().size(), 4950 );

    s_fcd.generateContacts(sccd.getPossibleContacts());

    WALBERLA_LOG_DEVEL( s_fcd.getContacts().size() );

    BodyID bd = (storage[0].begin() + 5).getBodyID();
    storage[0].remove( bd );

    sccd.generatePossibleContacts();
    hg.generatePossibleContacts();

    WALBERLA_CHECK_EQUAL( sccd.getPossibleContacts().size(), 4851 );

    s_fcd.generateContacts(sccd.getPossibleContacts());
    hg_fcd.generateContacts(hg.getPossibleContacts());
    WALBERLA_CHECK_EQUAL( s_fcd.getContacts().size(), hg_fcd.getContacts().size());

    bs.clear();

    bs.add( std::make_unique<Sphere>(UniqueID<Sphere>::createGlobal(), 0, Vec3( math::realRandom(real_c(0), real_c(10)), math::realRandom(real_c(0), real_c(10)), math::realRandom(real_c(0), real_c(10))), Vec3(0,0,0), Quat(), real_t(1), iron, false, false, false) );

    WcTimingPool pool;
    for (int runs = 0; runs < 10; ++runs)
    {
       auto oldSize = bs.size();
       for (uint_t i = 0; i < oldSize; ++i)
         bs.add( std::make_unique<Sphere>(UniqueID<Sphere>::createGlobal(), 0, Vec3( math::realRandom(real_c(0), real_c(10)), math::realRandom(real_c(0), real_c(10)), math::realRandom(real_c(0), real_c(10))), Vec3(0,0,0), Quat(), real_t(0.5), iron, false, false, false) );
       pool["SCCD"].start();
       sccd.generatePossibleContacts();
       pool["SCCD"].end();
       pool["HG"].start();
       hg.generatePossibleContacts();
       pool["HG"].end();
       WALBERLA_CHECK_GREATER_EQUAL(sccd.getPossibleContacts().size(), hg.getPossibleContacts().size());

       s_fcd.generateContacts(sccd.getPossibleContacts());
       hg_fcd.generateContacts(hg.getPossibleContacts());
       WALBERLA_CHECK_EQUAL( s_fcd.getContacts().size(), hg_fcd.getContacts().size());

       WALBERLA_LOG_DEVEL_ON_ROOT(bs.size() << "\t" << pool["SCCD"].last() << "\t" << pool["HG"].last() << "\t" << sccd.getPossibleContacts().size() << "\t" << hg.getPossibleContacts().size() << "\t" << hg.active());
       //std::cout << pool << std::endl;
    }

    return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}