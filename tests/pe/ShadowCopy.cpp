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
//! \file ShadowCopy.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"
#include "pe/utility/GetBody.h"
#include "pe/utility/DestroyBody.h"

#include "blockforest/Initialization.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "core/debug/TestSubsystem.h"

using namespace walberla;
using namespace walberla::pe;

typedef Union< boost::tuple<Sphere> > UnionT;
typedef boost::tuple<Sphere, UnionT> BodyTuple ;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            math::AABB(0,0,0,15,15,15),
            uint_c( 3), uint_c( 3), uint_c( 3), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            false,                               // max blocks per process
            true, true, true,                   // full periodicity
            false);

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
                              forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "HCCD");
                              forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");


//   logging::Logging::instance()->setStreamLogLevel(logging::Logging::DETAIL);
//   logging::Logging::instance()->setFileLogLevel(logging::Logging::DETAIL);
//   logging::Logging::instance()->includeLoggingToFile("ShadowCopy");

   bool syncShadowOwners = false;
   boost::function<void(void)> syncCall;
   if (!syncShadowOwners)
   {
      syncCall = boost::bind( pe::syncNextNeighbors<BodyTuple>, boost::ref(forest->getBlockForest()), storageID, static_cast<WcTimingTree*>(NULL), real_c(0.0), false );
   } else
   {
      syncCall = boost::bind( pe::syncShadowOwners<BodyTuple>, boost::ref(forest->getBlockForest()), storageID, static_cast<WcTimingTree*>(NULL), real_c(0.0), false );
   }

   WALBERLA_LOG_PROGRESS_ON_ROOT( " *** SPHERE *** ");
   SphereID sp = pe::createSphere(
            *globalBodyStorage,
            forest->getBlockStorage(),
            storageID,
            999999999,
            Vec3(real_t(4.9),2,2),
            real_c(1.2));
   auto sid = sp->getSystemID();
   sp->setLinearVel(1,2,3);
   sp->setAngularVel(1,2,3);
   syncCall();
   sp->setPosition(7,2,2);
   syncCall();
   sp = static_cast<SphereID> (getBody( *globalBodyStorage, forest->getBlockStorage(), storageID, sid, StorageSelect::LOCAL ));
   WALBERLA_CHECK_NOT_NULLPTR(sp);
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(1,2,3) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getAngularVel(), Vec3(1,2,3) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getRadius(), real_t(1.2) );
   destroyBodyBySID( *globalBodyStorage, forest->getBlockStorage(), storageID, sid );

   WALBERLA_LOG_PROGRESS_ON_ROOT( " *** UNION *** ");
   MaterialID iron = Material::find("iron");
   UnionT* un   = createUnion< boost::tuple<Sphere> >( *globalBodyStorage, forest->getBlockStorage(), storageID, 0, Vec3(2,2,2) );
   SphereID sp1 = new Sphere( 10, 0, Vec3(real_t(4.9),2,2), Vec3(0,0,0), Quat(), real_t(1)  , iron, false, true, false );
   SphereID sp2 = new Sphere( 11, 0, Vec3(3,2,2)          , Vec3(0,0,0), Quat(), real_t(1.5), iron, false, true, false );
   un->add(sp1);
   un->add(sp2);
   un->setPosition( Vec3( real_t(4.9), 2, 2) );
   auto relPosSp1 = sp1->getRelPosition();
   auto relPosSp2 = sp2->getRelPosition();
   sid = un->getSystemID();
   syncCall();

   un->setPosition( Vec3( real_t(5.9), 2, 2) );
   WALBERLA_LOG_PROGRESS_ON_ROOT( un->getPosition() );
   auto posUnion  = un->getPosition();
   syncCall();

   un  = static_cast<UnionT*> (getBody( *globalBodyStorage, forest->getBlockStorage(), storageID, sid, StorageSelect::LOCAL ));
   sp1 = static_cast<SphereID> (*(un->begin()));
   sp2 = static_cast<SphereID> (*(++(un->begin())));
   WALBERLA_CHECK_NOT_NULLPTR(sp1);
   WALBERLA_CHECK_NOT_NULLPTR(sp2);
   WALBERLA_CHECK_EQUAL( sp1->getTypeID(), Sphere::getStaticTypeID() );
   WALBERLA_CHECK_EQUAL( sp2->getTypeID(), Sphere::getStaticTypeID() );
   WALBERLA_CHECK_FLOAT_EQUAL( sp1->getRadius(), real_t(1.0) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp2->getRadius(), real_t(1.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( un->getPosition(), posUnion );
   WALBERLA_CHECK_FLOAT_EQUAL( relPosSp1, sp1->getRelPosition() );
   WALBERLA_CHECK_FLOAT_EQUAL( relPosSp2, sp2->getRelPosition() );

   un->setPosition(real_t(7.5),2,2);
   WALBERLA_LOG_PROGRESS_ON_ROOT( un->getPosition() );
   syncCall();

   un->setPosition(real_t(9.9),2,2);
   WALBERLA_LOG_PROGRESS_ON_ROOT( un->getPosition() );
   syncCall();

   un->setPosition(real_t(10.9),2,2);
   WALBERLA_LOG_PROGRESS_ON_ROOT( un->getPosition() );
   syncCall();

   un  = static_cast<UnionT*> (getBody( *globalBodyStorage, forest->getBlockStorage(), storageID, sid, StorageSelect::LOCAL ));
   un->setPosition(real_t(12.5),2,2);
   WALBERLA_LOG_PROGRESS_ON_ROOT( un->getPosition() );
   syncCall();

   un->setPosition(real_t(14.9),2,2);
   WALBERLA_LOG_PROGRESS_ON_ROOT( un->getPosition() );
   syncCall();

   un->setPosition(real_t(15.9),2,2);
   WALBERLA_LOG_PROGRESS_ON_ROOT( un->getPosition() );
   syncCall();

   posUnion = Vec3(real_t(0.9),2,2);
   un  = static_cast<UnionT*> (getBody( *globalBodyStorage, forest->getBlockStorage(), storageID, sid, StorageSelect::LOCAL ));
   sp1 = static_cast<SphereID> (*(un->begin()));
   sp2 = static_cast<SphereID> (*(++(un->begin())));
   WALBERLA_CHECK_NOT_NULLPTR(sp1);
   WALBERLA_CHECK_NOT_NULLPTR(sp2);
   WALBERLA_CHECK_EQUAL( sp1->getTypeID(), Sphere::getStaticTypeID() );
   WALBERLA_CHECK_EQUAL( sp2->getTypeID(), Sphere::getStaticTypeID() );
   WALBERLA_CHECK_FLOAT_EQUAL( sp1->getRadius(), real_t(1.0) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp2->getRadius(), real_t(1.5) );
   WALBERLA_CHECK_FLOAT_EQUAL( un->getPosition(), posUnion );
   WALBERLA_CHECK_FLOAT_EQUAL( relPosSp1, sp1->getRelPosition() );
   WALBERLA_CHECK_FLOAT_EQUAL( relPosSp2, sp2->getRelPosition() );

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt)
      {
         WALBERLA_LOG_DEVEL( blockIt->getAABB() );
         WALBERLA_LOG_DEVEL(*bodyIt );
      }
   }

   return EXIT_SUCCESS;
}
