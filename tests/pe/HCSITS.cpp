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
//! \file HCSITS.cpp
//! \brief checks equality of hash grids and simple ccd
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "core/debug/TestSubsystem.h"

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Sphere, Plane> BodyTuple ;

void normalReactionTest(cr::HCSITS& cr, SphereID sp)
{
   contactThreshold = Thresholds<real_t>::contactThreshold();
   // plane at 5,5,5
   // radius 1.1
   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setErrorReductionParameter( real_t(1.0) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(6.1)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(0.1)) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setErrorReductionParameter( real_t(0.5) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(6.05)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(0.05)) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setErrorReductionParameter( real_t(0.0) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,-1) );
   cr.setErrorReductionParameter( real_t(1.0) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(6.1)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(0.1)) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,-1) );
   cr.setErrorReductionParameter( real_t(0.5) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(6.05)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(0.05)) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,-1) );
   cr.setErrorReductionParameter( real_t(0.0) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0) );

   sp->setPosition(  Vec3(5,5,real_t(6.2)) );
   sp->setLinearVel( Vec3(0,0,real_t(-0.2)) );
   cr.setErrorReductionParameter( real_t(1.0) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(6.0)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(-0.2)) );

   contactThreshold = real_t(1.0);
   sp->setPosition(  Vec3(5,5,real_t(6.2)) );
   sp->setLinearVel( Vec3(0,0,real_t(-0.2)) );
   cr.setErrorReductionParameter( real_t(1.0) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(6.1)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(-0.1)) );

   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(6.1)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(0)) );
   contactThreshold = Thresholds<real_t>::contactThreshold();
}

void speedLimiterTest(cr::HCSITS& cr, SphereID sp)
{
   cr.setErrorReductionParameter( real_t(1.0) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setSpeedLimiter( true, real_t(0.2) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(6.1)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(0.1)) );

   sp->setPosition(  Vec3(5,5,5.5) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setSpeedLimiter( true, real_t(0.2) );
   cr.timestep( real_c( real_t(1.0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,real_t(5.94)) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,real_t(0.44)) );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            math::AABB(0,0,0,10,10,10),
            uint_c( 1), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            true,                               // max blocks per process
            true, true, true,                   // full periodicity
            false);

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto hccdID              = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "HCCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");
   cr::HCSITS cr(globalBodyStorage, forest->getBlockStoragePointer(), storageID, hccdID, fcdID);
   cr.setMaxIterations( 10 );
   cr.setRelaxationParameter    ( real_t(0.7) );
   cr.setErrorReductionParameter( real_t(1.0) );
   cr.setGlobalLinearAcceleration( Vec3(0,0,0) );

   pe::createPlane( *globalBodyStorage, 0, Vec3(0, 0, 1), Vec3(5, 5, 5) );

   SphereID sp = pe::createSphere(
            *globalBodyStorage,
            forest->getBlockStorage(),
            storageID,
            999999999,
            Vec3(5,5,6),
            real_c(1.1));

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::PROGRESS);

   WALBERLA_LOG_PROGRESS("Normal Reaction Test: InelasticFrictionlessContact");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticFrictionlessContact );
   normalReactionTest(cr, sp);
   WALBERLA_LOG_PROGRESS( "Normal Reaction Test: ApproximateInelasticCoulombContactByDecoupling");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling );
   normalReactionTest(cr, sp);
   //    WALBERLA_LOG_PROGRESS( "Normal Reaction Test: ApproximateInelasticCoulombContactByOrthogonalProjections");
   //    cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByOrthogonalProjections );
   //    normalReactionTest(cr, sp);
   WALBERLA_LOG_PROGRESS( "Normal Reaction Test: InelasticCoulombContactByDecoupling");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticCoulombContactByDecoupling );
   normalReactionTest(cr, sp);
   //    WALBERLA_LOG_PROGRESS( "Normal Reaction Test: InelasticCoulombContactByOrthogonalProjections");
   //    cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticCoulombContactByOrthogonalProjections );
   //    normalReactionTest(cr, sp);
   WALBERLA_LOG_PROGRESS( "Normal Reaction Test: InelasticGeneralizedMaximumDissipationContact");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticGeneralizedMaximumDissipationContact );
   normalReactionTest(cr, sp);

   WALBERLA_LOG_PROGRESS("SpeedLimiter Test: InelasticFrictionlessContact");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticFrictionlessContact );
   speedLimiterTest(cr, sp);

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}