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
//! \file ParallelEquivalence.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/basic.h"
#include "pe/rigidbody/Union.h"
#include "pe/rigidbody/UnionFactory.h"
#include "pe/ccd/SimpleCCDDataHandling.h"
#include "pe/synchronization/SyncNextNeighbors.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"
#include "vtk/VTKOutput.h"

#include <tuple>

#include <algorithm>
#include <vector>

namespace walberla {
using namespace walberla::pe;

using UnionType = Union<Sphere> ;
typedef std::tuple<Sphere, Plane, UnionType> BodyTuple ;

void SnowManFallingOnPlane()
{
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 1), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c(10),                         // dx: length of one cell in physical coordinates
            0,                                  // max blocks per process
            false, false,                       // include metis / force metis
            false, false, false );                 // full periodicity

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   cr::HCSITS cr(globalBodyStorage, forest->getBlockStoragePointer(), storageID, ccdID, fcdID);
   cr.setMaxIterations( uint_c(10) );
   cr.setRelaxationModel( cr::HCSITS::InelasticGeneralizedMaximumDissipationContact );
   cr.setRelaxationParameter( real_t(0.7) );
   cr.setErrorReductionParameter( real_t(0.8) );
   cr.setGlobalLinearAcceleration( Vec3(0,0,-1) );

   createPlane( *globalBodyStorage, 0, Vec3(0,0,1), Vec3(0,0,0) );

   UnionType* un  = createUnion<Sphere>( *globalBodyStorage, forest->getBlockStorage(), storageID, 0, Vec3(5,5,5) );
   auto sp1       = createSphere(un, 10, Vec3(5,5,1), real_t(1));
   auto sp2       = createSphere(un, 11, Vec3(real_t(6.7),5,real_t(1.2)), real_t(1.1));

   auto distance = (sp1->getPosition() - sp2->getPosition()).length();

   //auto vtkOutput   = make_shared<SphereVtkOutput>(storageID, forest->getBlockStorage()) ;
   //auto vtkWriter   = vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, "vtk_out", "simulation_step", false, false);

   for (unsigned int i = 0; i < 1000; ++i)
   {
      //vtkWriter->write( true );
      cr.timestep( real_t(0.1) );
   }

   //WALBERLA_CHECK_FLOAT_EQUAL( sp1->getLinearVel().length(), real_t(0) );
   //WALBERLA_CHECK_FLOAT_EQUAL( sp2->getLinearVel().length(), real_t(0) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( sp1->getPosition()[2], real_t(1)  , real_t(0.001) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( sp2->getPosition()[2], real_t(1.1), real_t(0.001) );
   WALBERLA_CHECK_FLOAT_EQUAL( (sp1->getPosition() - sp2->getPosition()).length(), distance );

   //WALBERLA_LOG_DEVEL(un);
}

void ImpulsCarryover()
{
   MaterialID iron = Material::find("iron");

   auto un  = std::make_unique<UnionType>(12, 0, Vec3(0,0,0), Vec3(0,0,0), Quat(), false, true, false);
   auto sp1 = std::make_unique<Sphere>( 10, 0, Vec3( 1,0,0), Vec3(0,0,0), Quat(), real_t(1), iron, false, true, false );
   auto sp2 = std::make_unique<Sphere>( 11, 0, Vec3(-1,0,0), Vec3(0,0,0), Quat(), real_t(1), iron, false, true, false );

   sp1->setLinearVel(Vec3(0,real_c(+1),0));
   sp2->setLinearVel(Vec3(0,real_c(-1),0));

   un->add( std::move(sp1) );
   un->add( std::move(sp2) );

   WALBERLA_CHECK_FLOAT_EQUAL( un->getPosition(),  Vec3(0,0,0) );
   WALBERLA_CHECK_FLOAT_EQUAL( un->getLinearVel(), Vec3(0,0,0) );
   WALBERLA_CHECK_FLOAT_EQUAL( un->getAngularVel() * Vec3(1,0,0), real_t(0) );
   WALBERLA_CHECK_FLOAT_EQUAL( un->getAngularVel() * Vec3(0,1,0), real_t(0) );
   WALBERLA_CHECK_GREATER( un->getAngularVel() * Vec3(0,0,1), real_t(0) );
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   SnowManFallingOnPlane();
   ImpulsCarryover();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}