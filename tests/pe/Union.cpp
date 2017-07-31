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
#include "pe/ccd/SimpleCCDDataHandling.h"
#include "pe/synchronization/SyncNextNeighbors.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"
#include "vtk/VTKOutput.h"

#include <boost/tuple/tuple.hpp>

#include <algorithm>
#include <vector>

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::blockforest;

typedef Union< boost::tuple<Sphere> >          UnionType ;
typedef boost::tuple<Sphere, Plane, UnionType> BodyTuple ;

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

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

   MaterialID iron = Material::find("iron");
   UnionType* un   = createUnion< boost::tuple<Sphere> >( *globalBodyStorage, forest->getBlockStorage(), storageID, 0, Vec3(5,5,5) );
   SphereID sp1 = new Sphere( 10, 0, Vec3(5,5,1), Vec3(0,0,0), Quat(), real_t(1)  , iron, false, true, false );
   SphereID sp2 = new Sphere( 11, 0, Vec3(6,5,2), Vec3(0,0,0), Quat(), real_t(1.5), iron, false, true, false );
   un->add(sp1);
   un->add(sp2);

   auto vtkOutput   = make_shared<SphereVtkOutput>(storageID, forest->getBlockStorage()) ;
   auto vtkWriter   = vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, "vtk_out", "simulation_step", false, false);

   for (unsigned int i = 0; i < 10000; ++i)
   {
      vtkWriter->write( true );
      cr.timestep( real_t(0.001) );
   }

   //WALBERLA_CHECK_FLOAT_EQUAL( sp1->getLinearVel().length(), real_t(0) );
   //WALBERLA_CHECK_FLOAT_EQUAL( sp2->getLinearVel().length(), real_t(0) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( sp1->getPosition()[2], real_t(1)  , real_t(0.001) );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( sp2->getPosition()[2], real_t(1.5), real_t(0.001) );
   WALBERLA_CHECK_FLOAT_EQUAL( (sp1->getPosition() - sp2->getPosition()).length(), real_t(sqrt(2)) );

   //WALBERLA_LOG_DEVEL(un);

   return EXIT_SUCCESS;
}
