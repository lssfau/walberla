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
//! \file Union.cpp
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
   WALBERLA_LOG_INFO("- Performing snowman on plane test");
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

void UnionConstruction()
{

   MaterialID iron = Material::find("iron");

   // Construct a union of two spheres in a non-rotated state and check all its parameters
   auto un  = std::make_unique<UnionType>(12, 0, Vec3(0,0,0), Quat(), false, true, false);
   auto sp1 = std::make_unique<Sphere>( 10, 0, Vec3( 1,0,0), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5)), real_t(1), iron, false, true, false );
   auto sp2 = std::make_unique<Sphere>( 11, 0, Vec3(-1,0,0), Quat(), real_t(1), iron, false, true, false );

   WALBERLA_LOG_INFO("- Adding first body.");
   // Add first body
   un->add( std::move(sp1) );
   const RigidBody& subb1 =  *un->begin();

   // Check relative and global position and rotation
   WALBERLA_CHECK_FLOAT_EQUAL(un->getPosition(),Vec3(1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(un->getRelPosition(), Vec3(0,0,0));

   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getRelPosition(), Vec3(0, 0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getPosition(), Vec3(1,0,0));
   // Global quaternion
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getQuaternion(), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5)));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getRelQuaternion(), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5)));

   // Check mass volume and inertia
   WALBERLA_CHECK_FLOAT_EQUAL(un->getVolume(), (real_t(4./3.)*real_t(math::pi)));
   WALBERLA_CHECK_FLOAT_EQUAL(un->getMass(), un->getVolume()*Material::getDensity(iron));
   real_t scalar_inertia = real_t(0.4)*un->getMass(); // for sphere: I = 2/5*m*r*r
   WALBERLA_CHECK_EQUAL(un->getInertia(), Mat3(scalar_inertia,0,0,0,scalar_inertia,0,0,0,scalar_inertia));

   WALBERLA_LOG_INFO("- Adding second body.");
   un->add( std::move(sp2) );
   const RigidBody & subb2 =  *(un->begin()+1);
   // Check relative and global position
   WALBERLA_CHECK_FLOAT_EQUAL(un->getPosition(),Vec3(0,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(un->getRelPosition(), Vec3(0,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getRelPosition(), Vec3(1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getPosition(), Vec3(1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb2.getRelPosition(), Vec3(-1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb2.getPosition(), Vec3(-1,0,0));

   // Check mass volume and inertia
   WALBERLA_CHECK_FLOAT_EQUAL(un->getVolume(), (real_t(8./3.)*real_t(math::pi)));
   WALBERLA_CHECK_FLOAT_EQUAL(un->getMass(), un->getVolume()*Material::getDensity(iron));
   // Mass of one sphere
   real_t masssphere = real_t(4./3.)*real_t(math::pi)*Material::getDensity(iron);
   Mat3 bodyinertia(real_t(2.0)*scalar_inertia, 0, 0, 0, real_t(2.0)*(scalar_inertia + masssphere),0, 0, 0, real_t(2.0)*(scalar_inertia + masssphere));
   WALBERLA_CHECK_FLOAT_EQUAL(un->getInertia(), bodyinertia);

   WALBERLA_CHECK_FLOAT_EQUAL( un->getLinearVel(), Vec3(0,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL( un->getAngularVel(), Vec3(0,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL( un->getAngularVel(), Vec3(0,0,0));

   WALBERLA_LOG_INFO("- Performing Rotation.");
   //Check values for rotated union
   Quat rotz30(Vec3(0,0,1), real_t(math::pi/6.0)); // rotate by 30 deg via z axis
   real_t sin30 = real_t(0.5);
   real_t cos30 = real_t(sqrt(3.0)/2.0);
   un->setOrientation(rotz30);
   WALBERLA_CHECK_FLOAT_EQUAL(un->getQuaternion(), rotz30);
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getRelPosition(), Vec3(1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getPosition(), Vec3(cos30,sin30,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb2.getRelPosition(), Vec3(-1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb2.getPosition(), Vec3(-cos30,-sin30,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getRelQuaternion(), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5)));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getQuaternion(), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5))*rotz30 );
   WALBERLA_CHECK_FLOAT_EQUAL(un->getInertia(), rotz30.toRotationMatrix()*bodyinertia*rotz30.toRotationMatrix().getTranspose());

   WALBERLA_LOG_INFO("- Applying velocities.");
   //Apply a linear velocity to the union
   un->setLinearVel(Vec3(1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getLinearVel(), Vec3(1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb2.getLinearVel(), Vec3(1,0,0));

   un->setAngularVel(Vec3(0,0,1));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getAngularVel(), Vec3(0,0,1));
   WALBERLA_CHECK_FLOAT_EQUAL(subb2.getAngularVel(), Vec3(0,0,1));
   WALBERLA_CHECK_FLOAT_EQUAL(subb1.getLinearVel(), Vec3(real_t(1.0-sin30),cos30,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb2.getLinearVel(), Vec3(real_t(1.0+sin30),-cos30,0));

   WALBERLA_LOG_INFO("- Constructing rotated union.");
   //Part 2: Construct exactly the same union, but now in a rotated state.
   auto un2  = std::make_unique<UnionType>(12, 0, Vec3(0,0,0), rotz30, false, true, false);
   auto sp12 = std::make_unique<Sphere>( 10, 0, Vec3( 1,0,0), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5)), real_t(1), iron, false, true, false );
   auto sp22 = std::make_unique<Sphere>( 11, 0, Vec3(-1,0,0), Quat(real_t(0.3), real_t(0.9),real_t(0.1),real_t(0.3)), real_t(1), iron, false, true, false );

   un2->add( std::move(sp12) );
   un2->add( std::move(sp22) );
   const RigidBody & subb12 = *un2->begin();
   const RigidBody & subb22 =  *(un2->begin()+1);
   // Check relative and global position and rotation
   WALBERLA_CHECK_FLOAT_EQUAL(un2->getPosition(),Vec3(0,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(un2->getRelPosition(), Vec3(0,0,0));

   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getRelPosition(), Vec3(cos30,-sin30, 0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getPosition(), Vec3(1,0,0));

   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getRelPosition(), Vec3(-cos30,sin30, 0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getPosition(), Vec3(-1,0,0));

   // Global quaternion and local
   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getQuaternion(), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5)));
   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getRelQuaternion(), subb12.getQuaternion()*rotz30.getInverse());
   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getQuaternion(), Quat(real_t(0.3), real_t(0.9),real_t(0.1),real_t(0.3)));
   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getRelQuaternion(), subb22.getQuaternion()*rotz30.getInverse());

   // Inertia tensor
   WALBERLA_CHECK_FLOAT_EQUAL(un2->getInertia(), bodyinertia);

   //Rotate again by 30 deg and perform the same checks as before
   un2->rotate(rotz30);
   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getPosition(), Vec3(cos30,sin30,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getPosition(), Vec3(-cos30,-sin30,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getQuaternion(), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5))*rotz30 );
   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getQuaternion(), Quat(real_t(0.3), real_t(0.9),real_t(0.1),real_t(0.3))*rotz30 );
   WALBERLA_CHECK_FLOAT_EQUAL(un->getInertia(), rotz30.toRotationMatrix()*bodyinertia*rotz30.toRotationMatrix().getTranspose());

   WALBERLA_LOG_INFO("- Applying velocities.");
   //Apply a linear velocity to the union
   un2->setLinearVel(Vec3(1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getLinearVel(), Vec3(1,0,0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getLinearVel(), Vec3(1,0,0));

   un2->setAngularVel(Vec3(0,0,1));
   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getAngularVel(), Vec3(0,0,1));
   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getAngularVel(), Vec3(0,0,1));
   WALBERLA_CHECK_FLOAT_EQUAL(subb12.getLinearVel(), Vec3(real_t(1.0-sin30), cos30, 0));
   WALBERLA_CHECK_FLOAT_EQUAL(subb22.getLinearVel(), Vec3(real_t(1.0+sin30), -cos30, 0));

}

void checkAABB(const AABB &a1, const AABB& a2){
   WALBERLA_CHECK_FLOAT_EQUAL(a1.xMin(), a2.xMin());
   WALBERLA_CHECK_FLOAT_EQUAL(a1.yMin(), a2.yMin());
   WALBERLA_CHECK_FLOAT_EQUAL(a1.zMin(), a2.zMin());
   WALBERLA_CHECK_FLOAT_EQUAL(a1.xMax(), a2.xMax());
   WALBERLA_CHECK_FLOAT_EQUAL(a1.yMax(), a2.yMax());
   WALBERLA_CHECK_FLOAT_EQUAL(a1.zMax(), a2.zMax());
}

void UnionAABB() {
   WALBERLA_LOG_INFO("- Performing AABB test");
   MaterialID iron = Material::find("iron");
   // Construct a union of two spheres in a non-rotated state and check all its parameters
   auto un  = std::make_unique<UnionType>(12, 0, Vec3(0,0,0), Quat(), false, true, false);
   auto sp1 = std::make_unique<Sphere>( 10, 0, Vec3( 1,0,0), Quat(real_t(0.5), real_t(0.5),real_t(0.5),real_t(0.5)), real_t(1), iron, false, true, false );
   auto sp2 = std::make_unique<Sphere>( 11, 0, Vec3(-1,0,0), Quat(), real_t(1), iron, false, true, false );
   un->add( std::move(sp1) );
   un->add( std::move(sp2) );
   AABB aabb = un->getAABB();
   WALBERLA_CHECK_FLOAT_EQUAL(aabb.minCorner(), Vec3(-2,-1,-1));
   WALBERLA_CHECK_FLOAT_EQUAL(aabb.maxCorner(), Vec3(2, 1, 1));

   Vec3 shift(10,10,10);
   un->setPosition(shift);
   aabb = un->getAABB();
   checkAABB(aabb, AABB(real_t(8),real_t(9),real_t(9),real_t(12),real_t(11),real_t(11)));

   Quat rotz30(Vec3(0,0,1), real_t(math::pi/6.0)); // rotate by 30 deg via z axis
   real_t sin30 = real_t(0.5);
   real_t cos30 = real_t(sqrt(3.0)/2.0);
   un->setOrientation(rotz30);
   aabb = un->getAABB();
   checkAABB(aabb, AABB(real_t(9.0-cos30),real_t(9.0-sin30),real_t(9), real_t(11.0+cos30), real_t(11.0+sin30), real_t(11)));

   aabb = un->begin()->getAABB();
   checkAABB(aabb, AABB(real_t(9.0+cos30),real_t(9.0+sin30),real_t(9), real_t(11.0+cos30), real_t(11.0+sin30), real_t(11)));

   aabb = (un->begin()+1)->getAABB();
   checkAABB(aabb, AABB(real_t(9.0-cos30),real_t(9.0-sin30),real_t(9), real_t(11.0-cos30), real_t(11.0-sin30), real_t(11)));

}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   UnionConstruction();
   UnionAABB();
   SnowManFallingOnPlane();



   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
