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
//! \file RigidBody.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/Materials.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/Types.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/cr/PlainIntegrator.h"
#include "core/DataTypes.h"
#include "core/math/Constants.h"

#include "core/debug/TestSubsystem.h"

namespace walberla {
using namespace walberla::pe;

void move( BodyStorage& storage, real_t dt )
{
   for (auto it = storage.begin(); it != storage.end(); ++it)
   {
      // Checking the state of the body
      WALBERLA_ASSERT( it->checkInvariants(), "Invalid capsule state detected" );
      WALBERLA_ASSERT( !it->hasSuperBody(), "Invalid superordinate body detected" );

      // Moving the capsule according to the acting forces (don't move a sleeping body)
      if( it->isAwake() ) {
         if( !it->hasInfiniteMass() ) {
            // Calculating the linear acceleration by the equation
            //   force * m^(-1) + gravity
            const Vec3 vdot( it->getForce() * it->getInvMass() );

            const Vec3 wdot( it->getInvInertia() * it->getTorque() );

            // Updating the linear velocity
            it->setLinearVel(it->getLinearVel() + vdot * dt);

            // Updating the angular velocity
            it->setAngularVel( it->getAngularVel() + wdot * dt );
         }

         // Calculating the translational displacement
         it->setPosition( it->getPosition() + it->getLinearVel() * dt );

         // Calculating the rotation angle
         const Vec3 phi( it->getAngularVel() * dt );

         // Calculating the new orientation
         it->rotate( Quat( phi, phi.length() ) );
         WALBERLA_ASSERT( realIsEqual( it->getRotation().getDeterminant(), real_c(1) ), "Corrupted rotation matrix determinant" );

         // Setting the axis-aligned bounding box
         it->calcBoundingBox();

         // Calculating the current motion
         it->calcMotion();
      }

      it->resetForceAndTorque();

      // Checking the state of the capsule
      WALBERLA_ASSERT( it->checkInvariants(), "Invalid capsule state detected" );
   }
}

void checkRotationFunctions()
{
   MaterialID iron = Material::find("iron");
   auto sp1 = std::make_shared<Sphere>( 0, 0, Vec3(0,0,0), Quat(), real_t(1), iron, false, true, false );
   auto sp2 = std::make_shared<Sphere>( 0, 0, Vec3(0,0,0), Quat(), real_t(1), iron, false, true, false );
   auto sp3 = std::make_shared<Sphere>( 0, 0, Vec3(0,0,0), Quat(), real_t(1), iron, false, true, false );
   auto sp4 = std::make_shared<Sphere>( 0, 0, Vec3(0,0,0), Quat(), real_t(1), iron, false, true, false );

   sp1->rotate( 1, 0, 0, math::pi * real_t(0.5));
   sp1->rotate( 0, 1, 0, math::pi * real_t(0.5));
   sp1->rotate( 0, 0, 1, math::pi * real_t(0.5));

   sp2->rotate( 1, 0, 0, math::pi * real_t(0.5));
   sp2->rotate( 0, 1, 0, math::pi * real_t(0.5));
   sp2->rotate( 0, 0, 1, math::pi * real_t(0.5));

   sp3->rotate( math::pi * real_t(0.5), math::pi * real_t(0.5), math::pi * real_t(0.5) );
   sp4->rotate( Vec3(math::pi * real_t(0.5), math::pi * real_t(0.5), math::pi * real_t(0.5)) );

   WALBERLA_CHECK_FLOAT_EQUAL( sp1->getQuaternion(), Quat(math::root_two * real_t(0.5), 0, math::root_two * real_t(0.5), 0) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp2->getQuaternion(), Quat(math::root_two * real_t(0.5), 0, math::root_two * real_t(0.5), 0) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp3->getQuaternion(), Quat(math::root_two * real_t(0.5), 0, math::root_two * real_t(0.5), 0) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp4->getQuaternion(), Quat(math::root_two * real_t(0.5), 0, math::root_two * real_t(0.5), 0) );

   WALBERLA_CHECK_FLOAT_EQUAL( sp1->getPosition(), Vec3(0, 0, 0) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp2->getPosition(), Vec3(0, 0, 0) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp3->getPosition(), Vec3(0, 0, 0) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp4->getPosition(), Vec3(0, 0, 0) );

   sp1->rotateAroundPoint( Vec3(-10, 0, 0), Vec3(0, 1, 0), math::pi );
   WALBERLA_CHECK_FLOAT_EQUAL( sp1->getPosition(), Vec3(-20, 0, 0) );

   sp2->rotateAroundPoint( Vec3(-10, 0, 0), Vec3(0, 0, 1), math::pi );
   WALBERLA_CHECK_FLOAT_EQUAL( sp2->getPosition(), Vec3(-20, 0, 0) );

   sp3->rotateAroundPoint( Vec3(0, 10, 0), Vec3(0, 0, 1), math::pi );
   WALBERLA_CHECK_FLOAT_EQUAL( sp3->getPosition(), Vec3(0, 20, 0) );
}

void checkPointFunctions()
{
   MaterialID iron = Material::find("iron");
   auto sp1 = std::make_shared<Sphere>( 0, 0, Vec3(10,10,10), Quat(), real_t(1), iron, false, true, false );

   WALBERLA_CHECK( sp1->containsPoint( 10, 10, 10 ) );
   WALBERLA_CHECK( sp1->containsPoint( real_c(10.9), 10, 10 ) );
   WALBERLA_CHECK( !sp1->containsPoint( real_c(11.1), 10, 10 ) );
   Vec3 d(1,1,1);
   d = d.getNormalized();
   WALBERLA_CHECK( sp1->containsRelPoint( d * real_c(0.9) ) );
   WALBERLA_CHECK( !sp1->containsRelPoint( d * real_c(1.1) ) );

   WALBERLA_CHECK( !sp1->isSurfaceRelPoint( d * real_c(0.9) ) );
   WALBERLA_CHECK( sp1->isSurfaceRelPoint( d * real_c(1.0) ) );
   WALBERLA_CHECK( !sp1->isSurfaceRelPoint( d * real_c(1.1) ) );

   WALBERLA_CHECK( !sp1->isSurfacePoint( real_c(10.9), real_c(10.0), real_c(10.0) ) );
   WALBERLA_CHECK( sp1->isSurfacePoint( real_c(11.0), real_c(10.0), real_c(10.0) ) );
   WALBERLA_CHECK( !sp1->isSurfacePoint( real_c(11.1), real_c(10.0), real_c(10.0) ) );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   MaterialID iron = Material::find("iron");
   BodyStorage storage;
   SpherePtr spPtr = std::make_unique<Sphere>(0, 0, Vec3(0,0,0), Quat(), real_t(1), iron, false, true, false);
   SphereID sphere = static_cast<SphereID>(&storage.add(std::move(spPtr)));

   Vec3 x0 = Vec3(-2,2,0);
   Vec3 v0 = Vec3(-1,-1,1);

   for (auto it = storage.begin(); it != storage.end(); ++it){
      it->rotateAroundPoint(Vec3(1,1,0), Vec3(0,0,1), math::pi);
      WALBERLA_CHECK_FLOAT_EQUAL(it->getPosition(), Vec3(2,2,0));
      it->rotateAroundOrigin(Vec3(0,0,1), real_c(0.5) * math::pi);
      WALBERLA_CHECK_FLOAT_EQUAL(it->getPosition(), x0);
      it->setLinearVel(v0);
   }
   real_t const dt = real_c(0.01);
   for (real_t t = dt; t < 10; t+=dt){
      move(storage, dt);
      WALBERLA_CHECK_FLOAT_EQUAL(sphere->getPosition(), x0 + v0*t);
      WALBERLA_CHECK_FLOAT_EQUAL(sphere->getLinearVel(), v0);
   }

   checkRotationFunctions();
   checkPointFunctions();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}