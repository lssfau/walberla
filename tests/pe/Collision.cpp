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
//! \file Collision.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/fcd/AnalyticCollisionDetection.h"
#include "pe/utility/BodyCast.h"

#include "pe/contact/Contact.h"
#include "pe/fcd/SimpleFCD.h"
#include "pe/Materials.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Union.h"
#include "pe/rigidbody/UnionFactory.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Types.h"

#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"
#include "core/math/Vector2.h"

namespace walberla {
using namespace walberla::pe;
using walberla::pe::fcd::analytic::collide;

void checkContact(const Contact& c1, const Contact& c2)
{
   WALBERLA_CHECK_EQUAL( c1.getBody1(), c2.getBody1() );
   WALBERLA_CHECK_EQUAL( c1.getBody2(), c2.getBody2() );
   WALBERLA_CHECK_FLOAT_EQUAL( c1.getPosition(), c2.getPosition() );
   WALBERLA_CHECK_FLOAT_EQUAL( c1.getNormal(), c2.getNormal() );
   WALBERLA_CHECK_FLOAT_EQUAL( c1.getDistance(), c2.getDistance() );
}

void SphereTest()
{
   MaterialID iron = Material::find("iron");
   Sphere sp1(123, 1, Vec3(0,0,0), Quat(), 1, iron, false, true, false);
   Sphere sp2(124, 2, Vec3(real_t(1.5),0,0), Quat(), 1, iron, false, true, false);
   Sphere sp3(125, 3, Vec3(real_t(3.0),0,0), Quat(), 1, iron, false, true, false);
   Sphere sp4(124, 2, Vec3(0,real_t(1.5),0), Quat(), 1, iron, false, true, false);
   Plane  pl1(223, 1, Vec3(0,0,0), Vec3(1,1,1).getNormalized(), 0, iron);
   CylindricalBoundary cb1(333, 0, Vec3(-100,0,0), 2, iron);

   std::vector<Contact> contacts;
   fcd::AnalyticCollideFunctor< std::vector<Contact> > collideFunc(contacts);

   // SPHERE <-> SPHERE
   WALBERLA_LOG_INFO("SPHERE <-> SPHERE");
   WALBERLA_CHECK( !collideFunc(&sp1, &sp3) );
   WALBERLA_CHECK(  collideFunc(&sp1, &sp2) );
   checkContact( contacts.at(0),
                 Contact( &sp1, &sp2, Vec3(real_t(0.75), 0, 0), Vec3(-1, 0, 0), real_t(-0.5)) );


   // SPHERE <-> PLANE
   WALBERLA_LOG_INFO("SPHERE <-> PLANE");
   contacts.clear();
   WALBERLA_CHECK(  collideFunc(&sp1, &pl1) );
   WALBERLA_CHECK(  collideFunc(&sp2, &pl1) );
   checkContact( contacts.at(1),
                 Contact( &sp2, &pl1, Vec3(1,real_t(-0.5),real_t(-0.5)), Vec3(1, 1, 1).getNormalized(), real_t(-0.133974596215561181)) );
   // ORDER INVARIANCE
   contacts.clear();
   WALBERLA_CHECK(  collideFunc(&pl1, &sp2) );
   checkContact( contacts.at(0),
                 Contact( &sp2, &pl1, Vec3(1,real_t(-0.5),real_t(-0.5)), Vec3(1, 1, 1).getNormalized(), real_t(-0.133974596215561181)) );

   WALBERLA_CHECK( !collideFunc(&sp3, &pl1) );

   // SPHERE <-> CYLINDRICAL BOUNDARY
   WALBERLA_LOG_INFO("SPHERE <-> CYLINDRICAL BOUNDARY");
   contacts.clear();
   WALBERLA_CHECK(  !collideFunc(&sp1, &cb1) );
   WALBERLA_CHECK(  !collideFunc(&sp2, &cb1) );
   WALBERLA_CHECK(  collideFunc(&sp4, &cb1) );
   checkContact( contacts.at(0),
                 Contact( &sp4, &cb1, Vec3(0,real_t(2),real_t(0)), Vec3(0, -1, 0).getNormalized(), real_t(-0.5)) );
   cb1.rotateAroundOrigin( Vec3( 0,0,1), math::pi * real_t(0.25) );
   WALBERLA_CHECK(  !collideFunc(&sp1, &cb1) );
   WALBERLA_CHECK(  collideFunc(&sp2, &cb1) );
   WALBERLA_CHECK(  collideFunc(&sp4, &cb1) );
   const real_t xPos = real_t(3) / real_t(4) + real_t(2) / real_c(sqrt(real_t(2)));
   const real_t yPos = xPos - real_t(4) / real_c(sqrt(real_t(2)));
   const real_t dist = real_c(sqrt((xPos - real_t(1.5)) * (xPos - real_t(1.5)) + yPos * yPos)) - sp4.getRadius();
   checkContact( contacts.at(1),
                 Contact( &sp2, &cb1, Vec3(xPos, yPos, 0), Vec3(-1, +1, 0).getNormalized(), dist) );
   checkContact( contacts.at(2),
                 Contact( &sp4, &cb1, Vec3(yPos, xPos, 0), Vec3(+1, -1, 0).getNormalized(), dist) );
}

void BoxTest()
{
   MaterialID iron = Material::find("iron");
   Box b1(123, 0, Vec3(0,0,0), Quat(), Vec3(2,2,2), iron, false, true, false);
   Box b2(124, 0, Vec3(real_t(1.5),0,0), Quat(), Vec3(2,2,2), iron, false, true, false);
   Box b3(125, 0, Vec3(real_t(3.0),0,0), Quat(), Vec3(2,2,2), iron, false, true, false);
   Box b4(123, 0, Vec3(0,0,0), Quat(), Vec3(2,2,2), iron, false, true, false);
   b4.rotate( Vec3(1,1,0), real_t(atan(sqrt(2))) );

   Box b5(123, 0, Vec3(0,0,0), Quat(), Vec3(2,2,2), iron, false, true, false);
   b5.rotate( Vec3(0,0,1), real_t(math::pi * 0.25) );
   b5.rotate( Vec3(1,0,0), real_t(math::pi * 0.25) );

   std::vector<Contact> contacts;
   fcd::AnalyticCollideFunctor< std::vector<Contact> > collideFunc(contacts);

//   std::vector<Contact> contacts;

   // BOX <-> BOX
   WALBERLA_LOG_INFO("BOX <-> BOX");
   WALBERLA_CHECK( !collideFunc(&b1, &b3) );
//   WALBERLA_LOG_WARNING("contactPoint    : " << contactPoint);
//   WALBERLA_LOG_WARNING("contactNormal   : " << contactNormal);
//   WALBERLA_LOG_WARNING("penetrationDepth: " << penetrationDepth);
   WALBERLA_CHECK(  collideFunc(&b1, &b2) );
//   WALBERLA_LOG_WARNING("contactPoint    : " << contactPoint);
//   WALBERLA_LOG_WARNING("contactNormal   : " << contactNormal);
//   WALBERLA_LOG_WARNING("penetrationDepth: " << penetrationDepth);


   b4.setPosition( (Vec3(0,0,1) * real_t(sqrt(3)) + Vec3(0,0,1)) * 0.999);
   WALBERLA_CHECK( collideFunc(&b1, &b4) );
//   WALBERLA_LOG_WARNING("contactPoint    : " << contacts.back().getPosition());
//   WALBERLA_LOG_WARNING("contactNormal   : " << contacts.back().getNormal());
//   WALBERLA_LOG_WARNING("penetrationDepth: " << contacts.back().getDistance());

   b4.setPosition( (Vec3(0,0,1) * real_t(sqrt(3)) + Vec3(0,0,1)) * 1.001);
   WALBERLA_CHECK( !collideFunc(&b1, &b4) );

   b5.setPosition( (Vec3(0,0,1) * real_t(sqrt(3)) + Vec3(0,0,1)) * 0.99);
   WALBERLA_CHECK( collideFunc(&b1, &b5) );
//   WALBERLA_LOG_WARNING("contactPoint    : " << contacts.back().getPosition());
//   WALBERLA_LOG_WARNING("contactNormal   : " << contacts.back().getNormal());
//   WALBERLA_LOG_WARNING("penetrationDepth: " << contacts.back().getDistance());

   b5.setPosition( (Vec3(0,0,1) * real_t(sqrt(3)) + Vec3(0,0,1)) * 1.01);
   WALBERLA_CHECK( !collideFunc(&b1, &b5) );

   Sphere s1(126, 0, Vec3(real_t(1.5), real_t(1.5), real_t(1.5)), Quat(), 1, iron, false, true, false);
   WALBERLA_CHECK( collideFunc(&b1, &s1) );
//   WALBERLA_LOG_WARNING("contactPoint    : " << contactPoint);
//   WALBERLA_LOG_WARNING("contactNormal   : " << contactNormal);
//   WALBERLA_LOG_WARNING("penetrationDepth: " << penetrationDepth);
}

void CapsuleTest()
{
   MaterialID iron = Material::find("iron");
   Capsule c1(100, 100, Vec3(0,0,0), Quat(), 1, 2, iron, false, true, false);
   Sphere sp1(123, 123, Vec3(0,0,0), Quat(), 1, iron, false, true, false);

   std::vector<Contact> contacts;
   fcd::AnalyticCollideFunctor< std::vector<Contact> > collideFunc(contacts);

   // CAPSULE <-> SPHERE
   WALBERLA_LOG_INFO("CAPSULE <-> SPHERE");

   sp1.setPosition(0, real_t(1.9), 0);
   WALBERLA_CHECK( collideFunc(&c1, &sp1) );
//   WALBERLA_LOG_WARNING("contactPoint    : " << contacts.at(1).getPosition());
//   WALBERLA_LOG_WARNING("contactNormal   : " << contacts.at(1).getNormal());
//   WALBERLA_LOG_WARNING("penetrationDepth: " << contacts.at(1).getDistance());
   sp1.setPosition(0, real_t(2.1), 0);
   WALBERLA_CHECK( !collideFunc(&c1, &sp1) );

   // CAPSULE <-> PLANE
   WALBERLA_LOG_INFO("CAPSULE <-> PLANE");
   Plane pl1(124, 124, Vec3(0,0,real_t(-0.9)), Vec3(0,0,1), 0, iron);
   WALBERLA_CHECK( collideFunc(&c1, &pl1) );
//   WALBERLA_LOG_WARNING("contactPoint    : " << contacts.at(2).getPosition());
//   WALBERLA_LOG_WARNING("contactNormal   : " << contacts.at(2).getNormal());
//   WALBERLA_LOG_WARNING("penetrationDepth: " << contacts.at(2).getDistance());
   pl1.setPosition(0,0, real_t(-1.1));
   WALBERLA_CHECK( !collideFunc(&c1, &pl1) );

   Plane pl2(124, 124, Vec3(real_t(1.9),0,0), Vec3(-1,0,0), 0, iron);
   WALBERLA_CHECK( collideFunc(&c1, &pl2) );
//   WALBERLA_LOG_WARNING("contactPoint    : " << contacts.at(3).getPosition());
//   WALBERLA_LOG_WARNING("contactNormal   : " << contacts.at(3).getNormal());
//   WALBERLA_LOG_WARNING("penetrationDepth: " << contacts.at(3).getDistance());
   pl2.setPosition(real_t(2.1),0, 0);
   WALBERLA_CHECK( !collideFunc(&c1, &pl2) );

}

void CapsuleTest2()
{
   const real_t   static_cof  ( real_t(0.1) / 2 );   // Coefficient of static friction. Roughly 0.85 with high variation depending on surface roughness for low stresses. Note: pe doubles the input coefficient of friction for material-material contacts.
   const real_t   dynamic_cof ( static_cof ); // Coefficient of dynamic friction. Similar to static friction for low speed friction.
   MaterialID     material = createMaterial( "granular", real_t( 1.0 ), 0, static_cof, dynamic_cof, real_t( 0.5 ), 1, 1, 0, 0 );
   //create obstacle
   Capsule c1(100, 100, Vec3(10,10,0), Quat(), 3, 40, material, false, true, false);
   c1.rotate( Vec3(0,1,0), math::pi * real_t(0.5) );
   Sphere sp1(123, 123, Vec3(real_t(6.5316496854295262864), real_t(10.099999999999999645), real_t(0.46999999991564372914) ), Quat(), real_t(0.47), material, false, true, false);

   std::vector<Contact> contacts;

   WALBERLA_LOG_DEVEL( c1 );
   WALBERLA_LOG_DEVEL( sp1 );

   WALBERLA_LOG_INFO("CAPSULE TEST");
   Vec2 distance;
   distance[0] = (sp1.getPosition() - c1.getPosition())[0];
   distance[1] = (sp1.getPosition() - c1.getPosition())[1];

   std::cout << std::setprecision(10);
   WALBERLA_LOG_DEVEL("DISTANCE: " << distance.length());
   WALBERLA_LOG_DEVEL(" SPHERE <-> CAPSULE (ANALYTICAL) ");
   WALBERLA_LOG_DEVEL( collide(&sp1, &c1, contacts) );
   WALBERLA_LOG_WARNING("contactPoint    : " << contacts.at(0).getPosition());
   WALBERLA_LOG_WARNING("contactNormal   : " << contacts.at(0).getNormal());
   WALBERLA_LOG_WARNING("penetrationDepth: " << contacts.at(0).getDistance());
}

void UnionTest()
{
   using UnionT = Union<Sphere>;
   UnionT  un1(120, 0, Vec3(0,0,0), Quat(), false, true, false);
   UnionT  un2(121, 0, Vec3(real_t(1.5),0,0), Quat(), false, true, false);
   auto sp1 = createSphere(&un1, 123, Vec3(0,0,0), 1);
   auto sp2 = createSphere(&un2, 124, Vec3(real_t(1.5),0,0), 1);

   std::vector<Contact> contacts;

   using namespace walberla::pe::fcd;
   // SPHERE <-> SPHERE
   WALBERLA_LOG_INFO("UNION <-> UNION");
   AnalyticCollideFunctor< std::vector<Contact> > func(contacts);
   func(&un1, &un2);

   checkContact( contacts.at(0),
                 Contact( sp1, sp2, Vec3(real_t(0.75), 0, 0), Vec3(-1, 0, 0), real_t(-0.5)) );
}

typedef std::tuple<Box, Capsule, Plane, Sphere> BodyTuple ;

int main( int argc, char** argv )
{
    walberla::debug::enterTestMode();

    walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

    SetBodyTypeIDs<BodyTuple>::execute();

    SphereTest();
    BoxTest();
    CapsuleTest();
//    CapsuleTest2();
    UnionTest();

    return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}