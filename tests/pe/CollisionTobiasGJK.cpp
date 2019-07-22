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
//! \file CollisionTobiasGJK.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================
#include "pe/Types.h"


#include "pe/contact/Contact.h"
//#include "pe/fcd/IFCD.h"
#include "pe/fcd/GenericFCD.h"
#include "pe/fcd/AnalyticCollisionDetection.h"
#include "pe/fcd/GJKEPACollideFunctor.h"
#include "pe/Materials.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Union.h"
#include "pe/rigidbody/UnionFactory.h"
#include "pe/rigidbody/Ellipsoid.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"


#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"
#include "core/math/Vector2.h"
#include "core/math/Constants.h"

#include "pe/collision/EPA.h"
#include "pe/collision/GJK.h"

namespace walberla {
using namespace walberla::pe;

typedef std::tuple<Box, Capsule, Plane, Sphere, Union<Sphere>, Union<Sphere, Union<Sphere>>, Ellipsoid> BodyTuple ;

bool gjkEPAcollideHybrid(GeomPrimitive &geom1, GeomPrimitive &geom2, Vec3& normal, Vec3& contactPoint, real_t& penetrationDepth)
{
   using namespace walberla::pe::fcd;
   // For more information on hybrid GJK/EPA see page 166 in "Collision Detecton in Interactive 3D
   // Environments" by Gino van den Bergen.

   //1. Run GJK with considerably enlarged objects.
   real_t margin = real_t(1e-4);
   GJK gjk;
   if(gjk.doGJKmargin(geom1, geom2, margin)){
      //2. If collision is possible perform EPA.
      //std::cerr << "Peforming EPA.";
      EPA epa;
      epa.useSphereOptimization( true );
      return epa.doEPAmargin(geom1, geom2, gjk, normal, contactPoint, penetrationDepth, margin);
   }else{
      return false;
   }
}

//Define Test values for different precision levels
#ifdef WALBERLA_DOUBLE_ACCURACY
static const int distancecount = 6;
static const real_t depth[distancecount] = {real_t(-1e-5), real_t(1e-5), real_t(1e-4), real_t(1e-2), real_t(0.1), real_t(1.0)};
static const real_t test_accuracy = real_t(1e-3);
#else
static const int distancecount = 3;
static const real_t depth[distancecount] = {real_t(1e-2), real_t(0.1), real_t(1.0)};
static const real_t test_accuracy = real_t(1e-2); //Single Precision is v. bad!
#endif


/** Compares Computed Contact c1 to analytical Contact c2,
 * and tests for equivalence.
 * The computed position must only be in the same plane, if planeNormal has not length 0. */
void checkContact(const Contact& c1, const Contact& c2, const Vec3& planeNormal, const real_t accuracy = test_accuracy )
{

   WALBERLA_CHECK_EQUAL( c1.getBody1(), c2.getBody1() );
   WALBERLA_CHECK_EQUAL( c1.getBody2(), c2.getBody2() );

   WALBERLA_CHECK_LESS( fabs((c1.getNormal() - c2.getNormal()).sqrLength()), accuracy*accuracy );
   WALBERLA_CHECK_LESS( fabs(c1.getDistance()- c2.getDistance()), accuracy );
   
   //Unfortunately position accuracy is one-two orders of magnitude lower...
   if(floatIsEqual(planeNormal.sqrLength(), real_t(0.0))){
      WALBERLA_CHECK_LESS( fabs((c1.getPosition()- c2.getPosition()).sqrLength()), real_t(1e4)*accuracy*accuracy  );
   }else{
      //check for containment in plane only.
      WALBERLA_CHECK_LESS( fabs(c1.getPosition()*planeNormal-c2.getPosition()*planeNormal), real_t(1e2)*accuracy );
   }
   
}

/** \brief Executes a test setup for collision data collection.
    * \param rb1 first rigid body
    * \param rb2 second rigid body
    * \param dir1 direction of rb2 moving towards rb1 (unit vector)
    * \param penetration_factor Increment of the penetration if rb2 is moved by dir1 (=1.0 in most cases)
    * \param real_axis Analytical collision normal (unit vector)
    * \param witnesspoint Analytical touching point of rb1 and rb2
    * \param witnessmove Movement of the touching point, if rb2 is moved by dir1
    * \param planeNormal The normal of the touching plane (if the touching point is unique,
    * a Vector of length 0.0 shall be passed)
    * \param accuracy Acceptance threshold
    * Before the test, rb1 and rb2 shall be in touching contact.
    * This function checks the collision data returned for different penetration depths and argument orders.
    */
void runCollisionDataTest(GeomPrimitive &rb1, GeomPrimitive &rb2, const Vec3& dir1, const real_t penetration_factor,
                          const Vec3& real_axis, const Vec3& witnesspoint, const Vec3& witnessmove, const Vec3& planeNormal, const real_t accuracy = test_accuracy){

   Vec3 org_pos = rb2.getPosition(); //Safe position

   Vec3 normal1, normal2;
   Vec3 pos1, pos2;
   real_t comp_pen_depth1, comp_pen_depth2;

   for(int j = 0; j < distancecount; j++){
      //move rb1.
      rb2.setPosition(org_pos + depth[j]*dir1);
      WALBERLA_LOG_INFO("Using depth: "+ std::to_string(depth[j]));
      //Compute collision between rb1 and rb2 and vice versa
      bool result1 = gjkEPAcollideHybrid(rb1, rb2, normal1, pos1, comp_pen_depth1);
      WALBERLA_LOG_DEVEL( normal1 << " " << pos1 << " " <<  comp_pen_depth1);
      bool result2 = gjkEPAcollideHybrid(rb2, rb1, normal2, pos2, comp_pen_depth2);
      WALBERLA_LOG_DEVEL( normal2 << " " << pos2 << " " <<  comp_pen_depth2);
      if(depth[j] > real_t(0.0)){
         WALBERLA_CHECK(result1);
         WALBERLA_CHECK(result2);
         //Check contact information
         checkContact( Contact( &rb1, &rb2, pos1, normal1, comp_pen_depth1),
                       Contact( &rb1, &rb2, witnesspoint + depth[j] * witnessmove, real_axis, -depth[j] * penetration_factor ), planeNormal, accuracy );
         checkContact( Contact( &rb2, &rb1, pos2, normal2, comp_pen_depth2),
                       Contact( &rb2, &rb1, witnesspoint + depth[j] * witnessmove, real_t(-1.0)*real_axis, -depth[j] * penetration_factor ), planeNormal, accuracy );
      }
      if(depth[j] < real_t(0.0)){
         WALBERLA_CHECK(!result1);
         WALBERLA_CHECK(!result2);
      }
   }
}

/** Test the GJK-EPA implementation on a variety of configuations 
 * and penetation depths */
void MainTest()
{
   MaterialID iron = Material::find("iron");

   // Original SPHERE <-> SPHERE
   Sphere sp1(123, 1, Vec3(0,0,0), Quat(), 1, iron, false, true, false);
   Sphere sp2(124, 2, Vec3(real_t(1.5),0,0), Quat(), 1, iron, false, true, false);
   Sphere sp3(125, 3, Vec3(real_t(3.0),0,0), Quat(), 1, iron, false, true, false);

   Vec3     normal;
   Vec3     contactPoint;
   real_t   penetrationDepth;


   WALBERLA_LOG_INFO("Original: SPHERE <-> SPHERE");
   WALBERLA_CHECK( !gjkEPAcollideHybrid(sp1, sp3, normal, contactPoint, penetrationDepth) );
   WALBERLA_CHECK(  gjkEPAcollideHybrid(sp1, sp2, normal, contactPoint, penetrationDepth) );
   checkContact( Contact( &sp1, &sp2, contactPoint,  normal, penetrationDepth),
                 Contact( &sp1, &sp2, Vec3(real_t(0.75), 0, 0), Vec3(real_t(-1.0), 0, 0), real_t(-0.5)), Vec3(0,0,0) );

   //Testcase 01 Box Sphere
   WALBERLA_LOG_INFO("Test 01: BOX <-> SPHERE");
   real_t sqr3_inv = real_t(1.0)/std::sqrt(real_t(3.0));
   real_t coordinate= real_t(5.0)* sqr3_inv + real_t(5.0); // 5*(1+ (1/sqrt(3)))
   Box box1_1(127, 5, Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   Sphere sphere1_2(130, 8, Vec3(coordinate, coordinate, coordinate), Quat(), 5, iron, false, true, false);
   Vec3 wp1(real_t(5.0), real_t(5.0), real_t(5.0));
   Vec3 wpm1(sqr3_inv*real_t(-0.5), sqr3_inv*real_t(-0.5), sqr3_inv*real_t(-0.5));
   Vec3 axis1(-sqr3_inv, -sqr3_inv, -sqr3_inv);
   runCollisionDataTest(box1_1, sphere1_2, axis1, real_t(1.0), axis1, wp1, wpm1, Vec3(0,0,0));

   //Testcase 02 Box LongBox (touching plane)
   //Reuse box1_1
   WALBERLA_LOG_INFO("Test 02: BOX <-> LONG BOX");
   Box box2_1(131, 9, Vec3(real_t(20.0),0,0),  Quat(), Vec3(real_t(30.0),1,1), iron, false, true, false);
   Vec3 wp2(5, 0, 0);
   Vec3 wpm2(real_t(-0.5),0,0);
   Vec3 axis2(-1,0,0);
   runCollisionDataTest(box1_1, box2_1, axis2, real_t(1.0), axis2, wp2, wpm2, axis2);

   //Testcase 03 Sphere Sphere
   WALBERLA_LOG_INFO("Test 03: SPHERE <-> SPHERE");
   Sphere sphere3_1(129, 7, Vec3(0,0,0), Quat(), 5, iron, false, true, false);
   Sphere sphere3_2(128, 6, Vec3(real_t(10.0),0,0), Quat(), 5, iron, false, true, false);
   Vec3 wp3(5, 0, 0);
   Vec3 wpm3(real_t(-0.5),0,0);
   Vec3 axis3(-1,0,0);
   runCollisionDataTest(sphere3_1, sphere3_2, axis3, real_t(1.0), axis3, wp3, wpm3, Vec3(0,0,0));

   //Testcase 04 Cube with turned Cube
   WALBERLA_LOG_INFO("Test 04: CUBE <-> TURNED CUBE");
   //compute rotation.
   real_t angle = walberla::math::pi/real_t(4.0);
   Vec3 zaxis(0, 0, 1);
   Quat q4(zaxis, angle);

   //create turned box
   real_t sqr2 = std::sqrt(real_t(2.0));
   Box box4_1(132, 10, Vec3(real_t(5.0)*(real_t(1.0)+sqr2), real_t(-5.0), 0), q4, Vec3(10, 10, 10), iron, false, true, false);
   Box box4_2(133, 11, Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   Vec3 wp4(5, -5, 0);
   Vec3 wpm4(real_t(-0.25),real_t(+0.25),0);
   Vec3 collision_axis4(-sqr2/real_t(2.0),+sqr2/real_t(2.0),0);
   Vec3 axis4(-1, 0, 0);

   runCollisionDataTest(box4_2, box4_1, axis4, sqr2/real_t(2.0), collision_axis4, wp4, wpm4, Vec3(0,real_t(1.0),0));

   //Testcase 05 Cube and Long Box non-centric (touching plane)
   WALBERLA_LOG_INFO("Test 05: CUBE <-> LONG BOX (NON_CENTRIC)");
   Box box5_1(133, 12, Vec3(0, 0, 0),  Quat(), Vec3(10, 10, 10), iron, false, true, false);
   Box box5_2(134, 13, Vec3(real_t(15.0),real_t(5.5), 0), Quat(), Vec3(real_t(30.0),1,1), iron, false, true, false);
   Vec3 wp5(real_t(3.75), 5, 0);
   Vec3 wpm5(0, real_t(-0.5), 0);
   Vec3 axis5(0, -1, 0);
   runCollisionDataTest(box5_1, box5_2, axis5, real_t(1.0), axis5, wp5, wpm5, axis5);  //check only for containment in plane.


   //Testcase 06:
   WALBERLA_LOG_INFO("Test 06: CUBE <-> TURNED CUBE 2");
   //compute rotation.

   real_t sqr6_2 = std::sqrt(real_t(2.0));
   real_t sqr6_3 = std::sqrt(real_t(3.0));
   real_t angle6 = std::acos(real_t(1.0)/sqr6_3); //acos(1/sqrt(3))
   Vec3 rot_axis6(0, real_t(1.0)/sqr6_2, -real_t(1.0)/sqr6_2);
   Quat q6(rot_axis6, angle6);

   //create turned box with pos = (5*(1+sqrt(3)), 0, 0)
   Box box6_1(136, 14, Vec3(real_t(5.0)*(real_t(1.0)+sqr6_3), 0, 0),  q6, Vec3(10, 10, 10), iron, false, true, false);
   Box box6_2(136, 15, Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   Vec3 wp6(5, 0, 0);
   Vec3 wpm6(real_t(-0.5), 0, 0);
   Vec3 axis6(-1, 0, 0);
   runCollisionDataTest(box6_2, box6_1, axis6, real_t(1.0), axis6, wp6, wpm6, Vec3(0,0,0));

   //Testcase 07:
   // BOX <-> SPHERE
   WALBERLA_LOG_INFO("Test 07: BOX <-> SPHERE");
   Sphere sphere7_1(137, 16, Vec3(0,0,0),  Quat(), 5, iron, false, true, false);
   Box box7_2(138, 17, Vec3(0, 0,real_t(7.5)), Quat(), Vec3(5, 5, 5), iron, false, true, false);
   Vec3 wpm7(0, 0, real_t(-0.5));
   Vec3 wp7(0, 0, real_t(5.0));
   Vec3 axis7(0, 0,  real_t(-1.0));
   runCollisionDataTest(sphere7_1, box7_2, axis7, real_t(1.0), axis7, wp7, wpm7, Vec3(0,0,0));

   //Testcase 08:
   // CAPSULE <-> CAPSULE
   WALBERLA_LOG_INFO("Test 08: CAPSULE <-> CAPSULE");
   Quat q8(Vec3(0,1,0), walberla::math::pi/real_t(2.0)); //creates a y-axis aligned capsule
   Capsule cap8_1(139, 18, Vec3(0,0,0), Quat(), real_t(4.0), real_t(10.0), iron, false, true, false);
   Capsule cap8_2(140, 19, Vec3(0,0, real_t(13.0)),  q8, real_t(4.0), real_t(10.0), iron, false, true, false);
   Vec3 wpm8(0, 0, real_t(-0.5));
   Vec3 wp8(0, 0, real_t(4.0));
   Vec3 axis8(0, 0,  real_t(-1.0));
   runCollisionDataTest(cap8_1, cap8_2, axis8, real_t(1.0), axis8, wp8, wpm8, Vec3(0,0,0));

   //Testcase 09:
   // ELLIPSOID <-> ELLIPSOID
   WALBERLA_LOG_INFO("Test 09: ELLIPSOID <-> ELLIPSOID");
   Ellipsoid ell9_1(141, 20, Vec3(0,0,0), Quat(), Vec3(10,5,5), iron, false, true, false);
   Ellipsoid ell9_2(142, 21, Vec3(15,0,0), Quat(), Vec3(5,10,5), iron, false, true, false);
   Vec3 wpm9(real_t(-0.5), 0, 0);
   Vec3 wp9(real_t(10), 0, 0);
   Vec3 axis9(real_t(-1.0), 0, 0);
   runCollisionDataTest(ell9_1, ell9_2, axis9, real_t(1.0), axis9, wp9, wpm9, Vec3(0,0,0));

}

/** Test the GJK-EPA implementation for a collision 
 *	of a plane and a body and test the interface calls. */
void PlaneTest()
{
   WALBERLA_LOG_INFO("PLANE AND INTERFACE TEST");
   MaterialID iron = Material::find("iron");
   fcd::GenericFCD<BodyTuple, fcd::GJKEPACollideFunctor> testFCD;

   Plane pl(1, 1, Vec3(0, 1, 0), Vec3(0, 1, 0), real_t(1.0), iron );
   Sphere sphere(2, 2, Vec3(0, real_t(1.9), 0), Quat(), 1, iron, false, true, false);
   Sphere sphere2(3, 3, Vec3(0, real_t(0.1), 0), Quat(), 1, iron, false, true, false);

   PossibleContacts pcs;

   pcs.push_back(std::pair<Sphere*, Sphere*>(&sphere, &sphere2));
   Contacts& container = testFCD.generateContacts(pcs);
   WALBERLA_CHECK(container.size() == 1);

   Contact &c = container.back();
   //
   WALBERLA_LOG_DEVEL( c.getDistance() << " " << c.getNormal() << " " << c.getPosition() );
   if(c.getBody1()->getID() == 2) {
      checkContact( c, Contact(&sphere, &sphere2,  Vec3(0, real_t(1), 0), Vec3(0, 1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else if (c.getBody1()->getID() == 3) {
      checkContact( c, Contact(&sphere2, &sphere, Vec3(0, real_t(1), 0), Vec3(0, -1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else {
      WALBERLA_ABORT("Unknown ID!");
   }
   pcs.clear();

   pcs.push_back(std::pair<Plane*, Sphere*>(&pl, &sphere));
   container = testFCD.generateContacts(pcs);
   WALBERLA_CHECK(container.size() == 1);

   c = container.back();
   //
   WALBERLA_LOG_DEVEL( c.getDistance() << " " << c.getNormal() << " " << c.getPosition() );
   if(c.getBody1()->getID() == 1) {
      checkContact( c, Contact(&pl, &sphere,  Vec3(0, real_t(0.95), 0), Vec3(0, -1, 0), real_t(-0.1)), Vec3(0,0,0));
   } else if (c.getBody1()->getID() == 2) {
      checkContact( c, Contact(&sphere, &pl, Vec3(0, real_t(0.95), 0), Vec3(0, 1, 0), real_t(-0.1)), Vec3(0,0,0));
   } else {
      WALBERLA_ABORT("Unknown ID!");
   }
   pcs.clear();

   pcs.push_back(std::pair<Sphere*, Plane*>(&sphere, &pl));

   container = testFCD.generateContacts(pcs);
   WALBERLA_CHECK(container.size() == 1);
   c = container.back();

   WALBERLA_LOG_DEVEL( c.getDistance() << " " << c.getNormal() << " " << c.getPosition() );
   if(c.getBody1()->getID() == 1) {
      checkContact( c, Contact(&pl, &sphere,  Vec3(0, real_t(0.95), 0), Vec3(0, -1, 0), real_t(-0.1)), Vec3(0,0,0));
   } else if (c.getBody1()->getID() == 2) {
      checkContact( c, Contact(&sphere, &pl, Vec3(0, real_t(0.95), 0), Vec3(0, 1, 0), real_t(-0.1)), Vec3(0,0,0));
   } else {
      WALBERLA_ABORT("Unknown ID!");
   }
}

/** Test the GJK-EPA implementation for a collision 
 *	of a union and a body and the interface calls. */
void UnionTest(){
   WALBERLA_LOG_INFO("UNION AND INTERFACE TEST");
   MaterialID iron = Material::find("iron");
   fcd::GenericFCD<BodyTuple, fcd::GJKEPACollideFunctor> testFCD;

   //A recursive union of three spheres is dropped on a box.
   Box box(179, 179, Vec3(0,0,0), Quat(), Vec3(real_t(10),real_t(2), real_t(10)), iron, false, true, false);


   using UnionT = Union<Sphere>;
   auto unsub = std::make_unique<UnionT>(192, 192, Vec3(0,real_t(3.8),0), Quat(), false, true, false);

   auto sp1 = createSphere(unsub.get(), 180, Vec3(-3,real_t(3.8),0), real_t(3.0));
   auto sp2 = createSphere(unsub.get(), 181, Vec3(3,real_t(3.8),0), real_t(3.0));

   //Create another union, and add sub union
   Union<Sphere, Union<Sphere>> un(193, 193, Vec3(0,0,0), Quat(), false, true, false);
   createSphere(&un, 182, Vec3(0,real_t(6),0), real_t(3.0));
   un.add(std::move(unsub));


   PossibleContacts pcs;
   pcs.push_back(std::pair<Union<Sphere,Union<Sphere>>*, Box*>(&un, &box));
   Contacts& container = testFCD.generateContacts(pcs);
   WALBERLA_CHECK(container.size() == 2);

   Contact &c = container.back();
   WALBERLA_LOG_DEVEL( c.getDistance() << " " << c.getNormal() << " " << c.getPosition() );
   if(c.getBody1()->getID() == 181) {
      checkContact( c, Contact(sp2, &box,  Vec3(real_t(3), real_t(0.9), 0), Vec3(0, 1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else if (c.getBody1()->getID() == 179) {
      checkContact( c, Contact(&box, sp2,  Vec3(real_t(3), real_t(0.9), 0), Vec3(0, -1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else {
      WALBERLA_ABORT("Unknown ID!");
   }
   container.pop_back();


   c = container.back();
   WALBERLA_LOG_DEVEL( c.getDistance() << " " << c.getNormal() << " " << c.getPosition() );
   if(c.getBody1()->getID() == 180) {
      checkContact( c, Contact(sp1, &box,  Vec3(real_t(-3), real_t(0.9), 0), Vec3(0, 1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else if (c.getBody1()->getID() == 179) {
      checkContact( c, Contact(&box, sp1,  Vec3(real_t(-3), real_t(0.9), 0), Vec3(0, -1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else {
      WALBERLA_ABORT("Unknown ID!");
   }
   pcs.clear();

   //Vice Versa
   pcs.push_back(std::pair<Box*, Union<Sphere, Union<Sphere>>* >(&box, &un));
   container = testFCD.generateContacts(pcs);
   WALBERLA_CHECK(container.size() == 2);

   c = container.back();
   WALBERLA_LOG_DEVEL( c.getDistance() << " " << c.getNormal() << " " << c.getPosition() );
   if(c.getBody1()->getID() == 181) {
      checkContact( c, Contact(sp2, &box,  Vec3(real_t(3), real_t(0.9), 0), Vec3(0, 1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else if (c.getBody1()->getID() == 179) {
      checkContact( c, Contact(&box, sp2,  Vec3(real_t(3), real_t(0.9), 0), Vec3(0, -1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else {
      WALBERLA_ABORT("Unknown ID!");
   }
   container.pop_back();

   c = container.back();
   WALBERLA_LOG_DEVEL( c.getDistance() << " " << c.getNormal() << " " << c.getPosition() );
   if(c.getBody1()->getID() == 180) {
      checkContact( c, Contact(sp1, &box,  Vec3(real_t(-3), real_t(0.9), 0), Vec3(0, 1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else if (c.getBody1()->getID() == 179) {
      checkContact( c, Contact(&box, sp1,  Vec3(real_t(-3), real_t(0.9), 0), Vec3(0, -1, 0), real_t(-0.2)), Vec3(0,0,0));
   } else {
      WALBERLA_ABORT("Unknown ID!");
   }
   pcs.clear();

}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   SetBodyTypeIDs<BodyTuple>::execute();
   MainTest();
   PlaneTest();
   UnionTest();
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}