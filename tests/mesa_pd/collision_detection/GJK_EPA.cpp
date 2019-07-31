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
//! \file GJK_EPA.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/AnalyticCollisionFunctions.h>
#include <mesa_pd/collision_detection/EPA.h>
#include <mesa_pd/collision_detection/GJK.h>
#include <mesa_pd/collision_detection/Support.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/Sphere.h>

#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/math/Vector2.h"

namespace walberla {
namespace mesa_pd {

using namespace walberla::mesa_pd::collision_detection;
using namespace walberla::mesa_pd::collision_detection::analytic;

bool gjkEPAcollideHybrid(Support &geom1, Support &geom2, Vec3& normal, Vec3& contactPoint, real_t& penetrationDepth)
{
   // For more information on hybrid GJK/EPA see page 166 in "Collision Detecton in Interactive 3D
   // Environments" by Gino van den Bergen.

   //1. Run GJK with considerably enlarged objects.
   real_t margin = real_t(1e-6);
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
void checkContact(const Vec3& contactPoint1,
                  const Vec3& normal1,
                  const real_t penetrationDepth1,
                  const Vec3& contactPoint2,
                  const Vec3& normal2,
                  const real_t penetrationDepth2,
                  const Vec3& planeNormal,
                  const real_t accuracy = test_accuracy )
{
   WALBERLA_CHECK_LESS( fabs((normal1 - normal2).sqrLength()), accuracy*accuracy );
   WALBERLA_CHECK_LESS( fabs(penetrationDepth1- penetrationDepth2), accuracy );
   
   //Unfortunately position accuracy is one-two orders of magnitude lower...
   if(floatIsEqual(planeNormal.sqrLength(), real_t(0.0))){
      WALBERLA_CHECK_LESS( fabs((contactPoint1- contactPoint2).sqrLength()), real_t(1e4)*accuracy*accuracy  );
   }else{
      //check for containment in plane only.
      WALBERLA_CHECK_LESS( fabs(contactPoint1*planeNormal-contactPoint2*planeNormal), real_t(1e2)*accuracy );
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
void runCollisionDataTest(Support &rb1, Support &rb2, const Vec3& dir1, const real_t penetration_factor,
                          const Vec3& real_axis, const Vec3& witnesspoint, const Vec3& witnessmove, const Vec3& planeNormal, const real_t accuracy = test_accuracy){

   Vec3 org_pos = rb2.pos_; //Safe position

   Vec3 normal1, normal2;
   Vec3 pos1, pos2;
   real_t comp_pen_depth1, comp_pen_depth2;

   for(int j = 0; j < distancecount; j++){
      //move rb1.
      rb2.pos_ = (org_pos + depth[j]*dir1);
//      WALBERLA_LOG_INFO("Using depth: "+ std::to_string(depth[j]));
      //Compute collision between rb1 and rb2 and vice versa
      bool result1 = gjkEPAcollideHybrid(rb1, rb2, normal1, pos1, comp_pen_depth1);
//      WALBERLA_LOG_DEVEL( normal1 << " " << pos1 << " " <<  comp_pen_depth1);
      bool result2 = gjkEPAcollideHybrid(rb2, rb1, normal2, pos2, comp_pen_depth2);
//      WALBERLA_LOG_DEVEL( normal2 << " " << pos2 << " " <<  comp_pen_depth2);
      if(depth[j] > real_t(0.0)){
         WALBERLA_CHECK(result1);
         WALBERLA_CHECK(result2);
         //Check contact information
         checkContact( pos1, normal1, comp_pen_depth1,
                       witnesspoint + depth[j] * witnessmove, real_axis, -depth[j] * penetration_factor, planeNormal, accuracy );
         checkContact( pos2, normal2, comp_pen_depth2,
                       witnesspoint + depth[j] * witnessmove, real_t(-1.0)*real_axis, -depth[j] * penetration_factor, planeNormal, accuracy );
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
   using namespace walberla::mesa_pd::data;
   // Original SPHERE <-> SPHERE
   auto sp = data::Sphere(real_t(1));
   auto sp1 = Support( Vec3(0,0,0), Rot3(), sp);
   auto sp2 = Support( Vec3(1.5,0,0), Rot3(), sp);
   auto sp3 = Support( Vec3(3.0,0,0), Rot3(), sp);

   Vec3     normal;
   Vec3     contactPoint;
   real_t   penetrationDepth;


   WALBERLA_LOG_INFO("Original: SPHERE <-> SPHERE");
   WALBERLA_CHECK( !gjkEPAcollideHybrid(sp1, sp3, normal, contactPoint, penetrationDepth) );
   WALBERLA_CHECK(  gjkEPAcollideHybrid(sp1, sp2, normal, contactPoint, penetrationDepth) );
   checkContact( contactPoint,  normal, penetrationDepth,
                 Vec3(real_t(0.75), 0, 0), Vec3(real_t(-1.0), 0, 0), real_t(-0.5), Vec3(0,0,0) );

   //Testcase 01 Box Sphere
   WALBERLA_LOG_INFO("Test 01: BOX <-> SPHERE");
   real_t sqr3_inv = real_t(1.0)/std::sqrt(real_t(3.0));
   real_t coordinate= real_t(5.0)* sqr3_inv + real_t(5.0); // 5*(1+ (1/sqrt(3)))
   Box bx_101010(Vec3(10, 10, 10));
   Support box1_1(Vec3(0, 0, 0), Rot3(), bx_101010);
   Sphere s_5(real_t(5));
   Support sphere1_2(Vec3(coordinate, coordinate, coordinate), Rot3(), s_5);
   Vec3 wp1(real_t(5.0), real_t(5.0), real_t(5.0));
   Vec3 wpm1(sqr3_inv*real_t(-0.5), sqr3_inv*real_t(-0.5), sqr3_inv*real_t(-0.5));
   Vec3 axis1(-sqr3_inv, -sqr3_inv, -sqr3_inv);
   runCollisionDataTest(box1_1, sphere1_2, axis1, real_t(1.0), axis1, wp1, wpm1, Vec3(0,0,0));

   //Testcase 02 Box LongBox (touching plane)
   //Reuse box1_1
   WALBERLA_LOG_INFO("Test 02: BOX <-> LONG BOX");
   Box bx_3011(Vec3(real_t(30.0),1,1));
   Support box2_1(Vec3(real_t(20.0),0,0), Rot3(), bx_3011);
   Vec3 wp2(5, 0, 0);
   Vec3 wpm2(real_t(-0.5),0,0);
   Vec3 axis2(-1,0,0);
   runCollisionDataTest(box1_1, box2_1, axis2, real_t(1.0), axis2, wp2, wpm2, axis2);

   //Testcase 03 Sphere Sphere
   WALBERLA_LOG_INFO("Test 03: SPHERE <-> SPHERE");
   Support sphere3_1(Vec3(0,0,0), Rot3(),    s_5);
   Support sphere3_2(Vec3(real_t(10.0),0,0), Rot3(), s_5);
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
   Support box4_1( Vec3(real_t(5.0)*(real_t(1.0)+sqr2), real_t(-5.0), 0), Rot3(q4), bx_101010);
   Support box4_2( Vec3(0, 0, 0),                       Rot3(), bx_101010);
   Vec3 wp4(5, -5, 0);
   Vec3 wpm4(real_t(-0.25),real_t(+0.25),0);
   Vec3 collision_axis4(-sqr2/real_t(2.0),+sqr2/real_t(2.0),0);
   Vec3 axis4(-1, 0, 0);

   runCollisionDataTest(box4_2, box4_1, axis4, sqr2/real_t(2.0), collision_axis4, wp4, wpm4, Vec3(0,real_t(1.0),0));

   //Testcase 05 Cube and Long Box non-centric (touching plane)
   WALBERLA_LOG_INFO("Test 05: CUBE <-> LONG BOX (NON_CENTRIC)");
   Support box5_1(Vec3(0, 0, 0),     Rot3(), bx_101010);
   Support box5_2(Vec3(real_t(15.0),real_t(5.5), 0), Rot3(), bx_3011);
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
   Support box6_1(Vec3(real_t(5.0)*(real_t(1.0)+sqr6_3), 0, 0), Rot3(q6), bx_101010);
   Support box6_2(Vec3(0, 0, 0), Rot3(), bx_101010);
   Vec3 wp6(5, 0, 0);
   Vec3 wpm6(real_t(-0.5), 0, 0);
   Vec3 axis6(-1, 0, 0);
   runCollisionDataTest(box6_2, box6_1, axis6, real_t(1.0), axis6, wp6, wpm6, Vec3(0,0,0));

   //Testcase 07:
   // BOX <-> SPHERE
   WALBERLA_LOG_INFO("Test 07: BOX <-> SPHERE");
   Support sphere7_1(Vec3(0,0,0), Rot3(), s_5);
   Box bx_555( Vec3(5, 5, 5));
   Support box7_2(Vec3(0, 0,real_t(7.5)), Rot3(), bx_555);
   Vec3 wpm7(0, 0, real_t(-0.5));
   Vec3 wp7(0, 0, real_t(5.0));
   Vec3 axis7(0, 0,  real_t(-1.0));
   runCollisionDataTest(sphere7_1, box7_2, axis7, real_t(1.0), axis7, wp7, wpm7, Vec3(0,0,0));

   //Testcase 09:
   // ELLIPSOID <-> ELLIPSOID
   WALBERLA_LOG_INFO("Test 09: ELLIPSOID <-> ELLIPSOID");
   Ellipsoid ell9_1_(Vec3(10,5,5));
   Support ell9_1(Vec3(0,0,0), Rot3(), ell9_1_);
   Ellipsoid ell9_2_(Vec3(5,10,5));
   Support ell9_2(Vec3(15,0,0), Rot3(), ell9_2_);
   Vec3 wpm9(real_t(-0.5), 0, 0);
   Vec3 wp9(real_t(10), 0, 0);
   Vec3 axis9(real_t(-1.0), 0, 0);
   runCollisionDataTest(ell9_1, ell9_2, axis9, real_t(1.0), axis9, wp9, wpm9, Vec3(0,0,0));

}

void FuzzyEllipsoidEllipsoid(uint_t iterations)
{
   using namespace walberla::mesa_pd::data;

   Vec3     normal;
   Vec3     contactPoint;
   real_t   penetrationDepth;

   Ellipsoid el(Vec3(2,2,2));
   Support e1(Vec3(0,0,0), Rot3(), el);
   Support e2(Vec3(real_t(3.999),0,0), Rot3(), el);

   for (uint_t i = 0; i < iterations; ++i)
   {
      e1.rot_.rotate( Vec3(math::realRandom(), math::realRandom(), math::realRandom()) );
      e2.rot_.rotate( Vec3(math::realRandom(), math::realRandom(), math::realRandom()) );
      WALBERLA_CHECK(gjkEPAcollideHybrid(e1, e2, normal, contactPoint, penetrationDepth));
      WALBERLA_CHECK_FLOAT_EQUAL(normal, Vec3(real_t(-1),0,0));
      WALBERLA_CHECK_FLOAT_EQUAL(contactPoint, Vec3(real_t(3.999) * real_t(0.5), 0, 0));
      WALBERLA_CHECK_FLOAT_EQUAL(penetrationDepth, real_t(3.999) - real_t(4));
   }
}

void FuzzyEllipsoidBox(uint_t iterations)
{
   using namespace walberla::mesa_pd::data;

   auto normal               = Vec3();
   auto contactPoint         = Vec3();
   auto penetrationDepth     = real_t(0);

   auto gjk_normal           = Vec3();
   auto gjk_contactPoint     = Vec3();
   auto gjk_penetrationDepth = real_t(0);

   Box       bx(Vec3(2,2,2));
   Ellipsoid el(Vec3(2,2,2));
   Support   p1(Vec3(0,0,0), Rot3(), el);
   Support   p2(Vec3(real_t(3.999),0,0), Rot3(), bx);

   for (uint_t i = 0; i < iterations; ++i)
   {
      //if(i%100==0) WALBERLA_LOG_INFO(i);
      p1.rot_.rotate( Vec3(math::realRandom(), math::realRandom(), math::realRandom()) );
      p2.rot_.rotate( Vec3(math::realRandom(), math::realRandom(), math::realRandom()) );
      p2.pos_ = Vec3(4,0,0);
      bool bCollision = false;
      while (!bCollision)
      {
         p2.pos_ -= Vec3(real_t(0.1),0,0);
         if (p2.pos_[0]<real_t(0)) WALBERLA_ABORT("ERROR");
         bCollision = detectSphereBoxCollision(p1.pos_,
                                               el.getSemiAxes()[0],
                                               p2.pos_,
                                               bx.getEdgeLength(),
                                               p2.rot_,
                                               contactPoint,
                                               normal,
                                               penetrationDepth,
                                               real_t(1e-4));
      }

      WALBERLA_CHECK(gjkEPAcollideHybrid(p1, p2, gjk_normal, gjk_contactPoint, gjk_penetrationDepth), penetrationDepth);
      //WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(gjk_normal, normal, real_t(1e-2));
      //WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(gjk_contactPoint, contactPoint, real_t(1e-2));
      //WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(gjk_penetrationDepth, penetrationDepth, real_t(1e-2));
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   if (std::is_same<walberla::real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   MainTest();
   FuzzyEllipsoidEllipsoid(1000);
   FuzzyEllipsoidBox(1000);

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::mesa_pd::main( argc, argv );
}
