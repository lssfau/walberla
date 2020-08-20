#include <pe/basic.h>
#include "pe/utility/BodyCast.h"

#include "pe/Materials.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Union.h"
#include "pe/rigidbody/Ellipsoid.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Types.h"

#include "blockforest/Initialization.h"
#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"
#include <core/math/Random.h>
#include "core/math/Vector3.h"

#include <pe/raytracing/Ray.h>
#include <pe/raytracing/Intersects.h>

#include "pe/rigidbody/BodyStorage.h"
#include <core/timing/TimingTree.h>

#include <pe/utility/GetBody.h>

#include <sstream>
#include <tuple>

namespace walberla {
using namespace walberla::pe;
using namespace walberla::pe::raytracing;

typedef std::tuple<Box, Plane, Sphere, Capsule, Ellipsoid> BodyTuple ;

void SphereIntersectsTest()
{
   MaterialID iron = Material::find("iron");
   Sphere sp1(123, 1, Vec3(3,3,3), Quat(), 2, iron, false, true, false);
   real_t t;
   Vec3 n;
   
   // ray through the center
   Ray ray1(Vec3(3,-5,3), Vec3(0,1,0));
   WALBERLA_LOG_INFO("RAY -> SPHERE");
   
   WALBERLA_CHECK(intersects(&sp1, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(6));
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(0));
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(-1));
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0));
   
   // ray tangential
   Ray ray2(Vec3(3,-5,3), Vec3(0,7.5,real_t(std::sqrt(real_t(15))/real_t(2))).getNormalized());
   WALBERLA_CHECK(intersects(&sp1, ray2, t, n));
   
   // sphere behind ray origin
   Sphere sp2(123, 1, Vec3(3,-8,3), Quat(), 2, iron, false, true, false);
   WALBERLA_CHECK(!intersects(&sp2, ray1, t, n));
   
   // sphere around ray origin
   Sphere sp3(123, 1, Vec3(3,-5,3), Quat(), 2, iron, false, true, false);
   WALBERLA_CHECK(intersects(&sp3, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(2));
}

void PlaneIntersectsTest() {
   MaterialID iron = Material::find("iron");
   // plane with center 3,3,3 and parallel to y-z plane
   Plane pl1(1, 1, Vec3(3, 3, 3), Vec3(1, 0, 0), real_t(1.0), iron);
   
   Ray ray1(Vec3(-5,3,3), Vec3(1,0,0));
   real_t t;
   Vec3 n;
   
   WALBERLA_LOG_INFO("RAY -> PLANE");
   WALBERLA_CHECK(intersects(&pl1, ray1, t, n), "ray through center did not hit");
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(8), "distance between ray and plane is incorrect");
   
   Ray ray2(Vec3(-5,3,3), Vec3(1,0,-1).getNormalized());
   WALBERLA_CHECK(intersects(&pl1, ray2, t, n), "ray towards random point on plane didn't hit");
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(sqrt(real_t(128))), "distance between ray and plane is incorrect");
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(-1), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0), "incorrect normal calculated");
   
   Plane pl1neg(1, 1, Vec3(3, 3, 3), Vec3(-1, 0, 0), real_t(1.0), iron);
   WALBERLA_CHECK(intersects(&pl1neg, ray2, t, n), "ray towards random point on plane didn't hit");
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(-1), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0), "incorrect normal calculated");
   
   Ray ray3(Vec3(-5,3,3), Vec3(-1,0,0).getNormalized());
   Plane pl5(1, 1, Vec3(-7, 3, 3), Vec3(1, 0, 0), real_t(1.0), iron);
   WALBERLA_CHECK(intersects(&pl5, ray3, t, n), "ray towards random point on plane didn't hit");
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(2), "distance between ray and plane is incorrect");
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(1), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0), "incorrect normal calculated");
   
   // plane with center 3,3,3 and parallel to x-z plane
   Plane pl2(1, 1, Vec3(3, 3, 3), Vec3(0, 1, 0), real_t(1.0), iron);
   WALBERLA_CHECK(!intersects(&pl2, ray1, t, n), "ray parallel to plane shouldnt hit");
   
   // plane with center -10,3,3 and parallel to y-z plane
   Plane pl4(1, 1, Vec3(-10, 3, 3), Vec3(1, 0, 0), real_t(1.0), iron);
   WALBERLA_CHECK(!intersects(&pl4, ray1, t, n), "ray hit plane behind origin");
   
   Plane pl6(1, 1, Vec3(3, 3, 0), Vec3(-1, 0, 0), real_t(1.0), iron);
   Ray ray4(Vec3(0,0,5), Vec3(1, 0, -1).getNormalized());
   WALBERLA_CHECK(intersects(&pl6, ray4, t, n), "ray didnt hit");
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(-1), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0), "incorrect normal calculated");
}

void BoxIntersectsTest() {
   WALBERLA_LOG_INFO("RAY -> BOX");
   
   MaterialID iron = Material::find("iron");
   real_t t;
   Vec3 n;
   
   Box box1(127, 5, Vec3(0, -15, 0),  Quat(), Vec3(10, 10, 10), iron, false, true, false);
   Ray ray1(Vec3(3,-5,3), Vec3(0,1,0));
   WALBERLA_CHECK(!intersects(&box1, ray1, t, n));
   
   Box box2(128, 5, Vec3(0, -2, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   WALBERLA_CHECK(intersects(&box2, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, real_t(8), real_t(1e-7));
   
   Box box3(128, 5, Vec3(0, 5, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   WALBERLA_CHECK(intersects(&box3, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(5));
   
   Ray ray6(Vec3(-8,5,0), Vec3(1,0,0));
   WALBERLA_CHECK(intersects(&box3, ray6, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(3));
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(-1), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0), "incorrect normal calculated");
   
   Ray ray7(Vec3(8,5,0), Vec3(-1,0,0));
   WALBERLA_CHECK(intersects(&box3, ray7, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(3));
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(1), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0), "incorrect normal calculated");
   
   // ray origin within box
   Ray ray2(Vec3(-2,0,0), Vec3(1,0,1).getNormalized());
   WALBERLA_CHECK(intersects(&box3, ray2, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, real_t(7.0710), real_t(1e-4));
   
   Ray ray3(Vec3(3,-5,3), Vec3(2, real_t(-1.5), real_t(0.5)).getNormalized());
   Box box4(128, 5, Vec3(0, 8, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   WALBERLA_CHECK(!intersects(&box4, ray3, t, n));
   
   Ray ray4(Vec3(3,-5,3), Vec3(-2, 3, real_t(0.5)).getNormalized());
   WALBERLA_CHECK(intersects(&box4, ray4, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, real_t(9.7068), real_t(1e-4));
   
   Box box5(128, 5, Vec3(4, 0, 0), Quat(), Vec3(4, 4, 4), iron, false, true, false);
   box5.rotate(0,0,math::pi/4);
   Ray ray5(Vec3(0,1.5,0), Vec3(1,0,0));
   WALBERLA_CHECK(intersects(&box5, ray5, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, real_t(2.67157), real_t(1e-4));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[0], real_t(-0.707107), real_t(1e-5), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[1], real_t(0.707107), real_t(1e-5), "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0), "incorrect normal calculated");
}

void AABBIntersectsTest() {
   WALBERLA_LOG_INFO("RAY -> AABB");
   
   Ray ray1(Vec3(-5,5,5), Vec3(1,0,0));
   real_t t;
   
   AABB aabb(0,0,0,
             10,10,10);
   
   WALBERLA_CHECK(intersects(aabb, ray1, t));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(5));
   
   WALBERLA_CHECK(intersects(aabb, ray1, t, 1.0));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(4));
   
   Ray ray2(Vec3(-5,5,10.5), Vec3(1,0,0)); // ray shooting over aabb, but within padding passed to intersects
   WALBERLA_CHECK(intersects(aabb, ray1, t, 1.0));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(4));
}

void CapsuleIntersectsTest() {
   MaterialID iron = Material::find("iron");
   real_t t;
   Vec3 n;
   
   Capsule cp1(0, 0, Vec3(2,3,3), Quat(), real_t(2), real_t(2), iron, false, true, false);
   
   // ray through the center
   Ray ray1(Vec3(3,-5,3), Vec3(0,1,0));
   WALBERLA_LOG_INFO("RAY -> CAPSULE");
   
   WALBERLA_CHECK(intersects(&cp1, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(6));
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(0));
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(-1));
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0));
   
   Ray ray2(Vec3(-5,3,3), Vec3(1,0,0));
   WALBERLA_CHECK(intersects(&cp1, ray2, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(4));
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(-1));
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0));
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0));
}

void EllipsoidTest() {
   WALBERLA_LOG_INFO("RAY -> ELLIPSOID");
   
   MaterialID iron = Material::find("iron");
   real_t t;
   Vec3 n;
   
   Ellipsoid el1(0,0, Vec3(2,3,3), Quat(), Vec3(2,3,1), iron, false, true, false);
   
   Ray ray1(Vec3(-2,3,3), Vec3(1,0,0).getNormalized());
   WALBERLA_CHECK(intersects(&el1, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(2));
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(-1));
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0));
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(0));

   Ray ray2(Vec3(2,3,0), Vec3(0,0,-1).getNormalized());
   WALBERLA_CHECK(!intersects(&el1, ray2, t, n));
   
   Ray ray3(Vec3(2,3,5), Vec3(0,0,-1).getNormalized());
   WALBERLA_CHECK(intersects(&el1, ray3, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(1));
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], real_t(0));
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], real_t(0));
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], real_t(1));

   Ray ray4(Vec3(-2,real_t(2),real_t(2)), Vec3(1,real_t(0),real_t(0.5)).getNormalized());
   WALBERLA_CHECK(intersects(&el1, ray4, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, real_t(2.36809), real_t(1e-5));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[0], real_t(-0.78193), real_t(1e-5));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[1], real_t(-0.62324), real_t(1e-5));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[2], real_t(0.012265), real_t(1e-5));
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   SetBodyTypeIDs<BodyTuple>::execute();
   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );
   
   SphereIntersectsTest();
   PlaneIntersectsTest();
   BoxIntersectsTest();
   AABBIntersectsTest();
   CapsuleIntersectsTest();
   EllipsoidTest();
   
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
