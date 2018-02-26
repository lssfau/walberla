#include <pe/basic.h>
#include "pe/utility/BodyCast.h"

#include "pe/Materials.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Union.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Types.h"

#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include <pe/raytracing/Ray.h>
#include <pe/raytracing/Intersects.h>
#include <pe/raytracing/Raytracer.h>
#include <pe/raytracing/Color.h>
#include <pe/raytracing/ShadingFunctions.h>

#include <pe/ccd/HashGrids.h>
#include "pe/rigidbody/BodyStorage.h"
#include <core/timing/TimingTree.h>

#include <pe/utility/GetBody.h>

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::pe::raytracing;

typedef boost::tuple<Box, Plane, Sphere, Capsule> BodyTuple ;

void SphereIntersectsTest()
{
   MaterialID iron = Material::find("iron");
   Sphere sp1(123, 1, Vec3(3,3,3), Vec3(0,0,0), Quat(), 2, iron, false, true, false);
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
   Sphere sp2(123, 1, Vec3(3,-8,3), Vec3(0,0,0), Quat(), 2, iron, false, true, false);
   WALBERLA_CHECK(!intersects(&sp2, ray1, t, n));
   
   // sphere around ray origin
   Sphere sp3(123, 1, Vec3(3,-5,3), Vec3(0,0,0), Quat(), 2, iron, false, true, false);
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
   
   Box box1(127, 5, Vec3(0, -15, 0), Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   Ray ray1(Vec3(3,-5,3), Vec3(0,1,0));
   WALBERLA_CHECK(!intersects(&box1, ray1, t, n));
   
   Box box2(128, 5, Vec3(0, -2, 0), Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   WALBERLA_CHECK(intersects(&box2, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, real_t(8), real_t(1e-7));
   
   Box box3(128, 5, Vec3(0, 5, 0), Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
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
   
   Ray ray3(Vec3(3,-5,3), Vec3(2, -1.5, 0.5).getNormalized());
   Box box4(128, 5, Vec3(0, 8, 0), Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   WALBERLA_CHECK(!intersects(&box4, ray3, t, n));
   
   Ray ray4(Vec3(3,-5,3), Vec3(-2, 3, 0.5).getNormalized());
   WALBERLA_CHECK(intersects(&box4, ray4, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, real_t(9.7068), real_t(1e-4));
   
   Box box5(128, 5, Vec3(4, 0, 0), Vec3(0, 0, 0), Quat(), Vec3(4, 4, 4), iron, false, true, false);
   box5.rotate(0,0,math::M_PI/4);
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
   
   Capsule cp1(0, 0, Vec3(2,3,3), Vec3(0,0,0), Quat(), real_t(2), real_t(2), iron, false, true, false);
   
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

ShadingParameters customBodyToShadingParams(const BodyID body) {
   if (body->getID() == 10) {
      return greenShadingParams(body).makeGlossy(30);
   } else if (body->getID() == 7) {
      return greenShadingParams(body).makeGlossy(10);
   } else if (body->getID() == 9) {
      return darkGreyShadingParams(body).makeGlossy(50);
   } else if (body->getID() == 3) {
      return redShadingParams(body).makeGlossy(200);
   } else {
      return defaultBodyTypeDependentShadingParams(body);
   }
}

void RaytracerTest() {
   WALBERLA_LOG_INFO("Raytracer");
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,10,10,10), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");

   Lighting lighting(Vec3(0, 5, 8), // 8, 5, 9.5 gut für ebenen, 0,5,8
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4, 0.4, 0.4)); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(640), size_t(480),
                       49.13,
                       Vec3(-5,5,5), Vec3(-1,5,5), Vec3(0,0,1), //-5,5,5; -1,5,5
                       lighting,
                       Color(0.2,0.2,0.2),
                       real_t(2),
                       customBodyToShadingParams);

   MaterialID iron = Material::find("iron");
   
   //PlaneID xNegPlane = createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(5,0,0), iron);
   // xNegPlane obstructs only the top left sphere and intersects some objects

   //PlaneID xNegPlaneClose = createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(1,0,0), iron);
   
   // Test Scene v1 - Spheres, (rotated) boxes, confining walls, tilted plane in right bottom back corner
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), Vec3(0,10,0), iron); // left wall
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), Vec3(0,0,0), iron); // right wall
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), Vec3(0,0,0), iron); // floor
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), Vec3(0,0,10), iron); // ceiling
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(10,0,0), iron); // back wall
   createPlane(*globalBodyStorage, 0, Vec3(1,0,0), Vec3(0,0,0), iron); // front wall, should not get rendered

   createPlane(*globalBodyStorage, 0, Vec3(-1,1,1), Vec3(8,2,2), iron); // tilted plane in right bottom back corner
   
   //createSphere(*globalBodyStorage, *forest, storageID, 2, Vec3(6,9.5,9.5), real_t(0.5));
   createSphere(*globalBodyStorage, *forest, storageID, 3, Vec3(4,5.5,5), real_t(1));
   //createSphere(*globalBodyStorage, *forest, storageID, 6, Vec3(3,8.5,5), real_t(1));
   BoxID box = createBox(*globalBodyStorage, *forest, storageID, 7, Vec3(5,6.5,5), Vec3(2,4,3));
   if (box != NULL) box->rotate(0,math::M_PI/4,math::M_PI/4);
   //createBox(*globalBodyStorage, *forest, storageID, 8, Vec3(5,1,8), Vec3(2,2,2));
   // Test scene v1 end
   
   // Test scene v2 additions start
   //createBox(*globalBodyStorage, *forest, storageID, 9, Vec3(9,9,5), Vec3(1,1,10));
   createCapsule(*globalBodyStorage, *forest, storageID, 10, Vec3(3, 9, 1), real_t(0.5), real_t(7), iron);
   CapsuleID capsule = createCapsule(*globalBodyStorage, *forest, storageID, 11, Vec3(7, 3.5, 7.5), real_t(1), real_t(2), iron);
   if (capsule != NULL) capsule->rotate(0,math::M_PI/3,math::M_PI/4-math::M_PI/8);
   // Test scene v2 end
   
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
      ccd::HashGrids* hashgrids = blockIt->getData<ccd::HashGrids>(ccdID);
      hashgrids->update();
   }

   WALBERLA_LOG_INFO("--- CAPSULE");
   WALBERLA_LOG_INFO(" hash: " << capsule->getHash());
   WALBERLA_LOG_INFO(" grid: " << capsule->getGrid());
   WALBERLA_LOG_INFO(" AABB: " << capsule->getAABB());
   WALBERLA_LOG_INFO("--- ROTATED BOX");
   WALBERLA_LOG_INFO(" hash: " << box->getHash());
   WALBERLA_LOG_INFO(" grid: " << box->getGrid());
   WALBERLA_LOG_INFO(" AABB: " << box->getAABB());
   
   Ray ray(Vec3(-5, 6.5, 5), Vec3(1, 0, 0));
   real_t t = INFINITY;
   Vec3 n;
   BodyID body = NULL;
   
   // output info about grids
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
      ccd::HashGrids* hashgrids = blockIt->getData<ccd::HashGrids>(ccdID);
      for (auto grid: hashgrids->gridList_) {
         WALBERLA_LOG_INFO("--- GRID " << grid << " ---");
         WALBERLA_LOG_INFO(" cellSpan: " << grid->getCellSpan());
         WALBERLA_LOG_INFO(" dims:     " << grid->xCellCount_ << "/" << grid->yCellCount_ << "/" << grid->zCellCount_);
         WALBERLA_LOG_INFO(" items:    " << grid->bodyCount_);
         WALBERLA_LOG_INFO(" enlargement threshold: " << grid->enlargementThreshold_);
      }
      body = hashgrids->getClosestBodyIntersectingWithRay<BodyTuple>(ray, blockIt->getAABB(), t, n);
      if (body != NULL) {
         WALBERLA_LOG_INFO("body found");
      } else {
         WALBERLA_LOG_INFO("body not found");
      }
   }
      
   //raytracer.setTBufferOutputDirectory("tbuffer");
   //raytracer.setTBufferOutputEnabled(true);
   raytracer.setImageOutputDirectory("image");
   raytracer.setImageOutputEnabled(true);
   //raytracer.setLocalImageOutputEnabled(true);
   
   raytracer.rayTrace<BodyTuple>(0);
}

ShadingParameters customSpheresBodyToShadingParams(const BodyID body) {
   if (body->getTypeID() == Plane::getStaticTypeID()) {
      return greyShadingParams(body);
   }
   
   switch (body->getID()) {
      case 0:
         return blueShadingParams(body).makeGlossy(1);
      case 1:
         return blueShadingParams(body).makeGlossy(10);
      case 2:
         return blueShadingParams(body).makeGlossy(30);
      case 3:
         return blueShadingParams(body).makeGlossy(80);
      case 4:
         return whiteShadingParams(body);
      case 5:
         return lightGreyShadingParams(body);
      case 6:
         return greyShadingParams(body);
      case 7:
         return darkGreyShadingParams(body);
      case 8:
         return blackShadingParams(body).makeGlossy(100);
      case 9:
         return redShadingParams(body);
      case 10:
         return blueShadingParams(body);
      case 11:
         return violetShadingParams(body);
      case 12:
         return greenShadingParams(body);
      case 13:
         return greenShadingParams(body).makeGlossy(30);
      case 14:
         return blueShadingParams(body).makeGlossy(1000);
      default:
         return lightGreyShadingParams(body);
   }
}

void RaytracerSpheresTest() {
   WALBERLA_LOG_INFO("Raytracer");
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,10,10,10), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");

   Lighting lighting(Vec3(0, 5, 8), // 8, 5, 9.5 gut für ebenen, 0,5,8
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4, 0.4, 0.4)); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(640), size_t(480),
                       49.13,
                       Vec3(-5,5,5), Vec3(-1,5,5), Vec3(0,0,1), //-5,5,5; -1,5,5
                       lighting,
                       Color(0.2,0.2,0.2),
                       real_t(2),
                       customSpheresBodyToShadingParams);
   
   MaterialID iron = Material::find("iron");
   
   //PlaneID xNegPlane = createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(5,0,0), iron);
   // xNegPlane obstructs only the top left sphere and intersects some objects
   
   //PlaneID xNegPlaneClose = createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(1,0,0), iron);
   
   // Test Scene v1 - Spheres, (rotated) boxes, confining walls, tilted plane in right bottom back corner
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), Vec3(0,10,0), iron); // left wall
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), Vec3(0,0,0), iron); // right wall
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), Vec3(0,0,0), iron); // floor
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), Vec3(0,0,10), iron); // ceiling
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(10,0,0), iron); // back wall
   createPlane(*globalBodyStorage, 0, Vec3(1,0,0), Vec3(0,0,0), iron); // front wall, should not get rendered
   
   walberla::id_t id=0;
   for (int j=0; j<4; j++) {
      for (int i=0; i<4; i++) {
         createSphere(*globalBodyStorage, *forest, storageID, id, Vec3(6,real_c(i+1)*real_t(2),real_c(j+1)*real_t(2)), real_t(0.9));
         id++;
      }
   }
   
   
   raytracer.setImageOutputDirectory("image");
   raytracer.setImageOutputEnabled(true);
   
   raytracer.rayTrace<BodyTuple>(1);
}

void hashgridsPlayground() {
   /*using namespace walberla::pe::ccd;
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(-2,-2,-2,2,2,2), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   HashGrids::HashGrid hashgrid(1.0); // initialize a 4x4x4 hashgrid with cellspan 1.0 (can hold a block with 4x4x4 in size.)
   
   std::vector<Vec3> centers = {
      Vec3(-0.25, -0.75, 0.5),
      Vec3(0.5, -0.5, 0.5),
      Vec3(-0.25, -0.75, 0.5),
      Vec3(1.25, -0.25, 0.5),
      Vec3(-1.25, 0.25, 0.5)
   };
   std::vector<Vec3> lengths = {
      Vec3(1, 1, 1),
      Vec3(0.5, 0.5, 0.5),
      Vec3(0.5, 0.5, 0.5),
      Vec3(1, 1, 1),
      Vec3(1, 1, 1)
   };

   for (size_t i = 0; i < centers.size(); i++) {
      BoxID box = createBox(*globalBodyStorage, *forest, storageID, walberla::id_t(i), centers[i], lengths[i]);
      hashgrid.add(box);
      WALBERLA_LOG_INFO("Box " << box->getID() << " at " << centers[i] << " hashed with " << box->getHash());
   }

   WALBERLA_LOG_INFO("domain: " << forest->getDomain());

   real_t inf = std::numeric_limits<real_t>::max();

   BodyID closestBody = NULL;
   real_t t_closest = inf;
   real_t t;
   Vec3 n;
   
   std::vector<Ray> rays = {
      Ray(Vec3(-2.5, -1.25, 0.5), Vec3(20, 7, 0).getNormalized()),
      Ray(Vec3(1.5, -0.5, 0.5), Vec3(-20, -7, 0).getNormalized()),
      Ray(Vec3(-2.25, 1.75, 0.5), Vec3(1, -1, 0).getNormalized())
   };

   int i = 0;
   for (const Ray& ray: rays) {
      closestBody = NULL;
      t_closest = inf;
      
      WALBERLA_LOG_INFO("--- ray " << i << ": " << ray.getOrigin() << ", " << ray.getDirection());
      
      std::vector<BodyID> bodiesContainer;
      hashgrid.possibleRayIntersectingBodies(ray, forest->getDomain(), bodiesContainer);
      
      WALBERLA_LOG_INFO("Resulting bodies:");
      IntersectsFunctor intersectsFunc(ray, t, n);
      for (BodyID const& body : bodiesContainer) {
         bool intersects = SingleCast<BodyTuple, IntersectsFunctor, bool>::execute(body, intersectsFunc);
         if (intersects && t < t_closest) {
            closestBody = body;
            t_closest = t;
         }
         std::cout << body->getID() << " ";
      }
      std::cout << std::endl;
      
      if (closestBody != NULL) {
         WALBERLA_LOG_INFO("closest: " << closestBody->getID() << " with t = " << t_closest);
      } else {
         WALBERLA_LOG_INFO("could not find closest body");
      }
      
      closestBody = NULL;
      t_closest = inf;
      
      closestBody = hashgrid.getRayIntersectingBody<BodyTuple>(ray, forest->getDomain(), t, n);
      if (closestBody == NULL) {
         WALBERLA_LOG_INFO("couldnt find closest body with direct resolution");
      } else {
         WALBERLA_LOG_INFO("closest with direct resolution: " << closestBody->getID() << " with t = " << t);
      }
      i++;
   }*/
}

ShadingParameters customHashgridsBodyToShadingParams(const BodyID body) {
   switch (body->getID()) {
      case 96:
         return blueShadingParams(body);
      /*case 203:
         return redShadingParams(body);*/
      case 140:
         return whiteShadingParams(body);
      case 50:
         return greyShadingParams(body);
   }
   if (body->getTypeID() == Sphere::getStaticTypeID()) {
      return yellowShadingParams(body).makeGlossy();
   }
   return defaultBodyTypeDependentShadingParams(body);

}

void HashGridsTest(size_t boxes, size_t capsules, size_t spheres) {
#if defined(USE_NAIVE_INTERSECTION_FINDING)
   WALBERLA_LOG_INFO("Using naive method for intersection testing");
#endif
#if !defined(USE_NAIVE_INTERSECTION_FINDING) && !defined(COMPARE_NAIVE_AND_HASHGRIDS_RAYTRACING)
   WALBERLA_LOG_INFO("Using hashgrids for intersection testing");
#endif
#if defined(COMPARE_NAIVE_AND_HASHGRIDS_RAYTRACING)
   WALBERLA_LOG_INFO("Comparing hashgrids and naive method for intersection testing");
#endif
   WALBERLA_LOG_INFO("Generating " << boxes << " boxes, " << capsules << " capsules and " << spheres << " spheres");
   
   using namespace walberla::pe::ccd;
   WcTimingTree tt;
   tt.start("Setup");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,4,4,4), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   std::vector<IBlock*> blocks;
   forest->getBlocks(blocks);
   //for (auto block: blocks) {
   //   WALBERLA_LOG_INFO("block AABB: " << block->getAABB());
   //}
   
   const AABB& forestAABB = forest->getDomain();
   
   bool removeUnproblematic = false;
   std::vector<walberla::id_t> problematicBodyIDs = {165, 5, 31}; //{50, 44, 66, 155, 170, 51};
   std::vector<walberla::id_t> bodySIDs;
   
   //WALBERLA_LOG_INFO("boxCount: " << boxCount);
   //WALBERLA_LOG_INFO("capsuleCount: " << capsuleCount);
   
   //BoxID box = createBox(*globalBodyStorage, *forest, storageID, 1000, Vec3(2.1, 0.5, 0.5), Vec3(0.25, 0.25, 0.25));
   
   // generate bodies for test
   std::vector<BodyID> bodies;
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(real_t(0.1), real_t(0.2)); //0.2 0.5
      real_t x_min = math::realRandom(forestAABB.xMin()+len, forestAABB.xMax());
      real_t y_min = math::realRandom(forestAABB.yMin()+len, forestAABB.yMax());
      real_t z_min = math::realRandom(forestAABB.zMin()+len, forestAABB.zMax());
      //real_t z_min = len+0.1;
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x_min, y_min, z_min), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != NULL);
      //box->rotate(0, math::realRandom(real_t(0), real_t(1))*math::M_PI, math::realRandom(real_t(0), real_t(1))*math::M_PI);
      bodies.push_back(box_);
      bodySIDs.push_back(box_->getSystemID());
   }
   
   for (size_t i = 0; i < capsules; i++) {
      real_t len = math::realRandom(real_t(0.1), real_t(0.3)); // 0.2 0.5
      real_t radius = real_t(0.1);
      real_t maxlen = len + 2*radius;
      real_t x = math::realRandom(forestAABB.xMin()+maxlen, forestAABB.xMax());
      real_t y = math::realRandom(forestAABB.yMin()+maxlen, forestAABB.yMax());
      real_t z = math::realRandom(forestAABB.zMin()+maxlen, forestAABB.zMax());
      walberla::id_t id = walberla::id_t(boxes+i);
      CapsuleID capsule = createCapsule(*globalBodyStorage, *forest, storageID, id, Vec3(x, y, z), radius, len);
      WALBERLA_CHECK(capsule != NULL);
      capsule->rotate(0, math::realRandom(real_t(0), real_t(1))*math::M_PI, math::realRandom(real_t(0), real_t(1))*math::M_PI);
      bodies.push_back(capsule);
      bodySIDs.push_back(capsule->getSystemID());
   }
   
   for (size_t i = 0; i < spheres; i++) {
      real_t radius = math::realRandom(real_t(0.1), real_t(0.2)); // 0.2 0.3
      // forestAABB.xMax()-radius gerechtfertigt?
      real_t x = math::realRandom(forestAABB.xMin()+radius, forestAABB.xMax());
      real_t y = math::realRandom(forestAABB.yMin()+radius, forestAABB.yMax());
      real_t z = math::realRandom(forestAABB.zMin()+radius, forestAABB.zMax());
      walberla::id_t id = walberla::id_t(boxes+capsules+i);
      SphereID sphere = createSphere(*globalBodyStorage, *forest, storageID, id, Vec3(x, y, z), radius);
      WALBERLA_CHECK(sphere != NULL);
      bodies.push_back(sphere);
      bodySIDs.push_back(sphere->getSystemID());
   }
   
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
      ccd::HashGrids* hashgrids = blockIt->getData<ccd::HashGrids>(ccdID);
      hashgrids->update();
      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt) {
         if (removeUnproblematic && std::find(problematicBodyIDs.begin(), problematicBodyIDs.end(), bodyIt->getID()) == problematicBodyIDs.end()) {
            bodyIt->setPosition(-100, -100, -100);
         }
      }
   }
   
   MaterialID iron = Material::find("iron");
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), Vec3(0,forestAABB.yMax(),0), iron); // left wall
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), Vec3(0,forestAABB.yMin(),0), iron); // right wall
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), Vec3(0,0,forestAABB.zMin()), iron); // floor
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), Vec3(0,0,forestAABB.zMax()), iron); // ceiling
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(forestAABB.xMax(),0,0), iron); // back wall
   createPlane(*globalBodyStorage, 0, Vec3(1,0,0), Vec3(forestAABB.xMin(),0,0), iron); // front wall, should not get rendered
   
   
   //createBox(*globalBodyStorage, *forest, storageID, boxCount+capsuleCount+1, Vec3(2, 2, 2), Vec3(1, 1, 1));
   
   
   
   Lighting lighting(Vec3(forestAABB.xSize()/real_t(2)+1, forestAABB.ySize()/real_t(2),
                          real_t(2)*forestAABB.zMax()+2), // 8, 5, 9.5 gut für ebenen, 0,5,8
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4, 0.4, 0.4)); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(640), size_t(480),
                       49.13,
                       Vec3(forestAABB.xSize()/real_t(2), forestAABB.ySize()/real_t(2),
                            real_t(2)*forestAABB.zMax()),
                       Vec3(forestAABB.xSize()/real_t(2), forestAABB.ySize()/real_t(2),
                            0),
                       Vec3(0,1,0), //-5,5,5; -1,5,5
                       lighting,
                       Color(0.2,0.2,0.2),
                       real_t(2),
                       customHashgridsBodyToShadingParams);
#if defined(USE_NAIVE_INTERSECTION_FINDING)
   raytracer.setImageOutputDirectory("image/naive");
#endif
#if !defined(USE_NAIVE_INTERSECTION_FINDING) && !defined(COMPARE_NAIVE_AND_HASHGRIDS_RAYTRACING)
   raytracer.setImageOutputDirectory("image/hashgrids");
#endif
#if defined(COMPARE_NAIVE_AND_HASHGRIDS_RAYTRACING)
   raytracer.setImageOutputDirectory("image/comparison");
#endif
   raytracer.setImageOutputEnabled(true);
   raytracer.setFilenameTimestepWidth(12);
   tt.stop("Setup");
   WALBERLA_LOG_INFO("output to: " << (boxes*100000000 + capsules*10000 + spheres));
   raytracer.rayTrace<BodyTuple>(boxes*100000000 + capsules*10000 + spheres, &tt);
   
   /*HashGrids::HashGrid* grid;
   for (auto bodyId: problematicBodyIDs) {
      BodyID body = getBody(*globalBodyStorage, *forest, storageID, bodySIDs[bodyId]);
      //WALBERLA_LOG_INFO(box);
      WALBERLA_LOG_INFO("--- BODY " << bodyId << " ---");
      WALBERLA_LOG_INFO("minCorner " << bodyId << ": " << body->getAABB().minCorner());
      WALBERLA_LOG_INFO("hash " << bodyId << ": " << body->getHash());
      WALBERLA_LOG_INFO("cell " << bodyId << ": " << body->getCellId());
      WALBERLA_LOG_INFO("grid " << bodyId << ": " << body->getGrid());
      grid = (HashGrids::HashGrid*)body->getGrid();
   }
   WALBERLA_CHECK(grid != NULL);
   WALBERLA_LOG_INFO("grid dims: " << grid->xCellCount_ << "/" << grid->yCellCount_ << "/" << grid->zCellCount_);*/
   
   auto temp = tt.getReduced();
   WALBERLA_ROOT_SECTION() {
      std::cout << temp;
   }
}

Vec3 minCornerToGpos(const Vec3& minCorner, real_t lengths) {
   return minCorner + Vec3(lengths/2, lengths/2, lengths/2);
}

void HashGridsTestScene() {
#if defined(USE_NAIVE_INTERSECTION_FINDING)
   WALBERLA_LOG_INFO("Using naive method for intersection testing");
#endif
#if !defined(USE_NAIVE_INTERSECTION_FINDING) && !defined(COMPARE_NAIVE_AND_HASHGRIDS_RAYTRACING)
   WALBERLA_LOG_INFO("Using hashgrids for intersection testing");
#endif
#if defined(COMPARE_NAIVE_AND_HASHGRIDS_RAYTRACING)
   WALBERLA_LOG_INFO("Comparing hashgrids and naive method for intersection testing");
#endif
   using namespace walberla::pe::ccd;
   WcTimingTree tt;
   tt.start("Setup");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,8,8,8), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   std::vector<BodyID> bodies;
   
   // create bodies
   size_t id = 0;
   real_t len = 2;
   //bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(1.1, 1.1, 1.1), len), Vec3(len, len, len)));
   len = 0.6;
   //bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(0.6, 0.5, 0.5), len), Vec3(len, len, len)));
   //bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(1.7, 0.5, 0.5), len), Vec3(len, len, len)));
   //bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(7.2, 0.5, 0.5), len), Vec3(len, len, len)));
   //bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(5, 0.5, 0.5), len), Vec3(len, len, len)));
   
   real_t x_min = 0, y_min = 0;
   len = 1.2;
   real_t gap = 0.4;
   for (int i = 0; ; ++i) {
      x_min = forestAABB.xMin() + i*(gap+len);
      if (x_min > forestAABB.max(0)) {
         break;
      }
      for (int j = 0; ; ++j) {
         y_min = forestAABB.yMin() + j*(gap+len);
         if (y_min > forestAABB.max(1)) {
            break;
         }
         
         bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(x_min, y_min, 0), len), Vec3(len, len, len)));
      }
   }
   
   // update hashgrids
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
      ccd::HashGrids* hashgrids = blockIt->getData<ccd::HashGrids>(ccdID);
      hashgrids->update();
   }
   
   // output info about bodies
   /*for (auto body: bodies) {
      WALBERLA_LOG_INFO("--- BODY " << body->getID() << " ---");
      WALBERLA_LOG_INFO(" corners: " << body->getAABB().minCorner() << ", " << body->getAABB().maxCorner());
      WALBERLA_LOG_INFO(" hash: " << body->getHash());
      //WALBERLA_LOG_INFO(" body index: " << body->getCellId());
      WALBERLA_LOG_INFO(" grid: " << body->getGrid());
   }*/
   
   // output info about grids
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
      ccd::HashGrids* hashgrids = blockIt->getData<ccd::HashGrids>(ccdID);
      for (auto grid: hashgrids->gridList_) {
         WALBERLA_LOG_INFO("--- GRID " << grid << " ---");
         WALBERLA_LOG_INFO(" cellSpan: " << grid->getCellSpan());
         WALBERLA_LOG_INFO(" dims:     " << grid->xCellCount_ << "/" << grid->yCellCount_ << "/" << grid->zCellCount_);
         WALBERLA_LOG_INFO(" items:    " << grid->bodyCount_);
         WALBERLA_LOG_INFO(" enlargement threshold: " << grid->enlargementThreshold_);
      }
   }
   
   
   /*Ray ray(Vec3(4,4,10), Vec3(0,0,-1));
   real_t t;
   Vec3 n;
   
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
      ccd::HashGrids* hashgrids = blockIt->getData<ccd::HashGrids>(ccdID);
      BodyID body = hashgrids->getClosestBodyIntersectingWithRay<BodyTuple>(ray, blockIt->getAABB(), t, n);
      if (body != NULL) {
         WALBERLA_LOG_INFO("found body " << body->getID() << ", t = " << t);
         WALBERLA_LOG_INFO(" hash: " << body->getHash());
         WALBERLA_LOG_INFO(" grid: " << body->getGrid());
      } else {
         WALBERLA_LOG_INFO("did not find body");
      }
   }*/
   
   
   // setup raytracer
   Lighting lighting(Vec3(1,2,15),
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4, 0.4, 0.4)); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(480), size_t(480),
                       49.13,
                       Vec3(4,4,18), //6,6,8
                       Vec3(4,4,10), //6,6,4
                       Vec3(0,1,0),
                       lighting,
                       Color(0.2,0.2,0.2),
                       real_t(2));
   
#if defined(USE_NAIVE_INTERSECTION_FINDING)
   raytracer.setImageOutputDirectory("image/naive");
#endif
#if !defined(USE_NAIVE_INTERSECTION_FINDING) && !defined(COMPARE_NAIVE_AND_HASHGRIDS_RAYTRACING)
   raytracer.setImageOutputDirectory("image/hashgrids");
#endif
#if defined(COMPARE_NAIVE_AND_HASHGRIDS_RAYTRACING)
   raytracer.setImageOutputDirectory("image/comparison");
#endif
   raytracer.setImageOutputEnabled(true);
   raytracer.setFilenameTimestepWidth(12);
   tt.stop("Setup");
   WALBERLA_LOG_INFO("output to: 1");
   raytracer.rayTrace<BodyTuple>(1, &tt);
}

void GridCellIntersectionPlayground() {
   real_t cellSpan = 2;
   
   AABB blockAABB(Vec3(-2.5,-2.5,0), Vec3(3.5,5.5,12));
   
   size_t blockXCellCount = uint_c(std::ceil(blockAABB.xSize() / cellSpan));
   size_t blockYCellCount = uint_c(std::ceil(blockAABB.ySize() / cellSpan));
   size_t blockZCellCount = uint_c(std::ceil(blockAABB.zSize() / cellSpan));

   //Ray ray(Vec3(2,-0.5,0), Vec3(2,13,0).getNormalized());
   //Ray ray(Vec3(1.5,6,0), Vec3(-4,-9,0).getNormalized());
   //Ray ray(Vec3(2,-2,0), Vec3(5,8,0).getNormalized());
   //Ray ray(Vec3(2,-2,0), Vec3(-5,8,0).getNormalized());
   Ray ray(Vec3(0.5,-3,0), Vec3(-4,7,0).getNormalized());

   const real_t inf = std::numeric_limits<real_t>::max();

   real_t invCellSpan = real_t(1)/cellSpan;

   std::vector<Vec3> intersectionPoints;
   std::vector<real_t> distances;
   std::vector<size_t> axii;
   
   Vec3 firstPoint;
   if (blockAABB.contains(ray.getOrigin())) {
      firstPoint = ray.getOrigin();
   } else {
      real_t t_start;
      if (intersects(blockAABB, ray, t_start)) {
         firstPoint = ray.getOrigin() + ray.getDirection()*t_start;
      } else {
         WALBERLA_LOG_INFO("Ray does not intersect block.");
         return;
      }
   }
   
   Vector3<int32_t> firstCell(int32_c(std::floor(std::max(firstPoint[0]-blockAABB.xMin()*invCellSpan, real_t(0)))),
                              int32_c(std::floor(std::max(firstPoint[1]-blockAABB.yMin()*invCellSpan, real_t(0)))),
                              int32_c(std::floor(std::max(firstPoint[2]-blockAABB.zMin()*invCellSpan, real_t(0)))));

   const int8_t stepX = ray.xDir() >= 0 ? 1 : -1;
   const int8_t stepY = ray.yDir() >= 0 ? 1 : -1;
   const int8_t stepZ = ray.zDir() >= 0 ? 1 : -1;
   
   Vec3 near((stepX >= 0) ? (firstCell[0]+1)*cellSpan-firstPoint[0] : firstPoint[0]-firstCell[0]*cellSpan,
             (stepY >= 0) ? (firstCell[1]+1)*cellSpan-firstPoint[1] : firstPoint[1]-firstCell[1]*cellSpan,
             (stepZ >= 0) ? (firstCell[2]+1)*cellSpan-firstPoint[2] : firstPoint[2]-firstCell[2]*cellSpan);
   
   real_t tMaxX = (ray.xDir() != 0) ? std::abs(near[0]*ray.xInvDir()) : inf;
   real_t tMaxY = (ray.yDir() != 0) ? std::abs(near[1]*ray.yInvDir()) : inf;
   real_t tMaxZ = (ray.zDir() != 0) ? std::abs(near[2]*ray.zInvDir()) : inf;
   
   real_t tDeltaX = (ray.xDir() != 0) ? std::abs(cellSpan*ray.xInvDir()) : inf;
   real_t tDeltaY = (ray.yDir() != 0) ? std::abs(cellSpan*ray.yInvDir()) : inf;
   real_t tDeltaZ = (ray.zDir() != 0) ? std::abs(cellSpan*ray.zInvDir()) : inf;
   
   Vector3<int32_t> currentCell = firstCell;

   // First cell needs extra treatment, as it might lay out of the blocks upper bounds
   // due to the nature of how it is calculated: If the first point lies on a upper border
   // it maps to the cell "above" the grid.
   if (currentCell[0] < blockXCellCount && currentCell[1] < blockYCellCount && currentCell[2] < blockZCellCount) {
      WALBERLA_LOG_INFO("found valid cell at " << currentCell);
   }
   
   while (true) {
      if (tMaxX < tMaxY) {
         if (tMaxX < tMaxZ) {
            tMaxX += tDeltaX;
            currentCell[0] += stepX;
            if (currentCell[0] == blockXCellCount || currentCell[0] == -1) {
               break;
            }
         } else {
            tMaxZ += tDeltaZ;
            currentCell[2] += stepZ;
            if (currentCell[2] == blockZCellCount || currentCell[2] == -1) {
               break;
            }
         }
      } else {
         if (tMaxY < tMaxZ) {
            tMaxY += tDeltaY;
            currentCell[1] += stepY;
            if (currentCell[1] == blockYCellCount || currentCell[1] == -1) {
               break;
            }
         } else {
            tMaxZ += tDeltaZ;
            currentCell[2] += stepZ;
            if (currentCell[2] == blockZCellCount || currentCell[2] == -1) {
               break;
            }
         }
      }
      
      WALBERLA_LOG_INFO("found cell at " << currentCell);
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   SetBodyTypeIDs<BodyTuple>::execute();
   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );
   
   //SphereIntersectsTest();
   //PlaneIntersectsTest();
   //BoxIntersectsTest();
   //AABBIntersectsTest();
   //CapsuleIntersectsTest();
   RaytracerTest();
   //RaytracerSpheresTest();
   

   std::vector<size_t> boxes = {127, 70, 20, 150};
   std::vector<size_t> capsules = {127, 60, 140, 100};
   std::vector<size_t> spheres = {0, 50, 40, 120};
   
   //for (size_t i = 0; i < boxes.size(); ++i) {
   //   HashGridsTest(boxes[i], capsules[i], spheres[i]);
   //}
   //HashGridsTest(boxes[1], boxes[1], boxes[1]);
   
   //HashGridsTestScene();
   
   //GridCellIntersectionPlayground();
   
   //HashGridsTest();
   
   //hashgridsPlayground();

   return EXIT_SUCCESS;
}
