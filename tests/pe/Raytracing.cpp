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
   Lighting lighting(Vec3(0, 5, 8), // 8, 5, 9.5 gut für ebenen, 0,5,8
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4, 0.4, 0.4)); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage,
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
   
   createSphere(*globalBodyStorage, *forest, storageID, 2, Vec3(6,9.5,9.5), real_t(0.5));
   createSphere(*globalBodyStorage, *forest, storageID, 3, Vec3(4,5.5,5), real_t(1));
   createSphere(*globalBodyStorage, *forest, storageID, 6, Vec3(3,8.5,5), real_t(1));
   BoxID box = createBox(*globalBodyStorage, *forest, storageID, 7, Vec3(5,6.5,5), Vec3(2,4,3));
   if (box != NULL) box->rotate(0,math::M_PI/4,math::M_PI/4);
   createBox(*globalBodyStorage, *forest, storageID, 8, Vec3(5,1,8), Vec3(2,2,2));
   // Test scene v1 end
   
   // Test scene v2 additions start
   createBox(*globalBodyStorage, *forest, storageID, 9, Vec3(9,9,5), Vec3(1,1,10));
   createCapsule(*globalBodyStorage, *forest, storageID, 10, Vec3(3, 9, 1), real_t(0.5), real_t(7), iron);
   CapsuleID capsule = createCapsule(*globalBodyStorage, *forest, storageID, 11, Vec3(7, 3.5, 7.5), real_t(1), real_t(2), iron);
   if (capsule != NULL) capsule->rotate(0,math::M_PI/3,math::M_PI/4-math::M_PI/8);
   // Test scene v2 end
   
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
   Lighting lighting(Vec3(0, 5, 8), // 8, 5, 9.5 gut für ebenen, 0,5,8
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4, 0.4, 0.4)); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage,
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
   
   raytracer.rayTrace<BodyTuple>(0);
}

void playground() {
   AABB gridAABB(0,0,0,3,3,3);
   AABB blockAABB(-1,-1,-1,2,2,2);
   real_t c = 1; // cell width / height
   real_t c_inv = real_t(1) / c;
   size_t n_x = 3;
   size_t n_y = 3;
   size_t n_z = 3;
   
   Ray ray(Vec3(-1, -2, 0), Vec3(2, 9, 0).getNormalized());
   
   for (int i = 0; i < n_x; i++) {
      real_t xValue = blockAABB.minCorner()[0] + i*c;
      real_t lambda = (xValue - ray.getOrigin()[0]) * ray.getInvDirection()[0];
      real_t yValue = ray.getOrigin()[1] + lambda * ray.getDirection()[1];
      real_t zValue = ray.getOrigin()[2] + lambda * ray.getDirection()[2];
      if (yValue > blockAABB.maxCorner()[1] || yValue < blockAABB.minCorner()[1] ||
          zValue > blockAABB.maxCorner()[2] || zValue < blockAABB.minCorner()[2] ||
          lambda != lambda) {
         WALBERLA_LOG_INFO("P_x" << i << " = (" << xValue << "/" << yValue << "/" << zValue << ") invalid");
      } else {
         size_t xIndex = size_t(xValue * c_inv) % n_x;
         size_t yIndex = size_t(yValue * c_inv) % n_y;
         size_t zIndex = size_t(zValue * c_inv) % n_z;
         if (xValue < 0) {
            xIndex = n_x - 1 - (size_t(-xValue * c_inv) % n_x);
         }
         if (yValue < 0) {
            yIndex = n_y - 1 - (size_t(-yValue * c_inv) % n_y);
         }
         if (zValue < 0) {
            zIndex = n_z - 1 - (size_t(-zValue * c_inv) % n_z);
         }
         size_t arrayIndex = xIndex + yIndex*n_x + zIndex*n_x*n_y;
         WALBERLA_LOG_INFO("P_x" << i << " = (" << xValue << "/" << yValue << "/" << zValue << ") maps to cell " << arrayIndex);
         WALBERLA_LOG_INFO("\t xIndex = " << xIndex << "; yIndex = " << yIndex << "; zIndex = " << zIndex);
      }
   }
   
   for (int i = 0; i < n_y; i++) {
      real_t yValue = blockAABB.minCorner()[1] + i*c;
      real_t lambda = (yValue - ray.getOrigin()[1]) * ray.getInvDirection()[1];
      real_t xValue = ray.getOrigin()[0] + lambda * ray.getDirection()[0];
      real_t zValue = ray.getOrigin()[2] + lambda * ray.getDirection()[2];
      if (xValue > blockAABB.maxCorner()[0] || xValue < blockAABB.minCorner()[0] ||
          zValue > blockAABB.maxCorner()[2] || zValue < blockAABB.minCorner()[2] ||
          lambda != lambda) {
         WALBERLA_LOG_INFO("P_y" << i << " = (" << xValue << "/" << yValue << "/" << zValue << ") invalid");
      } else {
         size_t xIndex = size_t(xValue * c_inv) % n_x;
         size_t yIndex = size_t(yValue * c_inv) % n_y;
         size_t zIndex = size_t(zValue * c_inv) % n_z;
         if (xValue < 0) {
            xIndex = n_x - 1 - (size_t(-xValue * c_inv) % n_x);
         }
         if (yValue < 0) {
            yIndex = n_y - 1 - (size_t(-yValue * c_inv) % n_y);
         }
         if (zValue < 0) {
            zIndex = n_z - 1 - (size_t(-zValue * c_inv) % n_z);
         }
         size_t arrayIndex = xIndex + yIndex*n_x + zIndex*n_x*n_y;
         WALBERLA_LOG_INFO("P_y" << i << " = (" << xValue << "/" << yValue << "/" << zValue << ") maps to cell " << arrayIndex);
         WALBERLA_LOG_INFO("\t xIndex = " << xIndex << "; yIndex = " << yIndex << "; zIndex = " << zIndex);
      }
   }
   
   for (int i = 0; i < n_z; i++) {
      real_t zValue = blockAABB.minCorner()[2] + i*c;
      real_t lambda = (zValue - ray.getOrigin()[2]) * ray.getInvDirection()[2];
      real_t xValue = ray.getOrigin()[0] + lambda * ray.getDirection()[0];
      real_t yValue = ray.getOrigin()[1] + lambda * ray.getDirection()[1];
      if (xValue > blockAABB.maxCorner()[0] || xValue < blockAABB.minCorner()[0] ||
          yValue > blockAABB.maxCorner()[1] || yValue < blockAABB.minCorner()[1] ||
          lambda != lambda) {
         WALBERLA_LOG_INFO("P_z" << i << " = (" << xValue << "/" << yValue << "/" << zValue << ") invalid");
      } else {
         size_t xIndex = size_t(xValue * c_inv) % n_x;
         size_t yIndex = size_t(yValue * c_inv) % n_y;
         size_t zIndex = size_t(zValue * c_inv) % n_z;
         if (xValue < 0) {
            xIndex = n_x - 1 - (size_t(-xValue * c_inv) % n_x);
         }
         if (yValue < 0) {
            yIndex = n_y - 1 - (size_t(-yValue * c_inv) % n_y);
         }
         if (zValue < 0) {
            zIndex = n_z - 1 - (size_t(-zValue * c_inv) % n_z);
         }
         size_t arrayIndex = xIndex + yIndex*n_x + zIndex*n_x*n_y;
         WALBERLA_LOG_INFO("P_z" << i << " = (" << xValue << "/" << yValue << "/" << zValue << ") maps to cell " << arrayIndex);
         WALBERLA_LOG_INFO("\t xIndex = " << xIndex << "; yIndex = " << yIndex << "; zIndex = " << zIndex);
      }
   }
}

void hashgridsPlayground() {
   using namespace walberla::pe::ccd;
   
   /*shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(-8,-8,-8,8,8,8), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");

   HashGrids::HashGrid hashgrid(1.0); // initialize a 16x16x16 hashgrid with cellspan 1.0 (can hold a block with 16x16x16 in size.)
   
   std::vector<Vec3> minCorners;
   minCorners.push_back(Vec3(0,0,0));
   minCorners.push_back(Vec3(0.1,0.1,0.1));
   minCorners.push_back(Vec3(0.2,0.2,0.2));
   minCorners.push_back(Vec3(0.3,0.3,0.3));
   minCorners.push_back(Vec3(0.9,0.9,0.9));
   minCorners.push_back(Vec3(1,1,1));
   minCorners.push_back(Vec3(1.1,1.1,1.1));
   minCorners.push_back(Vec3(1.4,1.4,1.4));
   minCorners.push_back(Vec3(1.6,1.6,1.6));

   minCorners.push_back(Vec3(-8+1e-5,-8+1e-5,-8+1e-5));
   

   Vec3 lengths(0.5,0.5,0.5);
   for (auto minCorner: minCorners) {
      BoxID box = createBox(*globalBodyStorage, *forest, storageID, 9, minCorner+lengths/2, lengths);
      if (box == NULL) {
         WALBERLA_LOG_INFO("could not create box at " << minCorner);
         continue;
      }
      hashgrid.add(box);
      WALBERLA_LOG_INFO("hash of box at " << minCorner << " (" << box->getAABB().minCorner() << "): " << box->getHash());
   }*/
   
   HashGrids::HashGrid hashgrid(1.0); // initialize a 4x4x4 hashgrid with cellspan 1.0 (can hold a block with 4x4x4 in size.)
   
   const AABB blockAABB(-2,-2,-2,2,2,2);
   
   WALBERLA_LOG_INFO("ray:");
   Ray ray(Vec3(-3, -1.6, 0), Vec3(4, 1, 0).getNormalized());
   hashgrid.possibleRayIntersectingBodies(ray, blockAABB);
   
   WALBERLA_LOG_INFO("ray2:");
   Ray ray2(Vec3(-2.1, -2.1, 0), Vec3(1, 1, 0).getNormalized());
   hashgrid.possibleRayIntersectingBodies(ray2, blockAABB);
   
   WALBERLA_LOG_INFO("ray3:");
   Ray ray3(Vec3(3, -1, 0), Vec3(-7, 2, 0).getNormalized());
   hashgrid.possibleRayIntersectingBodies(ray3, blockAABB);
   
   WALBERLA_LOG_INFO("");
   WALBERLA_LOG_INFO(hashgrid.hashPoint(-0.5, -0.5, 0));
   WALBERLA_LOG_INFO(hashgrid.hashPoint(-1.5, -1.5, 0));
   WALBERLA_LOG_INFO(hashgrid.hashPoint(-1.5, 1.5, 0));
   WALBERLA_LOG_INFO(hashgrid.hashPoint(1.5, -1.5, 0));
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   
   SetBodyTypeIDs<BodyTuple>::execute();
   
   //SphereIntersectsTest();
   //PlaneIntersectsTest();
   //BoxIntersectsTest();
   //AABBIntersectsTest();
   //CapsuleIntersectsTest();
   //RaytracerTest();
   //RaytracerSpheresTest();
   
   hashgridsPlayground();

   return EXIT_SUCCESS;
}
