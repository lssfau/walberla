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

void HashGridsTest() {
   using namespace walberla::pe::ccd;

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(-4,-4,-4,4,4,4), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   //HashGrids::HashGrid hashgrid(2.0); // initialize a 4x4x4 hashgrid with cellspan 2.0 (can hold a block with 8x8x8 in size.)
   
   //auto hccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "HCCD");
   
   std::vector<IBlock*> blocks;
   forest->getBlocks(blocks);
   WALBERLA_LOG_INFO(blocks.size());
   
   Storage* storage = blocks[0]->getData<Storage>(storageID);
   HashGrids hashGrids(*globalBodyStorage, (*storage)[0], (*storage)[1]);
   
   const AABB& blockAABB = forest->getDomain();
   real_t inf = std::numeric_limits<real_t>::max();
   
   size_t boxCount = 10000;
   size_t capsuleCount = 5000;
   
   // generate bodies for test
   std::vector<BodyID> bodies;
   for (size_t i = 0; i < boxCount; i++) {
      real_t x = math::realRandom(blockAABB.xMin()+0.5, blockAABB.xMax());
      real_t y = math::realRandom(blockAABB.yMin()+0.5, blockAABB.yMax());
      real_t z = math::realRandom(blockAABB.zMin()+0.5, blockAABB.zMax());
      real_t len = math::realRandom(real_t(0.2), real_t(0.5));
      BoxID box = createBox(*globalBodyStorage, *forest, storageID, walberla::id_t(i), Vec3(x, y, z), Vec3(len, len, len));
      //hashgrid.add(box);
      hashGrids.add(box);
      bodies.push_back(box);
   }
   for (size_t i = 0; i < capsuleCount; i++) {
      real_t x = math::realRandom(blockAABB.xMin()+0.5, blockAABB.xMax());
      real_t y = math::realRandom(blockAABB.yMin()+0.5, blockAABB.yMax());
      real_t z = math::realRandom(blockAABB.zMin()+0.5, blockAABB.zMax());
      real_t len = math::realRandom(real_t(0.2), real_t(0.5));
      real_t radius = real_t(0.1);
      CapsuleID capsule = createCapsule(*globalBodyStorage, *forest, storageID, walberla::id_t(boxCount+i), Vec3(x, y, z), radius, len);
      //hashgrid.add(capsule);
      hashGrids.add(capsule);
      bodies.push_back(capsule);
   }
   hashGrids.generatePossibleContacts();
   
   // rays used for testing
   std::vector<Ray> rays = {
      Ray(Vec3(-5, 0.2, 0.2), Vec3(1, 0, 0).getNormalized()),
      Ray(Vec3(5, 0.2, 0.2), Vec3(-1, 0, 0).getNormalized()),
      Ray(Vec3(5, -5, -5), Vec3(-20, -7, 2).getNormalized()),
      Ray(Vec3(-4.5, 3.5, 0.5), Vec3(1, -1, 0).getNormalized()),
      Ray(Vec3(4.5, 5, 0.5), Vec3(-1, -1, 1).getNormalized()),
      Ray(Vec3(0, 2, 6), Vec3(0, -0.1, -1).getNormalized()),
      Ray(Vec3(-5, -5, -5), Vec3(1, 1, 1).getNormalized())
   };
   
   real_t t;
   Vec3 n;
   
   BodyID closestBody_naive = NULL;
   real_t t_closest_naive = inf;
   
   BodyID closestBody_allpossible = NULL;
   real_t t_closest_allpossible = inf;
   
   BodyID closestBody_first = NULL;
   real_t t_closest_first = inf;
   
   BodyID closestBody_hashgrids_first = NULL;
   real_t t_closest_hashgrids_first = inf;
   
   int i = 0;
   for (const Ray& ray: rays) {
      WALBERLA_LOG_INFO("RAY " << i << ": " << ray.getOrigin() << ", " << ray.getDirection() << " ----");

      // reset variables for new ray
      closestBody_naive = NULL;
      t_closest_naive = inf;
      closestBody_allpossible = NULL;
      t_closest_allpossible = inf;
      closestBody_first = NULL;
      t_closest_first = inf;
      
      // -- naive method
      IntersectsFunctor naiveIntersectsFunc(ray, t, n);
      for (const BodyID& body: bodies) {
         bool intersects = SingleCast<BodyTuple, IntersectsFunctor, bool>::execute(body, naiveIntersectsFunc);
         if (intersects && t < t_closest_naive && t > 0) {
            closestBody_naive = body;
            t_closest_naive = t;
         }
      }
      WALBERLA_LOG_INFO(" Naive closest: \t\t\t" << (closestBody_naive != NULL ? std::to_string(closestBody_naive->getID()) : "none")
                        << " and t = " << t_closest_naive);
      if (closestBody_naive != NULL) {
         WALBERLA_LOG_INFO("  gpos: " << closestBody_naive->getPosition() << ", aabb: " << closestBody_naive->getAABB());
      }
      /*
      // -- hashgrid and collecting all possible bodies
      std::vector<BodyID> bodiesContainer;
      hashgrid.possibleRayIntersectingBodies(ray, forest->getDomain(), bodiesContainer);
      
      IntersectsFunctor intersectsFunc(ray, t, n);
      for (BodyID const& body : bodiesContainer) {
         bool intersects = SingleCast<BodyTuple, IntersectsFunctor, bool>::execute(body, intersectsFunc);
         if (intersects && t < t_closest_allpossible) {
            closestBody_allpossible = body;
            t_closest_allpossible = t;
         }
      }
      WALBERLA_LOG_INFO("All possible closest: \t" << (closestBody_allpossible != NULL ? std::to_string(closestBody_allpossible->getID()) : "none")
                        << " and t = " << t_closest_allpossible);
      
      // -- hashgrid and only until closest body found
      closestBody_first = hashgrid.getRayIntersectingBody<BodyTuple>(ray, forest->getDomain(), t_closest_first, n);
      WALBERLA_LOG_INFO("First closest: \t\t\t" << (closestBody_first != NULL ? std::to_string(closestBody_first->getID()) : "none")
                        << " and t = " << t_closest_first);
      
      //WALBERLA_CHECK(closestBody_naive == closestBody_allpossible && closestBody_naive == closestBody_first,
      //               "Different intersection methods dont match");
      */
      
      // -- using hashgrids and only until closest body found
      
      WALBERLA_LOG_INFO(" Hashgrids");
      closestBody_hashgrids_first = hashGrids.getClosestBodyIntersectingWithRay<BodyTuple>(ray, blockAABB, t_closest_hashgrids_first, n);
      
      WALBERLA_LOG_INFO(" Hashgrids closest: \t\t" << (closestBody_hashgrids_first != NULL ? std::to_string(closestBody_hashgrids_first->getID()) : "none")
                        << " and t = " << t_closest_hashgrids_first);
      if (closestBody_hashgrids_first != NULL) {
         WALBERLA_LOG_INFO("  gpos: " << closestBody_hashgrids_first->getPosition() << ", aabb: " << closestBody_hashgrids_first->getAABB());
      }
      
      if (closestBody_naive == closestBody_hashgrids_first) {
         WALBERLA_LOG_INFO("-> SAME")
      } else {
         WALBERLA_LOG_INFO("-> FAIL");
      }
      
      WALBERLA_CHECK(closestBody_naive == closestBody_hashgrids_first);

      WALBERLA_LOG_INFO("RAY " << i << " end. ----");
      i++;
   }

   WcTimingTree tt;
   tt.start("Hashgrids Loop");
   Ray ray(Vec3(0, 0, -7), Vec3(1, 0, 0));

   int rays_x = 1000;
   int rays_y = 1000;

   for (int x = 0; x < rays_x; ++x) {
      for (int y = 0; y < rays_y; ++y) {
         Vec3 dir = (Vec3(-4+8*real_t(x)/real_t(rays_x), -4+8*real_t(y)/real_t(rays_y), 4) - ray.getOrigin()).getNormalized();
         ray.setDirection(dir);
         closestBody_hashgrids_first = hashGrids.getClosestBodyIntersectingWithRay<BodyTuple>(ray, blockAABB, t_closest_hashgrids_first, n);
      }
   }
   tt.stop("Hashgrids Loop");

   /*tt.start("Naive Loop");
   IntersectsFunctor naiveIntersectsFunc(ray, t, n);
   for (int x = 0; x < rays_x; ++x) {
      for (int y = 0; y < rays_y; ++y) {
         Vec3 dir = (Vec3(-4+8*real_t(x)/real_t(rays_x), -4+8*real_t(y)/real_t(rays_y), 4) - ray.getOrigin()).getNormalized();
         ray.setDirection(dir);
         for (const BodyID& body: bodies) {
            SingleCast<BodyTuple, IntersectsFunctor, bool>::execute(body, naiveIntersectsFunc);
         }
      }
   }
   tt.stop("Naive Loop");*/
   
   auto temp = tt.getReduced( );
   WALBERLA_ROOT_SECTION()
   {
      std::cout << temp;
   }
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
   
   HashGridsTest();
   
   //hashgridsPlayground();

   return EXIT_SUCCESS;
}
