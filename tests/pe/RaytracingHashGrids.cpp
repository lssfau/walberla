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
#include <pe/raytracing/Raytracer.h>
#include <pe/raytracing/Color.h>
#include <pe/raytracing/ShadingFunctions.h>

#include <pe/ccd/HashGrids.h>
#include "pe/rigidbody/BodyStorage.h"
#include <core/timing/TimingTree.h>

#include <pe/utility/GetBody.h>

#include <sstream>
#include <tuple>

namespace walberla {
using namespace walberla::pe;
using namespace walberla::pe::raytracing;

typedef std::tuple<Box, Plane, Sphere, Capsule, Ellipsoid> BodyTuple ;

ShadingParameters customHashGridsBodyToShadingParams(const BodyID body) {
   if (body->getTypeID() == Sphere::getStaticTypeID()) {
      return yellowShadingParams(body);
   }
   
   return defaultBodyTypeDependentShadingParams(body);
}


void HashGridsTest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                   size_t boxes, size_t capsules, size_t spheres, size_t numberOfViews = 1,
                   real_t boxLenMin = real_t(0.1), real_t boxLenMax = real_t(0.2), bool boxRotation = false,
                   real_t capRadiusMin = real_t(0.1), real_t capRadiusMax = real_t(0.2), real_t capLenMin = real_t(0.1), real_t capLenMax = real_t(0.3),
                   real_t sphereRadiusMin = real_t(0.1), real_t sphereRadiusMax = real_t(0.3)) {
   WALBERLA_LOG_INFO("Generating " << boxes << " boxes, " << capsules << " capsules and " << spheres << " spheres");
   
   using namespace walberla::pe::ccd;
   WcTimingTree tt;
   tt.start("Setup");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   auto forest = blockforest::createBlockForest(AABB(0,0,0,4,4,4), Vector3<uint_t>(2,3,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   bool removeUnproblematic = false;
   std::vector<walberla::id_t> problematicBodyIDs = {165, 5, 31}; //{50, 44, 66, 155, 170, 51};
   std::vector<walberla::id_t> bodySIDs;
   
   // generate bodies for test
   std::vector<BodyID> bodies;
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(boxLenMin, boxLenMax); //0.2 0.5
      real_t x = math::realRandom(forestAABB.xMin(), forestAABB.xMax());
      real_t y = math::realRandom(forestAABB.yMin(), forestAABB.yMax());
      real_t z = math::realRandom(forestAABB.zMin(), forestAABB.zMax());
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x, y, z), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != nullptr);
      if (boxRotation) {
         box_->rotate(0, math::realRandom(real_t(0), real_t(1))*math::pi, math::realRandom(real_t(0), real_t(1))*math::pi);
      }
      bodies.push_back(box_);
      bodySIDs.push_back(box_->getSystemID());
   }
   
   for (size_t i = 0; i < capsules; i++) {
      real_t len = math::realRandom(capLenMin, capLenMax); // 0.2 0.5
      real_t radius = math::realRandom(capRadiusMin, capRadiusMax);
      real_t x = math::realRandom(forestAABB.xMin(), forestAABB.xMax());
      real_t y = math::realRandom(forestAABB.yMin(), forestAABB.yMax());
      real_t z = math::realRandom(forestAABB.zMin(), forestAABB.zMax());
      walberla::id_t id = walberla::id_t(boxes+i);
      CapsuleID capsule = createCapsule(*globalBodyStorage, *forest, storageID, id, Vec3(x, y, z), radius, len);
      WALBERLA_CHECK(capsule != nullptr);
      capsule->rotate(0, math::realRandom(real_t(0), real_t(1))*math::pi, math::realRandom(real_t(0), real_t(1))*math::pi);
      bodies.push_back(capsule);
      bodySIDs.push_back(capsule->getSystemID());
   }
   
   for (size_t i = 0; i < spheres; i++) {
      real_t radius = math::realRandom(sphereRadiusMin, sphereRadiusMax); // 0.2 0.3
      real_t x = math::realRandom(forestAABB.xMin(), forestAABB.xMax());
      real_t y = math::realRandom(forestAABB.yMin(), forestAABB.yMax());
      real_t z = math::realRandom(forestAABB.zMin(), forestAABB.zMax());
      walberla::id_t id = walberla::id_t(boxes+capsules+i);
      SphereID sphere = createSphere(*globalBodyStorage, *forest, storageID, id, Vec3(x, y, z), radius);
      WALBERLA_CHECK(sphere != nullptr);
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
   
   
   std::vector<std::tuple<Vec3, Vec3, Vec3>> viewVectors;
   
   // y up, in negative z direction
   viewVectors.emplace_back(Vec3(2, real_t(2.1), 7),
                                     Vec3(real_t(2.1), 2, 4),
                                     Vec3(0,1,0));
   // y up, in positive z direction
   viewVectors.emplace_back(Vec3(2, 2, -3),
                                     Vec3(2, real_t(2.1), real_t(0.1)),
                                     Vec3(0,1,0));
   // x up, in positive z direction
   viewVectors.emplace_back(Vec3(2, 2, -3),
                                     Vec3(2, real_t(2.1), real_t(0.1)),
                                     Vec3(1,0,0));
   // y and x up, in positive z direction
   viewVectors.emplace_back(Vec3(2, 2, -3),
                                     Vec3(2, real_t(2.1), real_t(0.1)),
                                     Vec3(1,1,0));
   // y and x up, in negative z direction
   viewVectors.emplace_back(Vec3(2, 2, 6.5),
                                     Vec3(real_t(2.1), real_t(2.1), 4),
                                     Vec3(real_t(0.5),1,0));
   // z up, in positive x direction
   viewVectors.emplace_back(Vec3(-3, 2, real_t(1.9)),
                                     Vec3(0, real_t(2.1), 2),
                                     Vec3(0,0,1));
   // z up, in negative x direction
   viewVectors.emplace_back(Vec3(7, 2, real_t(1.9)),
                                     Vec3(4, real_t(2.1), 2),
                                     Vec3(0,0,1));
   // z and y up, in negative x direction
   viewVectors.emplace_back(Vec3(7, 2, real_t(1.9)),
                                     Vec3(4, real_t(2.1), 2),
                                     Vec3(0,1,1));
   // z and x up, in negative y direction
   viewVectors.emplace_back(Vec3(2, 6, real_t(1.9)),
                                     Vec3(real_t(2.3), 4, 2),
                                     Vec3(1,0,1));
   // z up, in positive y direction
   viewVectors.emplace_back(Vec3(2, real_t(-3.6), real_t(1.9)),
                                     Vec3(real_t(2.3), 0, real_t(2.1)),
                                     Vec3(0,0,1));
   
   Lighting lighting0(Vec3(forestAABB.xSize()/real_t(2)+1, forestAABB.ySize()/real_t(2),
                           real_t(2)*forestAABB.zMax()+2), // 8, 5, 9.5 gut für ebenen, 0,5,8
                      Color(1, 1, 1), //diffuse
                      Color(1, 1, 1), //specular
                      Color(real_t(0.4), real_t(0.4), real_t(0.4))); //ambient
   tt.stop("Setup");

   size_t i = 0;
   for (auto& vector: viewVectors) {
      if (i == numberOfViews) {
         break;
      }
      
      Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                           size_t(640), size_t(480),
                           real_t(49.13), antiAliasFactor,
                           std::get<0>(vector),
                           std::get<1>(vector),
                           std::get<2>(vector),
                           lighting0,
                           Color(real_t(0.2),real_t(0.2),real_t(0.2)),
                           customHashGridsBodyToShadingParams);
      raytracer.setImageOutputEnabled(true);
      raytracer.setFilenameTimestepWidth(12);
      WALBERLA_LOG_INFO("output #" << i << " to: " << (boxes*100000000 + capsules*10000 + spheres) << " in " << raytracer.getImageOutputDirectory());
      raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
      raytracer.generateImage<BodyTuple>(boxes*100000000 + capsules*10000 + spheres, &tt);
      i++;
   }
   
   auto temp = tt.getReduced();
   WALBERLA_ROOT_SECTION() {
      std::cout << temp;
   }
}

ShadingParameters customArtifactsBodyToShadingParams(const BodyID body) {
   if (body->getTypeID() == Box::getStaticTypeID()) {
      return lightGreyShadingParams(body);
   }
   return defaultShadingParams(body);
}

void raytraceArtifactsForest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                             const shared_ptr<BlockStorage> forest, const BlockDataID storageID,
                             const shared_ptr<BodyStorage> globalBodyStorage,
                             const BlockDataID ccdID,
                             const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector,
                             size_t numberOfBoxes, walberla::uint8_t timestepWidth) {
   WcTimingTree tt;

   Lighting lighting(cameraPosition,
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(real_t(0.4), real_t(0.4), real_t(0.4))); //ambient
   
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(640), size_t(480),
                       real_t(49.13), antiAliasFactor,
                       cameraPosition,
                       lookAtPoint,
                       upVector,
                       lighting,
                       Color(real_t(0.2),real_t(0.2),real_t(0.2)),
                       customArtifactsBodyToShadingParams);
   raytracer.setImageOutputEnabled(true);
   raytracer.setFilenameTimestepWidth(timestepWidth);
   raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
   WALBERLA_LOG_INFO("output to: " << numberOfBoxes << " in " << raytracer.getImageOutputDirectory());
   raytracer.generateImage<BodyTuple>(numberOfBoxes, &tt);
   
   auto temp = tt.getReduced();
   WALBERLA_ROOT_SECTION() {
      std::cout << temp;
   }
}

void HashGridsArtifactsTest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                            size_t boxes, real_t boxLenMin = real_t(0.1), real_t boxLenMax = real_t(0.2)) {
   WALBERLA_LOG_INFO_ON_ROOT("HashGrids Artifacts Test - In negative Z direction");
   
   WALBERLA_LOG_INFO(" Generating " << boxes << " boxes");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   auto forest = blockforest::createBlockForest(AABB(0,0,0,4,4,4), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   // generate bodies for test
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(boxLenMin, boxLenMax); //0.2 0.5
      real_t x_min = math::realRandom(forestAABB.xMin()+len/real_t(2), forestAABB.xMax());
      real_t y_min = math::realRandom(forestAABB.yMin()+len/real_t(2), forestAABB.yMax());
      real_t z_min = math::realRandom(forestAABB.zMin()+len/real_t(2), forestAABB.zMax());
      if (i%5 == 0) {
         x_min = forestAABB.xMax() - math::realRandom(len/real_t(2), len);
      } else if (i%5 == 1){
         x_min = forestAABB.xMin() + math::realRandom(real_t(0), len/real_t(2));
      } else if (i%5 == 2){
         y_min = forestAABB.yMax() - math::realRandom(len/real_t(2), len);
      } else if (i%5 == 3){
         y_min = forestAABB.yMin() + math::realRandom(real_t(0), len/real_t(2));
      } else if (i%5 == 4){
         z_min = forestAABB.zMin() + math::realRandom(real_t(0), len/real_t(2));
      }
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x_min, y_min, z_min), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != nullptr);
   }
   
   raytraceArtifactsForest(raytracingAlgorithm, antiAliasFactor,
                           forest, storageID, globalBodyStorage, ccdID,
                           Vec3(2, 2, 7), Vec3(2, 2, 4), Vec3(0,1,0),
                           boxes, 3);
}

void HashGridsFromNegativeArtifactsTest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                                        size_t boxes, real_t boxLenMin = real_t(0.1), real_t boxLenMax = real_t(0.2)) {
   WALBERLA_LOG_INFO_ON_ROOT("HashGrids Artifacts Test - In positive Z direction");
   
   WALBERLA_LOG_INFO_ON_ROOT(" Generating " << boxes << " boxes");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   auto forest = blockforest::createBlockForest(AABB(0,0,0,4,4,4), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();

   // generate bodies for test
   std::vector<BodyID> bodies;
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(boxLenMin, boxLenMax); //0.2 0.5
      real_t x_min = math::realRandom(forestAABB.xMin()+len/real_t(2), forestAABB.xMax());
      real_t y_min = math::realRandom(forestAABB.yMin()+len/real_t(2), forestAABB.yMax());
      real_t z_min = math::realRandom(forestAABB.zMin()+len/real_t(2), forestAABB.zMax());
      
      if (i%5 == 0) {
         x_min = forestAABB.xMax() - math::realRandom(len/real_t(2), len);
      } else if (i%5 == 1){
         x_min = forestAABB.xMin() + math::realRandom(real_t(0), len/real_t(2));
      } else if (i%5 == 2){
         y_min = forestAABB.yMax() - math::realRandom(len/real_t(2), len);
      } else if (i%5 == 3){
         y_min = forestAABB.yMin() + math::realRandom(real_t(0), len/real_t(2));
      } else if (i%5 == 4){
         z_min = forestAABB.zMax() - math::realRandom(len/real_t(2), len);
      }
      
      //real_t z_min = len+0.1;
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x_min, y_min, z_min), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != nullptr);
   }
   
   raytraceArtifactsForest(raytracingAlgorithm, antiAliasFactor,
                           forest, storageID, globalBodyStorage, ccdID,
                           Vec3(2, 2, -3), Vec3(2, 2, 0), Vec3(0,1,0),
                           boxes, 4);
}

void HashGridsFromNegativeXArtifactsTest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                                         size_t boxes, real_t boxLenMin = real_t(0.1), real_t boxLenMax = real_t(0.2)) {
   WALBERLA_LOG_INFO_ON_ROOT("HashGrids Artifacts Test - In positive X direction");
   WALBERLA_LOG_INFO_ON_ROOT(" Generating " << boxes << " boxes");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   auto forest = blockforest::createBlockForest(AABB(0,0,0,4,4,4), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   // generate bodies for test
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(boxLenMin, boxLenMax); //0.2 0.5
      real_t x_min = math::realRandom(forestAABB.xMin()+len/real_t(2), forestAABB.xMax());
      real_t y_min = math::realRandom(forestAABB.yMin()+len/real_t(2), forestAABB.yMax());
      real_t z_min = math::realRandom(forestAABB.zMin()+len/real_t(2), forestAABB.zMax());
      
      if (i%5 == 0) {
         z_min = forestAABB.zMax() - math::realRandom(len/real_t(2), len);
      } else if (i%5 == 1){
         z_min = forestAABB.zMin() + math::realRandom(real_t(0), len/real_t(2));
      } else if (i%5 == 2){
         y_min = forestAABB.yMax() - math::realRandom(len/real_t(2), len);
      } else if (i%5 == 3){
         y_min = forestAABB.yMin() + math::realRandom(real_t(0), len/real_t(2));
      } else if (i%5 == 4){
         x_min = forestAABB.xMax() - math::realRandom(len/real_t(2), len);
      }
      
      //real_t z_min = len+0.1;
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x_min, y_min, z_min), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != nullptr);
   }
   
   raytraceArtifactsForest(raytracingAlgorithm, antiAliasFactor,
                           forest, storageID, globalBodyStorage, ccdID,
                           Vec3(-3, 2, 2), Vec3(0, 2, 2), Vec3(0,0,1),
                           boxes, 6);
}


Vec3 minCornerToGpos(const Vec3& minCorner, real_t lengths) {
   return minCorner + Vec3(lengths/2, lengths/2, lengths/2);
}

void HashGridsTestScene(Raytracer::Algorithm raytracingAlgorithm = Raytracer::RAYTRACE_HASHGRIDS, walberla::uint8_t antiAliasFactor = 1) {
   WALBERLA_LOG_INFO_ON_ROOT("HashGrids Test Scene");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   auto forest = blockforest::createBlockForest(AABB(0,0,0,8,8,8), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   std::vector<BodyID> bodies;
   
   // create bodies
   size_t id = 0;
   real_t len = real_t(0.6);
   
   real_t x_min = 0, y_min = 0;
   len = real_t(1.2);
   real_t gap = real_t(0.4);
   
   // cubes on z = 0 plane
   for (int i = 0; ; ++i) {
      x_min = forestAABB.xMin() + real_c(i)*(gap+len);
      if (x_min > forestAABB.max(0)) {
         break;
      }
      for (int j = 0; ; ++j) {
         y_min = forestAABB.yMin() + real_c(j)*(gap+len);
         if (y_min > forestAABB.max(1)) {
            break;
         }
         
         bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(x_min, y_min, 0), len), Vec3(len, len, len)));
      }
   }
   
   //cubes on z = max plane
   for (int i = 0; ; ++i) {
      x_min = forestAABB.xMin() + real_c(i)*(gap+len);
      if (x_min > forestAABB.max(0)) {
         break;
      }
      for (int j = 0; ; ++j) {
         y_min = forestAABB.yMin() + real_c(j)*(gap+len);
         if (y_min > forestAABB.max(1)) {
            break;
         }
         
         bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(x_min, y_min, forestAABB.zMax()-len), len), Vec3(len, len, len)));
      }
   }
   
   std::vector<std::tuple<Vec3, Vec3, Vec3>> viewVectors;
   
   // in negative x direction -> cubes to the right
   viewVectors.emplace_back(Vec3(15,4,4),
                                         Vec3(8,4,4),
                                         Vec3(0,1,0));
   // in negative x direction and negative z direction, up vector in y direction -> cubes from the right tilted
   viewVectors.emplace_back(Vec3(12,4,8),
                                         Vec3(6,4,2),
                                         Vec3(0,1,0));
   // in negative x direction and negative z direction, up vector in negative y direction
   viewVectors.emplace_back(Vec3(12,4,8),
                                         Vec3(6,4,2),
                                         Vec3(0,-1,0));
   // in positive x direction
   viewVectors.emplace_back(Vec3(-7,4,4),
                                         Vec3(0,4,4),
                                         Vec3(0,1,0));
   // in negative x direction
   viewVectors.emplace_back(Vec3(4,4,15),
                                         Vec3(4,4,8),
                                         Vec3(0,1,0));
   
   WcTimingTree tt;
   
   Lighting lighting(Vec3(1,2,15),
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(real_t(0.4), real_t(0.4), real_t(0.4))); //ambient
   
   int i = 0;
   for (auto& vector: viewVectors) {
      Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                          size_t(640), size_t(480),
                          real_t(49.13), antiAliasFactor,
                          std::get<0>(vector),
                          std::get<1>(vector),
                          std::get<2>(vector),
                          lighting,
                          Color(real_t(0.2),real_t(0.2),real_t(0.2)));
      
      raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
      raytracer.setImageOutputEnabled(true);
      raytracer.setFilenameTimestepWidth(1);
      WALBERLA_LOG_INFO("output to: " << i << " in " << raytracer.getImageOutputDirectory());
      raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
      raytracer.generateImage<BodyTuple>(size_t(i), &tt);
      
      auto temp = tt.getReduced();
      WALBERLA_ROOT_SECTION() {
         std::cout << temp;
      }
      
      i++;
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   SetBodyTypeIDs<BodyTuple>::execute();
   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   const Raytracer::Algorithm algorithm = Raytracer::RAYTRACE_COMPARE_BOTH_STRICTLY;
   const walberla::uint8_t antiAliasFactor = 1;
   HashGridsTestScene(algorithm, antiAliasFactor);
   
   if (argc >= 2 && strcmp(argv[1], "--longrun") == 0) {
      HashGridsTest(algorithm, antiAliasFactor,
                    50, 30, 130,
                    10);
      HashGridsTest(algorithm, antiAliasFactor,
                    60, 60, 3,
                    1,
                    real_t(0.1), real_t(0.3), true,
                    real_t(0.1), real_t(0.2), real_t(0.1), real_t(0.2),
                    real_t(0.5), real_t(0.6));
      HashGridsArtifactsTest(algorithm, antiAliasFactor, 750, real_t(0.2), real_t(0.3));
      HashGridsFromNegativeArtifactsTest(algorithm, antiAliasFactor, 750, real_t(0.2), real_t(0.3));
      HashGridsFromNegativeXArtifactsTest(algorithm, antiAliasFactor, 750, real_t(0.2), real_t(0.3));
   }
   
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
