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

void RaytracerSpheresTestScene(Raytracer::Algorithm raytracingAlgorithm = Raytracer::RAYTRACE_HASHGRIDS, walberla::uint8_t antiAliasFactor = 1) {
   WALBERLA_LOG_INFO("Raytracer Spheres Scene");
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   auto forest = blockforest::createBlockForest(AABB(0,0,0,10,10,10), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   
   Lighting lighting(Vec3(0, 5, 8), // 8, 5, 9.5 gut für ebenen, 0,5,8
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(real_t(0.4), real_t(0.4), real_t(0.4))); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(640), size_t(480),
                       real_t(49.13), antiAliasFactor,
                       Vec3(-5,5,5), Vec3(-1,5,5), Vec3(0,0,1), //-5,5,5; -1,5,5
                       lighting,
                       Color(real_t(0.2),real_t(0.2),real_t(0.2)),
                       customSpheresBodyToShadingParams);
   
   MaterialID iron = Material::find("iron");
   
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
   
   raytracer.setImageOutputEnabled(true);
   
   raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
   raytracer.generateImage<BodyTuple>(1);
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   SetBodyTypeIDs<BodyTuple>::execute();
   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   const Raytracer::Algorithm algorithm = Raytracer::RAYTRACE_COMPARE_BOTH_STRICTLY;
   const walberla::uint8_t antiAliasFactor = 1;
   RaytracerSpheresTestScene(algorithm, antiAliasFactor);
   
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
