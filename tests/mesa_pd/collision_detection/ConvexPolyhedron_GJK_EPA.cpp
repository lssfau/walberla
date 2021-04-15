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
//! \file
//! \author Lukas Werner
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"

#include <core/Environment.h>
#include <core/all.h>
#include <core/logging/Logging.h>
#include <mesa_pd/collision_detection/GeneralContactDetection.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/Sphere.h>

#include "mesa_pd/data/shape/ConvexPolyhedron.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/TriangleMeshes.h"

namespace walberla {
namespace mesa_pd {

using namespace walberla::mesa_pd::collision_detection;
using namespace walberla::mesa_pd::collision_detection::analytic;
using namespace walberla::mesa_pd::data;

bool gjkEPAcollideHybrid(Support &geom1, Support &geom2, Vec3& normal, Vec3& contactPoint, real_t& penetrationDepth) {
   // This function is stolen from the general mesapd GJK-EPA test.

   // For more information on hybrid GJK/EPA see page 166 in "Collision Detection in Interactive 3D
   // Environments" by Gino van den Bergen.

   //1. Run GJK with considerably enlarged objects.
   real_t margin = real_t(1e-6);
   GJK gjk;
   if(gjk.doGJKmargin(geom1, geom2, margin)){
      //2. If collision is possible perform EPA.
      //std::cerr << "Performing EPA.";
      EPA epa;
      epa.useSphereOptimization( true );
      return epa.doEPAmargin(geom1, geom2, gjk, normal, contactPoint, penetrationDepth, margin);
   } else {
      return false;
   }
}

/**
 * \brief Compare the collision results of a normal box with the results of a mesh resembling this box.
 * \param testSupport The particle with which the box and the mesh are colliding.
 * \param convPoly Convex polyhedron, contains a mesh resembling the box.
 * \param box Box, reference for convex polyhedron.
 */
void runComparisonTest(Support& testSupport, Support& convPoly, Support& box) {
   Vec3 polyNormal;
   Vec3 polyContactPoint;
   real_t polyPenetrationDepth;
   bool polyCollided = gjkEPAcollideHybrid(convPoly, testSupport, polyNormal, polyContactPoint, polyPenetrationDepth);
   //WALBERLA_LOG_INFO("conv poly: " << polyCollided << ", " << polyNormal << ", " << polyContactPoint << ", " << polyPenetrationDepth);

   Vec3 boxNormal;
   Vec3 boxContactPoint;
   real_t boxPenetrationDepth;
   bool boxCollided = gjkEPAcollideHybrid(box, testSupport, boxNormal, boxContactPoint, boxPenetrationDepth);
   //WALBERLA_LOG_INFO("box:       " << boxCollided << ", " << boxNormal << ", " << boxContactPoint << ", " << boxPenetrationDepth);

   WALBERLA_CHECK(polyCollided == boxCollided);
   if (polyCollided) {
      auto eps = real_t(1e-3);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(polyNormal, boxNormal, real_t(1e-2));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(polyContactPoint, boxContactPoint, eps);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(polyPenetrationDepth, boxPenetrationDepth, eps);
   }
}

/** Test the GJK-EPA implementation for collision detection with convex polyhedrons vs. multiple other shapes.
 */
int main(int argc, char** argv){
   walberla::debug::enterTestMode();
   walberla::mpi::Environment env(argc, argv);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   if (std::is_same<walberla::real_t, float>::value) {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   WALBERLA_LOG_INFO("ConvexPolyhedron <-> Sphere vs. Box <-> Sphere");
   mesh::TriangleMesh cubeMesh;
   std::string cubeInputFile = "../mesh/cube.obj";
   mesh::readFromFile< mesh::TriangleMesh >(cubeInputFile, cubeMesh);
   mesh::translate(cubeMesh, -mesh::toWalberla(mesh::computeCentroid(cubeMesh)));

   ConvexPolyhedron convPoly_(cubeMesh);
   Support convPoly(Vec3(0), Rot3(), convPoly_);

   Box box_(Vec3(1));
   Support box(Vec3(0), Rot3(), box_);

   Sphere sphere_(real_t(0.5));
   Support sphere(Vec3(real_t(0.5),real_t(0.5),real_t(0.9)), Rot3(), sphere_);
   runComparisonTest(sphere, convPoly, box);

   /// Fuzzy Testing

   AABB collisionTestAABB(mesa_pd::Vec3(real_t(-1)), mesa_pd::Vec3(real_t(1)));
   std::mt19937 rng;

   WALBERLA_LOG_INFO("Fuzzy ConvexPolyhedron <-> Sphere vs. Box <-> Sphere");
   for(uint_t i = 0; i < 500; ++i) {
      auto rndPoint = collisionTestAABB.randomPoint(rng);

      Sphere fuzzSphere_(real_t(0.5));
      Support fuzzSphere(rndPoint, Rot3(), sphere_);

      runComparisonTest(fuzzSphere, convPoly, box);
   }

   WALBERLA_LOG_INFO("Fuzzy ConvexPolyhedron <-> Ellipsoid vs. Box <-> Ellipsoid");
   for(uint_t i = 0; i < 500; ++i) {
      auto rndPoint = collisionTestAABB.randomPoint(rng);

      Ellipsoid el_(rndPoint);
      Support el(Vec3(real_t(1)), Rot3(), el_);

      el.rot_.rotate(Vec3(math::realRandom(), math::realRandom(), math::realRandom()));

      runComparisonTest(el, convPoly, box);
   }

   WALBERLA_LOG_INFO("Fuzzy ConvexPolyhedron <-> Turned Box vs. Box <-> Turned Box");
   for(uint_t i = 0; i < 500; ++i) {
      auto rndPoint = collisionTestAABB.randomPoint(rng);

      Box fuzzBox_(Vec3(real_t(1),real_t(0.5),real_t(0.5)));
      Support fuzzBox(rndPoint, Rot3(), fuzzBox_);

      fuzzBox.rot_.rotate(Vec3(math::realRandom(), math::realRandom(), math::realRandom()));

      runComparisonTest(fuzzBox, convPoly, box);
   }

   WALBERLA_LOG_INFO("Fuzzy Turned ConvexPolyhedron <-> Box vs. Turned Box <-> Box");
   for(uint_t i = 0; i < 500; ++i) {
      auto rndPoint = collisionTestAABB.randomPoint(rng);

      Box fuzzBox_(Vec3(real_t(1),real_t(0.5),real_t(0.5)));
      Support fuzzBox(rndPoint, Rot3(), fuzzBox_);

      // rotate all the things
      fuzzBox.rot_.rotate(Vec3(math::realRandom(), math::realRandom(), math::realRandom()));
      Vec3 rotPhi(math::realRandom(), math::realRandom(), math::realRandom());
      convPoly.rot_.rotate(rotPhi);
      box.rot_.rotate(rotPhi);

      runComparisonTest(fuzzBox, convPoly, box);
   }

   return EXIT_SUCCESS;
}

}
}

int main( int argc, char* argv[] ) {
   return walberla::mesa_pd::main( argc, argv );
}