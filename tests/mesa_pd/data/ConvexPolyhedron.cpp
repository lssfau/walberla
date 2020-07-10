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
//! \file MeshMesapdConvexPolyhedronTest.cpp
//! \author Lukas Werner
//
//======================================================================================================================

#include "mesa_pd/data/shape/ConvexPolyhedron.h"

#include "core/debug/TestSubsystem.h"
#include "core/mpi/all.h"

#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/mpi/ShapePackUnpack.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/PolyMeshes.h"

namespace walberla {

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   mpi::Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   /// Load Mesh File
   WALBERLA_LOG_INFO("LOAD MESH");
   // Ideally this should be done on only one process and
   // communicated to other processes, if or when required.
   mesh::TriangleMesh bunnyMesh;
   std::string meshInputFile = "../mesh/bunny.obj";
   mesh::readFromFile<mesh::TriangleMesh>(meshInputFile, bunnyMesh);
   WALBERLA_LOG_INFO("Read mesh file: " << meshInputFile << " (" << bunnyMesh.n_vertices() << " vertices, "
      << bunnyMesh.n_faces() << " faces, volume = " << mesh::computeVolume(bunnyMesh) << ")" );
   mesh::translate(bunnyMesh, -mesh::toWalberla(mesh::computeCentroid(bunnyMesh)));

   /// MESAPD Data
   WALBERLA_LOG_INFO("CREATE CONVEX POLYHEDRON");
   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(3);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   auto ac = walberla::make_shared<mesa_pd::data::ParticleAccessorWithShape>(ps, ss);
   auto bunnyShapeID = ss->create<mesa_pd::data::ConvexPolyhedron>(bunnyMesh);

   auto bunnyParticle = ps->create();
   bunnyParticle->setShapeID(bunnyShapeID);
   auto mesh1ParticleUid = bunnyParticle->getUid();
   auto mesh1ParticleIdx = ac->uidToIdx(mesh1ParticleUid);

   auto mesh1ParticleShape = dynamic_cast<mesa_pd::data::ConvexPolyhedron*>(ac->getShape(mesh1ParticleIdx));

   WALBERLA_CHECK_FLOAT_EQUAL(ac->getShape(mesh1ParticleIdx)->getVolume(), mesh::computeVolume(bunnyMesh));
   WALBERLA_CHECK_EQUAL(mesh1ParticleShape->getMesh().n_vertices(), bunnyMesh.n_vertices());
   WALBERLA_CHECK_EQUAL(mesh1ParticleShape->getMesh().n_faces(), bunnyMesh.n_faces());

   /// Check: Pack - Unpack
   WALBERLA_LOG_INFO("PACK - UNPACK");

   shared_ptr<mesa_pd::data::BaseShape> bs0 = make_shared<mesa_pd::data::ConvexPolyhedron>(bunnyMesh);
   mesh1ParticleShape->updateMassAndInertia(real_t(1));
   std::shared_ptr<mesa_pd::data::BaseShape> bs1 = nullptr;

   WALBERLA_LOG_INFO("Packing mesh shape");
   mpi::SendBuffer sb;
   sb << bs0;
   WALBERLA_LOG_INFO("Unpacking mesh shape");
   mpi::RecvBuffer rb(sb);
   rb >> bs1;

   //WALBERLA_CHECK_EQUAL(bs0->getShapeType(), mesa_pd::data::ConvexPolyhedron::SHAPE_TYPE);
   //WALBERLA_CHECK_EQUAL(bs1->getShapeType(), mesa_pd::data::ConvexPolyhedron::SHAPE_TYPE);
   WALBERLA_CHECK_IDENTICAL(bs0->getMass(), bs1->getMass());
   WALBERLA_CHECK_IDENTICAL(bs0->getInvMass(), bs1->getInvMass());
   WALBERLA_CHECK_IDENTICAL(bs0->getInertiaBF(), bs1->getInertiaBF());
   WALBERLA_CHECK_IDENTICAL(bs0->getInvInertiaBF(), bs1->getInvInertiaBF());

   auto bm0 = dynamic_cast<mesa_pd::data::ConvexPolyhedron*>(bs0.get());
   auto bm1 = dynamic_cast<mesa_pd::data::ConvexPolyhedron*>(bs1.get());
   WALBERLA_CHECK_EQUAL(bm0->getMesh().n_vertices(), bm1->getMesh().n_vertices());
   WALBERLA_CHECK_EQUAL(bm0->getMesh().n_faces(), bm1->getMesh().n_faces());

   /// Volume and Inertia
   WALBERLA_LOG_INFO("VOLUME AND INERTIA");

   mesh::TriangleMesh cubeMesh;
   std::string cubeMeshInputFile = "../mesh/cube.obj";
   mesh::readFromFile<mesh::TriangleMesh>(cubeMeshInputFile, cubeMesh);
   mesh::translate(cubeMesh, -mesh::toWalberla(mesh::computeCentroid(cubeMesh)));

   real_t cubeVolume = mesh::computeVolume(cubeMesh);
   real_t cubeSideLen = std::cbrt(cubeVolume);
   auto cubeDensity = real_t(123);
   mesa_pd::data::Box boxCubeShape((Vector3<real_t>(cubeSideLen)));
   boxCubeShape.updateMassAndInertia(cubeDensity);

   mesa_pd::data::ConvexPolyhedron meshCubeShape(cubeMesh);
   meshCubeShape.updateMassAndInertia(cubeDensity);

   // Test mass properties of MESAPD convexpolyhedron body

   WALBERLA_CHECK_FLOAT_EQUAL(boxCubeShape.getInertiaBF(), meshCubeShape.getInertiaBF());
   WALBERLA_CHECK_FLOAT_EQUAL(boxCubeShape.getVolume(), meshCubeShape.getVolume());
   WALBERLA_CHECK_FLOAT_EQUAL(boxCubeShape.getMass(), meshCubeShape.getMass());
   WALBERLA_CHECK_FLOAT_EQUAL(boxCubeShape.getInvMass(), meshCubeShape.getInvMass());

   // Test new mesh ops function to calculate all mass properties at once

   Matrix3<real_t> massPropInertia;
   Vector3<real_t> massPropCentroid;
   real_t massPropMass;

   mesh::computeMassProperties(cubeMesh, cubeDensity, massPropCentroid, massPropInertia, massPropMass);

   WALBERLA_CHECK_FLOAT_EQUAL(boxCubeShape.getInertiaBF(), massPropInertia);
   WALBERLA_CHECK_FLOAT_EQUAL(boxCubeShape.getMass(), massPropMass);
   WALBERLA_CHECK_FLOAT_EQUAL(massPropCentroid, Vector3<real_t>(real_t(0)))

   // Test new mesh ops inertia calculation against old one

   auto oldInertia = mesh::computeInertiaTensor(cubeMesh)*cubeDensity;
   WALBERLA_CHECK_FLOAT_EQUAL(massPropInertia, oldInertia);

   /// Interaction radius
   // for a cube: space diagonal length = side lengths * sqrt(3) -> bounding sphere radius = space diagonal / 2
   WALBERLA_CHECK_FLOAT_EQUAL(meshCubeShape.getBoundingSphereRadius(), real_t(1) * std::sqrt(real_t(3)) * real_t(0.5));

   auto bunnyShape = dynamic_cast<mesa_pd::data::ConvexPolyhedron*>(ss->shapes[bunnyShapeID].get());
   real_t bunnyRadius = bunnyShape->getBoundingSphereRadius();
   real_t maxSqRadius(0);
   for(auto vh : bunnyMesh.vertices()) {
      auto v = mesh::toWalberla(bunnyMesh.point(vh));
      auto centroidToVSqr = v.sqrLength();

      if (centroidToVSqr > maxSqRadius) {
         maxSqRadius = centroidToVSqr;
      }
   }
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(bunnyRadius, std::sqrt(maxSqRadius), real_t(1e-4));

   /// Support
   WALBERLA_LOG_INFO("SUPPORT");

   for (int x = -1; x <= 1; x+=2) {
      for (int y = -1; y <= 1; y+=2) {
         for (int z = -1; z <= 1; z+=2) {
            Vector3<real_t> d((real_t(x)), real_t(y), real_t(z));
            WALBERLA_CHECK_FLOAT_EQUAL(boxCubeShape.support(d), meshCubeShape.support(d));
         }
      }
   }

   AABB supportTestAABB(mesa_pd::Vec3(real_t(-1)), mesa_pd::Vec3(real_t(1)));
   std::mt19937 rng;
   for(uint_t i = 0; i < 500; ++i) {
      auto rndPoint = supportTestAABB.randomPoint(rng);
      WALBERLA_CHECK_FLOAT_EQUAL(boxCubeShape.support(rndPoint), meshCubeShape.support(rndPoint));
   }

   return EXIT_SUCCESS;
}
}

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}