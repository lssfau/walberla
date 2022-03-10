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

#include "blockforest/all.h"

#include "core/debug/TestSubsystem.h"
#include "core/math/Rot3.h"
#include "core/mpi/all.h"

#include "vtk/all.h"

#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/mpi/ShapePackUnpack.h>
#include <mesa_pd/vtk/ConvexPolyhedron/MeshParticleVTKOutput.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/SurfaceVelocityVertexDataSource.h>
#include <mesh_common/PolyMeshes.h>

#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"

namespace walberla {

using namespace walberla::mesa_pd;

int main(int argc, char** argv) {
   walberla::debug::enterTestMode();
   mpi::Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   std::string vtk_out = "vtk_out";

   /// BlockForest
   Vector3< real_t > minCorner(real_t(0), real_t(0), real_t(0));
   Vector3< real_t > maxCorner(real_t(4), real_t(4), real_t(4));
   math::AABB domainSize(minCorner, maxCorner);
   Vector3< bool > periodicity(true, true, true);
   Vector3< uint_t > numBlocks(2, 2, 1);
   shared_ptr< BlockForest > forest = blockforest::createBlockForest(domainSize, numBlocks, periodicity);
   auto domain = std::make_shared< mesa_pd::domain::BlockForestDomain >(forest);

   WALBERLA_LOG_INFO(domain->getNumLocalAABBs() << ": " << domain->getUnionOfLocalAABBs());

   /// MESAPD Data
   WALBERLA_LOG_INFO_ON_ROOT("CREATE MESAPD DATA STRUCTURES");
   auto ps = std::make_shared< data::ParticleStorage >(4);
   auto ss = std::make_shared< data::ShapeStorage >();
   auto ac = walberla::make_shared< data::ParticleAccessorWithShape >(ps, ss);

   /// Load Mesh File
   WALBERLA_LOG_INFO_ON_ROOT("LOAD MESH");
   // Ideally this should be done on only one process and
   // communicated to other processes, if or when required.
   mesh::TriangleMesh cubeMesh;
   std::string cubeMeshInputFile = "../mesh/cube.obj";
   mesh::readAndBroadcast< mesh::TriangleMesh >(cubeMeshInputFile, cubeMesh);
   WALBERLA_LOG_INFO_ON_ROOT("Read mesh file: " << cubeMeshInputFile << " (" << cubeMesh.n_vertices() << " vertices, "
                                                << cubeMesh.n_faces()
                                                << " faces, volume = " << mesh::computeVolume(cubeMesh) << ")");
   mesh::translate(cubeMesh, mesh::toWalberla(-mesh::computeCentroid(cubeMesh)));

   /// Mesh Shape
   auto cubeShapeID = ss->create< data::ConvexPolyhedron >(cubeMesh);

   /// Mesh Particles
   math::Rot3< real_t > cubeRotation(Vector3< real_t >(real_t(M_PI) / real_t(4), real_t(0), real_t(0)));
   for (uint_t x = 0; x <= 1; ++x) {
      for (uint_t y = 0; y <= 1; ++y) {
         for (uint_t z = 0; z <= 1; ++z) {
            // position of new mesh cube slightly shifted inwards
            Vector3< real_t > pos((real_t(0.05) + real_t(x) * maxCorner[0] - real_t(x) * real_t(0.1)),
                                  real_t(0.05) + real_t(y) * maxCorner[1] - real_t(y) * real_t(0.1),
                                  real_t(0.05) + real_t(z) * maxCorner[2] - real_t(z) * real_t(0.1));

            if (domain->isContainedInProcessSubdomain(uint_c(mpi::MPIManager::instance()->rank()), pos)) {
               auto rotatedCubeParticle = ps->create();
               rotatedCubeParticle->setShapeID(cubeShapeID);
               rotatedCubeParticle->setPosition(pos);
               rotatedCubeParticle->setRotation(cubeRotation);
               rotatedCubeParticle->setOwner(mpi::MPIManager::instance()->rank());
               rotatedCubeParticle->setLinearVelocity(pos);
               WALBERLA_LOG_INFO("Created cube particle");
            }
         }
      }
   }

   /// VTK Output

   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition(forest, "domain_decomposition", 1, vtk_out, "simulation_step");
   vtkDomainOutput->write();

   mesa_pd::MeshParticleVTKOutput< mesh::PolyMesh > meshParticleVTK(ps, "mesh", uint_t(1));
   meshParticleVTK.addFaceOutput< data::SelectParticleUid >("UID");
   meshParticleVTK.addVertexOutput< data::SelectParticleInteractionRadius >("InteractionRadius");
   meshParticleVTK.addFaceOutput< data::SelectParticleLinearVelocity >("LinearVelocity");
   meshParticleVTK.addVertexOutput< data::SelectParticlePosition >("Position");
   auto surfaceVelDataSource = make_shared<
      mesa_pd::SurfaceVelocityVertexDataSource< mesh::PolyMesh, mesa_pd::data::ParticleAccessorWithShape > >(
      "SurfaceVelocity", *ac);
   meshParticleVTK.addVertexDataSource(surfaceVelDataSource);
   meshParticleVTK(*ac);

   return EXIT_SUCCESS;
}
}

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}