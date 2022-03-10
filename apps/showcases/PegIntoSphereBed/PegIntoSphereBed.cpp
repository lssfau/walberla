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
//! \file   PegIntoSphereBed.cpp
//! \author Lukas Werner
//
//======================================================================================================================

/* This show case demonstrates a peg (consisting of a convex mesh) being rammed into
 * a bed of spheres. It is purely based on MESA-PD and OpenMesh.
 * The properties of the spheres and the peg, and additionally the height of the bed
 * can be controlled by adjusting the appropriate parameters in the config file.
 * The approximation of the peg's round body is controlled by the number of side
 * edges. Its initial position is given with respect to its tip, the velocity with
 * which it its afterwards moved is configurable as well.
 *
 *          peg radius
 *           <------>
 *            ______
 *           |      |            ^
 *           |      |            | peg
 *           |      |            | body height
 * oooooooooo|      |oooooooooo  v                   ^
 * ooooooooooo\    /ooooooooooo  ^                   |
 * oooooooooooo\  /oooooooooooo  | peg pike height   | sphere bed height
 * ooooooooooooo\/ooooooooooooo  v                   |
 * oooooooooooooooooooooooooooo                      |
 * oooooooooooooooooooooooooooo                      v
 * */

#include "blockforest/Initialization.h"
#include "blockforest/StructuredBlockForest.h"

#include "mesa_pd/vtk/ConvexPolyhedron/data_sources/SurfaceVelocityVertexDataSource.h"

#include "vtk/all.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <mesa_pd/collision_detection/GeneralContactDetection.h>
#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/SpringDashpot.h>
#include <mesa_pd/kernel/SpringDashpotSpring.h>
#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/ReduceContactHistory.h>
#include <mesa_pd/mpi/ReduceProperty.h>
#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>
#include <mesa_pd/vtk/ParticleVtkOutput.h>
#include <mesh_common/MeshOperations.h>
#include <mesh_common/vtk/DistributedVTKMeshWriter.h>

#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/shape/ConvexPolyhedron.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/vtk/ConvexPolyhedron/MeshParticleVTKOutput.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/PolyMeshes.h"
#include "mesh_common/TriangleMeshes.h"

namespace walberla {

void createPeg(mesh::TriangleMesh & mesh, mesa_pd::Vec3 & pegPikeTipPosition, real_t bodyHeight, real_t pikeHeight, real_t radius, uint_t numSideEdges) {
   real_t alpha = real_t(2) * math::pi / real_t(numSideEdges); // 360° / numSideEdges -> approximation of cylinder and cone
   real_t topCornerZ = pikeHeight + bodyHeight;
   real_t bottomCornerZ = pikeHeight;

   mesh::TriangleMesh::Point topCenterPoint(radius, radius, topCornerZ);
   mesh::TriangleMesh::VertexHandle topCenterVertex = mesh.add_vertex(topCenterPoint);

   mesh::TriangleMesh::Point bottomCenterPoint(radius, radius, real_t(0));
   mesh::TriangleMesh::VertexHandle bottomCenterVertex = mesh.add_vertex(bottomCenterPoint);

   mesh::TriangleMesh::VertexHandle firstTopVertex;
   mesh::TriangleMesh::VertexHandle firstBottomVertex;

   mesh::TriangleMesh::VertexHandle lastTopVertex;
   mesh::TriangleMesh::VertexHandle lastBottomVertex;

   for (uint_t e = 0; e < numSideEdges; ++e) {
      real_t x_corner = radius + radius * std::sin(real_t(e) * alpha);
      real_t y_corner = radius + radius * std::cos(real_t(e) * alpha);

      mesh::TriangleMesh::Point topPoint(real_t(x_corner), real_t(y_corner), topCornerZ);
      mesh::TriangleMesh::Point bottomPoint(real_t(x_corner), real_t(y_corner), bottomCornerZ);

      mesh::TriangleMesh::VertexHandle newTopVertex = mesh.add_vertex(topPoint);
      mesh::TriangleMesh::VertexHandle newBottomVertex = mesh.add_vertex(bottomPoint);

      if (e > 0) {
         // In the following, the order of vertices added to a face is very important!
         // The half edge data structure needs to be valid, normals need to show into the correct direction.

         // side faces (rectangle built up by two triangles)
         mesh.add_face(newTopVertex, newBottomVertex, lastBottomVertex);
         mesh.add_face(lastTopVertex, newTopVertex, lastBottomVertex);

         // bottom - "pike" ("pizza slices")
         mesh.add_face(newBottomVertex, bottomCenterVertex, lastBottomVertex);

         // top face ("pizza slices")
         mesh.add_face(topCenterVertex, newTopVertex, lastTopVertex);
      } else {
         firstTopVertex = newTopVertex;
         firstBottomVertex = newBottomVertex;
      }

      lastTopVertex = newTopVertex;
      lastBottomVertex = newBottomVertex;
   }

   // connect the first and the last sides
   // side faces (rectangle built up by two triangles)
   mesh.add_face(firstTopVertex, firstBottomVertex, lastBottomVertex);
   mesh.add_face(lastTopVertex, firstTopVertex, lastBottomVertex);
   // bottom - "pike"
   mesh.add_face(firstBottomVertex, bottomCenterVertex, lastBottomVertex);
   // top face
   mesh.add_face(topCenterVertex, firstTopVertex, lastTopVertex);

   // shift the mesh such that the centroid lies in 0,0,0 in its coordinate system (required)
   // the pike tip is at radius,radius,0 up to now, after shifting the centroid by centroidShift to 0,0,0,
   // the pike tip will be at -centroidShift.
   auto centroidShift = mesh::toWalberla(mesh::computeCentroid(mesh));
   mesh::translate(mesh, -centroidShift);

   pegPikeTipPosition = -centroidShift + mesa_pd::Vec3(radius, radius, real_t(0));
}

mesa_pd::data::ParticleStorage::iterator createPlane( mesa_pd::data::ParticleStorage& ps,
                                             mesa_pd::data::ShapeStorage& ss,
                                             const Vector3<real_t>& pos,
                                             const Vector3<real_t>& normal ) {
   auto p0              = ps.create(true);
   p0->getPositionRef() = pos;
   p0->getShapeIDRef()  = ss.create<mesa_pd::data::HalfSpace>( normal );
   p0->getOwnerRef()    = walberla::mpi::MPIManager::instance()->rank();
   p0->getTypeRef()     = 0;
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::NON_COMMUNICATING);
   return p0;
}

class DEM
{
public:
   DEM(const std::shared_ptr<mesa_pd::domain::BlockForestDomain>& domain, real_t dt, real_t mass)
         : dt_(dt)
         , domain_(domain)
   {
      sd_.setDampingT(0, 0, real_t(0));
      sd_.setFriction(0, 0, real_t(0));
      sd_.setParametersFromCOR(0, 0, real_t(0.9), dt*real_t(20), mass * real_t(0.5));

      sds_.setParametersFromCOR(0, 0, real_t(0.9), dt*real_t(20), mass * real_t(0.5));
      sds_.setCoefficientOfFriction(0,0,real_t(0.4));
      sds_.setStiffnessT(0,0,real_t(0.9) * sds_.getStiffnessN(0,0));
   }

   inline
   void operator()(const size_t idx1, const size_t idx2, mesa_pd::data::ParticleAccessorWithShape& ac)
   {

      ++contactsChecked_;
      if (double_cast_(idx1, idx2, ac, gcd_, ac)) {
         ++contactsDetected_;
         if (contact_filter_(gcd_.getIdx1(), gcd_.getIdx2(), ac, gcd_.getContactPoint(),
                             *domain_)) {
            ++contactsTreated_;
//            sd_(acd_.getIdx1(), acd_.getIdx2(), ac, acd_.getContactPoint(),
//                 acd_.getContactNormal(), acd_.getPenetrationDepth());
            sds_(gcd_.getIdx1(), gcd_.getIdx2(), ac, gcd_.getContactPoint(),
                 gcd_.getContactNormal(), gcd_.getPenetrationDepth(), dt_);
         }
      }
   }

   inline void resetCounters() {
      contactsChecked_ = 0;
      contactsDetected_ = 0;
      contactsTreated_ = 0;
   }

   inline int64_t getContactsChecked() const {
      return contactsChecked_;
   }

   inline int64_t getContactsDetected() const {
      return contactsDetected_;
   }

   inline int64_t getContactsTreated() const {
      return contactsTreated_;
   }

private:
   real_t dt_;
   mesa_pd::kernel::DoubleCast double_cast_;
   mesa_pd::mpi::ContactFilter contact_filter_;
   std::shared_ptr<mesa_pd::domain::BlockForestDomain> domain_;
   mesa_pd::kernel::SpringDashpot sd_ = mesa_pd::kernel::SpringDashpot(3);
   mesa_pd::kernel::SpringDashpotSpring sds_ = mesa_pd::kernel::SpringDashpotSpring(3);
   mesa_pd::collision_detection::GeneralContactDetection gcd_;
   int64_t contactsChecked_ = 0;
   int64_t contactsDetected_ = 0;
   int64_t contactsTreated_ = 0;
};

void updatePegPosition(mesa_pd::data::ParticleAccessorWithShape & accessor, walberla::id_t pegUID, real_t dt) {
   auto pegIdx = accessor.uidToIdx(pegUID);
   WALBERLA_CHECK(pegIdx != accessor.getInvalidIdx());
   auto newPegPosition = accessor.getPosition(pegIdx) + dt * accessor.getLinearVelocity(pegIdx);
   accessor.setPosition(pegIdx, newPegPosition);
}

int main( int argc, char ** argv ) {
   /// Setup
   Environment env(argc, argv);

   /// Config
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle mainConf = cfg->getBlock( "PegIntoSphereBed" );

   uint_t simulationSteps = mainConf.getParameter<uint_t>("simulationSteps", uint_t(1000));
   uint_t visSpacing = mainConf.getParameter<uint_t>("visSpacing", uint_t(100));
   real_t dt = mainConf.getParameter<real_t>("dt", real_t(0.0003));
   Vector3<real_t> shift = mainConf.getParameter<Vector3<real_t>>("shift", Vector3<real_t>(real_t(0.01)));

   real_t sphereBedHeight = mainConf.getParameter<real_t>("sphereBedHeight", real_t(1));
   real_t sphereRadius = mainConf.getParameter<real_t>("sphereRadius", real_t(0.5));
   real_t sphereSpacing = mainConf.getParameter<real_t>("sphereSpacing", real_t(0.5));
   real_t sphereDensity = mainConf.getParameter<real_t>("sphereDensity", real_t(1000));

   real_t pegBodyHeight = mainConf.getParameter<real_t>("pegBodyHeight", real_t(5));
   real_t pegPikeHeight = mainConf.getParameter<real_t>("pegPikeHeight", real_t(2));
   real_t pegRadius = mainConf.getParameter<real_t>("pegRadius", real_t(2));
   uint_t pegNumSideEdges = mainConf.getParameter<uint_t>("pegNumSideEdges", uint_t(4));
   Vector3<real_t> pegPikeTipPosition = mainConf.getParameter<Vector3<real_t>>("pegPikeTipPosition", Vector3<real_t>(real_t(1)));
   Vector3<real_t> pegVelocity = mainConf.getParameter<Vector3<real_t>>("pegVelocity", Vector3<real_t>(0, 0, real_t(-0.05)));

   std::string vtk_out = "vtk_out";

   /// BlockForest
   shared_ptr<BlockForest> forest = blockforest::createBlockForestFromConfig(mainConf);
   auto domainSize = forest->getDomain();

   /// MESAPD Domain
   auto domain = std::make_shared<mesa_pd::domain::BlockForestDomain>(forest);

   auto localDomain = forest->begin()->getAABB();
   for (auto& blk : *forest) {
      localDomain.merge(blk.getAABB());
   }

   /// MESAPD Data
   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(4);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   mesa_pd::data::ParticleAccessorWithShape ac(ps, ss);
   mesa_pd::data::LinkedCells lc(localDomain.getExtended(real_t(1)), real_t(2.1) * sphereRadius );
   mesa_pd::mpi::SyncNextNeighbors SNN;

   /// Mesh Peg
   shared_ptr<mesh::TriangleMesh> pegMesh = make_shared<mesh::TriangleMesh>();
   mesa_pd::Vec3 pegPikeTipOffset;
   createPeg(*pegMesh, pegPikeTipOffset, pegBodyHeight, pegPikeHeight, pegRadius, pegNumSideEdges);
   WALBERLA_LOG_INFO("Generated peg has " << pegMesh->n_vertices() << " vertices and " << pegMesh->n_faces() << " faces.");
   /*if (!OpenMesh::IO::write_mesh(pegMesh, "peg_mesh.ply")) {
      WALBERLA_ABORT("Error while writing peg mesh.");
   } else {
      WALBERLA_LOG_INFO("Wrote peg mesh to file.");
   }*/
   // We use simple distance calculation here due to a relatively low number of vertices. It it were much bigger, we would
   // certainly be better off using a DistanceOctree for optimization.
   auto pegTriDistance = make_shared<mesh::TriangleDistance<mesh::TriangleMesh>>(pegMesh);

   /// MESAPD Particles
   // Peg
   auto pegShapeID = ss->create<mesa_pd::data::ConvexPolyhedron>(*pegMesh);
   ss->shapes[pegShapeID]->updateMassAndInertia(real_t(1));
   auto pegShape = dynamic_cast<mesa_pd::data::ConvexPolyhedron*>(ss->shapes[pegShapeID].get());

   Vector3<real_t> pegPosition = pegPikeTipPosition - pegPikeTipOffset;
   //WALBERLA_CHECK(domainSize.contains(pegPosition), "The pegs position " << pegPosition << " (defined at the mesh centroid) is outside the domain AABB!");
   walberla::id_t pegUid = 0;

   auto pegParticle = ps->create();

   pegParticle->setShapeID(pegShapeID);
   pegParticle->setPosition(pegPosition);
   pegParticle->setOwner(walberla::MPIManager::instance()->rank());
   pegParticle->setLinearVelocity(pegVelocity);
   pegParticle->setInteractionRadius(pegShape->getBoundingSphereRadius());
   mesa_pd::data::particle_flags::set(pegParticle->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(pegParticle->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   mesa_pd::data::particle_flags::set(pegParticle->getFlagsRef(), mesa_pd::data::particle_flags::NON_COMMUNICATING);

   pegUid = pegParticle->getUid();

   // Spheres
   auto sphereShapeId = ss->create<mesa_pd::data::Sphere>(sphereRadius);
   ss->shapes[sphereShapeId]->updateMassAndInertia(real_t(sphereDensity));

   auto extendedPegAABB = mesh::computeAABB(*pegMesh).getExtended(sphereRadius);
   for (auto& iBlk : *forest) {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(),
                                            Vector3<real_t>(sphereSpacing) * real_c(0.5) + shift,
                                            sphereSpacing)) {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));

         if (pt[2] > sphereBedHeight) continue;
         // check if a sphere at position pt would protrude into the peg
         if (extendedPegAABB.contains(pt)) {
            auto dSq = pegTriDistance->sqSignedDistance(mesh::toOpenMesh(pt - pegPosition));
            if (dSq < sphereRadius*sphereRadius) {
               continue;
            }
         }

         auto sphereParticle = ps->create();

         sphereParticle->setShapeID(sphereShapeId);
         sphereParticle->setPosition(pt);
         sphereParticle->setOwner(walberla::MPIManager::instance()->rank());
         sphereParticle->setInteractionRadius(sphereRadius);
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   // Confining Planes
   auto planeShift = (sphereSpacing - sphereRadius - sphereRadius) * real_t(0.5);
   auto confiningDomain = domainSize.getExtended(planeShift);
   if (!forest->isPeriodic(0)) {
      createPlane(*ps, *ss, confiningDomain.minCorner()+ shift, Vector3<real_t>(+1,0,0));
      createPlane(*ps, *ss, confiningDomain.maxCorner()+ shift, Vector3<real_t>(-1,0,0));
   }
   if (!forest->isPeriodic(1)) {
      createPlane(*ps, *ss, confiningDomain.minCorner()+ shift, Vector3<real_t>(0,+1,0));
      createPlane(*ps, *ss, confiningDomain.maxCorner()+ shift, Vector3<real_t>(0,-1,0));
   }
   if (!forest->isPeriodic(2)) {
      createPlane(*ps, *ss, confiningDomain.minCorner()+ shift, Vector3<real_t>(0,0,+1));
      createPlane(*ps, *ss, confiningDomain.maxCorner()+ shift, Vector3<real_t>(0,0,-1));
   }

   SNN(*ps, *domain);

   /// VTK Output
   // domain output
   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition(forest, "domain_decomposition",
                                                                   uint_t(1), vtk_out, "simulation_step");
   vtkDomainOutput->write();
   // mesapd mesh output
   mesa_pd::MeshParticleVTKOutput<mesh::PolyMesh> meshParticleVTK(ps, "mesh", visSpacing);
   meshParticleVTK.addFaceOutput<mesa_pd::data::SelectParticleUid>("Uid");
   auto surfaceVelocityDataSource = make_shared<mesa_pd::SurfaceVelocityVertexDataSource<mesh::PolyMesh,
         mesa_pd::data::ParticleAccessorWithShape>>("SurfaceVelocity", ac);
   meshParticleVTK.addVertexDataSource(surfaceVelocityDataSource);
   // mesapd particle output
   auto particleVtkOutput = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
   particleVtkOutput->addOutput<mesa_pd::data::SelectParticleInteractionRadius>("interactionRadius");
   particleVtkOutput->setParticleSelector([sphereShapeId](const mesa_pd::data::ParticleStorage::iterator& pIt){
      return pIt->getShapeID() == sphereShapeId;
   });
   auto particleVtkWriter = walberla::vtk::createVTKOutput_PointData(particleVtkOutput, "particles", visSpacing,
                                                                     vtk_out, "simulation_step");


   /// MESAPD kernels
   mesa_pd::kernel::ExplicitEuler explicitEuler(dt);
   DEM dem(domain, dt, real_t(1));
   mesa_pd::kernel::InsertParticleIntoLinkedCells ipilc;
   mesa_pd::mpi::ReduceProperty RP;
   mesa_pd::mpi::ReduceContactHistory RCH;

   Vector3<real_t> globalAcceleration(real_t(0), real_t(0), real_t(-6));
   auto addGravitationalForce = [&globalAcceleration](const size_t idx, mesa_pd::data::ParticleAccessorWithShape& ac_) {
      auto mass = real_t(1) / ac_.getInvMass(idx);
      auto force = mass * globalAcceleration;
      mesa_pd::addForceAtomic(idx, ac_, force);
   };

   for (uint_t i = 0; i < simulationSteps; ++i) {
      if(i % visSpacing == 0){
         WALBERLA_LOG_INFO_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }

      // VTK
      meshParticleVTK(ac);
      particleVtkWriter->write();

      // Prepare Data Structures
      lc.clear();
      ps->forEachParticle(true, mesa_pd::kernel::SelectAll(), ac, ipilc, ac, lc);

      // Collision Resolution
      dem.resetCounters();
      lc.forEachParticlePairHalf(true, mesa_pd::kernel::SelectAll(), ac, dem, ac);
      RP.operator()<mesa_pd::ForceTorqueNotification>(*ps);
      RCH(*ps);

      // Force Application
      ps->forEachParticle(true, mesa_pd::kernel::SelectLocal(), ac, addGravitationalForce, ac);

      // Integration
      ps->forEachParticle(true, mesa_pd::kernel::SelectLocal(), ac, explicitEuler, ac);

      updatePegPosition(ac, pegUid, dt);

      SNN(*ps, *domain);
   }

   return EXIT_SUCCESS;
}

}

int main( int argc, char* argv[] ) {
   return walberla::main( argc, argv );
}