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

#pragma once

#include "core/Abort.h"

#include <waLBerlaDefinitions.h>

#include <mesa_pd/data/shape/BaseShape.h>


#ifdef WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE

/**
 * Full ConvexPolyhedron shape supporting all features, as the mesh_common module is available.
 */

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshOperations.h"

#include <core/math/Constants.h>

#include "core/mpi/BufferDataTypeExtensions.h"
#include "mesh_common/OpenMeshBufferTypeExtensions.h"

namespace walberla {
namespace mesa_pd {
namespace data {

using namespace walberla::mesa_pd;

/**
 * \brief Convex Polyhedron shape containing a mesh.
 *
 * Use ConvexPolyhedron::updateMeshQuantities() to initialize the particle properties (volume, mass, inertia).
 * This is done automatically upon initialization (and after shape unpacking), but if the mesh is manipulated
 * while the particle already exists, it has to be called manually.
 *
 * \attention The origin of the meshes coordinate system is required to be at the position of its centroid! As such, the particle and mesh COS also overlap.
 */
class ConvexPolyhedron : public walberla::mesa_pd::data::BaseShape {
public:
   explicit ConvexPolyhedron(const mesh::TriangleMesh& mesh = mesh::TriangleMesh())
         : BaseShape(ConvexPolyhedron::SHAPE_TYPE), mesh_(mesh)
   {
      if (mesh_.n_vertices() > 0) {
         WALBERLA_CHECK_GREATER_EQUAL(mesh_.n_vertices(), 4);
         updateMeshQuantities();
      }
   }

   constexpr static int IS_AVAILABLE = true;

   const mesh::TriangleMesh& getMesh() const { return mesh_; }
   
   mesh::TriangleMesh::VertexHandle supportVertex( const mesh::TriangleMesh::Normal & d,
                                                   const mesh::TriangleMesh::VertexHandle startVertex ) const;

   void updateMeshQuantities();

   real_t getBoundingSphereRadius() const;
   real_t getVolume() const override;
   void updateMassAndInertia(const real_t density) override;

   Vec3 support( const Vec3& d ) const override;

   void pack(walberla::mpi::SendBuffer& buf) override;
   void unpack(walberla::mpi::RecvBuffer& buf) override;

   constexpr static int SHAPE_TYPE = 5; ///< Unique shape type identifier for convex polyhedrons.\ingroup mesa_pd_shape

   bool meshQuantitiesAvailable = false;
   mesh::TriangleMesh::VertexHandle octandVertices_[8];
   
private:
   mesh::TriangleMesh mesh_;

   real_t unitMass_; // mass for density 1, equals volume
   Mat3 unitInertiaBF_; // inertia for density 1 relative to centroid
   real_t boundingSphereRadius_;
};

/**
 * \brief Get the bounding sphere radius around the centroid.
 * \return Interaction radius.
 */
inline real_t ConvexPolyhedron::getBoundingSphereRadius() const {
   return boundingSphereRadius_;
}

inline real_t ConvexPolyhedron::getVolume() const {
   WALBERLA_CHECK(meshQuantitiesAvailable,
         "unitMass and unitInertia first have to be initialized using updateMeshQuantities!");

   return unitMass_; // V = m*d = m*1 = m
}

inline void ConvexPolyhedron::updateMassAndInertia(const real_t density) {
   WALBERLA_CHECK(meshQuantitiesAvailable,
                  "unitMass and unitInertia first have to be initialized using updateMeshQuantities!");

   const real_t m = unitMass_ * density;
   const Mat3 I = unitInertiaBF_ * density;

   mass_ = m;
   invMass_ = real_t(1.) / m;

   inertiaBF_ = I;
   invInertiaBF_ = I.getInverse();
}

inline void ConvexPolyhedron::updateMeshQuantities() {
   WALBERLA_CHECK_GREATER(mesh_.n_vertices(), 0, "Cannot compute mesh quantities for an empty mesh!");

   mesh_.request_face_normals();
   mesh_.update_face_normals();

   octandVertices_[0] = supportVertex(mesh::TriangleMesh::Normal(real_t(1), real_t(1), real_t(1)), *mesh_.vertices_begin());
   octandVertices_[1] = supportVertex(mesh::TriangleMesh::Normal(real_t(1), real_t(1), real_t(-1)), *mesh_.vertices_begin());
   octandVertices_[2] = supportVertex(mesh::TriangleMesh::Normal(real_t(1), real_t(-1), real_t(1)), *mesh_.vertices_begin());
   octandVertices_[3] = supportVertex(mesh::TriangleMesh::Normal(real_t(1), real_t(-1), real_t(-1)), *mesh_.vertices_begin());
   octandVertices_[4] = supportVertex(mesh::TriangleMesh::Normal(real_t(-1), real_t(1), real_t(1)), *mesh_.vertices_begin());
   octandVertices_[5] = supportVertex(mesh::TriangleMesh::Normal(real_t(-1), real_t(1), real_t(-1)), *mesh_.vertices_begin());
   octandVertices_[6] = supportVertex(mesh::TriangleMesh::Normal(real_t(-1), real_t(-1), real_t(1)), *mesh_.vertices_begin());
   octandVertices_[7] = supportVertex(mesh::TriangleMesh::Normal(real_t(-1), real_t(-1), real_t(-1)), *mesh_.vertices_begin());

   Vec3 centroid;
   mesh::computeMassProperties(mesh_, real_t(1), centroid, unitInertiaBF_, unitMass_);
   WALBERLA_CHECK_FLOAT_EQUAL(centroid, Vector3<real_t>(0,0,0), "The mesh has to have its centroid at the origin of its coordinate system! Use mesh::computeCentroid and mesh::translate.");

   real_t maxSqRadius(0);
   for(auto vh : mesh_.vertices()) {
      auto v = mesh::toWalberla(mesh_.point(vh));
      auto centroidToVSqr = v.sqrLength();

      if (centroidToVSqr > maxSqRadius) {
         maxSqRadius = centroidToVSqr;
      }
   }
   boundingSphereRadius_ = std::sqrt(maxSqRadius);

   meshQuantitiesAvailable = true;
}

inline Vec3 ConvexPolyhedron::support( const Vec3& d ) const {
   WALBERLA_CHECK_UNEQUAL(mesh_.n_vertices(), 0, "Cannot compute support for a mesh with 0 vertices.");
   WALBERLA_CHECK(meshQuantitiesAvailable, "Octand vertices of this mesh first have to be initialized using updateMeshQuantities!");
   // taken from pe implementation

   if (math::equal(d.length(), real_t(0))) return Vec3(0,0,0);

   // d is already in the mesh coordinate system, as the origins of the particle and the mesh COS overlap
   mesh::TriangleMesh::Normal d_loc = mesh::toOpenMesh(d);

   mesh::TriangleMesh::VertexHandle startVertex;
   if(d_loc[0] >= real_t( 0 ))
   {
      if(d_loc[1] >= real_t( 0 ))
      {
         startVertex = d_loc[2] >= real_t( 0 ) ? octandVertices_[0] : octandVertices_[1];
      }
      else // d_loc[1] < 0
      {
         startVertex = d_loc[2] >= real_t( 0 ) ? octandVertices_[2] : octandVertices_[3];
      }
   }
   else // d_loc[0] < 0
   {
      if(d_loc[1] >= real_t( 0 ))
      {
         startVertex = d_loc[2] >= real_t( 0 ) ? octandVertices_[4] : octandVertices_[5];
      }
      else // d_loc[1] < 0
      {
         startVertex = d_loc[2] >= real_t( 0 ) ? octandVertices_[6] : octandVertices_[7];
      }
   }

   mesh::TriangleMesh::VertexHandle vh = supportVertex( d_loc, startVertex );

   // the resulting vertex has to be shifted back into the mesh particle coordinate system
   auto relativeSupport = mesh::toWalberla(mesh_.point( vh ));

   //WALBERLA_LOG_INFO("Conv poly support: " << relativeSupport << ", d: " << d);

   return relativeSupport;
}

inline mesh::TriangleMesh::VertexHandle ConvexPolyhedron::supportVertex( const mesh::TriangleMesh::Normal & d,
      const mesh::TriangleMesh::VertexHandle startVertex ) const {
   // taken from pe implementation

   mesh::TriangleMesh::VertexHandle maxScalarProductVertex = startVertex;
   real_t maxScalarProduct = mesh_.point(maxScalarProductVertex) | d;

   bool isExtremum = false;
   while( !isExtremum ) {
      isExtremum = true;
      for(auto vh : mesh_.vv_range( maxScalarProductVertex ))
      {
         const real_t sp = mesh_.point(vh) | d;
         if(sp > maxScalarProduct)
         {
            isExtremum = false;
            maxScalarProductVertex = vh;
            maxScalarProduct = sp;
            break;
         }
      }
   }

   return maxScalarProductVertex;
}

inline void ConvexPolyhedron::pack(walberla::mpi::SendBuffer& buf){
   // partially taken from pe implementation

   BaseShape::pack(buf);

   buf << mesh_.n_vertices(); // number of vertices

   int dbgIndex = 0;
   WALBERLA_UNUSED(dbgIndex);
   for(auto vh : mesh_.vertices()) {
      WALBERLA_ASSERT_EQUAL( vh.idx(), dbgIndex++ ); // assume vertices are compactly stored

      buf << mesh_.point(vh);
   }

   buf << mesh_.n_faces();

   for(auto fh: mesh_.faces()) {
      for (auto vhIt = mesh_.cfv_ccwbegin(fh); vhIt != mesh_.cfv_ccwend(fh); ++vhIt) {
         WALBERLA_ASSERT_GREATER_EQUAL(vhIt->idx(), 0);
         WALBERLA_ASSERT_LESS(vhIt->idx(), mesh_.n_vertices());

         buf << vhIt->idx();
      }
   }
}

inline void ConvexPolyhedron::unpack(walberla::mpi::RecvBuffer& buf){
   // partially taken from pe implementation

   BaseShape::unpack(buf);

   WALBERLA_CHECK_EQUAL(mesh_.n_vertices(), 0, "Mesh needs to be empty!");
   WALBERLA_CHECK_EQUAL(mesh_.n_faces(), 0, "Mesh needs to be empty!");

   size_t numVertices;
   buf >> numVertices;
   std::vector<mesh::TriangleMesh::VertexHandle> vertexHandles(numVertices);
   for(size_t i = 0; i < numVertices; ++i) {
      mesh::TriangleMesh::Point p;
      buf >> p;
      vertexHandles[size_t(i)] = mesh_.add_vertex( p );
   }

   size_t numFaces;
   buf >> numFaces;
   for(size_t i = 0; i < numFaces; ++i){
      int v0;
      int v1;
      int v2;
      buf >> v0 >> v1 >> v2;
      WALBERLA_ASSERT_GREATER_EQUAL( v0, 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( v1, 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( v2, 0 );
      WALBERLA_ASSERT_LESS( v0, numVertices );
      WALBERLA_ASSERT_LESS( v1, numVertices );
      WALBERLA_ASSERT_LESS( v2, numVertices );

      mesh_.add_face( vertexHandles[size_t(v0)], vertexHandles[size_t(v1)], vertexHandles[size_t(v2)] );
   }

   updateMeshQuantities();
}

}
}
}

#else

/**
 * Replacement for the "full" ConvexPolyhedron, which will throw compile time errors if any of its features are
 * actually used.
 */

namespace walberla {

namespace mesh {
// forward declaration failing if ConvexPolyhedron is actually used
class TriangleMesh;
}

namespace mesa_pd {
namespace data {

using namespace walberla::mesa_pd;

class ConvexPolyhedron : public walberla::mesa_pd::data::BaseShape {
public:
   explicit ConvexPolyhedron(const mesh::TriangleMesh&) : BaseShape(ConvexPolyhedron::SHAPE_TYPE) {
      WALBERLA_ABORT("Shape ConvexPolyhedron is not available! Ensure waLBerla is configured with OpenMesh support.");
   }

   ConvexPolyhedron() : BaseShape(ConvexPolyhedron::SHAPE_TYPE) {
      WALBERLA_ABORT("Shape ConvexPolyhedron is not available! Ensure waLBerla is configured with OpenMesh support.");
   }

   constexpr static int IS_AVAILABLE = false;
   constexpr static int SHAPE_TYPE = 5; ///< Unique shape type identifier for convex polyhedrons.\ingroup mesa_pd_shape

   void updateMassAndInertia(const real_t /*density*/) override {
      WALBERLA_ABORT("Shape ConvexPolyhedron is not available! Ensure waLBerla is configured with OpenMesh support.");
   };

   real_t getVolume() const override {
      WALBERLA_ABORT("Shape ConvexPolyhedron is not available! Ensure waLBerla is configured with OpenMesh support.");
   }

   Vec3 support( const Vec3& /*d*/ ) const override {
      WALBERLA_ABORT("Shape ConvexPolyhedron is not available! Ensure waLBerla is configured with OpenMesh support.");
   }

   void pack(walberla::mpi::SendBuffer& /*buf*/) override {
      WALBERLA_ABORT("Shape ConvexPolyhedron is not available! Ensure waLBerla is configured with OpenMesh support.");
   };

   void unpack(walberla::mpi::RecvBuffer& /*buf*/) override {
      WALBERLA_ABORT("Shape ConvexPolyhedron is not available! Ensure waLBerla is configured with OpenMesh support.");
   };
};

}
}
}


#endif