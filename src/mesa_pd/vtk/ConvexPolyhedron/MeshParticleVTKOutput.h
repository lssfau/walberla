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

#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/data/shape/ConvexPolyhedron.h>
#include <mesa_pd/vtk/ConvexPolyhedron/Types.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/DataSourceAdapters.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/FaceDataSource.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/OutputSelectorFaceDataSource.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/OutputSelectorVertexDataSource.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/VertexDataSource.h>
#include <mesa_pd/vtk/ConvexPolyhedron/tesselation/ConvexPolyhedronTesselation.h>
#include <mesh_common/vtk/DistributedVTKMeshWriter.h>
#include <utility>

namespace walberla {
namespace mesa_pd {

template<typename MeshType>
class MeshParticleVTKOutput {
   static_assert(MeshType::IsPolyMesh == 1, "We need polygonal meshes here!");

public:
   using ParticleSelectorFunc = std::function<bool (const walberla::mesa_pd::data::ParticleStorage::iterator& pIt)>;

   typedef typename mesh::DistributedVTKMeshWriter<MeshType>::Vertices Vertices;

   MeshParticleVTKOutput( shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps,
                          shared_ptr<walberla::mesa_pd::data::ShapeStorage> ss, const std::string & identifier,
                          const uint_t writeFrequency, const std::string & baseFolder = "vtk_out")
                          : ps_(std::move(ps)), ss_(std::move(ss)), mesh_(make_shared<MeshType>()),
                          faceToParticleIdxManager_(*mesh_, "particle"), vertexToParticleIdxManager_(*mesh_, "particle"),
                          meshWriter_(mesh_, identifier, writeFrequency, baseFolder) {

   }

   void setParticleSelector( const ParticleSelectorFunc& func) {particleSelector_ = func;}
   void assembleMesh();

   template <typename Selector>
   void addVertexOutput(const std::string& name);
   template<typename Selector>
   void addFaceOutput(const std::string &name);

   template <typename DataSourceType>
   void addVertexDataSource(const shared_ptr<DataSourceType> & dataSource);
   template <typename DataSourceType>
   void addFaceDataSource(const shared_ptr<DataSourceType> & dataSource);

   void operator()() {
      assembleMesh();
      // the mesh writer writes the mesh to vtk files and adds properties as defined by the data sources
      meshWriter_();
      mesh_->clean(); // the output mesh is no longer needed, thus discard its contents
   }

private:
   const shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps_;
   const shared_ptr<walberla::mesa_pd::data::ShapeStorage> ss_;
   shared_ptr<MeshType> mesh_; ///< the output mesh

   // these "managers" (which are just maps basically) map the faces and vertices of the output mesh to particles
   // such that we can assign those faces and vertices properties of the particle (i.e. ID or velocity)
   ParticleIdxFacePropertyManager<MeshType> faceToParticleIdxManager_;
   ParticleIdxVertexPropertyManager<MeshType> vertexToParticleIdxManager_;

   mesh::DistributedVTKMeshWriter<MeshType> meshWriter_;

   ParticleSelectorFunc particleSelector_ = [](const walberla::mesa_pd::data::ParticleStorage::iterator& /*pIt*/){
      return true;
   };
};

/**
 * \brief Adds a output selector for the vertices of the mesh
 * Similar to MESAPD's vtk output, one can add a property selector to select a property for a vertex
 * corresponding to the associated particle.
 *
 * \tparam MeshType
 * \tparam Selector Type of the selector
 * \param name Name of the output
 */
template<typename MeshType>
template<typename Selector>
void MeshParticleVTKOutput<MeshType>::addVertexOutput(const std::string &name) {
   typedef OutputSelectorVertexDataSource<MeshType, Selector, typename Selector::return_type> DataSourceType;
   auto ds = make_shared<DataSourceType>(name, Selector());
   addVertexDataSource(ds);
}

/**
 * \brief Adds a output selector for the faces of the mesh
 * Similar to MESAPD's vtk output, one can add a property selector to select a property for a face
 * corresponding to the associated particle.
 *
 * \tparam MeshType
 * \tparam Selector Type of the selector
 * \param name Name of the output
 */
template<typename MeshType>
template<typename Selector>
void MeshParticleVTKOutput<MeshType>::addFaceOutput(const std::string &name) {
   typedef OutputSelectorFaceDataSource<MeshType, Selector, typename Selector::return_type> DataSourceType;
   auto ds = make_shared<DataSourceType>(name, Selector());
   addFaceDataSource(ds);
}

/**
 * \brief Add a vertex data source assigning a piece of data to a vertex.
 * \tparam MeshType
 * \tparam DataSourceType Type of data source (has to be derived from VertexDataSource).
 * \param dataSource Data source responsible for picking data for a vertex.
 */
template<typename MeshType>
template<typename DataSourceType>
void MeshParticleVTKOutput<MeshType>::addVertexDataSource(const shared_ptr<DataSourceType> & dataSource) {
   typedef internal::VertexDataSourceAdapter<MeshType, typename DataSourceType::ComponentType> AdapterType;
   meshWriter_.addDataSource(make_shared<AdapterType>(dataSource, vertexToParticleIdxManager_, ps_));
}

/**
 * \brief Add a face data source assigning a piece of data to a face.
 * \tparam MeshType
 * \tparam DataSourceType Type of data source (has to be derived from FaceDataSource).
 * \param dataSource Data source responsible for picking data for a face.
 */
template<typename MeshType>
template<typename DataSourceType>
void MeshParticleVTKOutput<MeshType>::addFaceDataSource(const shared_ptr<DataSourceType> & dataSource) {
   typedef internal::FaceDataSourceAdapter<MeshType, typename DataSourceType::ComponentType> AdapterType;
   meshWriter_.addDataSource(make_shared<AdapterType>(dataSource, faceToParticleIdxManager_, ps_));
}


/**
 * \brief Creates the output mesh and writes it to mesh_.
 * \tparam MeshType
 */
template<typename MeshType>
void MeshParticleVTKOutput<MeshType>::assembleMesh() {
   // those will save the newly created vertices and faces for each mesh
   // to make it possible to map a vertex/face to the corresponding particle
   std::vector<typename MeshType::VertexHandle> newVertices;
   std::vector<typename MeshType::FaceHandle> newFaces;

   // ensure the mesh is empty, as this will contain the new output
   mesh_->clean();

   // then iterate over every particle and tessellate it to include it in the output mesh
   for (auto pIt = ps_->begin(); pIt != ps_->end(); ++pIt) {
      if (!particleSelector_(pIt)) continue;

      auto& shape = ss_->shapes[pIt->getShapeID()];

      newVertices.clear();
      newFaces.clear();

      if (shape->getShapeType() == walberla::mesa_pd::data::ConvexPolyhedron::SHAPE_TYPE) {
         const auto& convexPolyhedron = *static_cast<walberla::mesa_pd::data::ConvexPolyhedron*>(shape.get());
         const auto& particle = *pIt;

         // tessellate: add the shape at the particle's position into the output mesh
         tesselate(convexPolyhedron, particle, mesh_, newVertices, newFaces);
      }

      // save particle idx to managers
      for (const auto & vertex: newVertices) {
         vertexToParticleIdxManager_[vertex] = pIt->getIdx();
      }
      for (const auto & face: newFaces) {
         faceToParticleIdxManager_[face] = pIt->getIdx();
      }
   }

   //WALBERLA_LOG_INFO("MESA-PD VTK output mesh contains " << mesh_->n_vertices() << " vertices.")
}


}
}