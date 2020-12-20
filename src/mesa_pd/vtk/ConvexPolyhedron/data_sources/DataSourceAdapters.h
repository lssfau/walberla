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

#include <core/DataTypes.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/FaceDataSource.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/VertexDataSource.h>
#include <mesh_common/vtk/DistributedVTKMeshWriter.h>

namespace walberla {
namespace mesa_pd {
namespace internal {

/**
 * \brief Adapts a vertex data source for the MESAPD mesh output to the generic vertex data source class.
 * \tparam MeshType
 * \tparam T output type
 */
template<typename MeshType, typename T >
class VertexDataSourceAdapter : public mesh::DistributedVTKMeshWriter<MeshType>::template VertexDataSource<T> {
public:
   typedef typename mesh::DistributedVTKMeshWriter<MeshType>::template VertexDataSource<T>::Vertices Vertices;

   VertexDataSourceAdapter( const shared_ptr<VertexDataSource<MeshType, T>> & vertexDataSource,
         const ParticleIdxVertexPropertyManager<MeshType> & vertexToParticleIdxManager,
         const shared_ptr<walberla::mesa_pd::data::ParticleStorage> & ps)
         : mesh::DistributedVTKMeshWriter<MeshType>::template VertexDataSource<T>( vertexDataSource->name() ),
               vertexDataSource_(vertexDataSource), vertexToParticleIdxManager_(vertexToParticleIdxManager),
               ps_(ps) { }

   virtual void getData(const MeshType & mesh, const Vertices & vertices, std::vector<T> & data) {
      return vertexDataSource_->getData( mesh, vertices, data, vertexToParticleIdxManager_, ps_ );
   };

   virtual uint_t numComponents() { return vertexDataSource_->numComponents(); }

protected:
   shared_ptr<VertexDataSource<MeshType, T>> vertexDataSource_;
   const ParticleIdxVertexPropertyManager<MeshType> & vertexToParticleIdxManager_;
   const shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps_;
};

/**
 * \brief Adapts a face data source for the MESAPD mesh output to the generic face data source class.
 * \tparam MeshType
 * \tparam T output type
 */
template<typename MeshType, typename T >
class FaceDataSourceAdapter : public mesh::DistributedVTKMeshWriter<MeshType>::template FaceDataSource<T> {
public:
   typedef typename mesh::DistributedVTKMeshWriter<MeshType>::template FaceDataSource<T>::Faces Faces;

   FaceDataSourceAdapter(const shared_ptr<FaceDataSource<MeshType, T>> & faceDataSource,
         const ParticleIdxFacePropertyManager<MeshType> & faceToParticleIdxManager,
         shared_ptr<walberla::mesa_pd::data::ParticleStorage>  ps)
   : mesh::DistributedVTKMeshWriter<MeshType>::template FaceDataSource<T>(faceDataSource->name()),
         faceDataSource_(faceDataSource), faceToParticleIdxManager_(faceToParticleIdxManager),
         ps_(std::move(ps)) { }

   virtual void getData(const MeshType & mesh, const Faces & faces, std::vector<T> & data) {
      return faceDataSource_->getData( mesh, faces, data, faceToParticleIdxManager_, ps_ );
   };

   virtual uint_t numComponents() { return faceDataSource_->numComponents(); }

protected:
   shared_ptr<FaceDataSource<MeshType, T>> faceDataSource_;
   const ParticleIdxFacePropertyManager<MeshType> & faceToParticleIdxManager_;
   const shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps_;
};

}
}
}