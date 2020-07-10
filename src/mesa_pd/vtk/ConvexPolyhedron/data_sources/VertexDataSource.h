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
//! \file ParticleVertexDataSource.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <mesa_pd/vtk/ConvexPolyhedron/Types.h>
#include <mesh_common/vtk/DistributedVTKMeshWriter.h>

namespace walberla {
namespace mesa_pd {

template<typename MeshType, typename T>
class VertexDataSource {
public:
   typedef typename mesh::DistributedVTKMeshWriter<MeshType>::Vertices Vertices;

   explicit VertexDataSource(const std::string & name) : name_(name) {}

   const std::string & name() { return name_; }
   virtual uint_t numComponents() = 0;
   virtual void getData(const MeshType &, const Vertices &vertices, std::vector<T> &data,
         const ParticleIdxVertexPropertyManager<MeshType> & vertexToParticleIdxManager,
         shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps) = 0;

   virtual ~VertexDataSource() {}

protected:
   std::string name_;
};

}
}