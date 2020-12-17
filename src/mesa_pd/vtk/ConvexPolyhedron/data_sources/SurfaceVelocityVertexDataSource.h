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
#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/VertexDataSource.h>

namespace walberla {
namespace mesa_pd {

template<typename MeshType, typename Accessor, typename Type = real_t>
class SurfaceVelocityVertexDataSource : public VertexDataSource<MeshType, Type> {
public:
   typedef VertexDataSource<MeshType, Type> Base;
   typedef typename Base::Vertices Vertices;

   typedef Type ComponentType;

   SurfaceVelocityVertexDataSource(const std::string& name, const Accessor & ac) : Base(name), ac_(ac) { }

   virtual uint_t numComponents() {
      return uint_t(3);
   }

   using Base::getData;
   virtual void getData( const MeshType & mesh, const Vertices & vertices, std::vector<Type> & data,
                         const ParticleIdxVertexPropertyManager<MeshType> & vertexToParticleIdxManager,
                         shared_ptr<walberla::mesa_pd::data::ParticleStorage>) {
      for (const auto & vertex: vertices) {
         size_t particleIdx = vertexToParticleIdxManager[vertex];
         auto vertexPosition = mesh::toWalberlaNumericCast<real_t>(mesh.point(vertex));
         const Vector3<Type> d = walberla::mesa_pd::getVelocityAtWFPoint(particleIdx, ac_, vertexPosition);
         data.push_back(d[0]);
         data.push_back(d[1]);
         data.push_back(d[2]);
      }
   };
private:
   const Accessor & ac_;
};

}
}