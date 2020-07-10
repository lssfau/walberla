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
//! \file ParticleOutputSelectorVertexDataSource.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/VertexDataSource.h>

namespace walberla {
namespace mesa_pd {

template<typename MeshType, typename Selector, typename Type = typename Selector::return_type>
class OutputSelectorVertexDataSource : public VertexDataSource<MeshType, Type> {
public:
   typedef VertexDataSource<MeshType, Type> Base;
   typedef typename Base::Vertices Vertices;

   typedef Type ComponentType;

   OutputSelectorVertexDataSource(const std::string& name, Selector selector) : Base(name), selector_(selector) { }

   virtual uint_t numComponents() {
      return uint_t(1);
   }

   using Base::getData;
   virtual void getData( const MeshType &, const Vertices & vertices, std::vector<Type> & data,
                         const ParticleIdxVertexPropertyManager<MeshType> & vertexToParticleIdxManager,
                         shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps) {
      for (const auto & vertex: vertices) {
         size_t particleIdx = vertexToParticleIdxManager[vertex];
         auto p = (*ps)[particleIdx];
         data.push_back(selector_(p));
      }
   };

private:
   Selector selector_;
};

/**
 * \brief Data Source specialized for Vec3
 * \tparam MeshType
 * \tparam Selector
 * \tparam Type
 */
template<typename MeshType, typename Selector, typename Type>
class OutputSelectorVertexDataSource<MeshType, Selector, Vector3<Type>> : public VertexDataSource<MeshType, Type> {
public:
   typedef VertexDataSource<MeshType, Type> Base;
   typedef typename Base::Vertices Vertices;

   typedef Type ComponentType;

   OutputSelectorVertexDataSource(const std::string& name, Selector selector) : Base(name), selector_(selector) { }

   virtual uint_t numComponents() {
      return uint_t(3);
   }

   using Base::getData;
   virtual void getData( const MeshType &, const Vertices & vertices, std::vector<Type> & data,
                         const ParticleIdxVertexPropertyManager<MeshType> & vertexToParticleIdxManager,
                         shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps) {
      for (const auto & vertex: vertices) {
         size_t particleIdx = vertexToParticleIdxManager[vertex];
         auto p = (*ps)[particleIdx];
         const Vector3<Type> d = selector_(p);
         data.push_back(d[0]);
         data.push_back(d[1]);
         data.push_back(d[2]);
      }
   };

private:
   Selector selector_;
};

}
}