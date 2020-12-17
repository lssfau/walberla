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
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/vtk/ConvexPolyhedron/data_sources/FaceDataSource.h>

namespace walberla {
namespace mesa_pd {

/**
 * \brief Generic data source to pick data for a face using a particle selector.
 * \attention The underlying mesh data sources don't support Vec3's etc. natively, thus specializations need to be given.
 * \tparam MeshType
 * \tparam Selector Type of the selector.
 * \tparam Type Type of the data.
 */
template<typename MeshType, typename Selector, typename Type = typename Selector::return_type>
class OutputSelectorFaceDataSource : public FaceDataSource<MeshType, Type> {
public:
   typedef FaceDataSource<MeshType, Type> Base;
   typedef typename Base::Faces Faces;

   typedef Type ComponentType;

   OutputSelectorFaceDataSource(const std::string& name, Selector selector) : Base(name), selector_(selector) { }

   virtual uint_t numComponents() {
      return uint_t(1);
   }

   using Base::getData;
   virtual void getData( const MeshType &, const Faces & faces, std::vector<Type> & data,
                         const ParticleIdxFacePropertyManager<MeshType> & faceToParticleIdxManager,
                         shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps) {
      for (const auto & face: faces) {
         size_t particleIdx = faceToParticleIdxManager[face];
         auto p = (*ps)[particleIdx];
         data.push_back(selector_(p));
      }
   };

private:
   Selector selector_;
};

/**
 * \brief Data Source for particle selectors specialized for Vec3
 * \tparam MeshType
 * \tparam Selector
 * \tparam Type
 */
template<typename MeshType, typename Selector, typename Type>
class OutputSelectorFaceDataSource<MeshType, Selector, Vector3<Type>> : public FaceDataSource<MeshType, Type> {
public:
   typedef FaceDataSource<MeshType, Type> Base;
   typedef typename Base::Faces Faces;

   typedef Type ComponentType;

   OutputSelectorFaceDataSource(const std::string& name, Selector selector) : Base(name), selector_(selector) { }

   virtual uint_t numComponents() {
      return uint_t(3);
   }

   using Base::getData;
   virtual void getData( const MeshType &, const Faces & faces, std::vector<Type> & data,
                         const ParticleIdxFacePropertyManager<MeshType> & faceToParticleIdxManager,
                         shared_ptr<walberla::mesa_pd::data::ParticleStorage> ps) {
      for (const auto & face: faces) {
         size_t particleIdx = faceToParticleIdxManager[face];
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