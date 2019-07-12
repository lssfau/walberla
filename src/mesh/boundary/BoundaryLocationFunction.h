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
//! \file BoundaryLocationFunction.h
//! \ingroup mesh
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

namespace walberla {

template< typename MeshDistanceType, typename MeshType >
struct BoundaryLocationFunction
{
   BoundaryLocationFunction(const shared_ptr< MeshDistanceType >& meshDistanceObject,
                            const shared_ptr< mesh::BoundaryLocation< MeshType > >& boundaryLocation)
      : meshDistanceObject_(meshDistanceObject), boundaryLocation_(boundaryLocation)
   {}

   inline const mesh::BoundaryInfo& operator()(const Vector3< real_t >& p) const
   {
      typename MeshType::FaceHandle fh;
      meshDistanceObject_->sqSignedDistance(mesh::toOpenMesh(p), fh);
      return (*boundaryLocation_)[fh];
   }

   shared_ptr< MeshDistanceType > meshDistanceObject_;
   shared_ptr< mesh::BoundaryLocation< MeshType > > boundaryLocation_;
};

template< typename MeshDistanceType, typename MeshType >
inline BoundaryLocationFunction< MeshDistanceType, MeshType >
   makeBoundaryLocationFunction(const shared_ptr< MeshDistanceType >& meshDistanceObject,
                                const shared_ptr< mesh::BoundaryLocation< MeshType > >& boundaryLocation)
{
   return BoundaryLocationFunction< MeshDistanceType, MeshType >(meshDistanceObject, boundaryLocation);
}

} // namespace walberla