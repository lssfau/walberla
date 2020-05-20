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
//! \file DistanceFunction.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/math/Vector3.h"
#include "mesh_common/MatrixVectorOperations.h"

namespace walberla {

template<typename MeshDistanceType>
struct MeshDistanceFunction
{
   MeshDistanceFunction( const shared_ptr <MeshDistanceType> &meshDistanceObject ) : meshDistanceObject_(
           meshDistanceObject ) {}

   inline real_t operator()( const Vector3 <real_t> &p ) const
   {
      return real_c( meshDistanceObject_->sqSignedDistance( mesh::toOpenMesh( p )));
   }

   shared_ptr <MeshDistanceType> meshDistanceObject_;
};

template<typename MeshDistanceType>
inline MeshDistanceFunction<MeshDistanceType>
makeMeshDistanceFunction( const shared_ptr <MeshDistanceType> & meshDistanceObject )
{
   return MeshDistanceFunction<MeshDistanceType>( meshDistanceObject );
}

} // namespace walberla
