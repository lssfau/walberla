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
//! \file BoundarySetup.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "BoundaryLocation.h"

#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"

#include "mesh_common/vtk/VTKMeshWriter.h"

namespace walberla {
namespace mesh {

template< typename MeshType >
class BoundaryUIDFaceDataSource : public VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Faces Faces;
   typedef typename VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >::value_type value_type;

   BoundaryUIDFaceDataSource( const shared_ptr< BoundaryLocation< MeshType > > & boundaryLocation, const std::string & _name = "BoundaryUID" )
      : VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >( _name ), boundaryLocation_( boundaryLocation )
   {
   }

   virtual uint_t numComponents() { return uint_t(1); }
   virtual void   getData( const MeshType & /*mesh*/, const Faces & faces, std::vector<value_type> & data )
   {
      data.reserve( faces.size() );

      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         auto idx = (*boundaryLocation_)[ *it ].getUid().getUid();

         WALBERLA_CHECK_LESS_EQUAL( idx, std::numeric_limits<uint8_t>::max(), "You have to many BoundaryUIDs for BoundaryUIDFaceDataSource" );

         data.push_back( uint8_c( idx ) );
      }
   }

private:
   shared_ptr< BoundaryLocation< MeshType > > boundaryLocation_;
};

} // namespace mesh
} // namespace walberla