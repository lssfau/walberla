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
//! \file BoundaryInfo.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "BoundaryInfo.h"

#include <OpenMesh/Core/Utils/PropertyManager.hh>

namespace walberla {
namespace mesh {

template< typename MeshType >
class BoundaryLocation
{
public:
   typedef typename OpenMesh::FPropHandleT< BoundaryInfo > BoundaryInformation;
   typedef typename MeshType::FaceHandle FaceHandle;

   BoundaryLocation( MeshType & mesh ) : boundaryInfos_( mesh, "BoundaryInformation" ) {}

   const BoundaryInfo & operator[]( const FaceHandle & fh ) const { return boundaryInfos_[ fh ]; }
         BoundaryInfo & operator[]( const FaceHandle & fh )       { return boundaryInfos_[ fh ]; }

private:
   OpenMesh::PropertyManager< BoundaryInformation, MeshType > boundaryInfos_;
};

} // namespace mesh
} // namespace walberla