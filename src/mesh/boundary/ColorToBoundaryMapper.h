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
//! \file ColorToBoundaryMapper.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "BoundaryInfo.h"

#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"

#include <map>

namespace walberla {
namespace mesh {

template< typename MeshType >
class ColorToBoundaryMapper
{
public:
   typedef typename MeshType::Color Color;

   ColorToBoundaryMapper( const BoundaryInfo & defaultBoundaryInfo ) : defaultBoundaryInfo_(defaultBoundaryInfo) { }

   void set( const Color & c, const BoundaryInfo & bi )
   {
      boundaryInfoMap_[c] = bi;
   }

   const BoundaryInfo & get( const Color & c ) const
   {
      auto it = boundaryInfoMap_.find(c);
      return (it == boundaryInfoMap_.end()) ? defaultBoundaryInfo_ : it->second;
   }

   shared_ptr< BoundaryLocation< MeshType > > addBoundaryInfoToMesh( MeshType & mesh ) const
   {
      WALBERLA_CHECK(mesh.has_face_colors());

      auto boundaryLocations = make_shared< BoundaryLocation< MeshType > >(mesh);

      for (auto & faceHandle : mesh.faces())
      {
         (*boundaryLocations)[faceHandle] = get(mesh.color(faceHandle));
      }

      return boundaryLocations;
   }

private:
   BoundaryInfo                    defaultBoundaryInfo_;
   std::map< Color, BoundaryInfo > boundaryInfoMap_;
};

} // namespace walberla
} // namespace mesh