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
//! \file Types.h
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Forwards declarations for mesh::pe bodies
//
//======================================================================================================================

#include <memory>

namespace walberla {
namespace mesh {
namespace pe {

class ConvexPolyhedron;

typedef ConvexPolyhedron        ConvexPolyhedronType;    //!< Type of the convex polyhedron geometric primitive.
typedef ConvexPolyhedron*       ConvexPolyhedronID;      //!< Handle for a convexpolyhedron primitive.
typedef const ConvexPolyhedron* ConstConvexPolyhedronID; //!< Handle for a constant convex polyhedron primitive.
using   ConvexPolyhedronPtr   = std::unique_ptr<ConvexPolyhedron>;

} // namespace pe
} // namespace mesh
} // namespace walberla
