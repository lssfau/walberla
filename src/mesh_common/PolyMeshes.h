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
//! \file TriangleMeshes.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesh_common/TriangleMeshes.h"

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

namespace walberla {
namespace mesh {

using PythonPolyMesh = OpenMesh::PolyMesh_ArrayKernelT<OpenMesh::Python::MeshTraits>;
using PolyMesh = OpenMesh::PolyMesh_ArrayKernelT<RealTraits>;
using FloatPolyMesh = OpenMesh::PolyMesh_ArrayKernelT<FloatTraits>;

} // namespace mesh
} // namespace walberla