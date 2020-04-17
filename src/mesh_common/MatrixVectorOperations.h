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
//! \file MatrixVectorOperations.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Matrix3.h"
#include "core/math/Vector3.h"

#include <OpenMesh/Core/Geometry/VectorT.hh>

namespace walberla {
namespace mesh {

template< typename S >
inline OpenMesh::VectorT<S, 3> operator*( const math::Matrix3<S> & lhs, const OpenMesh::VectorT<S, 3> & rhs )
{
   return OpenMesh::VectorT<S, 3>( lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2],
                                   lhs[3]*rhs[0] + lhs[4]*rhs[1] + lhs[5]*rhs[2],
                                   lhs[6]*rhs[0] + lhs[7]*rhs[1] + lhs[8]*rhs[2] );
}


template< typename S >
inline math::Vector3<S> toWalberla( const OpenMesh::VectorT<S, 3> & v )
{
   return math::Vector3<S>( v[0], v[1], v[2] );
}


template< typename S >
inline OpenMesh::VectorT<S, 3> toOpenMesh( const math::Vector3<S> & v )
{
   return OpenMesh::VectorT<S, 3>( v[0], v[1], v[2] );
}

template< typename SWB, typename SOM >
inline math::Vector3<SWB> toWalberlaNumericCast( const OpenMesh::VectorT<SOM, 3> & v )
{
   return math::Vector3<SWB>( numeric_cast<SWB>( v[0] ), numeric_cast<SWB>( v[1] ), numeric_cast<SWB>( v[2] ) );
}


template< typename SOM, typename SWB >
inline OpenMesh::VectorT<SOM, 3> toOpenMeshNumericCast( const math::Vector3<SWB> & v )
{
   return OpenMesh::VectorT<SOM, 3>( numeric_cast<SOM>( v[0] ), numeric_cast<SOM>( v[1] ), numeric_cast<SOM>( v[2] ) );
}

} // namespace mesh
} // namespace walberla