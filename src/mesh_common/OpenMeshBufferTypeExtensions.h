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
//! \file OpenMeshBufferTypeExtensions.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================
#pragma once

#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh/Core/Geometry/Vector11T.hh>

namespace walberla {
namespace mpi {


// Handles

template< typename T,  // Element type of SendBuffer
          typename G > // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<(mpi::GenericSendBuffer<T,G> & buf, const OpenMesh::BaseHandle & handle )
{
   buf.addDebugMarker( "oh" );
   buf << handle.idx();
   return buf;
}

template< typename T,        // Element type  of RecvBuffer
          typename HandleT > // OpenMehs handle type 
mpi::GenericRecvBuffer<T>& unpackOpenMeshHandle(mpi::GenericRecvBuffer<T> & buf, HandleT & handle )
{
   buf.readDebugMarker( "oh" );
   int idx;
   buf >> idx;
   handle = HandleT( idx );
   return buf;
}

template< typename T> // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>(mpi::GenericRecvBuffer<T> & buf, OpenMesh::VertexHandle & handle )
{
   return unpackOpenMeshHandle( buf, handle );
}

template< typename T> // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>(mpi::GenericRecvBuffer<T> & buf, OpenMesh::FaceHandle & handle )
{
   return unpackOpenMeshHandle( buf, handle );
}

template< typename T> // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>(mpi::GenericRecvBuffer<T> & buf, OpenMesh::HalfedgeHandle & handle )
{
   return unpackOpenMeshHandle( buf, handle );
}

template< typename T> // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>(mpi::GenericRecvBuffer<T> & buf, OpenMesh::EdgeHandle & handle )
{
   return unpackOpenMeshHandle( buf, handle );
}

template<>
struct BufferSizeTrait< OpenMesh::VertexHandle > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait< int >::size + BUFFER_DEBUG_OVERHEAD;
};

template<>
struct BufferSizeTrait< OpenMesh::FaceHandle > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait< int >::size + BUFFER_DEBUG_OVERHEAD;
};

template<>
struct BufferSizeTrait< OpenMesh::HalfedgeHandle > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait< int >::size + BUFFER_DEBUG_OVERHEAD;
};

template<>
struct BufferSizeTrait< OpenMesh::EdgeHandle > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait< int >::size + BUFFER_DEBUG_OVERHEAD;
};

// Vectors

template< typename T,      // Element type of SendBuffer
          typename G,      // Growth policy of SendBuffer
          typename Scalar, // Vector scalar type
          int DIM >        // Vector dimension
mpi::GenericSendBuffer<T,G>& operator<<(mpi::GenericSendBuffer<T,G> & buf, const OpenMesh::VectorT<Scalar, DIM> & v )
{
   buf.addDebugMarker( "ov" );
   for( size_t i = 0; i < v.size(); ++i)
      buf << v[i];
   return buf;
}

template< typename T,      // Element type  of RecvBuffer
          typename Scalar, // Vector scalar type
          int DIM >        // Vector dimension
mpi::GenericRecvBuffer<T>& operator>>(mpi::GenericRecvBuffer<T> & buf, OpenMesh::VectorT<Scalar, DIM> & v )
{
   buf.readDebugMarker( "ov" );
   for( size_t i = 0; i < v.size(); ++i)
      buf >> v[i];
   return buf;
}

template< typename Scalar, // Vector scalar type
          int DIM >        // Vector dimension
struct BufferSizeTrait< OpenMesh::VectorT<Scalar, DIM> > {
   static const bool constantSize = true;
   static const uint_t size = DIM * BufferSizeTrait< Scalar >::size + BUFFER_DEBUG_OVERHEAD;
};

}

} // namespace walberla