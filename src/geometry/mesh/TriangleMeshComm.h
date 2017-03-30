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
//! \file TriangleMeshComm.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Functions for Packing/Unpacking TriangleMesh objects into Send/RecvBuffer objects.
//
//======================================================================================================================

#pragma once

#include "TriangleMesh.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/SendBuffer.h"


namespace walberla {
namespace geometry {


//**********************************************************************************************************************
/*! \brief Packs a TriangleMesh into a SendBuffer.
 *
 * \param buf  The SendBuffer object mesh is packed into.
 * \param mesh The TriangleMesh object to be packed into buf.
 * \returns    A reference to buf.
 */
//*********************************************************************************************************************/
template< typename T, typename G > // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<(mpi::GenericSendBuffer<T,G> & buf, const TriangleMesh & mesh )
{
   return buf << mesh.getVertexIndices() << mesh.getNormalIndices() << mesh.getVertices()
              << mesh.getVertexNormals() << mesh.getVertexColors();
}

//**********************************************************************************************************************
/*! \brief Unpacks a TriangleMesh from a RecvBuffer.
 *
 * \param buf  The RecvBuffer object mesh is unpacked from.
 * \param mesh The TriangleMesh object buf is unpacked into.
 * \returns    A reference to buf.
 */
//*********************************************************************************************************************/
template< typename T >
mpi::GenericRecvBuffer<T>& operator>>(mpi::GenericRecvBuffer<T> & buf, TriangleMesh & mesh )
{
   return buf >> mesh.getVertexIndices() >> mesh.getNormalIndices() >> mesh.getVertices()
              >> mesh.getVertexNormals() >> mesh.getVertexColors();
}

} // namespace geometry
} // namespace walberla
