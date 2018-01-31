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
//! \file MeshIO.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/mpi/Broadcast.h"
#include "core/mpi/BufferDataTypeExtensions.h"

#include <fstream>
#include <string>

#include "core/Filesystem.h"

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning( disable : 4456 )
#endif //_MSC_VER
#include <OpenMesh/Core/IO/MeshIO.hh>
#ifdef _MSC_VER
#  pragma warning(pop)
#endif //_MSC_VER

namespace walberla {
namespace mesh {


/**
* \brief Loads an OpenMesh in parallel
*
* The mesh is read from disk by a single process and then broadcasted. This ensures a scalable load process that does
* not put to much pressure on the file system at large scale.
* 
* \tparam MeshType The type of the OpenMesh
*
* \param filename filename of the mesh to be loaded
* \param mesh     The mesh data structure to be written to
*/
template< typename MeshType >
void readAndBroadcast( const std::string & filename, MeshType & mesh )
{
   if( !filesystem::exists( filename ) )
      WALBERLA_ABORT( "The mesh file \"" << filename << "\" does not exist!" );

   std::string extension = filesystem::path( filename ).extension().string();

   std::string str;

   WALBERLA_ROOT_SECTION()
   {
      std::ifstream t( filename.c_str() );
      if( !t )
         WALBERLA_ABORT( "Error while reading file \"" << filename << "\"!" );

      t.seekg( 0, std::ios::end );
      str.reserve( static_cast<std::string::size_type>( t.tellg() ) );
      t.seekg( 0, std::ios::beg );

      str.assign( ( std::istreambuf_iterator<char>(t) ),
                    std::istreambuf_iterator<char>() );
   }

   mpi::broadcastObject( str );

   std::istringstream iss( str );

   OpenMesh::IO::Options options;
   if( mesh.has_face_colors() )
      options += OpenMesh::IO::Options::FaceColor;
   if( mesh.has_vertex_colors() )
      options += OpenMesh::IO::Options::VertexColor;
   if( !OpenMesh::IO::read_mesh( mesh, iss, extension, options ) )
      WALBERLA_ABORT( "Error while reading file \"" << filename << "\"!" );
}

} // namespace mesh
} // namespace walberla