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
//! \file CommonFilters.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "VTKMeshWriter.h"

#include "core/math/AABB.h"

#include <OpenMesh/Core/Mesh/Status.hh>

#include <algorithm>
#include <set>
#include <vector>

namespace walberla {
namespace mesh {


template< typename MeshType >
struct SubsetFaceFilter
{
   template<typename InputIterator>
   SubsetFaceFilter( InputIterator facesBegin, InputIterator facesEnd ) : includedFaces_( facesBegin, facesEnd )
   {
      std::sort( includedFaces_.begin(), includedFaces_.end() );
      includedFaces_.erase( std::unique( includedFaces_.begin(), includedFaces_.end() ), includedFaces_.end() );
   }

   bool operator()( const MeshType & /*mesh*/, const typename MeshType::FaceHandle & fh ) { return std::binary_search( includedFaces_.begin(), includedFaces_.end(), fh ); }

private:
   std::vector< typename MeshType::FaceHandle > includedFaces_;
};


template< typename MeshType >
struct AABBFaceFilter
{
   AABBFaceFilter( const AABB & aabb ) : aabb_( aabb ) { }

   bool operator()( const MeshType & mesh, const typename MeshType::FaceHandle & fh )
   {
      for( auto it = mesh.cfv_begin( fh ); it != mesh.cfv_end( fh ); ++it )
         if( !aabb_.contains( toWalberla( mesh.point(*it) ) ) )
            return false;

      return true;
   }

private:
   AABB aabb_;
};


template< typename MeshType >
struct StatusFaceFilter
{
   StatusFaceFilter( const OpenMesh::Attributes::StatusBits & bit ) : bit_( bit ) { }

   bool operator()( const MeshType & mesh, const typename MeshType::FaceHandle & fh )
   {
      WALBERLA_ASSERT( mesh.has_face_status(), "The mesh's faces cannot be filtered by status bits, because the faces do not have them!" );

      return mesh.status( fh ).is_bit_set( bit_ );
   }

private:
   OpenMesh::Attributes::StatusBits bit_;
};



} // namespace mesh
} // namespace walberla