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
//! \file CommonDataSources.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "VTKMeshWriter.h"

#include "core/debug/CheckFunctions.h"
#include "core/mpi/MPIManager.h"

#include <OpenMesh/Core/Mesh/Status.hh>

#include <iostream>
#include <vector>
#include <string>

namespace walberla {
namespace mesh {

template< typename MeshType, typename OutputType = typename MeshType::Normal::value_type >
class NormalsVertexDataSource : public VTKMeshWriter<MeshType>::template VertexDataSource< OutputType >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Vertices Vertices;
   typedef typename VTKMeshWriter<MeshType>::template VertexDataSource< OutputType >::value_type value_type;

   NormalsVertexDataSource( const std::string & _name = "Normals" )
      : VTKMeshWriter<MeshType>::template VertexDataSource< OutputType >( _name ) { }
   virtual uint_t numComponents() { return uint_t(3); }
   virtual void   getData( const MeshType & mesh, const Vertices & vertices, std::vector<value_type> & data )
   {
      WALBERLA_CHECK( mesh.has_vertex_normals(), "You are trying to write vertex normals of a mesh which does not have any!" );

      data.reserve( size_t(3) * vertices.size() );

      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         const auto & n = mesh.normal( *it );

         data.push_back( numeric_cast<OutputType>( n[0] ) );
         data.push_back( numeric_cast<OutputType>( n[1] ) );
         data.push_back( numeric_cast<OutputType>( n[2] ) );
      }
   }
};



template< typename MeshType, typename OutputType = typename MeshType::Normal::value_type >
class NormalsFaceDataSource : public VTKMeshWriter<MeshType>::template FaceDataSource< OutputType >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Faces Faces;
   typedef typename VTKMeshWriter<MeshType>::template FaceDataSource< OutputType >::value_type value_type;

   NormalsFaceDataSource( const std::string & _name = "Normals" )
      : VTKMeshWriter<MeshType>::template FaceDataSource< OutputType >( _name ) { }
   virtual uint_t numComponents() { return uint_t(3); }
   virtual void   getData( const MeshType & mesh, const Faces & faces, std::vector<value_type> & data )
   {
      WALBERLA_CHECK( mesh.has_face_normals(), "You are trying to write face normals of a mesh which does not have any!" );

      data.reserve( size_t(3) * faces.size() );

      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         const auto & n = mesh.normal( *it );

         data.push_back( numeric_cast<OutputType>( n[0] ) );
         data.push_back( numeric_cast<OutputType>( n[1] ) );
         data.push_back( numeric_cast<OutputType>( n[2] ) );
      }
   }
};



template< typename MeshType, typename OutputType = typename MeshType::Normal::value_type >
class AreaFaceDataSource : public VTKMeshWriter<MeshType>::template FaceDataSource< OutputType >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Faces Faces;
   typedef typename VTKMeshWriter<MeshType>::template FaceDataSource< OutputType >::value_type value_type;

   AreaFaceDataSource( const std::string & _name = "Area" )
      : VTKMeshWriter<MeshType>::template FaceDataSource< OutputType >( _name ) { }
   virtual uint_t numComponents() { return uint_t(1); }
   virtual void   getData( const MeshType & mesh, const Faces & faces, std::vector<value_type> & data )
   {
      data.reserve( faces.size() );

      typename MeshType::Point v0, v1, v2;

      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         getVertexPositions( mesh, *it, v0, v1, v2 );
         data.push_back( ( ( v1 - v0 ) % ( v2 - v0 ) ).length() * typename MeshType::Normal::value_type( 0.5 ) );
      }
   }
};



template< typename MeshType >
class StatusBitFaceDataSource : public VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Faces Faces;
   typedef typename VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >::value_type value_type;

   StatusBitFaceDataSource( const OpenMesh::Attributes::StatusBits & bit, const std::string & _name )
      : VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >( _name ), bit_( bit ) {}
   virtual uint_t numComponents() { return uint_t(1); }
   virtual void   getData( const MeshType & mesh, const Faces & faces, std::vector<value_type> & data )
   {
      WALBERLA_CHECK( mesh.has_face_status(), "Cannot write face status bits, because the faces do not have them!" );

      data.reserve( faces.size() );
      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         data.push_back( mesh.status( *it ).is_bit_set( bit_ ) ? uint8_t(1) : uint8_t(0) );
      }
   }

private:
   OpenMesh::Attributes::StatusBits bit_;
};



template< typename MeshType >
class StatusBitVertexDataSource : public VTKMeshWriter<MeshType>::template VertexDataSource< uint8_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Vertices Vertices;
   typedef typename VTKMeshWriter<MeshType>::template VertexDataSource< uint8_t >::value_type value_type;

   StatusBitVertexDataSource( const OpenMesh::Attributes::StatusBits & bit, const std::string & _name )
      : VTKMeshWriter<MeshType>::template VertexDataSource< uint8_t >( _name ), bit_( bit ) {}
   virtual uint_t numComponents() { return uint_t(1); }
   virtual void   getData( const MeshType & mesh, const Vertices & vertices, std::vector<value_type> & data )
   {
      WALBERLA_CHECK( mesh.has_vertex_status(), "Cannot write vertex status bits, because the vertices do not have them!" );

      data.reserve( vertices.size() );
      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         data.push_back( mesh.status( *it ).is_bit_set( bit_ ) ? uint8_t(1) : uint8_t(0) );
      }
   }

private:
   OpenMesh::Attributes::StatusBits bit_;
};



template< typename MeshType >
class ColorFaceDataSource : public VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Faces Faces;
   typedef typename VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >::value_type value_type;

   ColorFaceDataSource( const std::string & _name = "color" )
      : VTKMeshWriter<MeshType>::template FaceDataSource< uint8_t >( _name ) {}

   virtual uint_t numComponents() { return uint_t(3); }
   virtual void   getData( const MeshType & mesh, const Faces & faces, std::vector<value_type> & data )
   {
      WALBERLA_CHECK( mesh.has_face_colors(), "Cannot write face colors, because the faces do not have them!" );

      data.reserve( faces.size() );
      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         const auto c = mesh.color( *it );
         data.push_back( c[0] );
         data.push_back( c[1] );
         data.push_back( c[2] );
      }
   }
};



template< typename MeshType >
class ColorVertexDataSource : public VTKMeshWriter<MeshType>::template VertexDataSource< uint8_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Vertices Vertices;
   typedef typename VTKMeshWriter<MeshType>::template VertexDataSource< uint8_t >::value_type value_type;

   ColorVertexDataSource( const std::string & _name = "color" )
      : VTKMeshWriter<MeshType>::template VertexDataSource< uint8_t >( _name ) {}

   virtual uint_t numComponents() { return uint_t(3); }
   virtual void   getData( const MeshType & mesh, const Vertices & vertices, std::vector<value_type> & data )
   {
      WALBERLA_CHECK( mesh.has_vertex_colors(), "Cannot write vertex colors, because the vertices do not have them!" );

      data.reserve( vertices.size() );
      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         const auto c = mesh.color( *it );
         data.push_back( c[0] );
         data.push_back( c[1] );
         data.push_back( c[2] );
      }
   }
};

template< typename MeshType >
class IndexFaceDataSource : public VTKMeshWriter<MeshType>::template FaceDataSource< int32_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Faces Faces;
   typedef typename VTKMeshWriter<MeshType>::template FaceDataSource< int32_t >::value_type value_type;

   IndexFaceDataSource( const std::string & _name = "index" )
      : VTKMeshWriter<MeshType>::template FaceDataSource< int32_t >( _name ) {}

   virtual uint_t numComponents() { return uint_t(1); }
   virtual void   getData( const MeshType & /*mesh*/, const Faces & faces, std::vector<value_type> & data )
   {
      data.reserve( faces.size() );
      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         data.push_back( it->idx() );
      }
   }
};


template< typename MeshType >
class IndexVertexDataSource : public VTKMeshWriter<MeshType>::template VertexDataSource< int32_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Vertices Vertices;
   typedef typename VTKMeshWriter<MeshType>::template VertexDataSource< int32_t >::value_type value_type;

   IndexVertexDataSource( const std::string & _name = "index" )
      : VTKMeshWriter<MeshType>::template VertexDataSource< int32_t >( _name ) {}

   virtual uint_t numComponents() { return uint_t(1); }
   virtual void   getData( const MeshType & /*mesh*/, const Vertices & vertices, std::vector<value_type> & data )
   {
      data.reserve( vertices.size() );
      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         data.push_back( it->idx() );
      }
   }
};


template< typename MeshType >
class RankFaceDataSource : public VTKMeshWriter<MeshType>::template FaceDataSource< int32_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Faces Faces;
   typedef typename VTKMeshWriter<MeshType>::template FaceDataSource< int32_t >::value_type value_type;

   RankFaceDataSource( const std::string & _name = "rank" )
      : VTKMeshWriter<MeshType>::template FaceDataSource< int32_t >( _name ) {}

   virtual uint_t numComponents() { return uint_t(1); }
   virtual void   getData( const MeshType & /*mesh*/, const Faces & faces, std::vector<value_type> & data )
   {
      int32_t rank = MPIManager::instance()->rank();
      data.assign( faces.size(), rank );
   }
};


template< typename MeshType >
class RankVertexDataSource : public VTKMeshWriter<MeshType>::template VertexDataSource< int32_t >
{
public:
   typedef typename VTKMeshWriter<MeshType>::Vertices Vertices;
   typedef typename VTKMeshWriter<MeshType>::template VertexDataSource< int32_t >::value_type value_type;

   RankVertexDataSource( const std::string & _name = "rank" )
      : VTKMeshWriter<MeshType>::template VertexDataSource< int32_t >( _name ) {}

   virtual uint_t numComponents() { return uint_t(1); }
   virtual void   getData( const MeshType & /*mesh*/, const Vertices & vertices, std::vector<value_type> & data )
   {
      int32_t rank = MPIManager::instance()->rank();
      data.assign( vertices.size(), rank );
   }
};


} // namespace mesh
} // namespace walberla