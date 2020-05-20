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

#include "PeVTKMeshWriter.h"

#include "core/debug/CheckFunctions.h"
#include "core/mpi/MPIManager.h"

#include "mesh_common/MatrixVectorOperations.h"

#include <OpenMesh/Core/Mesh/Status.hh>

#include <iostream>
#include <vector>
#include <string>

namespace walberla {
namespace mesh {
namespace pe {

template< typename MeshType, typename Tesselation, typename OutputType = uint64_t >
class SIDVertexDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType > Base;
   typedef typename Base::Vertices Vertices;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerVPropManager BodyPointerVPropManager;

   SIDVertexDataSource( const std::string & _name = "sid" )
      : Base( _name ) { }

   virtual uint_t numComponents() { return uint_t(1); }

   virtual void getData( const MeshType & /*mesh*/, const Vertices & vertices, std::vector<value_type> & data, const BodyPointerVPropManager & bodyPointer )
   {
      data.reserve( vertices.size() );

      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->getSystemID();
         data.push_back( numeric_cast<OutputType>( v ) );
      }
   }
};


template< typename MeshType, typename Tesselation, typename OutputType = uint64_t >
class SIDFaceDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template FaceDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template FaceDataSource< OutputType > Base;
   typedef typename Base::Faces Faces;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerFPropManager BodyPointerFPropManager;

   SIDFaceDataSource( const std::string & _name = "sid" )
      : Base( _name ) { }

   virtual uint_t numComponents() { return uint_t(1); }

   virtual void getData( const MeshType & /*mesh*/, const Faces & faces, std::vector<value_type> & data, const BodyPointerFPropManager & bodyPointer )
   {
      data.reserve( faces.size() );

      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->getSystemID();

         data.push_back( numeric_cast<OutputType>( v ) );
      }
   }
};

template< typename MeshType, typename Tesselation, typename OutputType = uint64_t >
class UIDVertexDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType > Base;
   typedef typename Base::Vertices Vertices;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerVPropManager BodyPointerVPropManager;

   UIDVertexDataSource( const std::string & _name = "uid" )
      : Base( _name ) { }

   virtual uint_t numComponents() { return uint_t(1); }

   virtual void getData( const MeshType & /*mesh*/, const Vertices & vertices, std::vector<value_type> & data, const BodyPointerVPropManager & bodyPointer )
   {
      data.reserve( vertices.size() );

      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->getID();
         data.push_back( numeric_cast<OutputType>( v ) );
      }
   }
};


template< typename MeshType, typename Tesselation, typename OutputType = uint64_t >
class UIDFaceDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template FaceDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template FaceDataSource< OutputType > Base;
   typedef typename Base::Faces Faces;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerFPropManager BodyPointerFPropManager;

   UIDFaceDataSource( const std::string & _name = "uid" )
      : Base( _name ) { }

   virtual uint_t numComponents() { return uint_t(1); }

   virtual void getData( const MeshType & /*mesh*/, const Faces & faces, std::vector<value_type> & data, const BodyPointerFPropManager & bodyPointer )
   {
      data.reserve( faces.size() );

      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->getID();

         data.push_back( numeric_cast<OutputType>( v ) );
      }
   }
};


template< typename MeshType, typename Tesselation, typename OutputType = real_t >
class LinearVelocityVertexDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType > Base;
   typedef typename Base::Vertices Vertices;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerVPropManager BodyPointerVPropManager;

   LinearVelocityVertexDataSource( const std::string & _name = "linearVelocity" )
      : Base( _name ) { }

   virtual uint_t numComponents() { return uint_t(3); }

   virtual void getData( const MeshType & /*mesh*/, const Vertices & vertices, std::vector<value_type> & data, const BodyPointerVPropManager & bodyPointer )
   {
      data.reserve( size_t(3) * vertices.size() );

      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->getLinearVel();

         data.push_back( numeric_cast<OutputType>( v[0] ) );
         data.push_back( numeric_cast<OutputType>( v[1] ) );
         data.push_back( numeric_cast<OutputType>( v[2] ) );
      }
   }
};


template< typename MeshType, typename Tesselation, typename OutputType = real_t >
class LinearVelocityFaceDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template FaceDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template FaceDataSource< OutputType > Base;
   typedef typename Base::Faces Faces;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerFPropManager BodyPointerFPropManager;

   LinearVelocityFaceDataSource( const std::string & _name = "linearVelocity" )
      : Base( _name ) { }

   virtual ~LinearVelocityFaceDataSource() = default;

   virtual uint_t numComponents() { return uint_t(3); }

   virtual void getData( const MeshType & /*mesh*/, const Faces & faces, std::vector<value_type> & data, const BodyPointerFPropManager & bodyPointer )
   {
      data.reserve( size_t(3) * faces.size() );

      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->getLinearVel();

         data.push_back( numeric_cast<OutputType>( v[0] ) );
         data.push_back( numeric_cast<OutputType>( v[1] ) );
         data.push_back( numeric_cast<OutputType>( v[2] ) );
      }
   }
};


template< typename MeshType, typename Tesselation, typename OutputType = real_t >
class AngularVelocityVertexDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType > Base;
   typedef typename Base::Vertices Vertices;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerVPropManager BodyPointerVPropManager;

   AngularVelocityVertexDataSource( const std::string & _name = "angularVelocity" )
      : Base( _name ) { }

   virtual uint_t numComponents() { return uint_t(3); }

   virtual void getData( const MeshType & /*mesh*/, const Vertices & vertices, std::vector<value_type> & data, const BodyPointerVPropManager & bodyPointer )
   {
      data.reserve( size_t(3) * vertices.size() );

      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->getAngularVel();

         data.push_back( numeric_cast<OutputType>( v[0] ) );
         data.push_back( numeric_cast<OutputType>( v[1] ) );
         data.push_back( numeric_cast<OutputType>( v[2] ) );
      }
   }
};


template< typename MeshType, typename Tesselation, typename OutputType = real_t >
class AngularVelocityFaceDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template FaceDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template FaceDataSource< OutputType > Base;
   typedef typename Base::Faces Faces;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerFPropManager BodyPointerFPropManager;

   AngularVelocityFaceDataSource( const std::string & _name = "angularVelocity" )
      : Base( _name ) { }

   virtual ~AngularVelocityFaceDataSource() = default;

   virtual uint_t numComponents() { return uint_t(3); }

   virtual void getData( const MeshType & /*mesh*/, const Faces & faces, std::vector<value_type> & data, const BodyPointerFPropManager & bodyPointer )
   {
      data.reserve( size_t(3) * faces.size() );

      for( auto it = faces.begin(); it != faces.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->getAngularVel();

         data.push_back( numeric_cast<OutputType>( v[0] ) );
         data.push_back( numeric_cast<OutputType>( v[1] ) );
         data.push_back( numeric_cast<OutputType>( v[2] ) );
      }
   }
};


template< typename MeshType, typename Tesselation, typename OutputType = real_t >
class SurfaceVelocityVertexDataSource : public PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType >
{
public:
   typedef typename PeVTKMeshWriter<MeshType, Tesselation>::template VertexDataSource< OutputType > Base;
   typedef typename Base::Vertices Vertices;
   typedef typename Base::value_type value_type;
   typedef typename Base::BodyPointerVPropManager BodyPointerVPropManager;

   SurfaceVelocityVertexDataSource( const std::string & _name = "surfaceVelocity" )
      : Base( _name ) { }

   virtual ~SurfaceVelocityVertexDataSource() = default;

   virtual uint_t numComponents() { return uint_t(3); }

   virtual void getData( const MeshType & mesh, const Vertices & vertices, std::vector<value_type> & data, const BodyPointerVPropManager & bodyPointer )
   {
      data.reserve( size_t(3) * vertices.size() );

      for( auto it = vertices.begin(); it != vertices.end(); ++it )
      {
         const auto & v = bodyPointer[*it]->velFromWF( toWalberlaNumericCast<real_t>( mesh.point( *it ) ) );

         data.push_back( numeric_cast<OutputType>( v[0] ) );
         data.push_back( numeric_cast<OutputType>( v[1] ) );
         data.push_back( numeric_cast<OutputType>( v[2] ) );
      }
   }
};


} // namespace pe
} // namespace mesh
} // namespace walberla
