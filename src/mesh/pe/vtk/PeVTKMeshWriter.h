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
//! \file PeVTKMeshWriter.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "domain_decomposition/BlockStorage.h"

#include "mesh_common/vtk/DistributedVTKMeshWriter.h"
#include "mesh_common/vtk/CommonDataSources.h"

#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/RigidBody.h"

#include <OpenMesh/Core/Utils/PropertyManager.hh>


namespace walberla {
namespace mesh {
namespace pe {

template< typename MeshType, typename Tesselation >
class PeVTKMeshWriter
{
public:
   static_assert( MeshType::IsPolyMesh == 1, "PeVTKMeshWriter only works with polygonal meshes!" );

   typedef typename MeshType::VertexHandle VertexHandle;
   typedef typename MeshType::FaceHandle   FaceHandle;
   typedef typename VTKMeshWriter<MeshType>::Faces    Faces;
   typedef typename VTKMeshWriter<MeshType>::Vertices Vertices;
   typedef typename OpenMesh::FPropHandleT< const walberla::pe::RigidBody * > BodyPointerFPropHandle;
   typedef typename OpenMesh::VPropHandleT< const walberla::pe::RigidBody * > BodyPointerVPropHandle;
   typedef OpenMesh::PropertyManager< BodyPointerFPropHandle, MeshType > BodyPointerFPropManager;
   typedef OpenMesh::PropertyManager< BodyPointerVPropHandle, MeshType > BodyPointerVPropManager;

   PeVTKMeshWriter( const shared_ptr<BlockStorage> & storage, const BlockDataID bodyStorageId, const Tesselation & tesselation,
                    const std::string & identifier, const uint_t writeFrequency, const std::string & baseFolder = "vtk_out" )
      : storage_(storage), bodyStorageId_(bodyStorageId), tesselation_(tesselation), 
        mesh_( make_shared<MeshType>() ), meshWriter_( mesh_, identifier, writeFrequency, baseFolder ),
        faceBodyPointer_( *mesh_, "bodyPointer" ), vertexBodyPointer_( *mesh_, "bodyPointer" )
   {
   }

   template< typename T >
   class DataSource
   {
   public:
      DataSource( const std::string & _name ) : name_( _name ) { }
      typedef T value_type;
      virtual uint_t      numComponents() = 0;
      const std::string & name() { return name_; }
   protected:
      std::string name_;
   };

   template< typename T >
   class VertexDataSource : public DataSource<T>
   {
   public:
      typedef typename DataSource<T>::value_type value_type;
      typedef typename PeVTKMeshWriter::Vertices Vertices;
      typedef typename PeVTKMeshWriter::BodyPointerVPropManager BodyPointerVPropManager;

      VertexDataSource( const std::string & _name ) : DataSource<T>( _name ) { }
      virtual void getData( const MeshType &, const Vertices &, std::vector<T> &, const BodyPointerVPropManager & ) = 0;
   };

   template< typename T >
   class FaceDataSource : public DataSource<T>
   {
   public:
      typedef typename DataSource<T>::value_type value_type;
      typedef typename PeVTKMeshWriter::Faces Faces;
      typedef typename PeVTKMeshWriter::BodyPointerFPropManager BodyPointerFPropManager;

      FaceDataSource( const std::string & _name ) : DataSource<T>( _name ) { }
      virtual void getData( const MeshType &, const Faces &, std::vector<T> &, const BodyPointerFPropManager & ) = 0;
   };

   void operator()()
   {
      if(meshWriter_.isWriteScheduled())
      {
         tesselateBodies();
      }
      meshWriter_();
   }

   template< typename T>
   inline void addDataSource( const shared_ptr< VertexDataSource<T> > & dataSource )
   {
      meshWriter_.addDataSource( make_shared<VertexDataSourceWrapper<T>>( dataSource, vertexBodyPointer_ ) ); 
   }

   template< typename T>
   inline void addDataSource( const shared_ptr< FaceDataSource<T> > & dataSource )
   {
      meshWriter_.addDataSource( make_shared<FaceDataSourceWrapper<T>>( dataSource, faceBodyPointer_ ) ); 
   }

   const BodyPointerFPropManager & getBodyPointerFPropManager() const { return faceBodyPointer_; }
   const BodyPointerVPropManager & getBodyPointerVPropManager() const { return vertexBodyPointer_; }

   void addVertexPropertyRank() { meshWriter_.addDataSource( make_shared<RankVertexDataSource<MeshType>>() ); }
   void addFacePropertyRank() { meshWriter_.addDataSource( make_shared<RankFaceDataSource<MeshType>>() ); }

   void setBodyFilter( const std::function<bool(const walberla::pe::RigidBody&)>& filter) { bodyFilter_ = filter; }
   const std::function<bool(const walberla::pe::RigidBody&)>& getBodyFilter() const { return bodyFilter_; }

protected:

   template< typename T >
   class VertexDataSourceWrapper : public DistributedVTKMeshWriter<MeshType>::template VertexDataSource<T>
   {
   public:
      typedef typename PeVTKMeshWriter::BodyPointerVPropManager BodyPointerVPropManager;
   
      VertexDataSourceWrapper( const shared_ptr<PeVTKMeshWriter::VertexDataSource<T>> & vertexDataSource, const BodyPointerVPropManager & bodyPointerProp )
         : DistributedVTKMeshWriter<MeshType>::template VertexDataSource<T>( vertexDataSource->name() ),
           vertexDataSource_(vertexDataSource), bodyPointerProp_(bodyPointerProp)
      { }
   
      virtual void getData( const MeshType & mesh, const Vertices & vertices, std::vector<T> & data )
      {
         return vertexDataSource_->getData( mesh, vertices, data, bodyPointerProp_ );
      };
   
      virtual uint_t numComponents() { return vertexDataSource_->numComponents(); }
   
   protected:
      shared_ptr<PeVTKMeshWriter::VertexDataSource<T>> vertexDataSource_;
      const BodyPointerVPropManager & bodyPointerProp_;
   };

   template< typename T >
   class FaceDataSourceWrapper : public DistributedVTKMeshWriter<MeshType>::template FaceDataSource<T>
   {
   public:
      typedef typename PeVTKMeshWriter::BodyPointerFPropManager BodyPointerFPropManager;

      FaceDataSourceWrapper( const shared_ptr<PeVTKMeshWriter::FaceDataSource<T>> & faceDataSource, const BodyPointerFPropManager & bodyPointerProp )
         : DistributedVTKMeshWriter<MeshType>::template FaceDataSource<T>( faceDataSource->name() ),
         faceDataSource_(faceDataSource), bodyPointerProp_(bodyPointerProp)
      { }

      virtual void getData( const MeshType & mesh, const Faces & faces, std::vector<T> & data )
      {
         return faceDataSource_->getData( mesh, faces, data, bodyPointerProp_ );
      };

      virtual uint_t numComponents() { return faceDataSource_->numComponents(); }

   protected:
      shared_ptr<PeVTKMeshWriter::FaceDataSource<T>> faceDataSource_;
      const BodyPointerFPropManager & bodyPointerProp_;
   };

   void tesselateBodies()
   {
      std::vector<VertexHandle> newVertices;
      std::vector<FaceHandle> newFaces;
      mesh_->clean();
      for(auto blockIt = storage_->begin(); blockIt != storage_->end(); ++blockIt)
      {
         walberla::pe::Storage * storage = blockIt->getData<walberla::pe::Storage>( bodyStorageId_ );
         const walberla::pe::BodyStorage & bodyStorage = (*storage)[0];
         for(auto bodyIt = bodyStorage.begin(); bodyIt != bodyStorage.end(); ++bodyIt)
         {
            if (!bodyFilter_(*bodyIt)) continue;
            newVertices.clear();
            newFaces.clear();
            tesselation_( *bodyIt, *mesh_, newVertices, newFaces );

            for( const auto & fh: newFaces )
               faceBodyPointer_[fh] = bodyIt.getBodyID();

            for( const auto & vh: newVertices )
               vertexBodyPointer_[vh] = bodyIt.getBodyID();

         }
      }
   }

   shared_ptr<BlockStorage> storage_;
   BlockDataID bodyStorageId_;
   Tesselation tesselation_;
   shared_ptr<MeshType> mesh_;
   DistributedVTKMeshWriter<MeshType> meshWriter_;
   BodyPointerFPropManager faceBodyPointer_;
   BodyPointerVPropManager vertexBodyPointer_;
   std::function<bool(const walberla::pe::RigidBody&)> bodyFilter_ = [](const walberla::pe::RigidBody&){ return true; };
};

} // namespace pe
} // namespace mesh
} // namespace walberla
