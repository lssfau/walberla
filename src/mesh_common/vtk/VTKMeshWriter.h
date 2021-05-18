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
//! \file VTKMeshWriter.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesh_common/MeshOperations.h"

#include "core/math/Vector3.h"
#include "core/mpi/MPIManager.h"

#include "vtk/Base64Writer.h"
#include "vtk/UtilityFunctions.h"

#include "core/Filesystem.h"

#include <iostream>
#include <vector>
#include <set>
#include <string>

namespace walberla {
namespace mesh {

template< typename MeshType >
class VTKMeshWriter {
public:
   typedef std::function<bool ( const MeshType &, const typename MeshType::FaceHandle & )> FaceFilterFunction;

   typedef std::vector< typename MeshType::FaceHandle >   Faces;
   typedef std::vector< typename MeshType::VertexHandle > Vertices;

   template< typename T >
   class DataSource
   {
   public:
      DataSource( const std::string & _name ) : name_( _name ) { }
      virtual ~DataSource() = default;
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
      typedef typename VTKMeshWriter::Vertices Vertices;

      VertexDataSource( const std::string & _name ) : DataSource<T>( _name ) { }
      virtual ~VertexDataSource() = default;
      virtual void getData( const MeshType &, const Vertices & vertices, std::vector<T> & ) = 0;
   };

   template< typename T >
   class FaceDataSource : public DataSource<T>
   {
   public:
      typedef typename DataSource<T>::value_type value_type;
      typedef typename VTKMeshWriter::Faces Faces;

      FaceDataSource( const std::string & _name ) : DataSource<T>( _name ) { }
      virtual ~FaceDataSource() = default;
      virtual void getData( const MeshType &, const Faces & faces, std::vector<T> & ) = 0;
   };

   VTKMeshWriter( const shared_ptr<const MeshType> & mesh, const std::string & identifier,
                  const uint_t writeFrequency, const std::string & baseFolder = "vtk_out" );

   void operator()();

   inline void addDataSource( const shared_ptr< VertexDataSource<float>   > & dataSource ) { floatVertexDataSources_.push_back( dataSource );  }
   inline void addDataSource( const shared_ptr< VertexDataSource<double>  > & dataSource ) { doubleVertexDataSources_.push_back( dataSource ); }
   inline void addDataSource( const shared_ptr< VertexDataSource<int32_t> > & dataSource ) { int32VertexDataSources_.push_back( dataSource );  }
   inline void addDataSource( const shared_ptr< VertexDataSource<uint8_t> > & dataSource ) { uint8VertexDataSources_.push_back( dataSource );  }
   inline void addDataSource( const shared_ptr< VertexDataSource<uint64_t> > & dataSource ) { uint64VertexDataSources_.push_back( dataSource );  }

   inline void addDataSource( const shared_ptr< FaceDataSource<float>   > & dataSource ) { floatFaceDataSources_.push_back( dataSource );  }
   inline void addDataSource( const shared_ptr< FaceDataSource<double>  > & dataSource ) { doubleFaceDataSources_.push_back( dataSource ); }
   inline void addDataSource( const shared_ptr< FaceDataSource<int32_t> > & dataSource ) { int32FaceDataSources_.push_back( dataSource );  }
   inline void addDataSource( const shared_ptr< FaceDataSource<uint8_t> > & dataSource ) { uint8FaceDataSources_.push_back( dataSource );  }
   inline void addDataSource( const shared_ptr< FaceDataSource<uint64_t> > & dataSource ) { uint64FaceDataSources_.push_back( dataSource );  }

   inline void setFaceFilter  ( const FaceFilterFunction & f ) { faceFilter_ = f;                    }
   inline void clearFaceFilter()                               { faceFilter_ = FaceFilterFunction(); }
   inline bool isFaceFilterSet()                         const { return static_cast<bool>(faceFilter_);        }

   void incrementTimeStep()      { ++timestep_; }
   bool isWriteScheduled() const { return writeFrequency_ > 0 && timestep_ % writeFrequency_ == 0; }

protected:
   void write( std::ostream & os ) const;
   void writePrefix( std::ostream & os ) const;
   void writePostfix( std::ostream & os ) const;
   void writePiece( std::ostream & os ) const;

   template<typename T>
   inline void writeVertexData( const T & vertexDataSources, const Vertices & vertices, std::ostream & os, vtk::Base64Writer & b64 ) const;

   template<typename T>
   inline void writeFaceData( const T & faceDataSources, const Faces & faces, std::ostream & os, vtk::Base64Writer & b64 ) const;

   void writePVD( std::ostream & os ) const;

   shared_ptr<const MeshType> mesh_;

   uint_t      writeFrequency_;
   std::string identifier_;
   std::string baseFolder_;

   uint_t      timestep_;

   FaceFilterFunction faceFilter_;

   std::vector< shared_ptr< VertexDataSource<float  > > >  floatVertexDataSources_;
   std::vector< shared_ptr< VertexDataSource<double > > >  doubleVertexDataSources_;
   std::vector< shared_ptr< VertexDataSource<int32_t> > >  int32VertexDataSources_;
   std::vector< shared_ptr< VertexDataSource<uint8_t> > >  uint8VertexDataSources_;
   std::vector< shared_ptr< VertexDataSource<uint64_t> > > uint64VertexDataSources_;

   std::vector< shared_ptr< FaceDataSource<float  > > >  floatFaceDataSources_;
   std::vector< shared_ptr< FaceDataSource<double > > >  doubleFaceDataSources_;
   std::vector< shared_ptr< FaceDataSource<int32_t> > >  int32FaceDataSources_;
   std::vector< shared_ptr< FaceDataSource<uint8_t> > >  uint8FaceDataSources_;
   std::vector< shared_ptr< FaceDataSource<uint64_t> > > uint64FaceDataSources_;
};


template< typename MeshType >
template< typename T >
void VTKMeshWriter<MeshType>::writeVertexData( const T & vertexDataSources, const Vertices & vertices, std::ostream & os, vtk::Base64Writer & b64 ) const
{
   typedef typename T::value_type::element_type::value_type value_type;

   for( auto dsIt = vertexDataSources.begin(); dsIt != vertexDataSources.end(); ++dsIt )
   {
      std::vector< value_type > data;
      (*dsIt)->getData( *mesh_, vertices, data );

      WALBERLA_CHECK_EQUAL( vertices.size() * (*dsIt)->numComponents(), data.size(), "The vertex data source \"" << (*dsIt)->name() << "\" wrote the wrong amount of data!" );

      os << "        <DataArray type=\"" << vtk::typeToString<value_type>() << "\" Name=\"" << (*dsIt)->name() << "\" NumberOfComponents=\"" << (*dsIt)->numComponents() << "\" format=\"binary\">\n";
      os << "          ";
      for( auto it = data.begin(); it != data.end(); ++it )
         b64 << *it;
      b64.toStream( os );
      os << "        </DataArray>\n";
   }
}


template< typename MeshType >
template< typename T >
void VTKMeshWriter<MeshType>::writeFaceData( const T & faceDataSources, const Faces & faces, std::ostream & os, vtk::Base64Writer & b64 ) const
{
   typedef typename T::value_type::element_type::value_type value_type;

   for( auto dsIt = faceDataSources.begin(); dsIt != faceDataSources.end(); ++dsIt )
   {
      std::vector< value_type > data;
      (*dsIt)->getData( *mesh_, faces, data );

      WALBERLA_CHECK_EQUAL( faces.size() * (*dsIt)->numComponents(), data.size(), "The face data source \"" << (*dsIt)->name() << "\" wrote the wrong amount of data!" );

      os << "        <DataArray type=\"" << vtk::typeToString<value_type>() << "\" Name=\"" << (*dsIt)->name() << "\" NumberOfComponents=\"" << (*dsIt)->numComponents() << "\" format=\"binary\">\n";
      os << "          ";
      for( auto it = data.begin(); it != data.end(); ++it )
         b64 << *it;
      b64.toStream( os );
      os << "        </DataArray>\n";
   }
}


template< typename MeshType >
VTKMeshWriter<MeshType>::VTKMeshWriter( const shared_ptr<const MeshType> & mesh, const std::string & identifier, const uint_t writeFrequency,
   const std::string & baseFolder )
   : mesh_( mesh ), writeFrequency_( writeFrequency ), identifier_( identifier ), baseFolder_( baseFolder ), timestep_( 0 )
{
   WALBERLA_ROOT_SECTION()
   {
      std::ostringstream folder;
      folder << baseFolder_ << '/' << identifier_;
      if( filesystem::exists( folder.str() ) )
         filesystem::remove_all( folder.str() );

      std::ostringstream pvdFile;
      pvdFile << baseFolder_ << '/' << identifier_ << ".pvd";
      if( filesystem::exists( pvdFile.str() ) )
         filesystem::remove( pvdFile.str() );

      filesystem::create_directories( folder.str() );
   }
   WALBERLA_MPI_BARRIER();
}


template< typename MeshType >
void VTKMeshWriter<MeshType>::writePrefix( std::ostream & os ) const
{
   os << "<VTKFile type=\"PolyData\" version=\"0.1\">\n"
      << "  <PolyData>\n";
}

template< typename MeshType >
void VTKMeshWriter<MeshType>::writePostfix( std::ostream & os ) const
{
   os << "  </PolyData>\n"
      << "</VTKFile>\n";
}

template< typename MeshType >
void VTKMeshWriter<MeshType>::writePiece( std::ostream & os ) const
{
   vtk::Base64Writer b64;

   Faces faces;
   Vertices vertices;

   if( isFaceFilterSet() )
   {
      for( auto it = mesh_->faces_sbegin(); it != mesh_->faces_end(); ++it )
         if( faceFilter_( *mesh_, *it ) )
            faces.push_back( *it );

      std::vector< typename MeshType::VertexHandle > tmpVertices = findConnectedVertices( *mesh_, faces.begin(), faces.end() );
      vertices.assign( tmpVertices.begin(), tmpVertices.end() );
   }
   else
   {
      std::copy( mesh_->faces_sbegin()   , mesh_->faces_end()   , std::back_inserter( faces )    );
      std::copy( mesh_->vertices_sbegin(), mesh_->vertices_end(), std::back_inserter( vertices ) );
   }

   const uint_t numVertices  = vertices.size();
   const uint_t numFaces     = faces.size();

   os << "    <Piece NumberOfPoints=\"" << numVertices << "\" NumberOfPolys=\"" << numFaces << "\">\n"
      << "      <Points>\n"
      << "        <DataArray type=\"" << vtk::typeToString<typename MeshType::Point::value_type>() << "\" NumberOfComponents=\"3\" format=\"binary\">\n";

   os << "          ";
   std::map<typename MeshType::VertexHandle, int > vertexIndizes;
   int32_t vertexCounter = 0;
   for( auto it = vertices.begin(); it != vertices.end(); ++it )
   {
      const auto & p = mesh_->point( *it );
      b64 << p[0] << p[1] << p[2];
      vertexIndizes[ *it ] = vertexCounter++;
   }
   b64.toStream( os );

   os << "        </DataArray>\n"
      << "      </Points>\n";

   os << "      <Polys>\n";

   os << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n";
   os << "          ";
   std::vector<int32_t> offsets;
   offsets.reserve( numFaces );
   int32_t currentOffset = 0;
   for( auto it = faces.begin(); it != faces.end(); ++it )
   {
      int32_t nv = 0;
      for( auto v_it = mesh_->cfv_ccwbegin( *it ); v_it != mesh_->cfv_ccwend( *it ); ++v_it, ++nv )
      {
         b64 << vertexIndizes[ *v_it ];
      }

      currentOffset += nv;
      offsets.push_back( currentOffset );

   }
   b64.toStream( os );

   os << "        </DataArray>\n";

   os << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n";
   os << "          ";
   for( auto it = offsets.begin(); it != offsets.end(); ++it )
   {
      b64 << *it;
   }
   b64.toStream( os );

   os << "        </DataArray>\n";

   os << "      </Polys>\n";

   os << "      <CellData>\n";

   writeFaceData( floatFaceDataSources_ ,  faces, os, b64 );
   writeFaceData( doubleFaceDataSources_,  faces, os, b64 );
   writeFaceData( int32FaceDataSources_ ,  faces, os, b64 );
   writeFaceData( uint8FaceDataSources_ ,  faces, os, b64 );
   writeFaceData( uint64FaceDataSources_ , faces, os, b64 );

   os << "      </CellData>\n";

   os << "      <PointData>\n";

   writeVertexData( floatVertexDataSources_ ,  vertices, os, b64 );
   writeVertexData( doubleVertexDataSources_,  vertices, os, b64 );
   writeVertexData( int32VertexDataSources_ ,  vertices, os, b64 );
   writeVertexData( uint8VertexDataSources_ ,  vertices, os, b64 );
   writeVertexData( uint64VertexDataSources_ , vertices, os, b64 );

   os << "      </PointData>\n";


   os << "    </Piece>\n";
}

template< typename MeshType >
void VTKMeshWriter<MeshType>::write( std::ostream & os ) const
{
   writePrefix( os );
   writePiece( os );
   writePostfix( os );
}


template< typename MeshType >
inline void VTKMeshWriter<MeshType>::writePVD( std::ostream & os ) const
{
   os << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
      << "  <Collection>\n";

   for( uint_t timestep = 0; timestep <= timestep_; timestep += writeFrequency_ )
   {
      os << "    <DataSet timestep=\"" << timestep << "\" "
         "file=\"" << identifier_ << '/' << identifier_ << '_' << timestep << ".vtp" << "\"/>\n";
   }

   os << "  </Collection>\n"
      << "</VTKFile>\n";
}


template< typename MeshType >
void VTKMeshWriter<MeshType>::operator()()
{
   if( isWriteScheduled() )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ostringstream filePathVtp;
         filePathVtp << baseFolder_ << '/' << identifier_ << '/' << identifier_ << '_' << timestep_ << ".vtp";
      
         std::ofstream ofsVtp;
         ofsVtp.open( filePathVtp.str().c_str() );
  
         write( ofsVtp );

         std::ostringstream filePathPvd;
         filePathPvd << baseFolder_ << '/' << identifier_ << ".pvd";
         std::ofstream ofsPvd( filePathPvd.str().c_str() );
         writePVD( ofsPvd );
      }
   }

   incrementTimeStep();
}


} // namespace mesh
} // namespace walberla
