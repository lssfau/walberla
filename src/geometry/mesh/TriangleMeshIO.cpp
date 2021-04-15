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
//! \file TriangleMeshIO.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Implementation of IO functions for Mesh data structure in OBJ and POV format
//
//======================================================================================================================

#include "TriangleMesh.h"
#include "TriangleMeshComm.h"
#include "TriangleMeshIO.h"

#include "core/Abort.h"
#include "core/logging/Logging.h"
#include "core/math/AABB.h"
#include "core/mpi/Broadcast.h"
#include "core/Regex.h"
#include "core/StringUtility.h"

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>


namespace walberla {
namespace geometry {

   void readMesh ( const std::string & meshFilename, TriangleMesh & mesh )
   {
      mesh.clear();

      WALBERLA_LOG_PROGRESS("Loading mesh " << meshFilename << "..." );

      std::ifstream is( meshFilename.c_str() );
      if( is.fail() )
         WALBERLA_ABORT( "Error while opening file " << meshFilename << "!" );

      if ( string_ends_with( meshFilename, ".obj")  )
         readMeshObj( is, mesh );
      else if ( string_ends_with( meshFilename, ".pov") )
         readMeshPov( is, mesh );
      else if ( string_ends_with( meshFilename, ".off") )
         readMeshOff( is, mesh );
      else
         WALBERLA_ABORT( "Unknown mesh file format when loading " << meshFilename << ". Supported formats are obj, pov and off." );

      WALBERLA_LOG_PROGRESS( "Loaded mesh " << meshFilename << ". It has " << mesh.getNumTriangles() << " triangles, "
         << mesh.getNumVertices() << " vertices, AABB " << mesh.getAABB() << " and a volume of " << mesh.volume() << "." );
   }



   void writeMesh ( const std::string & meshFilename, const TriangleMesh & mesh )
   {
      WALBERLA_LOG_PROGRESS("Writing mesh " << meshFilename << "..." );

      std::ofstream os( meshFilename.c_str() );
      if( os.fail() )
         WALBERLA_ABORT( "Error while opening file " << meshFilename << "!" );

      if ( string_ends_with( meshFilename, ".obj")  )
         writeMeshObj( os, mesh );
      else if ( string_ends_with( meshFilename, ".pov") )
         writeMeshPov( os, mesh );
      else if ( string_ends_with( meshFilename, ".off") )
         writeMeshOff( os, mesh );
      else if ( string_ends_with( meshFilename, ".vtp") )
         writeMeshVtp( os, mesh );
      else
         WALBERLA_ABORT( "Unknown mesh file format when writing " << meshFilename << ". Supported formats are: obj,pov,off and vtp.");
   }



   void readAndBroadcastMesh ( const std::string & meshFilename, TriangleMesh & mesh )
   {
      mesh.clear();

      WALBERLA_ROOT_SECTION()
      {
         readMesh( meshFilename, mesh );
      }

      WALBERLA_MPI_WORLD_BARRIER();

      mpi::broadcastObject( mesh );
      WALBERLA_MPI_WORLD_BARRIER();
   }



   void writeMeshOnRoot ( const std::string & meshFilename, const TriangleMesh & mesh )
   {
      WALBERLA_ROOT_SECTION()
      {
         writeMesh( meshFilename, mesh );
      }
   }



   //===================================================================================================================
   //
   //  OBJ FORMAT
   //
   //===================================================================================================================

   void readFaceVertex( std::istream & is, TriangleMesh::index_t & vertexIdx, TriangleMesh::index_t & textureIdx, TriangleMesh::index_t & normalIdx)
   {
      is >> vertexIdx;

      if( is.peek() == '/' )
      {
         is.ignore();
         if( is.peek() == '/' )
         {
            is.ignore();
            textureIdx = 0;
            is >> normalIdx;
         }
         else // TextureIdx specified
         {
            is >> textureIdx;
            if( is.peek() == '/' )
            {
               is.ignore();
               is >> normalIdx;
            }
            else
            {
               normalIdx = 0;
            }
         }
      }
      else // only vertexIdx specified
      {
         textureIdx = 0;
         normalIdx  = 0;
      }
   }

   void readMeshObj( std::istream & is, TriangleMesh & mesh)
   {
      using std::string;
      using std::stringstream;

      mesh.clear();

      string curLine;

      while (getline(is, curLine))
      {
         stringstream in(curLine);

         char type;
         in >> type;

         if( !in )
            continue;

         switch (type) {
            case '#': // Skip comments
               continue;
               break;
            case 'v': //Vertex
               if( in.peek() == 'n' )
               {
                  in.ignore();
                  real_t nx;
                  real_t ny;
                  real_t nz;
                  in >> nx >> ny >> nz;
                  mesh.addVertexNormal( TriangleMesh::normal_t(nx,ny,nz) );
                  continue;
               }
               real_t x;
               real_t y;
               real_t z;
               in >> x >> y >> z;
               float r;
               float g;
               float b;
               in >> r >> g >> b;
               if( in )
                  mesh.addVertex( TriangleMesh::vertex_t(x,y,z), TriangleMesh::color_t(r,g,b) );
               else
                  mesh.addVertex( TriangleMesh::vertex_t(x,y,z) );
               break;
            case 'f': //Face
               TriangleMesh::index_t i0, i1, i2, n0, n1, n2, dummy;
               readFaceVertex( in, i0, dummy, n0 );
               readFaceVertex( in, i1, dummy, n1 );
               readFaceVertex( in, i2, dummy, n2 );

               //In Obj-File indexing starts at 1, so subtract 1
               if( n0 > 0 && n1 > 0 && n2 > 0 )
                  mesh.addTriangle(--i0, --i1, --i2, --n0, --n1, --n2);
               else
                  mesh.addTriangle(--i0, --i1, --i2);
               break;
            default:
               //Skip other directives
               continue;
               break;
         }
      }
   }


   void writeMeshObj (std::ostream & os, const TriangleMesh & mesh)
   {
      os << "# walberla post processing Mesh\n\n";

      // Write vertices
      os << "# Vertices\n";
      for(TriangleMesh::index_t i = 0; i < mesh.getNumVertices(); ++i)
      {
         const TriangleMesh::vertex_t & v = mesh.getVertex(i);
         os << "v  " << v[0] << " " << v[1] << " " << v[2];
         if( mesh.hasVertexColors() )
         {
            const TriangleMesh::color_t & c = mesh.getVertexColor(i);
            os << " " << c[0] << " " << c[1] << " " << c[2];
         }
         os << '\n';
      }
      os << '\n';

      if( mesh.hasVertexNormals() )
      {
         // Write Normals
         os << "# Normals\n";
         for(TriangleMesh::index_t i = 0; i < mesh.getNumNormals(); ++i)
         {
            const TriangleMesh::normal_t & n = mesh.getVertexNormal(i);
            os << "vn " << n[0] << " " << n[1] << " " << n[2] << '\n';
         }

         // Write triangles
         os << "# Faces\n";
         for(size_t i = 0 ; i < mesh.getNumTriangles(); ++i)
         {
            os << "f  ";
            os << mesh.getVertexIndex(i,0)+1 << "//" << mesh.getNormalIndex(i,0)+1 << ' ';
            os << mesh.getVertexIndex(i,1)+1 << "//" << mesh.getNormalIndex(i,1)+1 << ' ';
            os << mesh.getVertexIndex(i,2)+1 << "//" << mesh.getNormalIndex(i,2)+1 << '\n';
         }
      }
      else // No vertex normals
      {
         os << "# Faces\n";
         for(size_t i = 0 ; i < mesh.getNumTriangles(); ++i)
         {
            os << "f  ";
            os << mesh.getVertexIndex(i,0) + 1 << ' ';
            os << mesh.getVertexIndex(i,1) + 1 << ' ';
            os << mesh.getVertexIndex(i,2) + 1 << '\n';
         }
      }
      os.flush();
   }



   //===================================================================================================================
   //
   //  POV FORMAT
   //
   //===================================================================================================================


   void removeComments( std::istream & is, std::ostream & os )
   {
      char c;
      enum State{ STATE_NORMAL, STATE_C, STATE_LC, STATE_MLC0, STATE_MLC1};
      State state = STATE_NORMAL;

      while( is.good() )          // loop while extraction from file is possible
      {
         is.get(c);              // get character from file

         switch(state)
         {
         case STATE_NORMAL:
            if( c == '/' )
               state = STATE_C;
            else if( c == ' ' || c == '<' || c == '>' || c == ',' )
               os.put(' ');
            else if( !isspace(c) )
               os.put(c);
            break;

         case STATE_C:
            if( c == '/' )
               state = STATE_LC;
            else if ( c == '*' )
               state = STATE_MLC0;
            else
               throw std::runtime_error("Invalid inputstream syntax");
            break;

         case STATE_LC:
            if( c == '\n' )
               state = STATE_NORMAL;
            break;

         case STATE_MLC0:
            if( c == '*' )
               state = STATE_MLC1;
            break;

         case STATE_MLC1:
            if( c == '/' )
               state = STATE_NORMAL;
            else if( c != '*' )
               state = STATE_MLC0;
            break;
         }
      }
   }

   enum State{ VERTEX, NORMAL, FACE };

   void readMeshPov  (std::istream & is, TriangleMesh & mesh, bool clear)
   {
      // TODO implement reading povray mesh2
      //throw std::runtime_error("not yet implemented");

      using std::string;
      using std::stringstream;

      std::string source( (std::istreambuf_iterator<char>( is ) ), std::istreambuf_iterator<char>( ) );

      // replace multiline comments: /\\*.*?\\*/
      // replace single line comments //.*?\n
      // replace chars: < > , \ n \t
      walberla::regex r( "/\\*.*?\\*/|//.*?\n|<|>|,|\n|\t" );
      std::string stripped = walberla::regex_replace( source , r , " " ) ;

      TriangleMesh::index_t faceOffset = 0u;
      if( clear )
         mesh.clear();
      else
         faceOffset = mesh.getNumVertices();

      std::vector< string > splitVec = string_split( stripped, "{}" );

      State state;

      std::map<State, std::string> stateMap;

      for( size_t i=1; i<splitVec.size(); ++i )
      {
         string_trim(splitVec[i]);
         if( splitVec[i] == "vertex_vectors" ) {
            state = VERTEX;
         } else if ( splitVec[i] == "normal_vectors" ) {
            state = NORMAL;
         } else if ( splitVec[i] == "face_indices" ) {
            state = FACE;
         } else {
            //std::cerr << "Unknown section in povray file: " << splitVec[i] << "\n";
            continue;
         }
         stateMap[state] = splitVec[++i];
      }

      if( stateMap.empty() )
         return;

      std::stringstream vin(stateMap[VERTEX]);
      std::stringstream nin(stateMap[NORMAL]);
      size_t vcount;
      size_t ncount;
      vin >> vcount;
      nin >> ncount;
      real_t x;
      real_t y;
      real_t z;
      real_t nx;
      real_t ny;
      real_t nz;
      for( size_t j=0; j<vcount; ++j )
      {
         vin >>  x >>  y >>  z;
         nin >> nx >> ny >> nz;
         mesh.addVertex      ( TriangleMesh::vertex_t(  x,  y,  z ) );
         mesh.addVertexNormal( TriangleMesh::normal_t( nx, ny, nz ) );
      }

      std::stringstream fin(stateMap[FACE]);
      size_t fcount;
      fin >> fcount;
      TriangleMesh::index_t ix;
      TriangleMesh::index_t iy;
      TriangleMesh::index_t iz;
      for( size_t j=0; j<fcount; ++j )
      {
         fin >> ix >> iy >> iz;
         mesh.addTriangle( ix+faceOffset, iy+faceOffset, iz+faceOffset );
      }
   }



   void writeMeshPov0 (std::ostream & os, const TriangleMesh & mesh, const size_t itemsPerLine)
   {
      os << "mesh2{\n\tvertex_vectors{\n\t\t" << mesh.getNumVertices() << ",\n\t\t";

      for( TriangleMesh::index_t i = 0; i < mesh.getNumVertices(); )
      {
         const TriangleMesh::vertex_t & v = mesh.getVertex(i);
         os << "<" << v[0] << "," << v[1] << "," << v[2] << ">";

         ++i;
         if( i != mesh.getNumVertices() ) os << ",";
         if( itemsPerLine > 0u && i % itemsPerLine == 0 ) os << "\n\t\t";
      }

      os << "\n\t}\n";

      if( mesh.hasVertexNormals() )
      {
         // Write normals
         os << "\tnormal_vectors {\n\t\t" << mesh.getNumVertices() << ",\n\t\t";

         const std::streamsize p = os.precision();
         os.precision(2);
         for(TriangleMesh::index_t i = 0; i < mesh.getNumVertices(); )
         {
            const TriangleMesh::normal_t & v = mesh.getVertexNormal(i);
            os << "<" << v[0] << "," << v[1] << "," << v[2] << ">";

            ++i;
            if( i != mesh.getNumVertices() ) os << ",";
            if( itemsPerLine > 0u && i % itemsPerLine == 0 ) os << "\n\t\t";
         }
         os.precision(p);

         os << "\n\t}\n";
      }

      // Write Face indices
      os << "\tface_indices {\n\t\t" << mesh.getNumTriangles() << ",\n\t\t";

      for( size_t i = 0 ; i < mesh.getNumTriangles(); )
      {
         os << "<" << mesh.getVertexIndex(i,0) << "," <<
                      mesh.getVertexIndex(i,1) << "," <<
                      mesh.getVertexIndex(i,2) << ">";

         ++i;
         if( i != mesh.getNumTriangles() ) os << ",";
         if( itemsPerLine > 0u && i % itemsPerLine == 0 ) os << "\n\t\t";
      }
      os << "\n\t}"; //end face_indices


      if( mesh.hasNormalIndices() )
      {
         // Write vertex normal indices
         os << "\tnormal_indices {\n\t\t" << mesh.getNumTriangles() << ",\n\t\t";

         for(size_t i = 0 ; i < mesh.getNumTriangles(); )
         {
            os << "<" << mesh.getNormalIndex(i,0) << "," <<
                         mesh.getNormalIndex(i,1) << "," <<
                         mesh.getNormalIndex(i,2) << ">";

            ++i;
            if( i != mesh.getNumTriangles() ) os << ",";
            if( itemsPerLine > 0u && i % itemsPerLine == 0 ) os << "\n\t\t";
         }
         os << "\n\t}"; //end normal_indices
      }

      os << "\n}\n"; //end mesh2
   }




   void writeMeshPov (std::ostream & os, const TriangleMesh & mesh, size_t itemsPerLine, size_t itemsPerMesh)
   {
      if( itemsPerMesh == 0u || mesh.getNumVertices() < itemsPerMesh ){
         writeMeshPov0(os, mesh, itemsPerLine);
         return;
      }

      std::vector<TriangleMesh> meshVec( 1u + mesh.getNumVertices() / itemsPerMesh );
      bool hasVertexNormals = mesh.hasVertexNormals();

      size_t j;

      for( TriangleMesh::index_t i = 0; i < mesh.getNumVertices(); ++i ){
         j = i / itemsPerMesh;
         meshVec[j].addVertex( mesh.getVertex(i) );
         if( hasVertexNormals )
            meshVec[j].addVertexNormal( mesh.getVertexNormal(i) );
      }

      for( size_t i = 0; i < mesh.getNumTriangles(); ++i ){
         TriangleMesh::index_t ix = mesh.getVertexIndex( i, 0 );
         TriangleMesh::index_t iy = mesh.getVertexIndex( i, 1 );
         TriangleMesh::index_t iz = mesh.getVertexIndex( i, 2 );

         size_t jx = ix / itemsPerMesh;
         size_t jy = iy / itemsPerMesh;
         size_t jz = iz / itemsPerMesh;

         TriangleMesh::index_t nix = uint32_c( ix % itemsPerMesh );
         TriangleMesh::index_t niy = uint32_c( iy % itemsPerMesh );
         TriangleMesh::index_t niz = uint32_c( iz % itemsPerMesh );

         if ( jx != jy || jy != jz ){
            j = std::max( std::max(jx, jy), jz );
            if( jx < j ){
               nix = meshVec[j].addVertex( mesh.getVertex(ix) );
               if( hasVertexNormals )
                  meshVec[j].addVertexNormal( mesh.getVertexNormal(ix) );
            }
            if( jy < j ){
               niy = meshVec[j].addVertex( mesh.getVertex(iy) );
               if( hasVertexNormals )
                  meshVec[j].addVertexNormal( mesh.getVertexNormal(iy) );
            }
            if( jz < j ){
               niz = meshVec[j].addVertex( mesh.getVertex(iz) );
               if( hasVertexNormals )
                  meshVec[j].addVertexNormal( mesh.getVertexNormal(iz) );
            }
         } else {
            j = jx;
         }
         meshVec[j].addTriangle( nix, niy, niz );
      }

      os << "union{\n";
      for( j = 0u; j < meshVec.size(); ++j ){
         meshVec[j].removeDuplicateVertices();
         writeMeshPov0(os, meshVec[j], itemsPerLine);
      }
      os << "}\n";
   }


   //===================================================================================================================
   //
   //  Geomview Object File Format (*.off)
   //
   //===================================================================================================================

   static void skipComments( std::istream & is )
   {
      while( is.peek() == '#' )
      {
         is.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
      }
   }

   void readMeshOff( std::istream & is, TriangleMesh & mesh)
   {
      mesh.clear();

      skipComments( is );

      std::string firstLine;
      std::getline( is, firstLine );
      if( firstLine != "OFF" )
         WALBERLA_ABORT( "Trying to read a mesh in Geomview Object File Format, but the first non-comment lin eis not \"OFF\"!" );

      skipComments( is );

      uint_t vertexCount;
      uint_t faceCount;
      uint_t edgeCount;
      is >> vertexCount >> faceCount >> edgeCount;

      skipComments( is );

      for( uint_t i = 0; i < vertexCount; ++i )
      {
         TriangleMesh::vertex_t vertex;
         is >> vertex[0] >> vertex[1] >> vertex[2];
         mesh.addVertex( vertex );
         skipComments( is );
      }

      for( uint_t i = 0; i < faceCount; ++i )
      {
         uint_t numVertices;
         TriangleMesh::index_t i0;
         TriangleMesh::index_t i1;
         TriangleMesh::index_t i2;
         is >> numVertices >> i0 >> i1 >> i2;
         if( numVertices != 3 )
            WALBERLA_ABORT( "Face with more or less than 3 vertices given while trying to read a mesh in Geomview Object File Format!" );
         mesh.addTriangle(i0, i1, i2);
         skipComments( is );
      }
   }

   void writeMeshOff( std::ostream & os, const TriangleMesh & mesh )
   {
      TriangleMesh::index_t numVertices = mesh.getNumVertices();
      size_t                numFaces    = mesh.getNumTriangles();

      os << "OFF\n";
      os << numVertices << ' ' << mesh.getNumTriangles() << " 0\n"; // Number of edges is unknown, we put 0

      for( TriangleMesh::index_t i = 0; i < numVertices; ++i )
      {
         TriangleMesh::vertex_t v = mesh.getVertex( i );
         os << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
      }

      for( size_t i = 0; i < numFaces; ++i )
      {
         os << "3 " << mesh.getVertexIndex(i, 0)
            << ' '  << mesh.getVertexIndex(i, 1)
            << ' '  << mesh.getVertexIndex(i, 2) << '\n';
      }
   }

   void writeMeshVtp  ( std::ostream & os, const TriangleMesh & mesh )
   {
      const TriangleMesh::index_t numVertices = mesh.getNumVertices();
      const size_t                numFaces    = mesh.getNumTriangles();

      os << "<VTKFile type=\"PolyData\" version=\"0.1\">\n"
         << "  <PolyData>\n"
         << "    <Piece NumberOfPoints=\"" << numVertices << "\" NumberOfPolys=\"" << numFaces << "\">\n";

      os << "      <Points>\n"
         << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

      os << "          ";
      for( auto it = mesh.getVertices().begin(); it != mesh.getVertices().end(); ++it )
      {
         os << (*it)[0] << ' ' << (*it)[1] << ' ' << (*it)[2] << ' ';
      }
      os  << '\n';

      os << "        </DataArray>\n"
         << "      </Points>\n";

      os << "      <Polys>\n";

      os << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
      os << "          ";
      for( auto it = mesh.getVertexIndices().begin(); it != mesh.getVertexIndices().end(); ++it )
      {
         os << *it << ' ';
      }
      os << '\n';
      os << "        </DataArray>\n";

      os << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
      os << "          ";
      for( size_t i = 1; i <= mesh.getNumTriangles(); ++i )
      {
         os << 3 * i << ' ';
      }
      os << '\n';
      os << "        </DataArray>\n";

      os << "      </Polys>\n";

      os << "      <CellData></CellData>\n";

      os << "      <PointData" << ( mesh.hasVertexColors() ? " Scalars=\"vertexColors\"" : "" )
                               << ( mesh.hasVertexColors() ? " Normals=\"vertexNormals\"" : "" )
                               << ">\n";

      if( mesh.hasVertexColors() )
      {
         os << "        <DataArray type=\"UInt8\" Name=\"vertexColors\" NumberOfComponents=\"3\" format=\"ascii\">\n";
         os << "          ";
         for( auto it = mesh.getVertexColors().begin(); it != mesh.getVertexColors().end(); ++it )
         {
            os << std::lround( (*it)[0] * 255.0f ) << ' '
               << std::lround( (*it)[1] * 255.0f ) << ' '
               << std::lround( (*it)[2] * 255.0f ) << ' ';
         }
         os << "        </DataArray>\n";
      }

      if( mesh.hasVertexNormals() )
      {
         os << "        <DataArray type=\"Float32\" Name=\"vertexNormals\" NumberOfComponents=\"3\" format=\"ascii\">\n";
         os << "          ";
         for( auto it = mesh.getVertexNormals().begin(); it != mesh.getVertexNormals().end(); ++it )
         {
            os << (*it)[0] << ' ' << (*it)[1] << ' ' << (*it)[2] << ' ';
         }
         os << "        </DataArray>\n";
      }

      os << "      </PointData>\n";


      os << "    </Piece>\n"
         << "  </PolyData>\n"
         << "</VTKFile>\n";
   }


} // namespace geometry
} // namespace walberla



