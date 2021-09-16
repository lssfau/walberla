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
//! \file TriangleMeshIO.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Input/Output functions for mesh data structure in OBJ and POV format
//
//======================================================================================================================

#pragma once

#include <iostream>


namespace walberla {
namespace geometry {

   class TriangleMesh;

   /**
    * \brief Reads mesh from file
    *
    * Mesh format is detected by file ending.
    *
    * \param meshFilename  Filename of mesh. Supported formats are obj, pov and off. Format is detected
    *                  by file ending
    * \param mesh      object where the mesh is stored, if the given mesh is not
    *                  empty, all contents are cleared
    */
   void readMesh ( const std::string & meshFilename, TriangleMesh & mesh );

   /**
    * \brief writes a mesh to file
    *
    * Mesh format is detected by file ending.
    *
    * \param meshFilename  Filename of mesh. Supported formats are obj, pov and off. Format is detected
    *                  by file ending
    * \param mesh      object where the mesh is read from
    */
   void writeMesh ( const std::string & meshFilename, const TriangleMesh & mesh );

   /**
    * \brief Reads mesh from file on root process and broadcasts it to all other processes
    *
    * Mesh format is detected by file ending.
    *
    * \param meshFilename  Filename of mesh. Supported formats are obj, pov and off. Format is detected
    *                  by file ending
    * \param mesh      object where the mesh is stored, if the given mesh is not
    *                  empty, all contents are cleared
    */
   void readAndBroadcastMesh ( const std::string & meshFilename, TriangleMesh & mesh );



   /**
    * \brief writes a mesh to file on root process
    *
    * Mesh format is detected by file ending.
    *
    * \param meshFilename  Filename of mesh. Supported formats are obj, pov and off. Format is detected
    *                  by file ending
    * \param mesh      object where the mesh is read from
    */
   void writeMeshOnRoot ( const std::string & meshFilename, const TriangleMesh & mesh );


   /**
    * \brief Reads mesh from input stream in obj file format
    *
    *  Reads "v" (vertex) and "f" (face) entries only.
    *
    * \param is    input stream, to read from file use std::fstream
    * \param mesh  object where the mesh is stored, if the given mesh is not
    *              empty, all contents are cleared
    */
   void readMeshObj  (std::istream & is, TriangleMesh & mesh);


   /**
    * \brief Writes mesh to output stream in obj format
    *
    * Writes obj file with vertex, normal and face information
    *
    * \param os   the output stream. To write to file call with ofstream("filename.obj"))
    * \param mesh the mesh to write
    */
   void writeMeshObj (std::ostream & os, const TriangleMesh & mesh);


   /**
    * \brief Reads mesh from input stream in povray's mesh2 format
    *
    * Vertices, normals and faces are read.
    *
    * \param is    input stream, to read from file use std::fstream
    * \param mesh  object where the mesh is stored, if the given mesh is not
    *              empty, all contents are cleared
    */
   void readMeshPov  (std::istream & is, TriangleMesh & mesh, bool clear = true);


   /**
    * \brief Writes mesh in povray's mesh2 format
    *
    * Vertices, normals and faces are written.
    *
    * \param os   the output stream. To write to file call with ofstream("filename.dat")
    * \param mesh the mesh to write
    */
   void writeMeshPov (std::ostream & os, const TriangleMesh & mesh, size_t itemsPerLine = 5u, size_t itemsPerMesh = 0u);


    /**
    * \brief Reads mesh from input stream in Geomview Object File Format
    *
    *  Reads "v" (vertex) and "f" (face) entries only. Only triangular faces are supported.
    *  Colors are currently also not supported.
    *  Format details: http://people.sc.fsu.edu/~jburkardt/data/off/off.html
    *
    * \param is    input stream, to read from file use std::fstream
    * \param mesh  object where the mesh is stored, if the given mesh is not
    *              empty, all contents are cleared
    */
   void readMeshOff  (std::istream & is, TriangleMesh & mesh);


    /**
    * \brief Writes a mesh to an output stream in Geomview Object File Format
    *
    *  Writes vertices and faces only. Colors are not (yet) supported.
    *  Format details: http://people.sc.fsu.edu/~jburkardt/data/off/off.html
    *
    * \param os    the output stream. To write to file call with ofstream("filename.off")
    * \param mesh  the mesh to write
    */
   void writeMeshOff  ( std::ostream & os, const TriangleMesh & mesh );


    /**
    * \brief Writes a mesh to an output stream in VTK Poly Data format
    *
    *  Writes colors and vertex normals if available
    *  Format details: http://www.vtk.org/VTK/img/file-formats.pdf
    *
    * \param os    the output stream. To write to file call with ofstream("filename.vtp")
    * \param mesh  the mesh to write
    */
   void writeMeshVtp  ( std::ostream & os, const TriangleMesh & mesh );


} // namespace geometry
} // namespace walberla


