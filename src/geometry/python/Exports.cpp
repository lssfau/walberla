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
//! \file PythonExports.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON
#include "python_coupling/helper/ModuleScope.h"

#include "core/math/AABB.h"
#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"



using namespace boost::python;


namespace walberla {
namespace geometry {

void triangleMesh_scaleReal( TriangleMesh & m, real_t scaleFactor )               { m.scale( scaleFactor ); }
void triangleMesh_scaleVec ( TriangleMesh & m, const Vector3<real_t> & scaleVec ) { m.scale( scaleVec );    }


shared_ptr<TriangleMesh> triangleMesh_load ( const std::string & filename, bool broadcast )
{
   auto mesh = make_shared<TriangleMesh>();

   try
   {
      if( broadcast )
         readAndBroadcastMesh( filename, *mesh );
      else
         readMesh( filename, *mesh );
   }
   catch( std::exception & e )
   {
      PyErr_SetString( PyExc_RuntimeError, e.what() );
      throw error_already_set();
   }
   return mesh;
}

void triangleMesh_save( TriangleMesh & mesh, const std::string & filename )
{
   try
   {
      writeMesh( filename, mesh );
   }
   catch( std::exception & e )
   {
      PyErr_SetString( PyExc_RuntimeError, e.what() );
      throw error_already_set();
   }
}



void exportModuleToPython()
{
   python_coupling::ModuleScope fieldModule( "geometry" );

   class_< TriangleMesh, shared_ptr<TriangleMesh>, boost::noncopyable > ( "TriangleMesh", no_init )
            .add_property( "numTriangles",             &TriangleMesh::getNumTriangles )
            .add_property( "numVertices" ,             &TriangleMesh::getNumVertices )
            .add_property( "numVertexNormals",         &TriangleMesh::getNumNormals )
            .def         ( "volume",                   &TriangleMesh::volume )
            .def         ( "scale",                    &triangleMesh_scaleReal, arg("factor") )
            .def         ( "exchangeAxes",             &TriangleMesh::exchangeAxes, ( arg("xAxisId"), arg("yAxisId"), arg("zAxisId") ) )
            .def         ( "scaleXYZ",                 &triangleMesh_scaleVec, arg("factors") )
            .def         ( "removeDuplicateVertices",  &TriangleMesh::removeDuplicateVertices, ( arg("tolerance") = 1e-4 ) )
            .def         ( "merge",                    &TriangleMesh::merge, ( arg("other"), arg("offset") = Vector3<real_t>(0) ) )
            .def         ( "getAABB",                  &TriangleMesh::getAABB )
            .def         ( "save",                     &triangleMesh_save, arg("filename") )
            .def         ( "load",                     &triangleMesh_load, ( arg("filename"), arg("broadcast") = true) ).staticmethod( "load" )
   ;
}


} // namespace geometry
} // namespace walberla


#endif //WALBERLA_BUILD_WITH_PYTHON
