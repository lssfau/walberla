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
//! \file Exports.impl.h
//! \ingroup mesh
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"
#include "python_coupling/helper/ModuleScope.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON
#ifdef WALBERLA_BUILD_WITH_OPENMESH

#include "python_coupling/Manager.h"
#include "python_coupling/helper/MplHelpers.h"
#include "python_coupling/helper/BlockStorageExportHelpers.h"
#include "python_coupling/helper/PythonIterableToStdVector.h"
#include "python_coupling/helper/SharedPtrDeleter.h"

#ifdef _MSC_VER
#  pragma warning( push )
#  pragma warning( disable : 4456 4459 )
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#ifdef _MSC_VER
#  pragma warning( pop )
#endif

#include "mesh/MeshIO.h"
#include "mesh/distance_octree/DistanceOctree.h"
#include "mesh/DistanceComputations.h"
#include "mesh/DistanceFunction.h"
#include "mesh/MeshOperations.h"
#include "mesh/TriangleMeshes.h"
#include "mesh/boundary/BoundarySetup.h"
#include "mesh/vtk/VTKMeshWriter.h"

using namespace boost::python;


namespace walberla {
namespace mesh {

typedef mesh::DistanceOctree<PythonTriangleMesh> Octree;
typedef mesh::VTKMeshWriter<PythonTriangleMesh> MeshWriter;

namespace internal
{
   void exp_readAndBroadcast(const std::string & fileName, PythonTriangleMesh & mesh )
   {
      mesh::readAndBroadcast(fileName, mesh);
   }


   shared_ptr<Octree> makeDistanceOctree(object meshObj, uint_t maxDepth, uint_t minNumTriangles)
   {
      auto mesh = python_coupling::createSharedPtrFromPythonObject<PythonTriangleMesh>(meshObj);
      auto triangleDistance = make_shared<TriangleDistance<PythonTriangleMesh> > (mesh);
      return make_shared< Octree >(triangleDistance, maxDepth, minNumTriangles );
   }

   shared_ptr<MeshWriter> makeVTKMeshWriter(object meshObj, const std::string & identifier, const uint_t writeFrequency, const std::string & baseFolder)
   {
      shared_ptr<PythonTriangleMesh> mesh = python_coupling::createSharedPtrFromPythonObject<PythonTriangleMesh>(meshObj);
      return walberla::make_shared< MeshWriter >( mesh, identifier, writeFrequency, baseFolder );
   }

   // ----------------------------------- markFlagField ----------------------------------------------------------------


   template<typename FlagField>
   void markFlagFieldTmpl(const shared_ptr<Octree> & octree, const shared_ptr<StructuredBlockStorage> & bs, BlockDataID flagFieldID,
                          const std::string & flagName, uint_t ghostLayers)
   {
      BoundarySetup setup( bs, makeMeshDistanceFunction( octree ), ghostLayers);
      setup.setFlag<FlagField>(flagFieldID, flagName, BoundarySetup::INSIDE);
   }

   FunctionExporterClass( markFlagFieldTmpl,
                          void ( const shared_ptr<Octree> &, const shared_ptr<StructuredBlockStorage> & ,
                                 BlockDataID ,const std::string &, uint_t ) );

   template<typename FlagFields>
   void exp_markFlagField(const shared_ptr<Octree> & octree, const shared_ptr<StructuredBlockStorage> & bs,
                          const std::string & blockDataStr, const std::string & flagName, uint_t ghostLayers)
   {
      if ( bs->begin() == bs->end() )
         return;
      IBlock * firstBlock =  & ( * bs->begin() );

      auto fieldID = python_coupling::blockDataIDFromString( *bs, blockDataStr );
      python_coupling::Dispatcher<FlagFields, Exporter_markFlagFieldTmpl > dispatcher( firstBlock );
      dispatcher( fieldID )(octree, bs, fieldID, flagName, ghostLayers);
   }
}


bool Octree_isAABBFullyInside(const Octree & octree, const AABB & aabb)
{
   for( auto corner: aabb.corners() )
   {
      const Octree::Point p ( numeric_cast<Octree::Scalar>(corner[0]),
                              numeric_cast<Octree::Scalar>(corner[1]),
                              numeric_cast<Octree::Scalar>(corner[2]) );
      if( octree.sqSignedDistance(p) > 0 )
         return false;
   }
   return true;
}


bool Octree_isAABBFullyOutside(const Octree & octree, const AABB & aabb)
{
   for( auto corner: aabb.corners() )
   {
      const Octree::Point p ( numeric_cast<Octree::Scalar>(corner[0]),
                              numeric_cast<Octree::Scalar>(corner[1]),
                              numeric_cast<Octree::Scalar>(corner[2]) );
      if( octree.sqSignedDistance(p) < 0 )
         return false;
   }
   return true;
}


template<typename FlagFields>
void exportModuleToPython()
{
   python_coupling::ModuleScope fieldModule( "mesh" );

   def( "readAndBroadcast", &internal::exp_readAndBroadcast, (arg("fileName"), arg("mesh")) );

   def( "aabb", &mesh::computeAABB<PythonTriangleMesh>, (arg("mesh")) );
   def( "translate", &mesh::translate<PythonTriangleMesh>, (arg("mesh"), arg("translationVector")) );
   def( "scale", &mesh::scale<PythonTriangleMesh>, (arg("mesh"), arg("scaleFactors")) );
   def( "volume", &mesh::computeVolume<PythonTriangleMesh>, (arg("mesh")) );
   def( "surfaceArea", &mesh::computeSurfaceArea<PythonTriangleMesh>, (arg("mesh")) );

   def( "markFlagField", &internal::exp_markFlagField<FlagFields> );

   Octree::Scalar (Octree::*sqSignedDistance1)(const Octree::Point&) const= &Octree::sqSignedDistance;
   Octree::Scalar (Octree::*sqDistance1)(const Octree::Point&) const= &Octree::sqDistance;

   class_<Octree, shared_ptr<Octree> >("DistanceOctree", no_init)
      .def("__init__", make_constructor(internal::makeDistanceOctree, default_call_policies() ,(arg("mesh"), arg("maxDepth")=20u, arg("minNumTriangles")=25u) ) )
      .def("sqSignedDistance", sqSignedDistance1, (arg("p")))
      .def("sqDistance", sqDistance1, (arg("p")))
      .def("height", &Octree::height)
      .def("writeVTKOutput", &Octree::writeVTKOutput, (arg("filestem")))
      .def("isAABBfullyOutside", &Octree_isAABBFullyOutside)
      .def("isAABBfullyInside", &Octree_isAABBFullyInside)
   ;

   class_<MeshWriter, shared_ptr<MeshWriter> >("VTKMeshWriter", no_init)
      .def("__init__", make_constructor(internal::makeVTKMeshWriter, default_call_policies() ,(arg("mesh"), arg("identifier"), arg("writeFrequency")=1u, arg("baseFolder")="vtk_out") ) )
      .def("__call__", &MeshWriter::operator())
   ;

}



} // namespace mesh
} // namespace walberla


#endif //WALBERLA_BUILD_WITH_PYTHON
#endif //WALBERLA_BUILD_WITH_OPENMESH