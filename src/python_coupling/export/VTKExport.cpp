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
//! \file VTKExport.cpp
//! \ingroup vtk
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "blockforest/StructuredBlockForest.h"
#include "vtk/VTKOutput.h"

namespace walberla {
namespace vtk {

namespace internal {

namespace py = pybind11;

void VTKOutput_write(const shared_ptr<VTKOutput> &vtkOut, int step)
{
    if (step < 0)
    {
       throw py::value_error("Step parameter has to be positive");
    }
    vtkOut->forceWrite(uint_c(step));
}
} // namespace internal

namespace py = pybind11;
using namespace pybind11::literals;
void exportModuleToPython(py::module_& m)
{
   py::module_ m2 = m.def_submodule("vtk", "VTK Extension of the waLBerla python bindings");


   void ( VTKOutput::*p_setSamplingResolution1) ( const real_t  ) = &VTKOutput::setSamplingResolution;
   void ( VTKOutput::*p_setSamplingResolution2) ( const real_t, const real_t, const real_t ) = &VTKOutput::setSamplingResolution;

   py::class_<BlockCellDataWriterInterface, shared_ptr<BlockCellDataWriterInterface> > (m2, "BlockCellDataWriterInterface");

   m2.def(
      "makeOutput",
      [](const shared_ptr<StructuredBlockForest> & sbf, const std::string & identifier,
         const std::string & baseFolder, const std::string & executionFolder,
         const bool binary, const bool littleEndian, const bool useMPIIO, uint_t ghostLayers=0)
      {
        return createVTKOutput_BlockData(*sbf, identifier, 1, ghostLayers, false, baseFolder, executionFolder,
                                         true, binary, littleEndian, useMPIIO, 0);
        },
        "blocks"_a, "name"_a, "baseFolder"_a=".", "executionFolder"_a="vtk", "binary"_a=true,
                              "littleEndian"_a=true, "useMPIIO"_a=true, "ghostLayers"_a=0);


   py::class_<VTKOutput, shared_ptr<VTKOutput> > (m2, "VTKOutput")
      .def( "addCellDataWriter"     , &VTKOutput::addCellDataWriter )
      .def( "write"                 , &internal::VTKOutput_write )
      .def( "__call__"              , &internal::VTKOutput_write )
      .def( "addAABBInclusionFilter", &VTKOutput::addAABBInclusionFilter )
      .def( "addAABBExclusionFilter", &VTKOutput::addAABBExclusionFilter )
      .def( "setSamplingResolution" , p_setSamplingResolution1 )
      .def( "setSamplingResolution" , p_setSamplingResolution2 );
}


} // namespace vtk
} // namespace walberla


#endif //WALBERLA_BUILD_WITH_PYTHON

