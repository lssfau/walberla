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
//! \file Exports.cpp
//! \ingroup timeloop
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "python_coupling/Manager.h"
#include "python_coupling/helper/ModuleScope.h"

#include "vtk/VTKOutput.h"

using namespace boost::python;


namespace walberla {
namespace vtk {

namespace internal {


   shared_ptr<VTKOutput> VTKOutput_create(const shared_ptr<StructuredBlockStorage> & sbs, const std::string & identifier,
                                          const std::string & baseFolder, const std::string & executionFolder,
                                          const bool binary, const bool littleEndian, const bool useMPIIO,
                                          uint_t ghostLayers=0)
   {
      return createVTKOutput_BlockData(*sbs, identifier, 1, ghostLayers, false, baseFolder, executionFolder,
                                       true, binary, littleEndian, useMPIIO, 0);
   }

   void VTKOutput_write(const shared_ptr<VTKOutput> &vtkOut, int step)
   {
       if (step < 0)
       {
          PyErr_SetString(PyExc_ValueError, "Step parameter has to be positive");
          throw boost::python::error_already_set();
       }
       vtkOut->forceWrite(uint_c(step));
   }
} // namespace internal


void exportModuleToPython()
{
   python_coupling::ModuleScope timeloopModule( "vtk" );


   void ( VTKOutput::*p_setSamplingResolution1) ( const real_t  ) = &VTKOutput::setSamplingResolution;
   void ( VTKOutput::*p_setSamplingResolution2) ( const real_t, const real_t, const real_t ) = &VTKOutput::setSamplingResolution;

   class_<BlockCellDataWriterInterface, //NOLINT
           boost::noncopyable,
           shared_ptr<BlockCellDataWriterInterface> > ("BlockCellDataWriterInterface", no_init)
       ;

   def("makeOutput", internal::VTKOutput_create, (arg("blocks"), arg("name"), arg("baseFolder")=".",
                                                  arg("executionFolder")="vtk", arg("binary")=true,
                                                  arg("littleEndian")=true, arg("useMPIIO")=true,
                                                  arg("ghostLayers")=0));

   class_<VTKOutput, shared_ptr<VTKOutput>, boost::noncopyable > ("VTKOutput", no_init)
      .def( "addCellDataWriter"     , &VTKOutput::addCellDataWriter )
      .def( "write"                 , &internal::VTKOutput_write )
      .def( "__call__"              , &internal::VTKOutput_write )
      .def( "addAABBInclusionFilter", &VTKOutput::addAABBInclusionFilter )
      .def( "addAABBExclusionFilter", &VTKOutput::addAABBExclusionFilter )
      .def( "setSamplingResolution" , p_setSamplingResolution1 )
      .def( "setSamplingResolution" , p_setSamplingResolution2 )
      ;
}


} // namespace vtk
} // namespace walberla


#endif //WALBERLA_BUILD_WITH_PYTHON

