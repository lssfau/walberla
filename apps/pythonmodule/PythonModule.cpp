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
//! \file PythonModule.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "field/GhostLayerField.h"

#include "python_coupling/Manager.h"
#include "python_coupling/export/BlockForestExport.h"
#include "python_coupling/export/FieldExports.h"
#include "python_coupling/export/VTKExport.h"
#include "python_coupling/helper/ModuleInit.h"

#include "stencil/all.h"

#ifdef WALBERLA_BUILD_WITH_CUDA
 #include "python_coupling/export/CUDAExport.h"
#endif


using namespace walberla;


#define FIELD_TYPES \
   Field<walberla::real_t,1>,\
   Field<walberla::real_t,2>,\
   Field<walberla::real_t,3>,\
   Field<walberla::real_t,9>,\
   Field<walberla::real_t,15>,\
   Field<walberla::real_t,19>,\
   Field<walberla::real_t,27>,\
   Field<walberla::int8_t,1>,\
   Field<walberla::int64_t,1>,\
   Field<walberla::int64_t,2>,\
   Field<walberla::int64_t,3>,\
   Field<walberla::uint8_t,1>,\
   Field<walberla::uint16_t,1>,\
   Field<walberla::uint32_t,1>

#define GPU_FIELD_TYPES \
   GPUField<real_t>,\
   GPUField<int8_t>,\
   GPUField<int32_t>,\
   GPUField<int64_t>,\
   GPUField<uint8_t>,\
   GPUField<uint64_t>

struct InitObject
{
   InitObject()
   {
      auto pythonManager = python_coupling::Manager::instance();
      // Field
      pythonManager->addExporterFunction( field::exportModuleToPython<FIELD_TYPES> );
      pythonManager->addExporterFunction( field::exportGatherFunctions<FIELD_TYPES> );
      pythonManager->addBlockDataConversion<FIELD_TYPES>();
      // Blockforest
      pythonManager->addExporterFunction(blockforest::exportModuleToPython<stencil::D2Q5, stencil::D2Q9, stencil::D3Q7, stencil::D3Q19, stencil::D3Q27>);
      // VTK
      pythonManager->addExporterFunction( vtk::exportModuleToPython );
      #ifdef WALBERLA_BUILD_WITH_CUDA
            using walberla::cuda::GPUField;

            pythonManager->addExporterFunction( cuda::exportModuleToPython<GPU_FIELD_TYPES> );
            pythonManager->addExporterFunction( cuda::exportCopyFunctionsToPython<FIELD_TYPES> );
            pythonManager->addBlockDataConversion<GPU_FIELD_TYPES>();
      #endif
      //
      python_coupling::initWalberlaForPythonModule();

   }
};
InitObject globalInitObject;