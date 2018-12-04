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
//
//======================================================================================================================
#include "python_coupling/PythonWrapper.h"
#include "waLBerlaDefinitions.h"
#include "blockforest/python/Exports.h"
#include "field/GhostLayerField.h"
#include "field/python/Exports.h"
#include "mesh/python/Exports.h"
#include "geometry/python/Exports.h"
#include "postprocessing/python/Exports.h"
#include "python_coupling/Manager.h"
#include "python_coupling/helper/ModuleInit.h"
#include "stencil/all.h"
#include "timeloop/python/Exports.h"
#include "vtk/python/Exports.h"

#ifdef WALBERLA_BUILD_WITH_CUDA
#include "cuda/python/Exports.h"
#endif

#include <boost/mpl/vector.hpp>
#include <boost/mpl/insert_range.hpp>



namespace bmpl = boost::mpl;
using namespace walberla;

typedef bmpl::vector<
            Field<walberla::real_t,1>,
            Field<walberla::real_t,2>,
            Field<walberla::real_t,3>,
            Field<walberla::real_t,4>,
            Field<walberla::real_t,5>,
            Field<walberla::real_t,6>,
            Field<walberla::real_t,9>,
            Field<walberla::real_t,15>,
            Field<walberla::real_t,19>,
            Field<walberla::real_t,27>,

            Field<walberla::int8_t,1>,
            Field<walberla::int16_t,1>,
            Field<walberla::int32_t,1>,

            Field<walberla::int64_t,1>,
            Field<walberla::int64_t,2>,
            Field<walberla::int64_t,3>,
            Field<walberla::int64_t,4>,

            Field<walberla::uint8_t,1>,
            Field<walberla::uint16_t,1>,
            Field<walberla::uint32_t,1>
      > FieldTypes;


typedef bmpl::vector<
                      GhostLayerField<walberla::real_t,1>,
                      GhostLayerField<walberla::real_t,3>
                    >  FieldTypesForMeshGeneration;


typedef bmpl::vector< FlagField<walberla::uint8_t>,
                      FlagField<walberla::uint16_t> > FlagFieldTypes;

typedef bmpl::vector< stencil::D2Q5,
                      stencil::D2Q9,
                      stencil::D3Q7,
                      stencil::D3Q19,
                      stencil::D3Q27 > Stencils;

typedef GhostLayerField<walberla::real_t,3> VecField_T;


using namespace walberla;

struct InitObject
{
   InitObject()
   {
      namespace bmpl = boost::mpl;

      auto pythonManager = python_coupling::Manager::instance();

      // Field
      pythonManager->addExporterFunction( field::exportModuleToPython<FieldTypes> );
      pythonManager->addExporterFunction( field::exportGatherFunctions<FieldTypes> );
      pythonManager->addBlockDataConversion<FieldTypes>() ;

      // Blockforest
      pythonManager->addExporterFunction( blockforest::exportModuleToPython<Stencils> );

      // Geometry
      pythonManager->addExporterFunction( geometry::exportModuleToPython );

      // VTK
      pythonManager->addExporterFunction( vtk::exportModuleToPython );

      // Postprocessing
      pythonManager->addExporterFunction( postprocessing::exportModuleToPython<FieldTypesForMeshGeneration, FlagFieldTypes> );

      // Timeloop
      pythonManager->addExporterFunction( timeloop::exportModuleToPython );

#ifdef WALBERLA_BUILD_WITH_OPENMESH
      pythonManager->addExporterFunction( mesh::exportModuleToPython<FlagFieldTypes> );
#endif

#ifdef WALBERLA_BUILD_WITH_CUDA
      using walberla::cuda::GPUField;
      typedef bmpl::vector<GPUField<double>, GPUField<float>,
                           GPUField<int8_t>,  GPUField<int16_t>,  GPUField<int32_t>, GPUField<int64_t>,
                           GPUField<uint8_t>, GPUField<uint16_t>, GPUField<uint32_t>,GPUField<uint64_t> > GPUFields;

      pythonManager->addExporterFunction( cuda::exportModuleToPython<GPUFields, FieldTypes> );
      pythonManager->addBlockDataConversion<GPUFields>();
#endif

      python_coupling::initWalberlaForPythonModule();
   }
};
InitObject globalInitObject;

