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
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "waLBerlaDefinitions.h"
#ifdef WALBERLA_BUILD_WITH_PYTHON

#   include "core/math/Constants.h"

#   include "field/communication/PackInfo.h"

#   include "python_coupling/Manager.h"
#   include "python_coupling/export/BlockForestExport.h"
#   include "python_coupling/export/FieldExports.h"

namespace walberla {
    using flag_t = uint8_t;
    // clang-format off
    void exportDataStructuresToPython() {

        auto pythonManager = python_coupling::Manager::instance();

        // Field
        pythonManager->addExporterFunction(field::exportModuleToPython<Field<walberla::real_t, 1>, Field<walberla::real_t, 2>, Field<walberla::real_t, 3>, Field<walberla::real_t, 9>, Field<walberla::real_t, 19>, Field<walberla::real_t, 27>>);
        pythonManager->addExporterFunction(field::exportGatherFunctions<Field<walberla::real_t, 1>, Field<walberla::real_t, 2>, Field<walberla::real_t, 3>, Field<walberla::real_t, 9>, Field<walberla::real_t, 19>, Field<walberla::real_t, 27>>);
        pythonManager->addBlockDataConversion<Field<walberla::real_t, 1>, Field<walberla::real_t, 2>, Field<walberla::real_t, 3>, Field<walberla::real_t, 9>, Field<walberla::real_t, 19>, Field<walberla::real_t, 27>>();

        // Blockforest
        pythonManager->addExporterFunction(blockforest::exportModuleToPython<stencil::D2Q9, stencil::D3Q19, stencil::D3Q27>);
    }
   // clang-format on
}

#else

namespace walberla {
   void exportDataStructuresToPython() {}
} // namespace walberla

#endif