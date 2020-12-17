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
//! \file CUDAExport.h
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#ifdef WALBERLA_BUILD_WITH_PYTHON


namespace walberla {
namespace cuda {


   template<typename... GpuFields>
   void exportModuleToPython(py::module_ &m);

   template<typename... CpuFields>
   void exportCopyFunctionsToPython(py::module_ &m);


} // namespace cuda
} // namespace walberla

#include "CUDAExport.impl.h"

#endif //WALBERLA_BUILD_WITH_PYTHON
