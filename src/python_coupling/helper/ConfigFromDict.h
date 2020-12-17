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
//! \file ConfigFromDict.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once


#ifdef WALBERLA_BUILD_WITH_PYTHON


#include "core/config/Config.h"
#include "python_coupling/DictWrapper.h"


namespace walberla {
namespace python_coupling {


   //*******************************************************************************************************************
   /*! Converts a python dictionary to a config::Block (recursively)
   */
   //*******************************************************************************************************************
   void configFromPythonDict( config::Config::Block & result, py::dict & pythonDict  );



   //*******************************************************************************************************************
   /*! Converts a python dictionary to a waLBerla config object
   */
   //*******************************************************************************************************************
   shared_ptr<Config> configFromPythonDict( py::dict & pythonDict );



} // namespace python_coupling
} // namespace walberla


#endif
