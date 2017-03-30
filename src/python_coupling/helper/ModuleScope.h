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
//! \file ModuleScope.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "python_coupling/PythonWrapper.h"

#include <string>

namespace walberla {
namespace python_coupling {


   class ModuleScope : public boost::python::scope
   {
   public:
      ModuleScope( const std::string & name)
         : boost::python::scope ( ModuleScope::createNew( name ) )
      {}

   private:
      static boost::python::object createNew ( const std::string & name )
      {
         using namespace boost::python;
         object module( handle<>( borrowed(PyImport_AddModule( name.c_str() ) ) ) );
         scope().attr( name.c_str() ) = module;
         return module;
      }

   };




} // namespace python_coupling
} // namespace walberla


