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
//! \file ExceptionDecode.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include <string>


namespace walberla {
namespace python_coupling {

   // Call this function when a boost::python::already_set exception was caught
   // returns formatted error string with traceback
   inline std::string decodeException()
   {
       namespace bp = boost::python;

       PyObject *exc,*val,*tb;

       bp::object formatted_list, formatted;
       PyErr_Fetch(&exc,&val,&tb);
       PyErr_NormalizeException(&exc, &val, &tb);
       bp::handle<> hexc( exc );
       bp::handle<> hval( bp::allow_null(val) );
       bp::handle<> htb ( bp::allow_null(tb)  );
       bp::object traceback( bp::import("traceback"));

       if (!tb) {
           bp::object format_exception_only( traceback.attr("format_exception_only"));
           formatted_list = format_exception_only(hexc,hval);
       } else {
          bp::object format_exception(traceback.attr("format_exception"));
          formatted_list = format_exception(hexc,hval,htb);
       }
       formatted = bp::str("").join(formatted_list);
       return bp::extract<std::string>(formatted);
   }


   inline void terminateOnPythonException( const std::string message )
   {
      if (PyErr_Occurred()) {
          std::string decodedException = decodeException();
          WALBERLA_ABORT_NO_DEBUG_INFO( message << "\n\n" << decodedException  );
      }
      WALBERLA_ABORT_NO_DEBUG_INFO( message << " (unable to decode Python exception) " );
   }


} // namespace python_coupling
} // namespace walberla


