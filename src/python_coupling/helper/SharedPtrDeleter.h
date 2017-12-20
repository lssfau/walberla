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
//! \file SharedPtrDeleter.h
//! \ingroup python_export
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "python_coupling/PythonWrapper.h"

namespace walberla {
namespace python_coupling {


namespace internal
{

template<typename T>
class SharedPtrDeleterTiedToPythonObject
{
public:
   SharedPtrDeleterTiedToPythonObject( PyObject *object ) : object_( object )
   {
   }

   void operator()( T * )
   {
       Py_DECREF( object_ );
   }

private:
   PyObject *object_;
};

} // namespace internal


template<typename T>
shared_ptr<T> createSharedPtrFromPythonObject(boost::python::object pythonObject) {
    T * ptr = boost::python::extract<T*>( pythonObject);
    auto deleter = internal::SharedPtrDeleterTiedToPythonObject<T>(pythonObject.ptr());
    Py_INCREF(pythonObject.ptr());
    return shared_ptr<T>(ptr, deleter);
}


} // namespace python_coupling
} // namespace walberla


