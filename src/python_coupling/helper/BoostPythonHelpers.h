
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
//! \file CppPythonTypeEquality.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "python_coupling/PythonWrapper.h"

#include <boost/python/converter/registry.hpp>


namespace walberla {
namespace python_coupling {

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

   // fallback - check boost python registry
   template<typename T>
   inline static bool isCppEqualToPythonType( PyTypeObject * obj)
   {
      boost::python::type_info info = boost::python::type_id<T>();
      const boost::python::converter::registration* reg = boost::python::converter::registry::query(info);
      if (reg == NULL)
         return false;

      try
      {
         reg->get_class_object();
         return ( reg->get_class_object() == obj  );
      }
      catch ( ... ) {
         return false;
      }
   }

   // native data types
   template<> inline bool isCppEqualToPythonType<bool>  ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.bool_"
                                                                                                          || n =="bool" );           }


   template<> inline bool isCppEqualToPythonType<float> ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.float32" ); }
   template<> bool inline isCppEqualToPythonType<double>( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.float64"
                                                                                                          || n == "numpy.float_"
                                                                                                          || n =="float");           }


   template<> inline bool isCppEqualToPythonType<uint8_t  > ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.uint8"  ); }
   template<> inline bool isCppEqualToPythonType<uint16_t > ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.uint16" ); }
   template<> inline bool isCppEqualToPythonType<uint32_t > ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.uint32" ); }
   template<> inline bool isCppEqualToPythonType<uint64_t > ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.uint64" ); }


   template<> inline bool isCppEqualToPythonType<int8_t  > ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.int8"  ); }
   template<> inline bool isCppEqualToPythonType<int16_t > ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.int16" ); }
   template<> inline bool isCppEqualToPythonType<int32_t > ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.int32" ); }
   template<> inline bool isCppEqualToPythonType<int64_t > ( PyTypeObject * o) { std::string n( o->tp_name ); return ( n == "numpy.int64"
                                                                                                             || n == "int" );         }

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif


   template< typename T >
   bool isTypeRegisteredInBoostPython( )
   {
      boost::python::type_info info = boost::python::type_id<T>();
      const boost::python::converter::registration* reg = boost::python::converter::registry::query(info);

      try {
         reg->get_class_object();
      }
      catch( ... ) {
         PyErr_Clear();
         return false;
      }
      PyErr_Clear();
      return (reg != NULL);
   }

} // namespace python_coupling
} // namespace walberla


