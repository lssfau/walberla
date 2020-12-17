
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
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "python_coupling/PythonWrapper.h"

namespace walberla {
namespace python_coupling {

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

// fallback - check for bool
template<typename T>
inline static bool isCppEqualToPythonType( std::string n )
{
   return ( n == "numpy.bool_" || n =="bool" );
}

// native data types
template<> inline bool isCppEqualToPythonType<float> ( std::string n) {return ( n == "numpy.float32" || n =="float32" );}
template<> inline bool isCppEqualToPythonType<double>( std::string n) {return ( n == "numpy.float64" || n == "numpy.float_" || n =="float64" || n=="float");}


template<> inline bool isCppEqualToPythonType<uint8_t  > ( std::string n) {return ( n == "numpy.uint8" || n =="uint8"  );}
template<> inline bool isCppEqualToPythonType<uint16_t > ( std::string n) {return ( n == "numpy.uint16"|| n =="uint16" );}
template<> inline bool isCppEqualToPythonType<uint32_t > ( std::string n) {return ( n == "numpy.uint32"|| n =="uint32" );}
template<> inline bool isCppEqualToPythonType<uint64_t > ( std::string n) {return ( n == "numpy.uint64"|| n =="uint64" );}


template<> inline bool isCppEqualToPythonType<int8_t  > ( std::string n) {return ( n == "numpy.int8"  || n =="int8");}
template<> inline bool isCppEqualToPythonType<int16_t > ( std::string n) {return ( n == "numpy.int16" || n =="int16");}
template<> inline bool isCppEqualToPythonType<int32_t > ( std::string n) {return ( n == "numpy.int32" || n =="int32");}
template<> inline bool isCppEqualToPythonType<int64_t > ( std::string n) {return ( n == "numpy.int64" || n =="int64" || n == "int" );}

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif

} // namespace python_coupling
} // namespace walberla
