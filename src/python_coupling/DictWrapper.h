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
//! \file DictWrapper.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Wrapper to store and extract values from pybind11
//
//! \warning: if you include this header you also have to include Python.h as first header in your
//!           cpp file
//
//======================================================================================================================

#pragma once

#include "PythonWrapper.h"
#include "core/DataTypes.h"
#include "core/Abort.h"

#include <functional>


namespace walberla {
namespace python_coupling {
   class DictWrapper
   {
   public:

      //** Expose Data *************************************************************************************************
      /*! \name Expose Data */
      //@{
      template<typename T>  inline void exposePtr(const char* name, T * var );
      template<typename T>  inline void exposePtr(const char* name, const shared_ptr<T> & var );
      template<typename T>  void exposeValue     ( const char* name, const T & var );
      //@}
      //****************************************************************************************************************


      //** Get Data  ***************************************************************************************************
      /*! \name Get Data */
      //@{
      template<typename T> inline T    get( const char* name );
      template<typename T> inline bool has( const char* name );
      template<typename T> inline bool checkedGet( const char* name, T output );
      //@}
      //****************************************************************************************************************


#ifdef WALBERLA_BUILD_WITH_PYTHON
            pybind11::dict & dict()        { return d_; }
      const pybind11::dict & dict() const  { return d_; }
   protected:
            pybind11::dict d_;
#endif
   };


} // namespace python_coupling
} // namespace walberla



#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "DictWrapper.impl.h"

#else

// Stubs when Python is not built
namespace walberla {
namespace python_coupling {

   template<typename T> void DictWrapper::exposePtr( const std::string & , T *  ) {}
   template<typename T> void DictWrapper::exposePtr( const std::string & , const shared_ptr<T> & ) {}
   template<typename T> void DictWrapper::exposeValue( const std::string & , const T &  ) {}

   template<typename T> bool DictWrapper::has( const std::string &  )                      { return false;  }
   template<typename T> bool DictWrapper::checkedGet( const std::string & name, T output ) { return false; }

   template<typename T> T DictWrapper::get( const std::string & ) {
      WALBERLA_ABORT("Not available - waLBerla was built without Python support");
#ifdef __IBMCPP__
      return *(reinterpret_cast< T * >( NULL )); // silencing incorrect IBM compiler warning
#endif
   }


} // namespace python_coupling
} // namespace walberla

#endif


