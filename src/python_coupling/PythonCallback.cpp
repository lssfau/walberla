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
//! \file PythonCallback.cpp
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "PythonCallback.h"
#include "DictWrapper.h"
#include "core/logging/all.h"
#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "Manager.h"
#include "core/Abort.h"
#include "core/Filesystem.h"

#include "pybind11/eval.h"

namespace walberla {
namespace python_coupling {

   static py::object importModuleOrFileInternal( const std::string & fileOrModuleName, const std::vector< std::string > & argv )
   {
      auto manager = python_coupling::Manager::instance();
      manager->triggerInitialization();

      namespace py = pybind11;

      std::string moduleName = fileOrModuleName;

      std::stringstream code;
      code << "import sys" << "\n" ;
      code << "sys.argv = [ ";
      for( auto argStrIt = argv.begin(); argStrIt != argv.end(); ++argStrIt )
         code << "'" << *argStrIt  << "',";
      code << "] \n";

      filesystem::path path ( fileOrModuleName );
      path = filesystem::absolute( path );

      if ( path.extension() == ".py" )
      {
         moduleName = path.stem().string();
         if ( ! path.parent_path().empty() ) {
            std::string p = filesystem::canonical(path.parent_path()).string();
            code << "sys.path.append( r'" << p << "')" << "\n";
         }
      }

      py::exec( code.str().c_str(), py::module::import("__main__").attr("__dict__") );

      try {
         return py::module::import( moduleName.c_str() );
      }
      catch ( py::error_already_set &e) {
         throw py::value_error(e.what());
      }
   }

   void importModuleOrFile( const std::string & fileOrModuleName, const std::vector< std::string > & argv )
   {
      importModuleOrFileInternal( fileOrModuleName, argv );
   }

   // initializes invalid/empty callback
   PythonCallback::PythonCallback()
      : exposedVars_( new DictWrapper() ), callbackDict_( new DictWrapper() )
   {
      Manager::instance()->triggerInitialization();
      callbackDict_->dict() = py::dict();
   }


   PythonCallback::PythonCallback( const std::string & fileOrModuleName, const std::string & functionName, const std::vector<std::string> & argv )
      : functionName_( functionName ), exposedVars_( new DictWrapper() ), callbackDict_( new DictWrapper() )
   {
      Manager::instance()->triggerInitialization();

      namespace py = pybind11;

      // Add empty callbacks module
      importModuleOrFileInternal( fileOrModuleName, argv );
      py::object callbackModule = py::module::import( "walberla_cpp.callbacks");

      callbackDict_->dict() = py::dict( callbackModule.attr( "__dict__" ) );
   }


   PythonCallback::PythonCallback( const std::string & functionName )
      : functionName_( functionName ), exposedVars_( new DictWrapper() ), callbackDict_( new DictWrapper() )
   {
      Manager::instance()->triggerInitialization();

      // Add empty callbacks module
      py::object callbackModule = py::module::import( "walberla_cpp.callbacks");

      callbackDict_->dict() = py::dict( callbackModule.attr( "__dict__" ) );
   }

   bool PythonCallback::isCallable() const
   {
      return callbackDict_->dict().contains( functionName_ );
   }

   void PythonCallback::operator() ()
   {
      if ( ! isCallable() )
         WALBERLA_ABORT_NO_DEBUG_INFO( "Could not call python function '" << functionName_ << "'. " <<
                                             "Did you forget to set the callback function?" )

      namespace py = pybind11;

      try
      {
         py::object function = callbackDict_->dict()[ functionName_.c_str() ];

         py::object returnVal;
         returnVal = function( *py::tuple(), **(exposedVars_->dict() ) );

         exposedVars_->dict()["returnValue"] = returnVal;
      }
      catch ( py::error_already_set &e ) {
         throw py::value_error(e.what());
      }
   }


} // namespace python_coupling
} // namespace walberla

#else

namespace walberla {
namespace python_coupling {

   void importModuleOrFile( const std::string &, const std::vector< std::string > & ) { }

   PythonCallback::PythonCallback( const std::string & functionName )
   : functionName_( functionName ), exposedVars_( new DictWrapper() ), callbackDict_( new DictWrapper() ) {}


   PythonCallback::PythonCallback( const std::string & functionName, const std::string & , const std::vector<std::string> &  )
   : functionName_( functionName ), exposedVars_( new DictWrapper() ), callbackDict_( new DictWrapper() ) {}

   bool PythonCallback::isCallable() const { return false; }
   void PythonCallback::operator() ()      {  }


} // namespace python_coupling
} // namespace walberla



#endif
