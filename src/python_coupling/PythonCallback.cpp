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
//
//======================================================================================================================

#include "PythonCallback.h"
#include "PythonWrapper.h"
#include "DictWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "Manager.h"
#include "core/Abort.h"
#include "core/logging/Logging.h"
#include "helper/ExceptionHandling.h"
#include "core/Filesystem.h"

#if defined(__GLIBCXX__) && __GLIBCXX__ <= 20160609
#define CURRENT_PATH_WORKAROUND
#include <unistd.h>
#include <errno.h>
#endif

namespace walberla {
namespace python_coupling {

   static boost::python::object importModuleOrFileInternal( const std::string & fileOrModuleName, const std::vector< std::string > & argv )
   {
      auto manager = python_coupling::Manager::instance();
      manager->triggerInitialization();

      namespace bp = boost::python;

      std::string moduleName = fileOrModuleName;

      std::stringstream code;
      code << "import sys" << "\n" ;
      code << "sys.argv = [ ";
      for( auto argStrIt = argv.begin(); argStrIt != argv.end(); ++argStrIt )
         code << "'" << *argStrIt  << "',";
      code << "] \n";

      filesystem::path path ( fileOrModuleName );
#ifdef CURRENT_PATH_WORKAROUND
      // workaround for double free in filesystem::current_path in libstdc++ 5.4 and lower
      size_t cwd_size = 16;
      char * cwd_buf = (char*) std::malloc(cwd_size * sizeof(char));

      while( getcwd( cwd_buf, cwd_size ) == NULL )
      {
         if (errno == ERANGE)
         {
            cwd_size *= 2;
            cwd_buf = (char*) std::realloc( cwd_buf, cwd_size * sizeof(char) );
         }
         else
         {
            python_coupling::terminateOnPythonException( std::string("Could not determine working directory") );
         }
      }

      std::string cwd(cwd_buf);
      std::free(cwd_buf);
      path = filesystem::absolute( path, cwd );
#else
      path = filesystem::absolute( path );
#endif
      if ( path.extension() == ".py" )
      {
         moduleName = path.stem().string();


         if ( ! path.parent_path().empty() )  {
#ifdef CURRENT_PATH_WORKAROUND
            std::string p = filesystem::canonical(path.parent_path(), cwd).string();
#else
            std::string p = filesystem::canonical(path.parent_path()).string();
#endif
            code << "sys.path.append( r'" << p << "')" << "\n";
         }
      }
      bp::exec( code.str().c_str(), bp::import("__main__").attr("__dict__") );

      try {
         return bp::import( moduleName.c_str() );
      }
      catch ( bp::error_already_set & ) {
         python_coupling::terminateOnPythonException( std::string("Python Error while loading ") + fileOrModuleName );
         return boost::python::object();
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
      callbackDict_->dict() = boost::python::dict();
   }


   PythonCallback::PythonCallback( const std::string & fileOrModuleName, const std::string & functionName, const std::vector<std::string> & argv )
      : functionName_( functionName ), exposedVars_( new DictWrapper() ), callbackDict_( new DictWrapper() )
   {
      Manager::instance()->triggerInitialization();

      using namespace boost::python;

      // Add empty callbacks module
      importModuleOrFileInternal( fileOrModuleName, argv );
      object callbackModule = import( "walberla_cpp.callbacks");

      callbackDict_->dict() = extract<dict>( callbackModule.attr( "__dict__" ) );
   }


   PythonCallback::PythonCallback( const std::string & functionName )
      : functionName_( functionName ), exposedVars_( new DictWrapper() ), callbackDict_( new DictWrapper() )
   {
      Manager::instance()->triggerInitialization();

      using namespace boost::python;

      // Add empty callbacks module
      object callbackModule = import( "walberla_cpp.callbacks");

      callbackDict_->dict() = extract<dict>( callbackModule.attr( "__dict__" ) );
   }

   bool PythonCallback::isCallable() const
   {
      return callbackDict_->dict().has_key( functionName_ );
   }

   void PythonCallback::operator() ()
   {
      if ( ! isCallable() )
         WALBERLA_ABORT_NO_DEBUG_INFO( "Could not call python function '" << functionName_ << "'. " <<
                                        "Did you forget to set the callback function?" );

      namespace bp = boost::python;

      try
      {
         if ( exposedVars_->dict().has_key("returnValue"))
            bp::api::delitem( exposedVars_->dict(), "returnValue" );

         bp::object function = callbackDict_->dict()[ functionName_ ];

         bp::object returnVal;
         returnVal = function( *bp::tuple(), **(exposedVars_->dict() ) );

         exposedVars_->dict()["returnValue"] = returnVal;
      }
      catch ( bp::error_already_set & ) {
         python_coupling::terminateOnPythonException( std::string("Error while running Python function ") + functionName_ );
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
