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
//! \file Manager.cpp
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================



// Do not reorder includes - the include order is important
#include "PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "core/waLBerlaBuildInfo.h"
#include "Manager.h"
#include "core/logging/Logging.h"
#include "python_coupling/basic_exports/BasicExports.h"

#include <cstdlib>


BOOST_PYTHON_MODULE( walberla_cpp )
{
   using namespace walberla::python_coupling;
   auto manager =Manager::instance();
   exportBasicWalberlaDatastructures();
   manager->exportAll();
}

using namespace boost::python;


namespace walberla {
namespace python_coupling {

Manager::Manager()
   : initialized_( false )
{
}

Manager::~Manager( ) //NOLINT
{
   // To work reliably this would have to be called at the end of the
   // main function. At this position this leads to a segfault in some cases
   // Py_Finalize();
}

void Manager::addEntryToPythonPath( const std::string & path )
{
   if ( initialized_ ){
      entriesForPythonPath_.push_back( path );
   }
   else
   {
      addPath ( path );
   }
}

void Manager::addPath( const std::string & path )
{
   WALBERLA_ASSERT( initialized_ );

   object sys = import("sys");
   list sys_path = extract<list>( sys.attr("path") );
   sys_path.append(path);
}

void Manager::triggerInitialization()
{
   using namespace boost::python;

   if ( initialized_ ){
      return;
   }
   initialized_ = true;

   try
   {
#if PY_MAJOR_VERSION >= 3
      PyImport_AppendInittab( (char*)"walberla_cpp", PyInit_walberla_cpp );
#else
      PyImport_AppendInittab( (char*)"walberla_cpp", initwalberla_cpp );
#endif

      Py_Initialize();
      import("__main__");
      import("walberla_cpp");

      // Setup python path
      addPath( std::string(WALBERLA_SOURCE_DIR) + "/python" );

      const char * pPath = std::getenv( "WALBERLAPYTHONPATH" );
      if ( pPath )
         addPath( pPath );

      for ( auto it = entriesForPythonPath_.begin(); it != entriesForPythonPath_.end(); ++it )
         addPath( *it );
      entriesForPythonPath_.clear();

   }
   catch ( boost::python::error_already_set & ) {
      PyObject *type_ptr = nullptr;

      PyObject *value_ptr = nullptr;

      PyObject *traceback_ptr = nullptr;
      PyErr_Fetch(&type_ptr, &value_ptr, &traceback_ptr);

      if( type_ptr )
      {
         extract<std::string> type_str(( str( handle<>( type_ptr ) ) ));
         if( type_str.check() )
            WALBERLA_LOG_DEVEL( type_str() );
      }
      if(value_ptr)
      {
         extract<std::string> value_str(( str( handle<>( value_ptr ) ) ));
         if ( value_str.check() )
            WALBERLA_LOG_DEVEL( value_str() );
      }

      WALBERLA_ABORT( "Error while initializing Python" );
   }
}



void Manager::exportAll()
{
   for( auto it = exporterFunctions_.begin(); it != exporterFunctions_.end(); ++it ) {
      (*it)();
   }
}

boost::python::object Manager::pythonObjectFromBlockData( IBlock & block, BlockDataID id )
{
   if( block.isDataOfType< boost::python::object > ( id )  )
      return *block.getData< boost::python::object > ( id );


   for( auto it = blockDataToObjectFunctions_.begin(); it != blockDataToObjectFunctions_.end(); ++it )
   {
      auto res = (*it)( block, id );
      if ( res != boost::python::object() )
         return res;
   }

   return boost::python::object();
}


} // namespace python_coupling
} // namespace walberla

#else

int someSymbolSoThatLinkerDoesNotComplain=0;

#endif



