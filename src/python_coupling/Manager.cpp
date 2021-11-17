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
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================



// Do not reorder includes - the include order is important
#include "PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#   include "core/logging/Logging.h"
#   include "core/waLBerlaBuildInfo.h"

#   include "python_coupling/export/BasicExport.h"

#   include <pybind11/embed.h>

#   include "Manager.h"

PYBIND11_MODULE( walberla_cpp, m)
{
   using namespace walberla::python_coupling;
   auto manager = Manager::instance();
   exportBasicWalberlaDatastructures(m);
   manager->exportAll(m);
}

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
   // py::finalize_interpreter();
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

   py::object sys = py::module::import("sys");
   py::list sys_path = py::list( sys.attr("path") );
   sys_path.append(path);
}

void Manager::triggerInitialization()
{
   if ( initialized_ ){
      return;
   }
   initialized_ = true;

   try
   {
      // The python module is used as embedded module here. There is a pybind11 macro for that called
      // PYBIND11_EMBEDDED_MODULE. However it can not be used here since we want a shared lib for the python coupling.
      // With the C-Call the so is embedded here.
      PyImport_AppendInittab( (char*)"walberla_cpp", PyInit_walberla_cpp );

      py::initialize_interpreter();
      py::module::import("__main__");
      py::module::import("walberla_cpp");


      // Setup python path
      addPath( std::string(WALBERLA_SOURCE_DIR) + "/python" );

      const char * pPath = std::getenv( "WALBERLAPYTHONPATH" );
      if ( pPath )
         addPath( pPath );

      for ( auto it = entriesForPythonPath_.begin(); it != entriesForPythonPath_.end(); ++it )
         addPath( *it );
      entriesForPythonPath_.clear();

   }
   catch ( py::error_already_set & e ) {
      if (e.matches(PyExc_ModuleNotFoundError)){
         py::print("The module walberla_cpp could not be found");
      }
      else {
         py::print("Unexpected Exception");
         throw;
      }
      WALBERLA_ABORT( "Error while initializing Python" );
   }
}



void Manager::exportAll(py::module_ &m)
{
   for( auto it = exporterFunctions_.begin(); it != exporterFunctions_.end(); ++it ) {
      (*it)(m);
   }
}

py::object Manager::pythonObjectFromBlockData( IBlock & block, BlockDataID id )
{
   if( block.isDataOfType< py::object > ( id )  ){
      return *block.getData< py::object > ( id );}


   for( auto it = blockDataToObjectFunctions_.begin(); it != blockDataToObjectFunctions_.end(); ++it )
   {
      auto res = (*it)( block, id );
      if ( !res.is(py::object()) ){
         return res;}
   }

   return py::object();
}


} // namespace python_coupling
} // namespace walberla

#else

int someSymbolSoThatLinkerDoesNotComplain=0;

#endif



