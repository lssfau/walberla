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
//! \file ConfigFromDict.cpp
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "waLBerlaDefinitions.h"
#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "ConfigFromDict.h"
#include "python_coupling/PythonWrapper.h"

#include "core/logging/Logging.h"
#include "core/config/Config.h"

namespace walberla {
namespace python_coupling {


void handlePythonBooleans( std::string & value) {
   if ( value=="True")  value="1";
   if ( value=="False") value="0";
}

namespace py = pybind11;

void configFromPythonDict( config::Config::Block & block, py::dict & pythonDict )
{
   for (auto item : pythonDict)
   {
      // py::print(item);
      if( py::isinstance<std::string>(item.first) ) {
         WALBERLA_LOG_WARNING( "Detected non-string key in waLBerla configuration" );
         continue;
      }

      std::string key = py::str(item.first);

      try
      {
         if( py::isinstance<std::string>(item.second) ){
            std::string value = py::str(item.second);
            handlePythonBooleans( value );
            block.addParameter( key, value );
         }
         else if ( py::isinstance<py::dict>(item.second) )
         {
            walberla::config::Config::Block & childBlock = block.createBlock( key );
            py::dict childDict = py::dict(pythonDict[key.c_str()]);
            configFromPythonDict( childBlock, childDict);
         }
         else if ( py::isinstance<py::list>(item.second) )
         {
            py::list childList = py::list(pythonDict[key.c_str()]);
            for(py::size_t i = 0; i < childList.size(); ++i){
               walberla::config::Config::Block & childBlock = block.createBlock( key );
               py::dict d = py::dict(childList[i]);
               configFromPythonDict( childBlock, d );
            }
         }
         else if (  py::isinstance<py::tuple>(item.second) )
         {
            std::stringstream ss;
            py::tuple childTuple = py::tuple(pythonDict[key.c_str()]);

            WALBERLA_ASSERT(len(childTuple) == 2 || len(childTuple) == 3,
                            "Config problem: " << key << ": Python tuples are mapped to walberla::Vector2 or Vector3. \n" <<
                               "So only tuples of size 2 or 3 are supported! Option " << key << "  is ignored ")

            if ( len(childTuple) == 2 )
            {
               std::string e0 = py::str( childTuple[0].attr("__str__" )() );
               std::string e1 = py::str( childTuple[1].attr("__str__" )() );
               handlePythonBooleans( e0 );
               handlePythonBooleans( e1 );
               ss << "< " << e0 << " , " << e1 << " > ";
               block.addParameter( key, ss.str() );
            }
            else if ( len(childTuple) == 3)
            {
               std::string e0 = py::str( childTuple[0].attr("__str__" )() );
               std::string e1 = py::str( childTuple[1].attr("__str__" )() );
               std::string e2 = py::str( childTuple[2].attr("__str__" )() );
               handlePythonBooleans( e0 );
               handlePythonBooleans( e1 );
               handlePythonBooleans( e2 );
               ss << "< " << e0 << ", " << e1 << ", " << e2 << " > ";
               block.addParameter( key, ss.str() );
            }
         }
         else
         {
            // if value is not a string try to convert it
            std::string value = py::str( pythonDict[key.c_str()].attr("__str__" )() );
            block.addParameter ( key, value );
         }
      }
      catch ( py::error_already_set & ) {
         WALBERLA_LOG_WARNING ( "Error when reading configuration option " << key << ". Could not be converted to string.");
      }
   }
}


shared_ptr<Config> configFromPythonDict( py::dict & pythonDict )
{
   shared_ptr<Config> config = make_shared<Config>();
   configFromPythonDict( config->getWritableGlobalBlock(), pythonDict );
   return config;
}



} // namespace python_coupling
} // namespace walberla

#endif
