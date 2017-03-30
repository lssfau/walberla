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

void configFromPythonDict( config::Config::Block & block, boost::python::dict & pythonDict )
{
   using namespace boost::python;

   boost::python::list keys = pythonDict.keys();
   keys.sort();

   for (int i = 0; i < boost::python::len( keys ); ++i)
   {
      // Extract key
      boost::python::extract<std::string> extracted_key( keys[i] );
      if( !extracted_key.check() ) {
         WALBERLA_LOG_WARNING( "Detected non-string key in waLBerla configuration" );
         continue;
      }
      std::string key = extracted_key;

      // Extract value
      extract<std::string>  extracted_str_val   ( pythonDict[key] );
      extract<dict>         extracted_dict_val  ( pythonDict[key] );
      extract<list>         extracted_list_val  ( pythonDict[key] );
      extract<tuple>        extracted_tuple_val ( pythonDict[key] );

      try
      {
         if( extracted_str_val.check() ){
            std::string value = extracted_str_val;
            handlePythonBooleans( value );
            block.addParameter( key, value );
         }
         else if ( extracted_dict_val.check() )
         {
            walberla::config::Config::Block & childBlock = block.createBlock( key );
            dict childDict = extracted_dict_val;
            configFromPythonDict( childBlock, childDict );
         }
         else if ( extracted_list_val.check() )
         {
            list childList = extracted_list_val;
            for( int l=0; l < len( childList ); ++l ) {
               walberla::config::Config::Block & childBlock = block.createBlock( key );
               dict d = extract<dict>( childList[l] );
               configFromPythonDict( childBlock, d );
            }
         }
         else if (  extracted_tuple_val.check() )
         {
            std::stringstream ss;
            tuple childTuple = extracted_tuple_val;
            if ( len(childTuple) == 2 )
            {
               std::string e0 = extract<std::string>( childTuple[0].attr("__str__" )() );
               std::string e1 = extract<std::string>( childTuple[1].attr("__str__" )() );
               handlePythonBooleans( e0 );
               handlePythonBooleans( e1 );
               ss << "< " << e0 << " , " << e1 << " > ";
               block.addParameter( key, ss.str() );
            }
            else if ( len(childTuple) == 3)
            {
               std::string e0 = extract<std::string>( childTuple[0].attr("__str__" )() );
               std::string e1 = extract<std::string>( childTuple[1].attr("__str__" )() );
               std::string e2 = extract<std::string>( childTuple[2].attr("__str__" )() );
               handlePythonBooleans( e0 );
               handlePythonBooleans( e1 );
               handlePythonBooleans( e2 );
               ss << "< " << e0 << " , " << e1 << ", " << e2 << " > ";
               block.addParameter( key, ss.str() );
            }
            else
            {
               WALBERLA_LOG_WARNING( "Config problem: " << key << ": Python tuples are mapped to walberla::Vector2 or Vector3. \n" <<
                                     "So only tuples of size 2 or 3 are supported! Option " << key << "  is ignored ");

            }
         }
         else
         {
            // if value is not a string try to convert it
            std::string value = extract<std::string>( pythonDict[key].attr("__str__" )() );
            block.addParameter ( key, value );
         }
      }
      catch ( error_already_set & ) {
         WALBERLA_LOG_WARNING ( "Error when reading configuration option " << key << ". Could not be converted to string.");
      }
   }
}


shared_ptr<Config> configFromPythonDict( boost::python::dict & pythonDict )
{
   using namespace boost::python;

   shared_ptr<Config> config = make_shared<Config>();
   configFromPythonDict( config->getWritableGlobalBlock(), pythonDict );
   return config;
}



} // namespace python_coupling
} // namespace walberla

#endif
