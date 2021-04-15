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
//! \file PythonCallback.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include <string>
#include <vector>


namespace walberla {
namespace python_coupling {

   class DictWrapper;


   void importModuleOrFile( const std::string & moduleOrFile,
                            const std::vector< std::string > & argv = std::vector<std::string>() );

   //*******************************************************************************************************************
   /*! Run a Python function from C++ code
   *
   * \ingroup python_coupling
   *
   * Example:
   *  \code
        python_coupling::PythonCallback callback ( "atEndOfTimestep" );
        callback.data().exposePtr("blockStorage", blocks ); // blockStorage is available as global variable in script,
                                                            // the blockStorage can be modified from the python script

        callback.data().exposeValue("someInteger", 42);   // expose copy of variable to script
        callback();                                       // run 'someFunction(data)' in script
                                                          // where data is a python dict that contains your exposed data


        // Access data that was modified in script
        int i = 0;
        if ( callback.data().has<int>("integerSetInScript") ) {
           i = callback.data().get<int>( "integerSetInScript" );
        }


   *  \endcode
   *
   *  This example calls the python function "walberla_cpp.callbacks.atEndOfTimestep", which can be either
   *  set directly or using a decorator:
   *  \code
   *     @waLBerla.callback( "atEndOfTimestep" )
   *  \endcode
   *
   *
   *
   *  The result of the function call will be stored in the data dictionary with the key "returnValue".
   *
   *  Depending on how many parameters the python function expects it is called in different ways:
   *     - no parameters: the data is not passed
   *     - multiple parameters: the data dictionary is unpacked i.e. yourFunction( **data )
   */
   //*******************************************************************************************************************
   class PythonCallback
   {
   public:
      PythonCallback();
      PythonCallback( const std::string & functionName );
      PythonCallback( const std::string & moduleOrFile,
                             const std::string & functionName,
                             const std::vector<std::string> & argv = std::vector<std::string>() );

            DictWrapper & data()       { return *exposedVars_; }
      const DictWrapper & data() const { return *exposedVars_; }

      bool isCallable() const;
      void operator() ();

   protected:
      const std::string functionName_;
      shared_ptr<DictWrapper> exposedVars_;
      shared_ptr<DictWrapper> callbackDict_;
   };



} // namespace python_coupling
} // namespace walberla


