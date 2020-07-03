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
//! \file ScriptRunner.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Runs a python script that can access waLBerla variables
//
//======================================================================================================================

#pragma once

#include <string>

#include "waLBerlaDefinitions.h"
#include "core/DataTypes.h"



namespace walberla {
namespace python_coupling {

   class DictWrapper;

   class Shell
   {
   public:
      Shell( const std::string & prompt = "waLBerla");
      ~Shell();

      void operator() ();

      inline void run()           { (*this)();           }
      inline DictWrapper & data() { return *exposedVars_; }

   protected:
      bool getCompleteCommand( std::string & result );
      bool isCompleteCommand ( const std::string & code );

      shared_ptr<DictWrapper> exposedVars_;
      std::string prompt1_;
      std::string prompt2_;
   };


} // namespace python_coupling
} // namespace walberla




