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
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

namespace walberla {
namespace python_coupling {

   class DictWrapper;
   class PythonCallback;


   class TimeloopIntercept
   {
   public:
      TimeloopIntercept( const std::string & functionName,
                         uint_t interval = 1, bool enableSignalInterrupt=true );

      TimeloopIntercept( const std::string & pythonFile, const std::string & functionName,
                         uint_t interval = 1, bool enableSignalInterrupt=true  );

      inline shared_ptr<PythonCallback> callback() { return callback_;         }

      void operator() ();

   protected:
      shared_ptr<PythonCallback> callback_;
      uint_t timestep_;
      uint_t interval_;
      bool   enableSignalInterrupt_;
   };


} // namespace python_coupling
} // namespace walberla





