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
//! \file Tracing.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Logging.h"
#include "core/logging/CMakeDefs.h"

#include <sstream>
#include <string>


//**********************************************************************************************************************
/*!
*   \brief Macro enabling tracing for a specific function
*
*   For details see documentation of class Tracer.
*/
//**********************************************************************************************************************
#ifdef WALBERLA_LOGLEVEL_TRACING
#  define WALBERLA_TRACE_IN walberla::logging::Tracer walberlaTracingObject( __FUNCTION__, __FILE__, __LINE__ )
#else
#  define WALBERLA_TRACE_IN (void(0))
#endif


namespace walberla {
namespace logging {


//**********************************************************************************************************************
/*!
*   \brief Tracer printing a log message upon creation and destruction
*
*   This class is intended to be used via the WALBERLA_TRACE_IN macro. Place "WALBERLA_TRACE_IN;" at the beginning of a
*   function to have a log message printed when control flow enters the function and another one when it leaves the
*   function's scope.
*   Note that messages are printed only if "WALBERLA_LOGLEVEL" is set to "TRACING" in CMake cache and the runtime log
*   level of Logging is set to TRACING!
*/
//**********************************************************************************************************************
class Tracer
{
public:
   Tracer( const std::string & functionName, const std::string & fileName, int lineNumber ) :
      functionName_(functionName), fileName_(fileName), lineNumber_(lineNumber)
   {
      Logging::instance()->logTracing( std::string("Entering ") + functionName_, fileName_, lineNumber_ );
   }

   ~Tracer()
   {
      Logging::instance()->logTracing( std::string("Leaving ") + functionName_, fileName_, lineNumber_ );
   }
private:
   std::string functionName_;
   std::string fileName_;
   int         lineNumber_;
};


} // namespace logging
} // namespace walberla
