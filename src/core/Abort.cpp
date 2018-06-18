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
//! \file Abort.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Abort.h"
#include "core/debug/PrintStacktrace.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"


#ifdef _MSC_VER
#  include <windows.h>
#endif



namespace walberla {



void Abort::abort( const std::string & logMessage, const std::string & callerPath, const int line )
{
   if( abortFunction_ )
      abortFunction_( logMessage, callerPath, line, true );

   defaultAbort( logMessage, callerPath, line, true ); // only called if no abort function is set or if the abort function that was provided returns (which is not supposed to happen!)
}


void Abort::abortNoDebugInfo( const std::string & logMessage,  const std::string & callerPath, const int line )
{
   if( abortFunction_ )
      abortFunction_( logMessage, callerPath, line, false );

   defaultAbort( logMessage, callerPath, line, false );// only called if no abort function is set or if the abort function that was provided returns (which is not supposed to happen!)
}


void Abort::defaultAbort( const std::string & logMessage, const std::string & callerPath, const int line, bool withDebugInfo )
{
   std::ostringstream ss;
   ss << logMessage << "\n\nFatal error came from " << callerPath << ":" << line << "\nAborting now ...\n\nStack backtrace:\n";
   if ( withDebugInfo )
      debug::printStacktrace( ss );
   ss << std::endl;


#ifdef _OPENMP
   #pragma omp critical (Abort_abort)
   {
#endif

   if( logging::Logging::isInstantiated() )
   {
      logging::Logging::instance()->logCallerPath(true);
      logging::Logging::instance()->logError( ss.str(), callerPath, line );
   }
   else
      std::cerr << ss.str() << std::endl;

#ifdef _OPENMP
   }
#endif

#ifdef _MSC_VER
   if( IsDebuggerPresent() )
      __debugbreak(); // Make the MSVC++ Debugger stop here if it is attached
#endif

   MPIManager::instance()->abort();
}



void Abort::exceptionAbort( const std::string & message, const std::string & /*callerPath*/, const int /*line*/, bool /*withDebugInfo*/ )
{
   throw std::runtime_error( message );
}



} // namespace walberla
