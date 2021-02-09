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
//! \file Abort.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/singleton/Singleton.h"

#include <functional>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>


namespace walberla {



class Abort : public singleton::Singleton<Abort>
{
   WALBERLA_BEFRIEND_SINGLETON;

public:

   typedef std::function<void ( const std::string & message, const std::string & callerPath, const int line, bool withDebugInfo  )> AbortFunction;

   void resetAbortFunction( const AbortFunction & function = AbortFunction() ) { abortFunction_ = function; }

   void abort( const std::string & logMessage, const std::string & callerPath, const int line );
   void abortNoDebugInfo( const std::string & logMessage, const std::string & callerPath, const int line );

   static void   defaultAbort( const std::string & message, const std::string & callerPath, const int line, bool withDebugInfo );
   static void exceptionAbort( const std::string & message, const std::string & callerPath, const int line, bool withDebugInfo );

private:

   Abort() = default;

   AbortFunction abortFunction_;
};



#define WALBERLA_ABORT(msg) {\
   std::ostringstream WALBERLA__ABORT__OSS;\
   WALBERLA__ABORT__OSS << msg;\
   walberla::Abort::instance()->abort( WALBERLA__ABORT__OSS.str(), __FILE__, __LINE__ ); /* aborts execution, will never return */ \
   std::exit( EXIT_FAILURE ); /* must not be deleted, prevents "not all paths return a value" warnings! */ \
}

#define WALBERLA_ABORT_NO_DEBUG_INFO(msg) {\
   std::ostringstream WALBERLA__ABORT__OSS;\
   WALBERLA__ABORT__OSS << msg;\
   walberla::Abort::instance()->abortNoDebugInfo( WALBERLA__ABORT__OSS.str(), __FILE__, __LINE__ ); /* aborts execution, will never return */ \
   std::exit( EXIT_FAILURE ); /* must not be deleted, prevents "not all paths return a value" warnings! */ \
}



} // namespace walberla
