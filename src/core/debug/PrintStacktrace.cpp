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
//! \file PrintStacktrace.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "PrintStacktrace.h"
#include "demangle.h"
#include "core/DataTypes.h"

#include <iostream>


namespace walberla {
namespace debug {

void printStacktrace()
{
   printStacktrace( std::cerr );
}

} // namespace debug
} // namespace walberla

#ifdef WALBERLA_BUILD_WITH_BACKTRACE

#include WALBERLA_BACKTRACE_HEADER
#include <cstdlib>
#include <string>

namespace walberla {
namespace debug {

   void printStacktrace( std::ostream & os )
   {
      using std::string;
      using std::cerr;
      using std::endl;

      const int BACKTRACE_LENGTH = 20;
      void * array[BACKTRACE_LENGTH];
      size_t size;
      char **strings;

      size = numeric_cast< size_t >( backtrace (array, BACKTRACE_LENGTH) );
      strings = backtrace_symbols (array, int_c(size) );

      os << "Backtrace: " << std::endl;

      for (size_t i = 0; i < size; i++)
      {
         std::string line ( strings[i] );
#ifdef __APPLE__
         // one line might look like this:
         //0   PrintStacktraceTest                 0x0000000000408c6b _ZN8walberla4core15printStacktraceEv + 75

         size_t plusPos       = line.find_last_of('+');
         size_t funcPos       = line.find_last_of(' ', plusPos-2)+1;

         string functionName = line.substr( funcPos, plusPos-funcPos-1 );

         size_t addrPos       = line.find_last_of('x', funcPos-2)-2;
         size_t afterLevelPos = line.find_first_of(' ');
         size_t beforeAppPos  = line.find_first_not_of(" ", afterLevelPos);

         string appName      = line.substr( beforeAppPos, addrPos-beforeAppPos );

         size_t afterAppPos   = appName.find_last_not_of(" ")+1;
         appName             = appName.substr( 0, afterAppPos );
#else
         // one line might look like this:
         //./PrintStacktraceTest(_ZN8walberla4core15printStacktraceEv+0x4b) [0x408c6b]

         // extract the portion in the brackets () and demangle it
         size_t leftBracket  = line.find_first_of( '('  );
         size_t rightBracket = line.find_first_of( ')', leftBracket );

         string appName     = line.substr( 0, leftBracket );
         string bracketPart = line.substr( leftBracket+1, rightBracket - leftBracket -1 );
         string rest        = line.substr( rightBracket +1 );

         // split the bracketPart on plus sign
         size_t plusPos = bracketPart.find_first_of('+');
         string functionName = bracketPart.substr(0, plusPos );
         string offset       = bracketPart.substr( plusPos+1 );
#endif

         // try to demangle -> no return code if successful
         // but the returned string starts with "demangle" if demangling failed
         // really hacky :(
         string demangled    = demangle( functionName );
         if ( demangled.substr(0, 8) == "demangle")
            demangled = functionName;

         os << "\t" << appName << " \t " << demangled << endl;
      }
      free (strings);
   }


} // namespace debug
} // namespace walberla

#elif defined(WALBERLA_CXX_COMPILER_IS_MSVC)

#pragma warning( push )
#pragma warning( disable : 4091 )
#include "extern/StackWalker.h"

namespace walberla {
namespace debug {

class WalberlaStackWalker : public stack_walker::StackWalker
{
public:
   WalberlaStackWalker( std::ostream & os ) : StackWalker(), os_( os ) { }

protected:
   virtual void OnOutput( LPCSTR szText ) { os_ << szText; }

   virtual void OnSymInit ( LPCSTR /*szSearchPath*/, DWORD /*symOptions*/, LPCSTR /*szUserName*/ ) { }
   virtual void OnLoadModule ( LPCSTR /*img*/, LPCSTR /*mod*/, DWORD64 /*baseAddr*/, DWORD /*size*/, DWORD /*result*/,
                               LPCSTR /*symType*/, LPCSTR /*pdbName*/, ULONGLONG /*fileVersion*/ ) { }
   virtual void OnDbgHelpErr ( LPCSTR /*szFuncName*/, DWORD /*gle*/, DWORD64 /*addr*/ ) { }

   std::ostream & os_;
};

void printStacktrace( std::ostream & os )
{
   WalberlaStackWalker sw( os );
   sw.ShowCallstack();
}

} // namespace debug
} // namespace walberla

#pragma warning( pop )

#else

namespace walberla {
namespace debug {

   void printStacktrace( std::ostream &  )
   {}

} // namespace debug
} // namespace walberla



#endif
