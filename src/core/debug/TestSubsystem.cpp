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
//! \file TestSubsystem.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Implementation of common Test functions
//
//======================================================================================================================

#include "TestSubsystem.h"

#ifdef _MSC_VER
#  include <Windows.h>
#  include <crtdbg.h>
#endif // _MSC_VER



namespace walberla {
namespace debug {



/*******************************************************************************************************************//**
 *
 * \brief   Enters test mode. Call this at the beginning of every tests's main() function!
 *
 * This function may be used to set up error handling for tests. (Or anything else for that
 * matter...)
 *
 **********************************************************************************************************************/
void enterTestMode()
{
#  ifdef _MSC_VER
      if( !IsDebuggerPresent() )
      {
         // Send all reports to STDERR instead of opening an Dialog Box on Windows
         _CrtSetReportMode( _CRT_WARN,   _CRTDBG_MODE_FILE   );
         _CrtSetReportFile( _CRT_WARN,   _CRTDBG_FILE_STDERR );
         _CrtSetReportMode( _CRT_ERROR,  _CRTDBG_MODE_FILE   );
         _CrtSetReportFile( _CRT_ERROR,  _CRTDBG_FILE_STDERR );
         _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE   );
         _CrtSetReportFile( _CRT_ASSERT, _CRTDBG_FILE_STDERR );
      }
#  endif // _MSC_VER
}



} // namespace debug
} // namespace walberla
