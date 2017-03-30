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
//! \file Version.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "waLBerlaDefinitions.h"
#include "core/debug/TestSubsystem.h"

#include <sstream>

using namespace walberla;

int main( int /*argc*/, char** /*argv*/ )
{
   debug::enterTestMode();   
   
   int major       = WALBERLA_MAJOR_VERSION;
   int patch_level = WALBERLA_PATCH_LEVEL;
   
   WALBERLA_CHECK_GREATER_EQUAL( major,       0 );  
   WALBERLA_CHECK_GREATER_EQUAL( patch_level, 0 );  
   
   std::ostringstream oss;
   oss << major  << '.' << patch_level;
   WALBERLA_CHECK_EQUAL( oss.str(), WALBERLA_VERSION_STRING );
   
   int version = major * 100 + patch_level;   
   int calced_version = WALBERLA_VERSION_CALC( major, patch_level );
   
   WALBERLA_CHECK_EQUAL( version, calced_version );
   WALBERLA_CHECK_EQUAL( version, WALBERLA_VERSION );
   
   bool result = WALBERLA_VERSION_COMPARE( ==, major, patch_level );
   WALBERLA_CHECK( result );
   
   result = WALBERLA_VERSION_COMPARE( <, major + 1, patch_level );
   WALBERLA_CHECK( result );
   
   result = WALBERLA_VERSION_COMPARE( <, major, patch_level + 1 );
   WALBERLA_CHECK( result );
   
   result = WALBERLA_VERSION_COMPARE( <, major - 1, patch_level );
   WALBERLA_CHECK( !result );
   
   result = WALBERLA_VERSION_COMPARE( <, major, patch_level - 1 );
   WALBERLA_CHECK( !result );
}