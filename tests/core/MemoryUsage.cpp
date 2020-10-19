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
//! \file MemoryUsage.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/MemoryUsage.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"

#include <vector>

using namespace walberla;

int main( int argc, char** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   debug::enterTestMode();

   std::vector<int> stuff(200*1024*1024/sizeof(int),1);

   auto memSize = getResidentMemorySize();
#ifdef __linux__
   WALBERLA_CHECK_GREATER(memSize, 200000);
   WALBERLA_CHECK_LESS   (memSize, 250000);
#else
   WALBERLA_CHECK_EQUAL(memSize, 0);
#endif
   printResidentMemoryStatistics();
   WALBERLA_LOG_DEVEL(stuff[0]); //make sure stuff is not optimized
   return EXIT_SUCCESS;
}
