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

#include "MemoryUsage.h"

#ifdef __linux__
#include <sys/resource.h>
#endif
#include "waLBerlaDefinitions.h"
#include "core/logging/Logging.h"
#include "core/math/Sample.h"

namespace walberla {

long getResidentMemorySize()
{
#ifdef __linux__
   struct rusage usage;
   int gru = getrusage(RUSAGE_SELF, &usage);

   if (gru!=0)
      WALBERLA_LOG_WARNING("There was an error getting memory statistics!");

   return usage.ru_maxrss;
#else
   WALBERLA_LOG_WARNING("Getting memory statistics is currently not supported on non-Linux systems!");
   return 0;
#endif
}

void printResidentMemoryStatistics()
{
   math::Sample memory;
   memory.castToRealAndInsert(getResidentMemorySize());
   memory.mpiGatherRoot();
   WALBERLA_LOG_INFO_ON_ROOT("resident memory: " << memory.format());
}

} // namespace walberla
