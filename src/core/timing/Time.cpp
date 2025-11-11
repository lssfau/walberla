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
//! \file Time.h
//! \ingroup core
//! \author Klaus Iglberger
//! \brief Header file for time functions
//!
//! Copyright (C) 2009 Klaus Iglberger
//! Taken from "pe Physics Engine" with small changes
//
//======================================================================================================================

#include "Time.h"

#include <chrono>

#  include <sys/resource.h>
#  include <sys/time.h>
#  include <sys/types.h>


namespace walberla {
namespace timing {

double getCpuTime()
{
   struct rusage ruse;
   getrusage( RUSAGE_SELF, &ruse );
   return ( static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6 );
}

} // namespace timing
} // namespace walberla