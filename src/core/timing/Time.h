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

#pragma once

#include "core/DataTypes.h"

#include <ctime>
#include <sstream>
#include <string>

#include <boost/chrono/chrono.hpp>

#if defined(_MSC_VER)
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#  include <winsock.h>
#  include <time.h>
#  include <sys/timeb.h>
#else
#  include <sys/resource.h>
#  include <sys/time.h>
#  include <sys/types.h>
#endif


namespace walberla {
namespace timing {


//======================================================================================================================
//
//  TIME FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\name Time functions */
//@{
inline double      getWcTime();
inline double      getCpuTime();
//@}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the current wall clock time in seconds.
// \ingroup util
//
// \return The current wall clock time in seconds.
*/
inline double getWcTime()
{
#ifdef WIN32
   LARGE_INTEGER perfCounter, perfFrequency;
   QueryPerformanceCounter( &perfCounter );
   QueryPerformanceFrequency( &perfFrequency );
   return static_cast<double>(perfCounter.QuadPart) / static_cast<double>(perfFrequency.QuadPart);
#elif defined __bg__ // boost::chrono seems broken on BG/Q
   struct timeval tp;
   gettimeofday( &tp, NULL );
   return ( static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6 );
#else
   return static_cast<double>(boost::chrono::duration_cast<boost::chrono::nanoseconds>(boost::chrono::high_resolution_clock::now().time_since_epoch()).count()) * 1e-9;
#endif
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the current CPU time in seconds.
// \ingroup util
//
// \return The current CPU time in seconds.
*/
inline double getCpuTime()
{
#ifdef WIN32
   FILETIME CreateTime, ExitTime, KernelTime, UserTime;
   SYSTEMTIME SysTime;

   if( GetProcessTimes( GetCurrentProcess(), &CreateTime, &ExitTime, &KernelTime, &UserTime ) != 1 ) {
      return 0.0;
   }
   else {
      FileTimeToSystemTime( &UserTime, &SysTime );
      return ( static_cast<double>( SysTime.wSecond ) + static_cast<double>( SysTime.wMilliseconds )/1E3 );
   }
#else
   struct rusage ruse;
   getrusage( RUSAGE_SELF, &ruse );
   return ( static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6 );
#endif
}
//**********************************************************************************************************************



/*******************************************************************************************************************//**
 * \brief Converts a timespan to a string.
 *
 * \param time  The timespan in seconds.
 *
 * \return   A string of the form "3d 2h 12m 9s". Zero values are omitted.
 **********************************************************************************************************************/
inline std::string timeToString(real_t time)
{
   const unsigned int sPerMin = 60; // Seconds per Minute
   const unsigned int sPerH   = 60 * 60; // Seconds per Hour
   const unsigned int sPerD   = 60 * 60 * 24; // Seconds per Day

   std::stringstream result;

   real_t tmp;

   if(time > sPerD)
   {
      tmp = std::floor(time / sPerD);
      result << tmp << "d ";
      time -= tmp * sPerD;
   }

   if(time > sPerH)
   {
      tmp = std::floor(time / sPerH);
      result << tmp << "h ";
      time -= tmp * sPerH;
   }

   if(time > sPerMin)
   {
      tmp = std::floor(time / sPerMin);
      result << tmp << "m ";
      time -= tmp * sPerMin;
   }

   result << std::floor(time + real_t(0.5)) << "s ";

   return result.str();
}



} // namespace timing
} // namespace walberla
