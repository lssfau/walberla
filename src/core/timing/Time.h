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

#include <chrono>


namespace walberla {
namespace timing {


//======================================================================================================================
//
//  TIME FUNCTIONS
//
//======================================================================================================================


//**********************************************************************************************************************
/*!\brief Returns the current wall clock time in seconds.
// \ingroup util
//
// \return The current wall clock time in seconds.
*/
inline double getWcTime()
{
   return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count()) * 1e-9;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the current CPU time in seconds.
// \ingroup util
//
// \return The current CPU time in seconds.
*/
double getCpuTime();
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
