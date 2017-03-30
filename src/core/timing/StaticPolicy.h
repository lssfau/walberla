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
//! \file StaticPolicy.h
//! \ingroup core
//! \author Sebastian Eibl
//! \brief Static Clock Timing Policy
//
//======================================================================================================================

#pragma once

#include "Time.h"


namespace walberla {
namespace timing {

//======================================================================================================================
//
//  CLASS DEFINITION
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Timing policy for fake measurement with a static settable clock. This clock is always stopped!
// \ingroup timing
*/
struct StaticPolicy
{
private:
   static double time_;
public:
   /// Sets the clock to a certain point in time
   /// \param time Time in seconds
   static inline void setTime(double time);

   /// Adds a certain amount of seconds to the clock
   /// \param time Time in seconds
   static inline void addTime(double time);


   //**Timing functions****************************************************************************
   /*!\name Timing functions */
   //@{
   static inline double getTimestamp();
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************

inline void StaticPolicy::setTime(double time)
{
   time_ = time;
}

inline void StaticPolicy::addTime(double time)
{
   time_ += time;
}

//======================================================================================================================
//
//  TIMING FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Returns a timestamp of the current wall clock time in seconds.
//
// \return Wall clock timestamp in seconds.
*/
inline double StaticPolicy::getTimestamp()
{
   return time_;
}
//**********************************************************************************************************************

} // namespace timing
} // namespace walberla
