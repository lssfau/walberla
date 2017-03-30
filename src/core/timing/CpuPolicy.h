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
//! \file CpuPolicy.h
//! \ingroup core
//! \author Klaus Iglberger
//! \brief CPU Timing Policy
//!
//! Copyright (C) 2009 Klaus Iglberger
//! Taken from "pe Physics Engine" with small changes
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
/*!\brief Timing policy for the measurement of the CPU time.
// \ingroup timing
//
// The CpuPolicy class represents the timing policy for CPU time measurements that can be used
// in combination with the Timer class template. This combination is realized with the CpuTimer
// type definition.
*/
struct CpuPolicy
{
 public:
   //**Timing functions****************************************************************************
   /*!\name Timing functions */
   //@{
   static inline double getTimestamp();
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************




//======================================================================================================================
//
//  TIMING FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Returns a timestamp of the current CPU time in seconds.
//
// \return CPU timestamp in seconds.
*/
inline double CpuPolicy::getTimestamp()
{
   return getCpuTime();
}
//**********************************************************************************************************************

} // namespace timing
} // namespace walberla
