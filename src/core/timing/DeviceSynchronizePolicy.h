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
//! \file DeviceSynchronizePolicy.h
//! \ingroup core
//! \author Richard Angersbach
//! \brief Gpu Timing Policy
//
//======================================================================================================================

#pragma once

#include "gpu/DeviceWrapper.h"

#include "Time.h"

namespace walberla
{
namespace timing
{

//======================================================================================================================
//
//  CLASS DEFINITION
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Timing policy for the measurement of the GPU time.
// \ingroup timing
//
// The DeviceSynchronizePolicy class represents the timing policy for GPU time measurements that can be used
// in combination with the Timer class template. This combination is realized with the DeviceSynchronizePolicy
// type definition.
// This class uses device synchronization internally and is therefore not suited for CUDA
// applications with overlapping kernels.
*/
struct DeviceSynchronizePolicy
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
/*!\brief Returns a timestamp of the current GPU time in seconds. Uses wall clock time and device synchronization
internally.
//
// \return GPU timestamp in seconds.
*/
inline double DeviceSynchronizePolicy::getTimestamp()
{
   // synchronize device before getting timestamp
   WALBERLA_DEVICE_SECTION() { gpuDeviceSynchronize(); }

   return getWcTime();
}
//**********************************************************************************************************************

} // namespace timing
} // namespace walberla
