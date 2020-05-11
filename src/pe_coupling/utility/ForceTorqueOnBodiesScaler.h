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
//! \file ForceTorqueOnBodiesScaler.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"
#include "BodySelectorFunctions.h"

namespace walberla {
namespace pe_coupling {

// scales force/torque on all bodies (local and remote) by a constant scalar value
// can e.g. be used to average the force/torque over two time steps
class ForceTorqueOnBodiesScaler
{  
public:

   ForceTorqueOnBodiesScaler( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & bodyStorageID,
                              const real_t & scalingFactor, const std::function<bool(
            pe::BodyID)> &bodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ), scalingFactor_( scalingFactor ), bodySelectorFct_( bodySelectorFct )
     { }

   // resets forces and torques on all (local and remote) bodies
   void operator()();

   void resetScalingFactor( const real_t newScalingFactor );

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_;
   real_t scalingFactor_;
   const std::function<bool(pe::BodyID)> bodySelectorFct_;
};

} // namespace pe_coupling
} // namespace walberla
