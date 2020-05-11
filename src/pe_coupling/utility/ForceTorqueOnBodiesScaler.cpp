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
//! \file ForceTorqueOnBodiesScaler.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "ForceTorqueOnBodiesScaler.h"

#include "core/math/Vector3.h"
#include "pe/rigidbody/BodyIterators.h"

namespace walberla {
namespace pe_coupling {

void ForceTorqueOnBodiesScaler::operator()()
{
   Vector3<real_t> force(real_t(0));
   Vector3<real_t> torque(real_t(0));
   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {
      for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( !bodySelectorFct_(bodyIt.getBodyID()) ) continue;
         force  = scalingFactor_ * bodyIt->getForce();
         torque = scalingFactor_ * bodyIt->getTorque();

         bodyIt->resetForceAndTorque();

         bodyIt->setForce ( force );
         bodyIt->setTorque( torque );
      }
   }
}

void ForceTorqueOnBodiesScaler::resetScalingFactor( const real_t newScalingFactor )
{
   scalingFactor_ = newScalingFactor;
}

} // namespace pe_coupling
} // namespace walberla
