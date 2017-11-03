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

#include "core/math/Vector3.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "pe/rigidbody/BodyIterators.h"


namespace walberla {
namespace pe_coupling {

// scales force/torquew on all bodies (local and remote) by a constant scalar value
// can e.g. be used to average the force/torque over two time steps
class ForceTorqueOnBodiesScaler
{  
public:

   ForceTorqueOnBodiesScaler( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & bodyStorageID,
                              const real_t & scalingFactor )
   : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ), scalingFactor_( scalingFactor )
     { }

   // resets forces and torques on all (local and remote) bodies
   void operator()()
   {
      Vector3<real_t> force(real_t(0));
      Vector3<real_t> torque(real_t(0));
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            force  = scalingFactor_ * bodyIt->getForce();
            torque = scalingFactor_ * bodyIt->getTorque();

            bodyIt->resetForceAndTorque();

            bodyIt->setForce ( force );
            bodyIt->setTorque( torque );
         }
      }
   }


private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_;
   const real_t scalingFactor_;
};

} // namespace pe_coupling
} // namespace walberla
