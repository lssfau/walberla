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
//! \file TimeStep.cpp
//! \ingroup pe_coupling
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "TimeStep.h"

#include "pe/rigidbody/BodyIterators.h"

#include <map>
#include <array>

namespace walberla {
namespace pe_coupling {

void TimeStep::operator()()
{
   if( numberOfSubIterations_ == 1 )
   {
      forceEvaluationFunc_();

      collisionResponse_.timestep( timeStepSize_ );
      synchronizeFunc_();
   }
   else
   {
      // during the intermediate time steps of the collision response, the currently acting forces
      // (interaction forces, gravitational force, ...) have to remain constant.
      // Since they are reset after the call to collision response, they have to be stored explicitly before.
      // Then they are set again after each intermediate step.

      // generate map from all known bodies (process local) to total forces/torques
      // this has to be done on a block-local basis, since the same body could reside on several blocks from this process
      using BlockID_T = domain_decomposition::IBlockID::IDType;
      std::map< BlockID_T, std::map< walberla::id_t, std::array< real_t, 6 > > > forceTorqueMap;

      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         BlockID_T blockID = blockIt->getId().getID();
         auto& blockLocalForceTorqueMap = forceTorqueMap[blockID];

         // iterate over local and remote bodies and store force/torque in map
         for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            auto & f = blockLocalForceTorqueMap[ bodyIt->getSystemID() ];

            const auto & force = bodyIt->getForce();
            const auto & torque = bodyIt->getTorque();

            f = {{force[0], force[1], force[2], torque[0], torque[1], torque[2] }};
         }
      }

      // perform pe time steps
      const real_t subTimeStepSize = timeStepSize_ / real_c( numberOfSubIterations_ );
      for( uint_t i = 0; i != numberOfSubIterations_; ++i )
      {

         // in the first iteration, forces are already set
         if( i != 0 )
         {
            for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
            {
               BlockID_T blockID = blockIt->getId().getID();
               auto& blockLocalForceTorqueMap = forceTorqueMap[blockID];

               // re-set stored force/torque on bodies
               for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
               {

                  const auto f = blockLocalForceTorqueMap.find( bodyIt->getSystemID() );

                  if( f != blockLocalForceTorqueMap.end() )
                  {
                     const auto & ftValues = f->second;
                     bodyIt->addForce ( ftValues[0], ftValues[1], ftValues[2] );
                     bodyIt->addTorque( ftValues[3], ftValues[4], ftValues[5] );
                  }
               }
            }
         }

         // evaluate forces (e.g. lubrication forces)
         forceEvaluationFunc_();

         collisionResponse_.timestep( subTimeStepSize );
         synchronizeFunc_();
      }
   }
}


} // namespace pe_coupling
} // namespace walberla
