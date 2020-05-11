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
//! \file BodiesForceTorqueContainer.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "BodiesForceTorqueContainer.h"

#include "blockforest/StructuredBlockForest.h"
#include "core/math/Vector3.h"
#include "pe/rigidbody/BodyIterators.h"

namespace walberla {
namespace pe_coupling {



void BodiesForceTorqueContainer::store()
{
   // clear map
   clear();

   // (re-)build map
   for( auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt )
   {
      auto bodyForceTorqueStorage = blockIt->getData<ForceTorqueStorage_T>(bodyForceTorqueStorageID_);

      for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( !bodySelectorFct_(bodyIt.getBodyID()) ) continue;
         auto & f = (*bodyForceTorqueStorage)[ bodyIt->getSystemID() ];

         const auto & force = bodyIt->getForce();
         const auto & torque = bodyIt->getTorque();
         f = {{force[0], force[1], force[2], torque[0], torque[1], torque[2] }};
      }
   }
}

void BodiesForceTorqueContainer::setOnBodies()
{
   // set the force/torque stored in the block-local map onto all bodies
   for( auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt )
   {
      auto bodyForceTorqueStorage = blockIt->getData<ForceTorqueStorage_T>(bodyForceTorqueStorageID_);

      for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( !bodySelectorFct_(bodyIt.getBodyID()) ) continue;
         const auto f = bodyForceTorqueStorage->find( bodyIt->getSystemID() );

         if( f != bodyForceTorqueStorage->end() )
         {
            const auto & ftValues = f->second;
            bodyIt->addForce ( ftValues[0], ftValues[1], ftValues[2] );
            bodyIt->addTorque( ftValues[3], ftValues[4], ftValues[5] );
         }
         // else: new body has arrived that was not known before
      }
   }
}

void BodiesForceTorqueContainer::clear()
{
   for( auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt )
   {
      auto bodyForceTorqueStorage = blockIt->getData<ForceTorqueStorage_T>(bodyForceTorqueStorageID_);
      bodyForceTorqueStorage->clear();
   }
}

void BodiesForceTorqueContainer::swap( BodiesForceTorqueContainer & other )
{
   std::swap( bodyForceTorqueStorageID_, other.bodyForceTorqueStorageID_);
}







} // namespace pe_coupling
} // namespace walberla
