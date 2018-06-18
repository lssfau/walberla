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
//! \file BodiesForceTorqueContainer.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "core/math/Vector3.h"
#include "pe/rigidbody/BodyIterators.h"
#include "pe/synchronization/SyncForces.h"

#include <map>
#include <array>

namespace walberla {
namespace pe_coupling {

class BodiesForceTorqueContainer
{  
public:

   typedef std::map< walberla::id_t, std::array<real_t,6> > ForceTorqueStorage_T;

   BodiesForceTorqueContainer( const shared_ptr<StructuredBlockForest> & blockForest, const BlockDataID & bodyStorageID )
   : blockForest_( blockForest ), bodyStorageID_( bodyStorageID )
   {
      // has to be added to the forest (not the storage) to register correctly
      bodyForceTorqueStorageID_ = blockForest->addBlockData(make_shared<blockforest::AlwaysCreateBlockDataHandling<ForceTorqueStorage_T> >(), "BodiesForceTorqueContainer");
   }

   void operator()()
   {
      store();
   }

   void store()
   {
      // clear map
      clear();

      // (re-)build map
      for( auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt )
      {
         auto bodyForceTorqueStorage = blockIt->getData<ForceTorqueStorage_T>(bodyForceTorqueStorageID_);

         for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            auto & f = (*bodyForceTorqueStorage)[ bodyIt->getSystemID() ];

            const auto & force = bodyIt->getForce();
            const auto & torque = bodyIt->getTorque();
            f = {{force[0], force[1], force[2], torque[0], torque[1], torque[2] }};
         }
      }

   }

   void setOnBodies()
   {
      // set the force/torque stored in the block-local map onto all bodies
      for( auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt )
      {
         auto bodyForceTorqueStorage = blockIt->getData<ForceTorqueStorage_T>(bodyForceTorqueStorageID_);

         for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
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

   void clear()
   {
      for( auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt )
      {
         auto bodyForceTorqueStorage = blockIt->getData<ForceTorqueStorage_T>(bodyForceTorqueStorageID_);
         bodyForceTorqueStorage->clear();
      }
   }

   void swap( BodiesForceTorqueContainer & other )
   {
      std::swap( bodyForceTorqueStorageID_, other.bodyForceTorqueStorageID_);
   }

private:

   shared_ptr<StructuredBlockStorage> blockForest_;
   const BlockDataID bodyStorageID_;
   BlockDataID bodyForceTorqueStorageID_;
};


class BodyContainerSwapper
{
public:
   BodyContainerSwapper( const shared_ptr<BodiesForceTorqueContainer> & cont1, const shared_ptr<BodiesForceTorqueContainer> & cont2 )
   : cont1_( cont1 ), cont2_( cont2 )
   { }

   void operator()()
   {
      cont1_->swap( *cont2_ );
   }

private:
   shared_ptr<BodiesForceTorqueContainer> cont1_;
   shared_ptr<BodiesForceTorqueContainer> cont2_;
};

} // namespace pe_coupling
} // namespace walberla
