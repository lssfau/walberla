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
//! \file BodiesForceAndTorqueContainer.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "pe/rigidbody/BodyIterators.h"
#include "pe/synchronization/SyncForces.h"

#include <map>
#include <vector>

namespace walberla {
namespace pe_coupling {

class BodiesForceAndTorqueContainer
{  
public:

   BodiesForceAndTorqueContainer( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & bodyStorageID )
   : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID )
     { }

   void operator()()
   {
      store();
   }

   void store()
   {
      // clear map
      clear();

      // sum up all forces/torques from shadow copies on local body (owner)
      pe::reduceForces( blockStorage_->getBlockStorage(), bodyStorageID_ );
      // send total forces/torques to shadow owners
      pe::distributeForces( blockStorage_->getBlockStorage(), bodyStorageID_ );

      // (re-)build map
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            auto & f = bodyForceAndTorqueMap_[ bodyIt->getSystemID() ];

            // only add if body has not been added already before (from another block)
            if( f.empty() )
            {
               const auto & force = bodyIt->getForce();
               f.push_back( force[0] );
               f.push_back( force[1] );
               f.push_back( force[2] );

               const auto & torque = bodyIt->getTorque();
               f.push_back( torque[0] );
               f.push_back( torque[1] );
               f.push_back( torque[2] );
            }

            // reset of force/torque on remote bodies necessary to erase the multiple occurrences of forces/torques on bodies
            // (due to call to distributeForces() before)
            if ( bodyIt->isRemote() )
            {
               bodyIt->resetForceAndTorque();
            }
         }
      }

   }

   void setOnBodies()
   {
      // owning process sets the force/torque on the bodies
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            const auto &f = bodyForceAndTorqueMap_[bodyIt->getSystemID()];
            bodyIt->addForce ( f[0], f[1], f[2] );
            bodyIt->addTorque( f[3], f[4], f[5] );
         }
      }
   }

   void clear()
   {
      bodyForceAndTorqueMap_.clear();
   }

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_;
   std::map< walberla::id_t, std::vector<real_t> > bodyForceAndTorqueMap_;
};

} // namespace pe_coupling
} // namespace walberla
