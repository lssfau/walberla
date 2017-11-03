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
//! \file TimeStep.h
//! \ingroup pe_coupling
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/timing/TimingTree.h"

#include "domain_decomposition/BlockStorage.h"

#include "pe/cr/ICR.h"
#include "pe/rigidbody/BodyIterators.h"
#include "pe/synchronization/SyncForces.h"

#include <boost/function.hpp>

#include <map>


namespace walberla {
namespace pe_coupling {

class TimeStep
{
public:

   explicit TimeStep( const shared_ptr<StructuredBlockStorage> & blockStorage,
                      const BlockDataID & bodyStorageID,
                      pe::cr::ICR & collisionResponse,
                      boost::function<void (void)> synchronizeFunc,
                      const real_t timeStep = real_t(1), const uint_t intermediateSteps = uint_c(1) )
         : timeStep_( timeStep )
         , intermediateSteps_( ( intermediateSteps == 0 ) ? uint_c(1) : intermediateSteps )
         , blockStorage_( blockStorage )
         , bodyStorageID_(bodyStorageID)
         , collisionResponse_(collisionResponse)
         , synchronizeFunc_(synchronizeFunc)
   {}

   void operator()()
   {
      if( intermediateSteps_ == 1 )
      {
         collisionResponse_.timestep( timeStep_ );
         synchronizeFunc_( );
      }
      else
      {
         // during the intermediate time steps of the collision response, the currently acting forces
         // (interaction forces, gravitational force, ...) have to remain constant.
         // Since they are reset after the call to collision response, they have to be stored explicitly before.
         // Then they are set again after each intermediate step.

         // sum up all forces/torques from shadow copies on local body (owner)
         pe::reduceForces( blockStorage_->getBlockStorage(), bodyStorageID_ );
         // send total forces/torques to shadow owners
         pe::distributeForces( blockStorage_->getBlockStorage(), bodyStorageID_ );

         // generate map from all known bodies (process local) to total forces
         std::map< walberla::id_t, std::vector< real_t > > forceMap;
         for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
         {
            // iterate over local and remote bodies and store force/torque in map
            // Remote bodies are required since during the then following collision response time steps,
            // bodies can change owner ship (i.e. a now remote body could become a local body ).
            // Since only the owning process re-sets the force/torque later, remote body values have to be stores as well.
            for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
            {
               auto & f = forceMap[ bodyIt->getSystemID() ];

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

         // perform pe time steps
         const real_t subTimestep = timeStep_ / real_c( intermediateSteps_ );
         for( uint_t i = 0; i != intermediateSteps_; ++i )
         {
            // in the first set forces on local bodies are already set by force synchronization
            if( i != 0 )
            {
               for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
               {
                  // only owning process sets force/torque on bodies
                  for( auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
                  {
                     const auto & f = forceMap[ bodyIt->getSystemID() ];
                     WALBERLA_ASSERT( !f.empty(), "When attempting to set force/torque on local body " << bodyIt->getSystemID() << " at position " << bodyIt->getPosition() << ", body was not found in map!");
                     bodyIt->addForce ( f[0], f[1], f[2] );
                     bodyIt->addTorque( f[3], f[4], f[5] );
                  }
               }
            }

            collisionResponse_.timestep( subTimestep );
            synchronizeFunc_( );
         }
      }
   }

protected:

   const real_t timeStep_;

   const uint_t intermediateSteps_;

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID &  bodyStorageID_;

   pe::cr::ICR & collisionResponse_;
   boost::function<void (void)> synchronizeFunc_;

}; // class TimeStep



} // namespace pe_coupling
} // namespace walberla
