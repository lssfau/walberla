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
//! \file SubgridTimeStep.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/timing/TimingTree.h"

#include "domain_decomposition/BlockStorage.h"

#include "field/GhostLayerField.h"

#include "pe/cr/ICR.h"
#include "pe/rigidbody/BodyIterators.h"
#include "pe/synchronization/SyncForces.h"

#include <functional>

#include <map>

/*!\brief Carries out the the PE time steps, including sub iteration functionality.
 *
 * This class is similar to src/pe_coupling/utility/TimeStep.
 *
 * It executes 'intermediateSteps' PE steps within one timestep of size 'timeStepSize'.
 *
 * These PE sub iterations require, that the current external (e.g. hydrodynamic, gravitational,..) forces and torques
 * acting on each particle remains unchanged. Thus, a map is set up internally that stores and re-sets these forces
 * and torques in each PE sub iteration
 * .
 * Additionally, a function 'forceEvaluationFunc' can be given that allows to evaluate different forces before the PE
 * step is carried out. An example are particle-particle lubrication forces that have to be updated in each sub iteration.
 *
 */
namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

class SubgridTimeStep
{
public:

   typedef field::GhostLayerField< Vector3<real_t>, 1 > ForceField_T;

   explicit SubgridTimeStep( const shared_ptr<StructuredBlockStorage> & blockStorage,
                             const BlockDataID & bodyStorageID,
                             pe::cr::ICR & collisionResponse,
                             const std::function<void (void)> & synchronizeFunc,
                             const std::function<void ()> & forceEvaluationFunc,
                             const real_t timeStepSize = real_t(1), const uint_t intermediateSteps = uint_t(1) )
   : timeStepSize_( timeStepSize )
   , intermediateSteps_( ( intermediateSteps == 0 ) ? uint_t(1) : intermediateSteps )
   , blockStorage_( blockStorage )
   , bodyStorageID_(bodyStorageID)
   , collisionResponse_(collisionResponse)
   , synchronizeFunc_(synchronizeFunc)
   , forceEvaluationFunc_(forceEvaluationFunc)
   {}

   void operator()()
   {
      if( intermediateSteps_ == 1 )
      {
         forceEvaluationFunc_();

         collisionResponse_.timestep( timeStepSize_ );
         synchronizeFunc_();
      }
      else
      {
         // sum up all forces from shadow copies on local body (owner)
         pe::reduceForces(blockStorage_->getBlockStorage(), bodyStorageID_);
         // send total forces to shadow owners
         pe::distributeForces(blockStorage_->getBlockStorage(), bodyStorageID_);

         // store force/torques. Those should only contain external forces (gravity, buoyancy,..)
         // generate map from all known bodies (process local) to total forces
         std::map< walberla::id_t, std::vector< real_t > > forceMap;
         for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
         {
            for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
            {
               auto & f = forceMap[ bodyIt->getSystemID() ];

               const auto & force = bodyIt->getForce();
               f.push_back( force[0] );
               f.push_back( force[1] );
               f.push_back( force[2] );

               const auto & torque = bodyIt->getTorque();
               f.push_back( torque[0] );
               f.push_back( torque[1] );
               f.push_back( torque[2] );

               if ( bodyIt->isRemote() )
               {
                  bodyIt->resetForceAndTorque();
               }
            }
         }

         // perform pe sub time steps
         const real_t subTimeStepSize = timeStepSize_ / real_c( intermediateSteps_ );
         for( uint_t i = 0; i != intermediateSteps_; ++i )
         {
            // evaluate forces (e.g. lubrication forces)
            forceEvaluationFunc_();

            // in the first set forces on local bodies are already set by force synchronization
            if( i != 0 ) {
               for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
               {
                  for( auto body = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); body != pe::LocalBodyIterator::end(); ++body )
                  {
                     const auto & f = forceMap[ body->getSystemID() ];
                     body->addForce ( f[0], f[1], f[2] );
                     body->addTorque( f[3], f[4], f[5] );
                  }
               }
            }

            collisionResponse_.timestep( subTimeStepSize );
            synchronizeFunc_();
         }
      }
   }

protected:

   const real_t timeStepSize_;

   const uint_t intermediateSteps_;

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID &  bodyStorageID_;

   pe::cr::ICR & collisionResponse_;
   std::function<void (void)> synchronizeFunc_;
   std::function<void ()> forceEvaluationFunc_;

};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
