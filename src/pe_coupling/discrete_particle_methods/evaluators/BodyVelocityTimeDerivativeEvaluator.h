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
//! \file BodyVelocityTimeDerivativeEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "pe/rigidbody/BodyIterators.h"

#include "pe_coupling/utility/BodySelectorFunctions.h"

#include <functional>
#include <map>

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Evaluates the time derivative of the body velocity
 *
 * Using the operator()(), an internal map is filled with pairs of the body ID and its current velocity,
 * after clearing the old map.
 * Calling get(..), the body's former velocity is fetched from the map and used to calculate a simple approximation of the
 * body velocity time derivative (du/dt = ( u_new - u_old ) / deltaT )
 *
 * Whether or not a body gets treated by the evaluator depends on the return value of 'dpmBodySelectorFct'.
 *
 */
class BodyVelocityTimeDerivativeEvaluator
{  
public:

   BodyVelocityTimeDerivativeEvaluator( const shared_ptr<StructuredBlockStorage> & blockStorage,
                                        const BlockDataID & bodyStorageID, const real_t & deltaT = real_t(1),
                                        const std::function<bool(pe::BodyID)> & dpmBodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ), deltaTinv_( real_t(1) / deltaT ),
     dpmBodySelectorFct_( dpmBodySelectorFct)
     { }

   void operator()()
   {
      // rebuild velocity map for all known bodies (process local)
      bodyVelocityMap_.clear();
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            if(!dpmBodySelectorFct_(bodyIt.getBodyID())) continue;

            bodyVelocityMap_.insert( std::pair<walberla::id_t, Vector3< real_t > >( bodyIt->getSystemID(), bodyIt->getLinearVel() ) );
         }
      }
   }

   void resetDeltaT( const real_t & deltaT )
   {
      deltaTinv_ = real_t(1) / deltaT;
   }

   void get( Vector3<real_t> & particleVelocityTimeDerivative, const Vector3<real_t> currentParticleVelocity, const walberla::id_t & bodySystemID )
   {
      auto it = bodyVelocityMap_.find( bodySystemID );
      WALBERLA_ASSERT( it != bodyVelocityMap_.end(), "body with ID " << bodySystemID << " not found in body velocity map!" );
      particleVelocityTimeDerivative = ( currentParticleVelocity - it->second ) * deltaTinv_;
   }

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_;
   std::map< walberla::id_t, Vector3< real_t > > bodyVelocityMap_;
   real_t deltaTinv_;
   std::function<bool(pe::BodyID)> dpmBodySelectorFct_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
