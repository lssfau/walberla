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
//! \file BodyVelocityInitializer.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/interpolators/all.h"
#include "field/distributors/all.h"
#include "field/GhostLayerField.h"

#include "pe/rigidbody/BodyIterators.h"
#include "pe/rigidbody/RigidBody.h"
#include "pe/Types.h"

#include "pe_coupling/utility/BodySelectorFunctions.h"

#include <functional>

#include <cmath>

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Initializes the bodies with the velocity of the surrounding fluid.
 *
 * The class uses an interpolator to obtain the approximate value of the fluid velocity at the bodies' position,
 * and then sets the bodies' velocity accordingly.
 *
 * Whether or not a body gets treated by this initializer depends on the return value of 'dpmBodySelectorFct'.
 *
 * For more infos on interpolators, see field interpolators in src/field/interpolators.
 */

template< typename FlagField_T, template<typename,typename> class FieldInterpolator_T >
class BodyVelocityInitializer
{  

public:

   typedef GhostLayerField< Vector3<real_t>, 1>          Vec3Field_T;
   typedef FieldInterpolator_T<Vec3Field_T, FlagField_T> Vec3FieldInterpolator_T;

   BodyVelocityInitializer( const shared_ptr<StructuredBlockStorage> & blockStorage,
                            const BlockDataID & bodyStorageID, const BlockDataID & flagFieldID, const Set< FlagUID > & domainFlags,
                            const BlockDataID & velocityFieldID,
                            const std::function<bool(pe::BodyID)> & dpmBodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ),
     dpmBodySelectorFct_( dpmBodySelectorFct)
   {
      velocityFieldInterpolatorID_ = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blockStorage, velocityFieldID, flagFieldID, domainFlags );
   }

   void operator()()
   {
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {

         Vec3FieldInterpolator_T * velocityInterpolator = blockIt->getData< Vec3FieldInterpolator_T >( velocityFieldInterpolatorID_ );

         for( auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            if(!dpmBodySelectorFct_(bodyIt.getBodyID())) continue;

            Vector3<real_t> bodyPosition = bodyIt->getPosition();

            // interpolate fluid velocity to body position
            Vector3<real_t> fluidVelocity( real_t(0) );
            velocityInterpolator->get( bodyPosition, &fluidVelocity );

            WALBERLA_ASSERT( !math::isnan(fluidVelocity[0]) && !math::isnan(fluidVelocity[1]) && !math::isnan(fluidVelocity[2]),
                             "NaN found in interpolated fluid velocity at position " << bodyPosition );

            bodyIt->setLinearVel( fluidVelocity );
         }
      }
   }

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   BlockDataID bodyStorageID_;
   BlockDataID velocityFieldInterpolatorID_;
   std::function<bool(pe::BodyID)> dpmBodySelectorFct_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
