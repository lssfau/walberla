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
//! \file LiftForceEvaluator.h
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
#include "pe/Materials.h"

#include "pe_coupling/geometry/SphereEquivalentDiameter.h"
#include "pe_coupling/utility/BodySelectorFunctions.h"

#include <functional>

#include <cmath>

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Evaluator of the lift force based on the given lift correlation.
 *
 * The lift force is evaluated for each body (requires interpolation of several quantities from the fluid fields)
 * and then applied onto them.
 * The corresponding reaction force is added to the fluid via the given distributor.
 *
 * Note that the fluid velocity, contained in the velocityField, has to be the fluid-phase velocity (and not the volume-averaged velocity).
 *
 * Whether or not a body gets treated by the evaluator depends on the return value of 'dpmBodySelectorFct'.
 *
 * For more infos on interpolators, see field interpolators in src/field/interpolators.
 * For more infos on distributors, see src/field/distributors.
 */
template< typename FlagField_T, template<typename,typename> class FieldInterpolator_T, template<typename,typename> class Distributor_T >
class LiftForceEvaluator
{  

public:

   using Vec3Field_T = GhostLayerField<Vector3<real_t>, 1>;
   using Vec3FieldInterpolator_T = FieldInterpolator_T<Vec3Field_T, FlagField_T>;
   using ForceDistributor_T = Distributor_T<Vec3Field_T, FlagField_T>;

   LiftForceEvaluator( const shared_ptr<StructuredBlockStorage> & blockStorage,
                       const BlockDataID & forceFieldID, const BlockDataID & bodyStorageID, const BlockDataID & flagFieldID, const Set< FlagUID > & domainFlags,
                       const BlockDataID & velocityFieldID, const BlockDataID & velocityCurlFieldID,
                       const std::function<Vector3<real_t> ( const Vector3<real_t> &, const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t, real_t )> & liftForceCorrelationFunction,
                       real_t fluidDynamicViscosity,
                       const std::function<bool(pe::BodyID)> & dpmBodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ), liftForceCorrelationFunction_( liftForceCorrelationFunction ),
     fluidDynamicViscosity_( fluidDynamicViscosity ), dpmBodySelectorFct_( dpmBodySelectorFct)
   {
      velocityFieldInterpolatorID_     = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blockStorage, velocityFieldID, flagFieldID, domainFlags );
      velocityCurlFieldInterpolatorID_ = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blockStorage, velocityCurlFieldID, flagFieldID, domainFlags );
      forceDistributorID_              = field::addDistributor< ForceDistributor_T, FlagField_T >( blockStorage, forceFieldID, flagFieldID, domainFlags );
   }

   void operator()();

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   BlockDataID bodyStorageID_;

   std::function<Vector3<real_t> ( const Vector3<real_t> &, const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t, real_t )> liftForceCorrelationFunction_;

   real_t fluidDynamicViscosity_;

   std::function<bool(pe::BodyID)> dpmBodySelectorFct_;

   BlockDataID velocityFieldInterpolatorID_;
   BlockDataID velocityCurlFieldInterpolatorID_;
   BlockDataID forceDistributorID_;
};

template< typename FlagField_T, template<typename,typename> class FieldInterpolator_T, template<typename,typename> class Distributor_T >
void LiftForceEvaluator< FlagField_T, FieldInterpolator_T, Distributor_T >
::operator()()
{

   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {

      Vec3FieldInterpolator_T * velocityInterpolator     = blockIt->getData< Vec3FieldInterpolator_T >( velocityFieldInterpolatorID_ );
      Vec3FieldInterpolator_T * velocityCurlInterpolator = blockIt->getData< Vec3FieldInterpolator_T >( velocityCurlFieldInterpolatorID_ );
      ForceDistributor_T * forceDistributor              = blockIt->getData< ForceDistributor_T >( forceDistributorID_ );

      for( auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
      {
         if(!dpmBodySelectorFct_(bodyIt.getBodyID())) continue;

         Vector3<real_t> forceOnFluid( real_t(0) );

         Vector3<real_t> bodyPosition = bodyIt->getPosition();
         Vector3<real_t> bodyVelocity = bodyIt->getLinearVel();

         real_t fluidDensity( real_t(1) );
         real_t bodyDiameter = getSphereEquivalentDiameter( *bodyIt );

         // interpolate fluid velocity and fluid curl to body position
         Vector3<real_t> fluidVelocity( real_t(0) );
         velocityInterpolator->get( bodyPosition, &fluidVelocity );

         Vector3<real_t> velocityCurl( real_t(0) );
         velocityCurlInterpolator->get( bodyPosition, &velocityCurl );

         // evaluate lift force according to empirical model
         Vector3<real_t> liftForce = liftForceCorrelationFunction_(fluidVelocity, velocityCurl, bodyVelocity, bodyDiameter, fluidDynamicViscosity_, fluidDensity );
         WALBERLA_ASSERT( !math::isnan(liftForce[0]) && !math::isnan(liftForce[1]) && !math::isnan(liftForce[2]),
                          "NaN found in lift force for body at position " << bodyPosition );
         bodyIt->addForce( liftForce );
         forceOnFluid += ( -liftForce );

         // set/distribute force on fluid
         forceDistributor->distribute( bodyPosition, &forceOnFluid );
      }
   }

}


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
