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
//! \file InteractionForceEvaluator.h
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

#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/ForceModel.h"

#include "pe/rigidbody/BodyIterators.h"
#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/Types.h"
#include "pe/Materials.h"

#include "pe_coupling/geometry/SphereEquivalentDiameter.h"
#include "pe_coupling/utility/BodySelectorFunctions.h"

#include "stencil/Directions.h"

#include <functional>

#include <vector>

#include <cmath>

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Evaluator of the two most important interaction forces: drag and pressure gradient
 *
 * The evaluation of the drag force is based on the given drag correlation.
 * The pressure gradient force requires a valid pressure gradient field.
 * These two forces are evaluated for each body (requires interpolation of several quantities from the fluid fields)
 * and then applied onto them.
 * However, only the drag force (with negative sign) is applied (i.e. distributed) onto the fluid since it is assumed that
 * the pressure force is already implicitly included in the formulation of the equations that describe the fluid flow.
 *
 * Note that the fluid velocity, contained in the velocityField, has to be the fluid-phase velocity (and not the volume-averaged velocity).
 *
 * Whether or not a body gets treated by the evaluator depends on the return value of 'dpmBodySelectorFct'.
 *
 * For more infos on interpolators, see field interpolators in src/field/interpolators.
 * For more infos on distributors, see src/field/distributors.
 */

template< typename FlagField_T, template<typename,typename> class FieldInterpolator_T, template<typename,typename> class Distributor_T >
class InteractionForceEvaluator
{  

public:

   using Vec3Field_T = GhostLayerField<Vector3<real_t>, 1>;
   using ScalarField_T = GhostLayerField<real_t, 1>;
   using Vec3FieldInterpolator_T = FieldInterpolator_T<Vec3Field_T, FlagField_T>;
   using ScalarFieldInterpolator_T = FieldInterpolator_T<ScalarField_T, FlagField_T>;
   using ForceDistributor_T = Distributor_T<Vec3Field_T, FlagField_T>;

   InteractionForceEvaluator( const shared_ptr<StructuredBlockStorage> & blockStorage,
                              const BlockDataID & forceFieldID, const BlockDataID & bodyStorageID,
                              const BlockDataID & flagFieldID, const Set< FlagUID > & domainFlags,
                              const BlockDataID & velocityFieldID, const BlockDataID & solidVolumeFractionFieldID, const BlockDataID & pressureGradientFieldID,
                              const std::function<Vector3<real_t> ( const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t, real_t, real_t )> & dragForceCorrelationFunction,
                              real_t fluidDynamicViscosity,
                              const std::function<bool(pe::BodyID)> & dpmBodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ),
     dragForceCorrelationFunction_( dragForceCorrelationFunction ), fluidDynamicViscosity_( fluidDynamicViscosity ),
     dpmBodySelectorFct_( dpmBodySelectorFct)
   {
      velocityFieldInterpolatorID_            = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blockStorage, velocityFieldID, flagFieldID, domainFlags );
      solidVolumeFractionFieldInterpolatorID_ = field::addFieldInterpolator< ScalarFieldInterpolator_T, FlagField_T >( blockStorage, solidVolumeFractionFieldID, flagFieldID, domainFlags );
      pressureGradientFieldInterpolatorID_    = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blockStorage, pressureGradientFieldID, flagFieldID, domainFlags );
      forceDistributorID_                     = field::addDistributor< ForceDistributor_T, FlagField_T >( blockStorage, forceFieldID, flagFieldID, domainFlags );
   }

   void operator()();

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_ {};
   const BlockDataID pdfFieldID_ {};

   std::function<Vector3<real_t> ( const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t, real_t, real_t )> dragForceCorrelationFunction_;

   real_t fluidDynamicViscosity_;

   std::function<bool(pe::BodyID)> dpmBodySelectorFct_;

   BlockDataID velocityFieldInterpolatorID_;
   BlockDataID solidVolumeFractionFieldInterpolatorID_;
   BlockDataID pressureGradientFieldInterpolatorID_;
   BlockDataID forceDistributorID_;
};

template< typename FlagField_T, template<typename,typename> class FieldInterpolator_T, template<typename,typename> class Distributor_T >
void InteractionForceEvaluator< FlagField_T, FieldInterpolator_T, Distributor_T >
::operator()()
{

   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {

      Vec3FieldInterpolator_T * velocityInterpolator               = blockIt->getData< Vec3FieldInterpolator_T >( velocityFieldInterpolatorID_ );
      ScalarFieldInterpolator_T * solidVolumeFractionInterpolator  = blockIt->getData< ScalarFieldInterpolator_T >( solidVolumeFractionFieldInterpolatorID_ );
      Vec3FieldInterpolator_T * pressureGradientInterpolator       = blockIt->getData< Vec3FieldInterpolator_T >( pressureGradientFieldInterpolatorID_ );
      ForceDistributor_T * forceDistributor                        = blockIt->getData< ForceDistributor_T >( forceDistributorID_ );

      for( auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
      {
         if(!dpmBodySelectorFct_(bodyIt.getBodyID())) continue;

         Vector3<real_t> forceOnFluid( real_t(0) );

         Vector3<real_t> bodyPosition = bodyIt->getPosition();

         // interpolate fluid velocity to body position
         Vector3<real_t> fluidVelocity( real_t(0) );
         velocityInterpolator->get( bodyPosition, &fluidVelocity );

         // interpolate solid volume fraction to body position
         real_t solidVolumeFraction( real_t(0) );
         solidVolumeFractionInterpolator->get( bodyPosition, &solidVolumeFraction );

         WALBERLA_ASSERT_GREATER( solidVolumeFraction, real_t(0) );

         // evaluate drag force
         Vector3<real_t> bodyVelocity = bodyIt->getLinearVel();
         real_t bodyDiameter = getSphereEquivalentDiameter( *bodyIt );
         real_t bodyVolume = bodyIt->getVolume();
         real_t fluidDensity( real_t(1) );

         Vector3<real_t> dragForce = dragForceCorrelationFunction_( fluidVelocity, bodyVelocity, solidVolumeFraction, bodyDiameter, fluidDynamicViscosity_, fluidDensity );

         WALBERLA_ASSERT( !math::isnan(dragForce[0]) && !math::isnan(dragForce[1]) && !math::isnan(dragForce[2]),
                          "NaN found in drag force " << dragForce << " for body at position " << bodyPosition );

         bodyIt->addForce( dragForce );

         forceOnFluid += ( -dragForce );

         // evaluate pressure gradient force = - V * grad(p)
         Vector3<real_t> pressureGradient( real_t(0) );
         pressureGradientInterpolator->get( bodyPosition, &pressureGradient );
         Vector3<real_t> pressureGradientForce = -bodyVolume * pressureGradient;
         WALBERLA_ASSERT( !math::isnan(pressureGradientForce[0]) && !math::isnan(pressureGradientForce[1]) && !math::isnan(pressureGradientForce[2]),
                          "NaN found in pressure gradient force " << pressureGradientForce << " for body at position " << bodyPosition );
         bodyIt->addForce( pressureGradientForce );

         // set/distribute force on fluid
         forceDistributor->distribute( bodyPosition, &forceOnFluid );
      }
   }

}

} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
